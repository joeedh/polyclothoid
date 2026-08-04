/**
 * The Gauss-Newton solver for {@link SPowerClothoid} — `docs/plans/spower-solver.md` §4–§6.
 *
 * The unknowns are vertex-owned curvature blocks (§3), so `Gᵖ` continuity is structural and
 * the only thing left to solve for is tangent continuity. That makes each step a small
 * equality-constrained quadratic program,
 *
 * ```
 *   minimize  ½ xᵀ H x     subject to  c(x) = 0
 * ```
 *
 * with `H` the bending-plus-foundation energy of §6 and `c` one wrapped angle gap per
 * interior joint. Written as a KKT system and ordered so the multipliers interleave with the
 * blocks, it is banded, and `math/banded.ts` factors it in `O(V·(p+1)³)`.
 *
 * ## Two Jacobians
 *
 * `KTH_e` and `L_e` are smooth functions of the DOF, and §5 is emphatic that *not*
 * differentiating through them is the wrong sparsity pattern rather than a small lag. Both
 * versions are here, selected by {@link SPowerSolverOptions.jacobian}: §12 phased the frozen
 * one first because it is the same assembly minus the derivative accumulators, so an
 * assembly bug and a Jacobian bug could not be confused for one another. It is kept because
 * "the wrong sparsity pattern" is a claim worth being able to re-measure.
 *
 * How wrong, measured: finite-differencing the true residual against a frozen row (§14) gives
 * `−0.53` where the truth is `−0.16`, and `0` where the truth is `−0.15` — off by a factor of
 * three on the entries it has, and missing the leaving edge's far block outright. Its visible
 * symptom is {@link SPowerSolverReport.unenforced}: with `KTH` frozen an edge reaches a joint
 * only through its *total turning*, hence only at the end the turning integral arrives at, so
 * one arrangement of edge directions leaves a joint with no row at all and it is dropped
 * rather than given a multiplier that runs off to `c/δ`. Its invisible symptom is that a
 * chain with more than one interior joint does not converge at all — the coupling between
 * neighbouring joints is exactly what the missing blocks carried.
 *
 * The exact row reaches every block both edges read, so no joint is ever unenforced, and the
 * step count stops depending on how the edges happened to be built: four, eight and six on
 * the three chains §14 tabulates, identical across all three orientations of each, against
 * frozen counts ranging from six to a hundred-and-stall.
 *
 * `run` still backtracks on §5's `ℓ1` merit rather than taking the full step. That is not
 * only insurance for the frozen mode — a Gauss-Newton step is exact only to first order, and
 * without a merit test the coefficients ratchet, reaching `1e+70` by the sixty-fourth pass
 * on a four-edge zigzag before the factorization hands back `NaN`.
 */
import { BandedSymmetric, type KKTSolveOptions, defaultKKTSolveOptions, rankRatio, solveKKT } from "../math/index.js";
import { type EdgeFrame, type ProfileDOF, referenceLength, sPowerDOF } from "./blocks.js";
import {
  type Diagnostic,
  type DiagnosticCondition,
  RATE_WINDOW,
  type SolveReport,
  type TraceStep,
  diagnosticThresholds,
  geometricRate,
  stabilityThresholds,
} from "./diagnostics.js";
import { type CurveSolver, type SolvableEdge, type SolvableMesh, type SolvableVertex } from "./mesh_types.js";
import { sharing } from "./pairing.js";
import { SPOWER_ORDER, type SPowerClothoid } from "./spower_clothoid.js";
import { MIN_CANONICAL_CHORD, type TransformJacobian, transformJacobian } from "./quadrature.js";
import { type Chain, chains } from "./topology.js";

const TAU = Math.PI * 2.0;

/** Backtracking line search on the merit — see {@link ChainSystem.run}. */
const SHRINK = 0.5;
const MIN_RELAXATION = 1e-4;

/** Signed angle in `(−π, π]`. §4: the residual has to carry a sign and stay differentiable. */
function wrapAngle(a: number) {
  return a - TAU * Math.round(a / TAU);
}

function infinityNorm(x: Float64Array) {
  let worst = 0.0;

  for (const v of x) {
    worst = Math.max(worst, Math.abs(v));
  }

  return worst;
}

export interface SPowerSolverOptions {
  /** The s-power order `p`. Curvature degree is `2p + 1` and each vertex owns `p + 1` scalars. */
  order: number;

  /**
   * Energy locality, §6's `ε` made dimensionless: `E_e = (K + α²M)/C_e³`.
   *
   * A bending energy on an elastic foundation, so an edit at one vertex decays like
   * `exp(−√(α/2)·u)` in units of the edge — `α = 0.1` is about four and a half edges. It is
   * also what makes `H` strictly positive definite rather than merely semi-definite: the
   * mass term removes the constant-curvature null space of `∫q′²`. Zero is therefore not a
   * legal value, and larger values trade shape fidelity for locality.
   */
  alpha: number;

  /** Outer Gauss-Newton steps. */
  iterations: number;

  /**
   * Longest step allowed, as a fraction of the Newton step. `1.0` is undamped.
   *
   * This is only a ceiling — {@link ChainSystem.run} backtracks below it whenever a step
   * fails to improve the measured angle gap, so there is rarely a reason to lower it.
   */
  relaxation: number;

  /** Stop once every joint's angle gap is below this, in radians. */
  tolerance: number;

  /** How DOF become coefficients — §10. Swap this for §12's control experiment. */
  dof: ProfileDOF;

  /**
   * Which constraint Jacobian to assemble.
   *
   * `"exact"` differentiates through `KTH_e` as §5 specifies. `"frozen"` is Phase 2's
   * version, kept because §12 asks for the iteration counts side by side and because it is
   * the only thing that makes the *cost* of the derivative accumulators measurable. It does
   * not solve chains with more than one interior joint — see §14.
   */
  jacobian: "exact" | "frozen";

  /**
   * Record a {@link SolveTrace} per chain — §8's per-iteration history.
   *
   * Off by default because it is the one part of the diagnostics that allocates per step.
   * Everything in {@link SPowerSolverReport.diagnostics} is measured either way.
   */
  trace: boolean;

  /**
   * What to do about a joint the solve cannot hold — §9's split, as a mode on the solve
   * rather than a decision per joint.
   *
   * `"artistic"` degrades and reports: the joint loses a continuity level and the report says
   * which, where and why. A stroke that renders slightly wrong beats one that does not render.
   * `"engineering"` refuses and reports: detection is identical, the geometry is left where
   * the last attempt put it, and {@link SolveReport.ok} goes false. Silently delivering less
   * than was asked for is what neither mode is allowed to do.
   */
  mode: "artistic" | "engineering";

  /**
   * How many levels one chain may lose in a single solve, across all its joints.
   *
   * Three walks a `p = 1` joint the whole way from `G3` to a corner, so the default cannot
   * run out of room on a chain; it is a bound on thrash, not on depth.
   */
  breaks: number;

  kkt: Partial<KKTSolveOptions>;
}

export const defaultSPowerSolverOptions: SPowerSolverOptions = {
  order     : SPOWER_ORDER,
  alpha     : 0.1,
  iterations: 100,
  relaxation: 1.0,
  tolerance : 1e-10,
  dof       : sPowerDOF,
  jacobian  : "exact",
  trace     : false,
  mode      : "artistic",
  breaks    : 3,
  kkt       : {},
};

export interface SPowerSolverReport extends SolveReport {
  /** Chains found by {@link chains}, including any that were skipped. */
  chains: number;

  /** Chains skipped for closing on themselves — §5's cut point, deferred to Phase 8. */
  skippedClosed: number;

  /**
   * Interior joints whose G1 row vanished, and which were therefore left unconstrained.
   *
   * Only reachable with the Jacobian frozen; the exact row never vanishes.
   */
  unenforced: number;

  /** Largest `|wrapped angle gap|` over every enforced joint, measured after the last step. */
  maxResidual: number;

  /** Gauss-Newton steps actually taken, summed over chains. */
  steps: number;

  /**
   * Joints delivered below the level they were authored with — §8's breaks, counted.
   *
   * Zero in engineering mode by construction: it refuses rather than degrades. Every one of
   * these has a `"degraded"` record in {@link SolveReport.diagnostics} saying where and why.
   */
  degraded: number;

  /** Chain solves run, including the attempts thrown away. `chains − skippedClosed` when nothing broke. */
  attempts: number;
}

function emptyReport(): SPowerSolverReport {
  return {
    chains       : 0,
    skippedClosed: 0,
    unenforced   : 0,
    maxResidual  : 0.0,
    steps        : 0,
    degraded     : 0,
    attempts     : 0,
    ok           : true,
    diagnostics  : [],
  };
}

/**
 * What one chain's {@link ChainSystem.run} measured on its way through, for {@link
 * ChainSystem.diagnose} to turn into §8 records.
 *
 * All of it is a byproduct: the history is the convergence test, the backtracks and the
 * relaxation floor are the line search, the refinement residual comes back from `solveKKT`,
 * and the multiplier is read off the solution vector. §8's second rule is what makes the
 * unconditional measurement of the first rule affordable.
 */
export interface ChainRun {
  steps: number;
  residual: number;

  /** False if a KKT factorization underflowed. */
  factored: boolean;

  /**
   * False if the chain produced nothing usable — a failed factorization, or geometry that
   * went non-finite. Not a convergence test: a chain that stalls at a visible angle gap has
   * still produced an answer, and says so through an `error` diagnostic instead.
   */
  ok: boolean;

  /** The last {@link RATE_WINDOW} residuals, oldest first — a window, not a history. */
  history: number[];

  /** Most backtracks any one line search needed. */
  backtracks: number;

  /** True if a line search ran out of room without improving the merit. */
  starved: boolean;

  /** Worst `‖Ax − b‖∞` after refinement, over the steps taken. `Infinity` if one failed. */
  refinement: number;

  /** Largest `|λ_v|` seen. §8's corner evidence — measured here, acted on by nobody. */
  multiplier: number;

  /** Present only with {@link SPowerSolverOptions.trace} on. */
  trace?: TraceStep[];
}

/**
 * A stability criterion of §8, crossed, and pinned to the joint that will pay for it.
 *
 * Localization is the part that is not in §8, because §8 describes the criteria and not what
 * they are evidence *about*. Two of the three measure something wider than one joint — a
 * degenerate chord belongs to an edge, a divergent iteration to a whole chain — so each has
 * to name a joint before the response can be per-joint, and the choice is documented at the
 * point it is made in {@link ChainSystem.faults}.
 */
export interface Fault {
  condition: DiagnosticCondition;

  /** Vertex index within the chain. Always an interior one: a chain end has no level to lose. */
  at: number;

  vertex: SolvableVertex;
  measured: number;
  threshold: number;
}

/** `Rᵥ` for every vertex, from chord lengths only — §3, and §5's note that `Rᵥ` is shared. */
export function referenceLengths(mesh: SolvableMesh) {
  const refs = new Map<SolvableVertex, number>();
  const chords: number[] = [];

  for (const v of mesh.verts) {
    chords.length = 0;

    for (const e of v.edges) {
      chords.push(e.v1.vectorDistance(e.v2));
    }

    refs.set(v, referenceLength(chords));
  }

  return refs;
}

/**
 * Which of an edge's two ends the frozen G1 row can see, and with what sign.
 *
 * The chain-direction tangent angle at an end is `θ_e(u) + (edge runs backwards ? π : 0)`,
 * and `θ_e(u) = ∫₀^u q + KTH_e`. Freezing `KTH_e` leaves the turning integral, which is only
 * nonzero at `u = 1` — the edge's own far end. So an edge enters the row exactly when the
 * joint sits at its `v2`-side in chain order, which for the arriving edge means it runs
 * forwards and for the leaving edge means it runs backwards.
 */
function seesTurning(frame: EdgeFrame, arriving: boolean) {
  return arriving === frame.forward;
}

/**
 * One chain's linear system: DOF layout, the banded KKT matrix, and the Newton loop.
 *
 * Exported for tests, not for the barrel. Assembly is the part of Phase 2 with no
 * independent oracle — a wrong ordering and a wrong sign both just look like slow
 * convergence — so `tests/spower_solver.test.ts` reaches in and checks the assembled matrix
 * against quantities the curves compute for themselves.
 */
export class ChainSystem {
  p: number;
  block: number;

  frames: EdgeFrame[] = [];
  curves: SPowerClothoid[] = [];

  /** Index of the *shared* part of `block(verts[i])`. All of it wherever the level is at ceiling. */
  blockAt: Int32Array;

  /** Where the arriving and leaving edge-ends keep the entries they do not share. */
  arrivingAt: Int32Array;
  leavingAt: Int32Array;

  /** Entries shared at `verts[i]`, i.e. `delivered − 1` clamped — see {@link slot}. */
  shared: Int32Array;

  /** Authored continuity level at each vertex, ceiling-clamped. `-1` at a chain end. */
  requested: Int32Array;

  /** What the solve actually holds there, after any §8 lowering. Never above {@link requested}. */
  delivered: Int32Array;

  /** Index of the multiplier at `verts[i]`, or `-1` where there is no enforced joint. */
  lambdaAt: Int32Array;

  rows: number[] = [];
  n = 0;

  z: Float64Array;
  kkt: BandedSymmetric;

  /** Per-edge geometry read back from the curves each step: turning, placement, arclength. */
  turning: Float64Array;
  angle: Float64Array;
  arclength: Float64Array;

  /** `∂KTH_e/∂a` and `∂L_e/∂a` per edge — §5, and empty with the Jacobian frozen. */
  jacobians: TransformJacobian[] = [];

  /** Scratch for one coefficient-space constraint row. */
  weights: Float64Array;

  /** Scratch for one assembled G1 row, column-indexed — see {@link g1Row}. */
  g1 = new Map<number, number>();

  /** Element Hessians as of the last {@link assemble}, kept so the merit can reuse them. */
  elements: Float64Array[] = [];

  residuals: Float64Array;

  constructor(
    public chain: Chain,
    refs: Map<SolvableVertex, number>,
    public options: SPowerSolverOptions,
    lowered?: Map<SolvableVertex, number>
  ) {
    const { verts, edges } = chain;
    const dof = options.dof;

    this.p = options.order;
    this.block = dof.blockLength(this.p);

    for (let i = 0; i < edges.length; i++) {
      const e = edges[i] as SolvableEdge;
      const curve = e.curve as SPowerClothoid;

      curve.update(e);
      curve.setProfile(dof.profile, dof.quadrature);

      const count = dof.coefficientLength(this.p);

      if (curve.order !== this.p || curve.klen !== count) {
        curve.setOrder(this.p, count);
      }

      this.curves.push(curve);

      const chord = verts[i].vectorDistance(verts[i + 1]);

      this.frames.push({
        chord,
        arclength: chord,
        rEarlier : refs.get(verts[i]) ?? 1.0,
        rLater   : refs.get(verts[i + 1]) ?? 1.0,
        forward  : e.v1 === verts[i],
      });
    }

    this.blockAt = new Int32Array(verts.length);
    this.arrivingAt = new Int32Array(verts.length).fill(-1);
    this.leavingAt = new Int32Array(verts.length).fill(-1);
    this.lambdaAt = new Int32Array(verts.length).fill(-1);

    this.shared = new Int32Array(verts.length).fill(this.block);
    this.requested = new Int32Array(verts.length).fill(-1);
    this.delivered = new Int32Array(verts.length).fill(-1);

    for (let i = 1; i < edges.length; i++) {
      const { entries, tangent } = sharing(verts[i], edges[i - 1], edges[i], this.p);

      // Coarsening makes the shared orders a prefix, so the level is recoverable from a count.
      const level = tangent ? entries + 1 : 0;
      const cap = lowered?.get(verts[i]);

      this.requested[i] = level;
      this.delivered[i] = cap === undefined ? level : Math.min(level, cap);
      this.shared[i] = Math.max(0, Math.min(this.delivered[i] - 1, this.block));
    }

    this.layout();

    this.z = new Float64Array(this.n);
    this.kkt = new BandedSymmetric(this.n, this.bandwidth());

    this.turning = new Float64Array(edges.length);
    this.angle = new Float64Array(edges.length);
    this.arclength = new Float64Array(edges.length);
    this.residuals = new Float64Array(verts.length);
    this.weights = new Float64Array(dof.coefficientLength(this.p));
  }

  /**
   * True if the row at interior vertex `i` has any nonzero entry at all.
   *
   * Always, with the exact Jacobian: both edges reach the joint through `∂KTH/∂x`, which is
   * nonzero even on a straight edge. Only the frozen row can vanish.
   */
  enforced(i: number) {
    const f = this.frames;

    if (this.delivered[i] < 1) {
      return false;
    }

    return this.options.jacobian === "exact" || seesTurning(f[i - 1], true) || seesTurning(f[i], false);
  }

  /**
   * Interleave multipliers with blocks: `[b₀] [λ₁] [b₁] … [λ_{m−1}] [b_m]`.
   *
   * Placing `λ_i` immediately before `block(v_i)` is what keeps the system banded — its row
   * reaches the two adjacent blocks in either direction, and no further.
   *
   * A vertex below the ceiling level splits its block in three: the entries both ends share,
   * then a private tail for the arriving end and one for the leaving end. At the ceiling the
   * tails are empty and this is the single block of §3.
   */
  layout() {
    const m = this.chain.edges.length;
    let at = 0;

    for (let i = 0; i < this.chain.verts.length; i++) {
      const interior = i > 0 && i < m;

      if (interior && this.enforced(i)) {
        this.lambdaAt[i] = at++;
        this.rows.push(this.lambdaAt[i]);
      }

      this.blockAt[i] = at;
      at += this.shared[i];

      if (interior) {
        const rest = this.block - this.shared[i];

        this.arrivingAt[i] = at;
        at += rest;

        this.leavingAt[i] = at;
        at += rest;
      }
    }

    this.n = at;
  }

  /**
   * Widest offset any assembled entry spans, measured off the layout rather than predicted
   * from it — the split blocks make the old closed form wrong in both directions.
   */
  bandwidth() {
    const m = this.chain.edges.length;
    const d = 2 * this.block;
    let w = 0;

    for (let i = 0; i < m; i++) {
      let lo = Infinity;
      let hi = -Infinity;

      for (let k = 0; k < d; k++) {
        const at = this.at(i, k);

        lo = Math.min(lo, at);
        hi = Math.max(hi, at);
      }

      w = Math.max(w, hi - lo);
    }

    for (let i = 1; i < m; i++) {
      const at = this.lambdaAt[i];

      if (at < 0) {
        continue;
      }

      for (const edge of [i - 1, i]) {
        for (let k = 0; k < d; k++) {
          w = Math.max(w, Math.abs(this.at(edge, k) - at));
        }
      }
    }

    return Math.min(w, Math.max(0, this.n - 1));
  }

  /**
   * Global index of entry `n` of `verts[i]`, as read by the end arriving at it or leaving it.
   *
   * The two agree on the first {@link shared} entries and diverge after — which is the whole
   * of what a lowered level does to the state space. It is a change of *structure*: there is
   * no residual anywhere that measures the resulting `Gⁿ` discontinuity, because it is not a
   * quantity the solve can see.
   */
  slot(i: number, n: number, arriving: boolean) {
    if (n < this.shared[i]) {
      return this.blockAt[i] + n;
    }

    return (arriving ? this.arrivingAt[i] : this.leavingAt[i]) + n - this.shared[i];
  }

  /** Global index of local DOF `k` of edge `i`, whose DOF are its two blocks in chain order. */
  at(i: number, k: number) {
    return k < this.block ? this.slot(i, k, false) : this.slot(i + 1, k - this.block, true);
  }

  /** Push the current DOF through `M_e` into each curve's coefficients. */
  write() {
    const { p, block, frames, curves } = this;
    const dof = this.options.dof;
    const scratch = new Float64Array(2 * block);

    for (let i = 0; i < curves.length; i++) {
      for (let k = 0; k < 2 * block; k++) {
        scratch[k] = this.z[this.at(i, k)];
      }

      dof.coefficients(p, frames[i], scratch, curves[i].ks);
      curves[i].recalc = 1;
    }
  }

  /**
   * Read the realized geometry back and measure every joint.
   *
   * Reading `length` is what forces `_update`, so this is also where the next step's `L_e`
   * estimate comes from. The residual is measured on the actual curves rather than on
   * `J·x`, which is what lets it converge despite the frozen Jacobian.
   */
  measure() {
    const { curves, frames } = this;
    const dof = this.options.dof;
    const exact = this.options.jacobian === "exact";

    for (let i = 0; i < curves.length; i++) {
      this.arclength[i] = curves[i].length;
      this.angle[i] = curves[i].th;
      this.turning[i] = curves[i].turning;

      if (exact) {
        this.jacobians[i] = transformJacobian(
          dof.profile,
          curves[i].ks,
          this.weights.length,
          frames[i].chord,
          curves[i].quadrature,
          this.jacobians[i]
        );
      }
    }

    let worst = 0.0;

    for (let i = 1; i < curves.length; i++) {
      this.residuals[i] = 0.0;

      if (this.lambdaAt[i] < 0) {
        continue;
      }

      const incoming = frames[i - 1].forward ? this.turning[i - 1] + this.angle[i - 1] : this.angle[i - 1] + Math.PI;

      const outgoing = frames[i].forward ? this.angle[i] : this.turning[i] + this.angle[i] + Math.PI;

      this.residuals[i] = wrapAngle(outgoing - incoming);
      worst = Math.max(worst, Math.abs(this.residuals[i]));
    }

    return worst;
  }

  /**
   * Scatter `Σ M_eᵀ E_e M_e` and the G1 rows into the band.
   *
   * {@link measure} has to have run since the last {@link write}: the exact rows read the
   * transform derivatives it caches.
   */
  assemble() {
    const { p, block, frames, kkt } = this;
    const dof = this.options.dof;
    const d = 2 * block;

    kkt.zero();

    for (let i = 0; i < frames.length; i++) {
      const h = dof.hessian(p, frames[i], this.options.alpha, this.elements[i]);

      this.elements[i] = h;

      // Local and global index order agree, so the local lower triangle is the global one.
      for (let r = 0; r < d; r++) {
        for (let c = 0; c <= r; c++) {
          kkt.add(this.at(i, r), this.at(i, c), h[r * d + c]);
        }
      }
    }

    for (let i = 1; i < frames.length; i++) {
      if (this.lambdaAt[i] >= 0) {
        this.constraintRow(i);
      }
    }
  }

  /**
   * The G1 row at interior vertex `i`: `∂(outgoing − incoming)/∂x`.
   *
   * The tangent angles {@link measure} compares are built from `KTH_e` and, at an edge's own
   * far end, its total turning. Both are differentiated in coefficient space and pulled back
   * through `M_e` together, so the row touches all four blocks the two edges read. `∂KTH/∂a`
   * is what the frozen version threw away, and it is the larger of the two terms.
   *
   * `L_e` still lags: it enters only through `M_e`, so the transform used here is the one the
   * coefficients were built with, refreshed between passes. See §14.
   */
  constraintRow(i: number) {
    const at = this.lambdaAt[i];

    for (const [column, value] of this.g1Row(i)) {
      this.kkt.add(at, column, value);
    }
  }

  /**
   * The same row as a sparse column-to-coefficient map, which is what the rank test needs.
   *
   * Assembly could scatter as it goes, but a joint whose two edges share a block writes the
   * same column twice and {@link BandedSymmetric.add} accumulates — so the *matrix* is right
   * either way while the *row* is not, and a rank test on a row with duplicate columns is
   * measuring the wrong thing. One code path, summed first, keeps them the same object.
   */
  g1Row(i: number) {
    const { p, block, frames } = this;
    const dof = this.options.dof;
    const frozen = this.options.jacobian === "frozen";
    const row = this.g1;

    row.clear();

    // The arriving edge enters negatively, the leaving edge positively — as in `measure`.
    for (const arriving of [true, false]) {
      const edge = arriving ? i - 1 : i;
      const frame = frames[edge];
      const sign = arriving ? -1.0 : 1.0;
      const far = seesTurning(frame, arriving);

      let entries: Float64Array;

      if (frozen) {
        if (!far) {
          continue;
        }

        entries = dof.turningRow(p, frame);
      } else {
        const w = this.weights;

        w.set(this.jacobians[edge].dKth);

        if (far) {
          const turning = dof.turningWeights(p);

          for (let c = 0; c < w.length; c++) {
            w[c] += turning[c];
          }
        }

        entries = dof.pullback(p, frame, w);
      }

      for (let k = 0; k < 2 * block; k++) {
        const column = this.at(edge, k);

        row.set(column, (row.get(column) ?? 0.0) + sign * entries[k]);
      }
    }

    return row;
  }

  /**
   * Run the outer iteration. Returns the steps taken, the final residual, and whether every
   * factorization held.
   *
   * Each pass writes the DOF out, measures the geometry they produced, then linearizes about
   * *that* state — so the transform used to build the coefficients and the transform used to
   * build `H` are the same one. The arclength estimate only moves at the end of a pass,
   * which is what makes `L_e` an outer-iteration quantity rather than a circular one.
   *
   * Steps are accepted by a backtracking search on §5's `ℓ1` merit
   *
   * ```
   *   ϕ(z) = ½ zᵀHz + μ·Σ|c_i(z)| ,       μ > max|λ|
   * ```
   *
   * which is the standard exact penalty for an equality-constrained minimization: with `μ`
   * above the multipliers the Newton direction is a descent direction for it, and a solution
   * of the original problem is a local minimum of it. `H` is held at the linearization point
   * for the duration of one search — only `c` is remeasured, on the actual curves — so the
   * comparison is between two states of the same model.
   */
  run(): ChainRun {
    const opts = this.options;
    const grad = new Float64Array(this.n);
    const rhs = new Float64Array(this.n);
    const base = new Float64Array(this.n);

    const history: number[] = [];
    const trace = opts.trace ? ([] as TraceStep[]) : undefined;

    let steps = 0;
    let factored = true;
    let mu = 1.0;

    let backtracks = 0;
    let starved = false;
    let refinement = 0.0;
    let multiplier = 0.0;

    this.write();

    let residual = this.measure();

    for (let iter = 0; ; iter++) {
      if (iter >= opts.iterations || residual < opts.tolerance || this.rows.length === 0) {
        break;
      }

      this.assemble();
      this.kkt.apply(this.z, grad);

      for (let i = 0; i < this.n; i++) {
        rhs[i] = -grad[i];
      }

      for (let i = 1; i < this.frames.length; i++) {
        if (this.lambdaAt[i] >= 0) {
          rhs[this.lambdaAt[i]] = -this.residuals[i];
        }
      }

      const step = solveKKT(this.kkt, rhs, this.rows, opts.kkt);

      refinement = Math.max(refinement, step.residual);

      if (!step.ok) {
        factored = false;
        break;
      }

      for (const row of this.rows) {
        mu = Math.max(mu, 2.0 * Math.abs(this.z[row] + step.x[row]));
      }

      base.set(this.z);

      const start = this.merit(mu);
      let t = opts.relaxation;
      let cuts = 0;

      for (;;) {
        for (let i = 0; i < this.n; i++) {
          this.z[i] = base[i] + t * step.x[i];
        }

        this.write();
        this.measure();

        if (this.merit(mu) < start) {
          break;
        }

        // Out of room rather than out of patience: the direction was not a descent one.
        if (t < MIN_RELAXATION) {
          starved = true;
          break;
        }

        t *= SHRINK;
        cuts++;
      }

      backtracks = Math.max(backtracks, cuts);

      // `L_e` only moves here, between passes — the coefficients are rebuilt through the
      // updated `M_e` so the next linearization point is a self-consistent one.
      for (let i = 0; i < this.frames.length; i++) {
        this.frames[i].arclength = this.arclength[i];
      }

      this.write();
      residual = this.measure();

      steps++;

      const largest = this.largestMultiplier();

      multiplier = Math.max(multiplier, largest);

      if (history.length >= RATE_WINDOW) {
        history.shift();
      }

      history.push(residual);

      trace?.push({
        merit: this.merit(mu),
        residual,
        stepNorm     : t * infinityNorm(step.x),
        maxMultiplier: largest,
        relaxation   : t,
        refinement   : step.residual,
      });
    }

    return {
      steps,
      residual,
      factored,
      ok: factored && Number.isFinite(residual),
      history,
      backtracks,
      starved,
      refinement,
      multiplier,
      trace,
    };
  }

  /** `max|λ_v|` over the chain's enforced joints, read off the current solution vector. */
  largestMultiplier() {
    let worst = 0.0;

    for (const row of this.rows) {
      worst = Math.max(worst, Math.abs(this.z[row]));
    }

    return worst;
  }

  /**
   * `σ_min/σ_max` of the G1 rows at `verts[i]`, by column-pivoted QR — §8's first criterion.
   *
   * §8 is emphatic that this is measured on the *local* block and not by watching the pivots
   * of the global factorization: pivot magnitude is not a rank test, Kahan's matrix being the
   * standard counterexample. So the block is rebuilt here, in its own column order, where
   * pivoting is free to choose.
   *
   * At valence 2 the block is a single row and the ratio is therefore `1` unless the row
   * vanishes outright — which only the frozen Jacobian can arrange. The criterion is a
   * junction phenomenon; the chain case is the degenerate one it has to survive, and does.
   */
  rankAt(i: number) {
    if (this.lambdaAt[i] < 0) {
      return 0.0;
    }

    const row = this.g1Row(i);
    const cols = row.size;

    if (cols === 0) {
      return 0.0;
    }

    const a = new Float64Array(cols);
    let c = 0;

    for (const value of row.values()) {
      a[c++] = value;
    }

    return rankRatio(a, 1, cols);
  }

  /**
   * Every §8 criterion that came out past its threshold, worst structural cause first.
   *
   * `broken` is the hysteresis state: a joint already carrying a break is judged against the
   * *restore* thresholds instead, so it has to become comfortably healthy to keep the level
   * this attempt handed back rather than merely stop being unhealthy. That is the whole of the
   * hysteresis, and it is per-joint rather than per-chain because the caps are.
   *
   * Ordering is rank, then chord, then divergence — cause before symptom. A rank-deficient
   * joint makes the iteration diverge, and lowering the joint the *divergence* pointed at
   * would be treating the second reading.
   */
  faults(run: ChainRun, broken?: Map<SolvableVertex, number>): Fault[] {
    const limits = stabilityThresholds;
    const { verts } = this.chain;
    const found: Fault[] = [];
    const healing = (v: SolvableVertex) => broken?.has(v) ?? false;

    for (let i = 1; i < this.frames.length; i++) {
      if (this.lambdaAt[i] < 0) {
        continue;
      }

      const measured = this.rankAt(i);
      const threshold = healing(verts[i]) ? limits.rankRestore : limits.rankBreak;

      if (!(measured > threshold)) {
        found.push({ condition: "g1-rank-deficiency", at: i, vertex: verts[i], measured, threshold });
      }
    }

    for (let e = 0; e < this.curves.length; e++) {
      const measured = this.curves[e].canonical;

      // The edge is degenerate, but an edge has no level; its two ends do. Whichever end the
      // solve is working hardest at is the one whose continuity is buying the degeneracy.
      const at = this.strainedEnd(e);

      if (at < 0) {
        continue;
      }

      const threshold = healing(verts[at]) ? limits.chordRestore : limits.chordBreak;

      if (!(measured <= 1.0) || measured < threshold) {
        found.push({ condition: "chord-degeneracy", at, vertex: verts[at], measured, threshold });
      }
    }

    const rate = geometricRate(run.history);
    const stalled = run.steps >= this.options.iterations && !(run.residual < this.options.tolerance);
    const at = this.worstJoint();

    if (at >= 0) {
      const threshold = healing(verts[at]) ? limits.rateRestore : limits.rateBreak;

      // All three of §8's signals, and the fit is required rather than assumed: a run too
      // short to fit a rate over three iterations is a budget, not a divergence, and a
      // solve cut off at two steps must not cost anyone a continuity level.
      const diverging = Number.isFinite(rate) && rate >= threshold && (stalled || run.starved);

      if (!run.ok || diverging) {
        const measured = Number.isFinite(rate) ? rate : Infinity;

        found.push({ condition: "newton-not-converging", at, vertex: verts[at], measured, threshold });
      }
    }

    return found;
  }

  /**
   * Which end of edge `e` a break would have to act on: the interior one, or the busier of two.
   *
   * `-1` when neither end is an interior joint, which is a two-vertex chain — nothing to
   * lower, so nothing to report as a fault either.
   */
  strainedEnd(e: number) {
    const interior = (i: number) => i > 0 && i < this.frames.length && this.lambdaAt[i] >= 0;
    const a = interior(e) ? e : -1;
    const b = interior(e + 1) ? e + 1 : -1;

    if (a < 0 || b < 0) {
      return Math.max(a, b);
    }

    return Math.abs(this.residuals[a]) >= Math.abs(this.residuals[b]) ? a : b;
  }

  /** The enforced joint with the largest angle gap — where a non-converging chain is stuck. */
  worstJoint() {
    let at = -1;
    let worst = -1.0;

    for (let i = 1; i < this.frames.length; i++) {
      if (this.lambdaAt[i] >= 0 && Math.abs(this.residuals[i]) > worst) {
        worst = Math.abs(this.residuals[i]);
        at = i;
      }
    }

    return at;
  }

  /**
   * Turn one run into §8's records: located, typed, and carrying the measurement that
   * produced them alongside the threshold it crossed.
   *
   * Severity says which side of the threshold the *solve* came out on, not which side the
   * measurement did — a chain can converge slowly (`warning`) or not at all (`error`) on the
   * same fitted rate, and only the second is a failure. {@link Diagnostic.action} is `"none"`
   * throughout: these are the records of the attempt that was *kept*, so by construction
   * nothing was done about them. The records for what the solver did are written by
   * {@link SPowerSolver.solve}, which is the only thing that knows an attempt was one of
   * several.
   */
  diagnose(chain: number, run: ChainRun, into: Diagnostic[]) {
    const limits = diagnosticThresholds;
    const opts = this.options;

    for (let i = 0; i < this.curves.length; i++) {
      const measured = this.curves[i].canonical;

      // `|C(1)| ≤ 1` exactly, the canonical curve having unit arclength. Above it — or not a
      // number at all — is the quadrature having lost the curve, which is the same fault
      // arriving from the other side and is worth catching in the same place.
      const lost = !(measured <= 1.0);

      if (!lost && measured >= limits.chordDegeneracy) {
        continue;
      }

      into.push({
        condition: "chord-degeneracy",
        severity : lost || measured <= MIN_CANONICAL_CHORD ? "error" : "warning",
        action   : "none",
        chain,
        at  : i,
        edge: this.chain.edges[i],
        measured,
        threshold: limits.chordDegeneracy,
      });
    }

    if (!run.factored) {
      const delta = opts.kkt.delta ?? defaultKKTSolveOptions.delta;

      into.push({
        condition: "factorization-underflow",
        severity : "error",
        action   : "none",
        chain,
        at       : -1,
        measured : delta,
        threshold: delta,
      });
    }

    const rate = geometricRate(run.history);

    // Negated rather than `>=`, so a residual that went non-finite counts as not converged.
    const stalled = run.steps >= opts.iterations && !(run.residual < opts.tolerance);

    if (stalled || rate >= limits.newtonRate) {
      into.push({
        condition: "newton-not-converging",
        severity : stalled ? "error" : "warning",
        action   : "none",
        chain,
        at       : -1,
        measured : rate,
        threshold: limits.newtonRate,
      });
    }

    if (run.starved || run.backtracks >= limits.lineSearch) {
      into.push({
        condition: "line-search-failing",
        severity : run.starved ? "error" : "warning",
        action   : "none",
        chain,
        at       : -1,
        measured : run.backtracks,
        threshold: limits.lineSearch,
      });
    }

    // `Infinity` here is the failed factorization, which has its own record above.
    if (run.refinement > limits.refinement && Number.isFinite(run.refinement)) {
      into.push({
        condition: "refinement-not-converging",
        severity : "warning",
        action   : "none",
        chain,
        at       : -1,
        measured : run.refinement,
        threshold: limits.refinement,
      });
    }

    this.multipliers(chain, into);

    return into;
  }

  /**
   * Joints whose multiplier stands out from the chain's — §8's corner evidence, for a client.
   *
   * `|λ_v|` is the cost of holding the joint tangent, in the units of the energy, so what it
   * says is "a corner may have been wanted here". §8 is explicit that this is *modelling* and
   * not stability, and that the solver may not act on it: hence `info`, `action: "none"`, and
   * no entry in {@link faults}. It is measured against the chain's own median rather than an
   * absolute, because the energy scale is set by `Rᵥ` and means nothing across strokes.
   */
  multipliers(chain: number, into: Diagnostic[]) {
    const ratio = stabilityThresholds.multiplierRatio;
    const magnitudes = this.rows.map((row) => Math.abs(this.z[row])).sort((a, b) => a - b);

    if (magnitudes.length === 0) {
      return;
    }

    const median = magnitudes[magnitudes.length >> 1];

    if (!(median > 0.0)) {
      return;
    }

    for (let i = 1; i < this.frames.length; i++) {
      const row = this.lambdaAt[i];

      if (row < 0) {
        continue;
      }

      const measured = Math.abs(this.z[row]);

      if (measured < ratio * median) {
        continue;
      }

      into.push({
        condition: "large-multiplier",
        severity : "info",
        action   : "none",
        chain,
        at    : i,
        vertex: this.chain.verts[i],
        measured,
        threshold: ratio * median,
      });
    }
  }

  /**
   * `½zᵀHz + μ·Σ|c_i|`, using the element Hessians as of the last {@link assemble}.
   *
   * Reading `H` from the cached elements rather than from the band keeps the multiplier slots
   * out of it — `kkt.apply` would fold `Jᵀλ` into the same numbers.
   */
  merit(mu: number) {
    const d = 2 * this.block;
    let energy = 0.0;

    for (let e = 0; e < this.elements.length; e++) {
      const h = this.elements[e];

      for (let r = 0; r < d; r++) {
        const zr = this.z[this.at(e, r)];

        for (let c = 0; c < d; c++) {
          energy += h[r * d + c] * zr * this.z[this.at(e, c)];
        }
      }
    }

    let infeasible = 0.0;

    for (let i = 1; i < this.frames.length; i++) {
      infeasible += Math.abs(this.residuals[i]);
    }

    return 0.5 * energy + mu * infeasible;
  }

  /**
   * Interior joints this system chose not to constrain. See the module doc comment.
   *
   * A joint delivered at level 0 has no row *by construction* and is not counted — that is an
   * authored corner, or a break the report already accounts for elsewhere.
   */
  unenforced() {
    let count = 0;

    for (let i = 1; i < this.frames.length; i++) {
      if (this.delivered[i] >= 1 && this.lambdaAt[i] < 0) {
        count++;
      }
    }

    return count;
  }
}

/**
 * Fits every {@link SPowerClothoid} in a mesh so the segments meet tangentially.
 *
 * Chains are solved independently, which §5 licenses as long as no vertex block is shared
 * across a chain boundary — true here, since Phase 2 predates the pairing levels of §3 and
 * every node is fully split. `Rᵥ` still couples them, benignly, because it is built from
 * chord lengths and so does not move during the solve.
 */
export class SPowerSolver implements CurveSolver {
  options: SPowerSolverOptions;
  report = emptyReport();

  /**
   * Levels currently withheld, per vertex — the whole of §8's hysteresis state.
   *
   * A solve always *tries* the authored levels first, so what is in here never lowers anything
   * by itself; it only says which joints must clear the stricter restore thresholds rather
   * than the break ones to keep what the attempt handed them. That keeps the delivered
   * geometry a function of the input alone, while still making a joint that sits on a
   * threshold stay where it is instead of blinking between a corner and a curve as the input
   * jitters — which is the failure §8 asks hysteresis to prevent.
   *
   * Survives across solves on the solver instance, so a caller that rebuilds the solver each
   * frame gets correct but memoryless behaviour. `Mesh` keeps one per spline type.
   */
  broken = new Map<SolvableVertex, number>();

  constructor(
    public mesh: SolvableMesh,
    options: Partial<SPowerSolverOptions> = {}
  ) {
    this.options = { ...defaultSPowerSolverOptions, ...options };

    if (!(this.options.alpha > 0.0)) {
      throw new Error("alpha must be positive — see SPowerSolverOptions.alpha");
    }
  }

  solve() {
    const report = emptyReport();
    const refs = referenceLengths(this.mesh);

    if (this.options.trace) {
      report.traces = [];
    }

    // Skipped chains take an index too, so the diagnostic's `chain` is a position and not a rank.
    let index = -1;

    for (const chain of chains(this.mesh)) {
      index++;
      report.chains++;

      if (chain.closed) {
        report.skippedClosed++;
        continue;
      }

      this.solveChain(chain, refs, index, report);
    }

    this.report = report;

    return report;
  }

  /**
   * Solve one chain, lowering a joint and retrying for as long as §8 says the answer cannot
   * be trusted.
   *
   * The first attempt always asks for the authored levels, whatever is in {@link broken} —
   * hysteresis lives in the thresholds {@link ChainSystem.faults} compares against, not in the
   * levels it starts from. So the delivered geometry depends on the input and not on the order
   * of solves that reached it, which is the stronger half of §8's determinism requirement, and
   * a joint still has to earn its level back against the restore thresholds.
   *
   * Each pass lowers exactly one joint by exactly one level — the shallowest response that
   * addresses the fault — and rebuilds the system, because a level is a change of *layout* and
   * not a coefficient that could be edited in place.
   */
  solveChain(chain: Chain, refs: Map<SolvableVertex, number>, index: number, report: SPowerSolverReport) {
    const opts = this.options;
    const caps = new Map<SolvableVertex, number>();
    const actions: Diagnostic[] = [];

    let system = new ChainSystem(chain, refs, opts, caps);
    let run = system.run();
    let refused = false;

    report.attempts++;

    for (let attempt = 0; attempt < opts.breaks; attempt++) {
      const faults = system.faults(run, this.broken);

      if (faults.length === 0) {
        break;
      }

      if (opts.mode === "engineering") {
        refused = true;

        for (const fault of faults) {
          actions.push(this.record(fault, "refused", index, system, system.delivered[fault.at]));
        }

        break;
      }

      // A joint already at 0 has nothing left to give; look past it rather than giving up,
      // since the fault it still reports may be the second-worst one and fixable.
      const fault = faults.find((f) => system.delivered[f.at] > 0);

      if (!fault) {
        break;
      }

      const to = system.delivered[fault.at] - 1;

      actions.push(this.record(fault, "degraded", index, system, to));
      caps.set(fault.vertex, to);

      system = new ChainSystem(chain, refs, opts, caps);
      run = system.run();

      report.attempts++;
    }

    for (let i = 1; i < chain.edges.length; i++) {
      const v = chain.verts[i];
      const cap = caps.get(v);

      if (cap === undefined) {
        this.broken.delete(v);
      } else {
        this.broken.set(v, cap);
        report.degraded++;
      }
    }

    report.steps += run.steps;
    report.unenforced += system.unenforced();
    report.maxResidual = Math.max(report.maxResidual, run.residual);
    report.ok &&= run.ok && !refused;

    system.diagnose(index, run, report.diagnostics);
    report.diagnostics.push(...actions);

    if (run.trace) {
      report.traces?.push({ chain: index, steps: run.trace });
    }
  }

  /**
   * One §8 record for something the solver *did*: the fault that caused it, and the levels
   * either side of it.
   *
   * A joint that is only being put back where it already was gets `level-lowered` at `info`
   * instead of its fault condition at `warning`. Both are true; the distinction is the one a
   * caller reading the report actually wants, which is whether this is news.
   */
  record(fault: Fault, action: "degraded" | "refused", chain: number, system: ChainSystem, delivered: number) {
    const carried = this.broken.get(fault.vertex);
    const retained = action === "degraded" && carried !== undefined && delivered >= carried;

    return {
      condition: retained ? "level-lowered" : fault.condition,
      severity : retained ? "info" : action === "refused" ? "error" : "warning",
      action,
      chain,
      at       : fault.at,
      vertex   : fault.vertex,
      measured : fault.measured,
      threshold: fault.threshold,
      requested: system.requested[fault.at],
      delivered,
    } satisfies Diagnostic;
  }
}
