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
import {
  BandedSymmetric,
  DenseLU,
  KKTFactorization,
  type KKTSolveOptions,
  defaultKKTSolveOptions,
  rankRatio,
} from "../math/index.js";
import { type EdgeFrame, type ProfileDOF, continuationEntry, referenceLength, sPowerDOF } from "./blocks.js";
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
import { type ChainEnd, type Component, type GammaSlot, components, interfaceSlots } from "./junctions.js";
import { type CurveSolver, type SolvableEdge, type SolvableMesh, type SolvableVertex } from "./mesh_types.js";
import { sharing } from "./pairing.js";
import { SPOWER_ORDER, type SPowerClothoid } from "./spower_clothoid.js";
import { MIN_CANONICAL_CHORD, type TransformJacobian, transformJacobian } from "./quadrature.js";
import { type Chain, chains, cutOpen } from "./topology.js";

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
   * Climb `p = 0 … order`, warm-starting each rung from the one below — §9.
   *
   * A warm start and not an embedding: the low-order entries of a block mean the same thing at
   * every order, so they carry over unchanged, but the new top entry cannot preserve the curve
   * — §9's argument is that the two edges at a joint want different values for it and the
   * ceiling makes them share one. So the rung below supplies a *starting point*, not a
   * solution, and the claim being tested is "fewer Newton steps" and nothing stronger. §14 has
   * what it measured.
   *
   * Off by default: it is only worth its extra assemblies where the top rung is struggling,
   * which at `p = 1` it generally is not.
   */
  continuation: boolean;

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
  order       : SPOWER_ORDER,
  alpha       : 0.1,
  iterations  : 100,
  relaxation  : 1.0,
  tolerance   : 1e-10,
  dof         : sPowerDOF,
  jacobian    : "exact",
  trace       : false,
  mode        : "artistic",
  continuation: false,
  breaks      : 3,
  kkt         : {},
};

export interface SPowerSolverReport extends SolveReport {
  /** Chains found by {@link chains}. */
  chains: number;

  /**
   * Chains that closed on themselves and were cut open at `verts[0]` — §5's cut point.
   *
   * The cut vertex becomes a junction and its two halves are reconciled through the same
   * interface machinery as any other node, so a cut chain is solved rather than skipped. What
   * this counts is how many of {@link chains} had to be treated that way.
   */
  cutClosed: number;

  /**
   * Interior joints whose G1 row vanished, and which were therefore left unconstrained.
   *
   * Only reachable with the Jacobian frozen; the exact row never vanishes.
   */
  unenforced: number;

  /** Largest `|wrapped angle gap|` over every enforced joint, measured after the last step. */
  maxResidual: number;

  /** Gauss-Newton steps actually taken, summed over components. */
  steps: number;

  /**
   * Joints delivered below the level they were authored with — §8's breaks, counted.
   *
   * Zero in engineering mode by construction: it refuses rather than degrades. Every one of
   * these has a `"degraded"` record in {@link SolveReport.diagnostics} saying where and why.
   */
  degraded: number;

  /** Component solves run, including the attempts thrown away. One per component when nothing broke. */
  attempts: number;

  /**
   * Steps spent on continuation rungs below `order`, and not counted in {@link steps} — §9.
   *
   * The price of the warm start, kept separate from what it bought so the two can be compared
   * without arithmetic. Zero unless {@link SPowerSolverOptions.continuation} is on.
   */
  seedSteps: number;
}

function emptyReport(): SPowerSolverReport {
  return {
    chains     : 0,
    cutClosed  : 0,
    unenforced : 0,
    maxResidual: 0.0,
    steps      : 0,
    degraded   : 0,
    attempts   : 0,
    seedSteps  : 0,
    ok         : true,
    diagnostics: [],
  };
}

/**
 * What one {@link ComponentSystem.run} measured on its way through, for {@link
 * ChainSystem.diagnose} to turn into §8 records.
 *
 * All of it is a byproduct: the history is the convergence test, the backtracks and the
 * relaxation floor are the line search, the refinement residual comes back from the solve,
 * and the multiplier is read off the solution vector. §8's second rule is what makes the
 * unconditional measurement of the first rule affordable.
 *
 * It is per *component* and the chains inside one share it, because the quantities are
 * global to the Newton loop — one line search, one residual history, one factorization
 * verdict. §8's localization runs the other way: each chain reads the same run and reports
 * the joints of its own that are answerable for what it says.
 */
export interface ComponentRun {
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

/**
 * §8's precedence among the localizable criteria: cause before symptom.
 *
 * {@link ChainSystem.faults} emits in this order already; the constant exists so the same
 * order can be re-imposed after faults from several chains of a component are merged, where
 * concatenation would otherwise let one chain's symptom outrank another chain's cause.
 */
const FAULT_ORDER: DiagnosticCondition[] = ["g1-rank-deficiency", "chord-degeneracy", "newton-not-converging"];

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
 * Whether the given end of an edge is the edge's *own* far end, `u = 1`.
 *
 * The chain-direction tangent angle at an end is `θ_e(u) + (edge runs backwards ? π : 0)`,
 * and `θ_e(u) = ∫₀^u q + KTH_e`. The turning integral is nonzero only at `u = 1`, so this is
 * exactly the test for whether an end's angle carries the edge's total turning — and, with
 * the Jacobian frozen, whether it enters the G1 row at all.
 *
 * `later` means the end at the higher chain index: the `v2`-side when the edge runs forwards
 * along the chain, the `v1`-side when it runs backwards.
 */
function seesTurning(frame: EdgeFrame, later: boolean) {
  return later === frame.forward;
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

  /**
   * Interface index of each local unknown, or `-1` where the unknown is private — §5.
   *
   * Set by {@link ComponentSystem}, and all `-1` until one claims this chain. Nothing else in
   * this class reads it: the layout, the band and the assembly are the same whether or not a
   * chain is coupled to another, and the split into `A_i`, `B_i` and `D` happens where the
   * matrix is *read*. That is what makes general valence a decomposition rather than a
   * change of representation.
   */
  localToGlobal: Int32Array;

  /**
   * `±1` relating each mirrored unknown to the interface value it mirrors — see
   * {@link ComponentSystem.orientation}.
   *
   * A block is written in the direction the *chain* travels through its vertex, and two chains
   * meeting at a node need not travel the same way through it. Identification across the
   * interface is therefore `z_local = σ · γ`, not `z_local = γ`, and `1` everywhere is what a
   * chain-interior joint gets by construction.
   */
  localSign: Float64Array;

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
    this.localToGlobal = new Int32Array(this.n).fill(-1);
    this.localSign = new Float64Array(this.n).fill(1.0);
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

  /**
   * The vertex, edge and block slot at one of the chain's two ends — `0` for `verts[0]`, `1`
   * for the last vertex, matching {@link ChainEnd.end}.
   *
   * A chain end is never a split vertex, since only interior joints have two edge-ends to
   * split between, so {@link slot} agrees with itself there whichever way it is asked.
   */
  endVertex(end: 0 | 1) {
    return end === 1 ? this.frames.length : 0;
  }

  endEdge(end: 0 | 1) {
    return end === 1 ? this.frames.length - 1 : 0;
  }

  endSlot(end: 0 | 1, n: number) {
    return this.slot(this.endVertex(end), n, end === 1);
  }

  /**
   * Chain-direction tangent angle at one end of edge `e`.
   *
   * `θ_e(u) + (edge runs backwards ? π : 0)`, evaluated at whichever of `u = 0, 1` this end
   * is — see {@link seesTurning} for which that is.
   */
  endAngle(e: number, later: boolean) {
    const base = this.angle[e] + (seesTurning(this.frames[e], later) ? this.turning[e] : 0.0);

    return this.frames[e].forward ? base : base + Math.PI;
  }

  /**
   * The direction the curve leaves a vertex in, along edge `e`.
   *
   * {@link endAngle} points along the chain, so at the `later` end it points *into* the
   * vertex and has to be reversed. Two edge-ends are G1-continuous exactly when their
   * outgoing directions are antiparallel, which is the one form that says the same thing at a
   * chain joint and at a junction of any valence: there is no "arriving" edge at a node.
   */
  outgoing(e: number, later: boolean) {
    return this.endAngle(e, later) + (later ? Math.PI : 0.0);
  }

  /**
   * Take this system's starting point from a converged solve one order below — §9's rung.
   *
   * Three things move across. The entries both orders have mean the same thing at every order
   * — `(Rᵥⁿ/n!)·dⁿκ/dsⁿ` is not a function of `p` — so they carry over untouched. The
   * arclengths carry over too, which replaces the `L_e = C_e` seed with a measured one and is
   * the part that costs nothing and is never wrong. And the new top entry is read off the
   * curve the rung below actually produced, via {@link continuationEntry}.
   *
   * That last one is where §9's "warm start, not an embedding" bites: the two edges at a joint
   * disagree about the new derivative, and at the ceiling they must share one number, so what
   * lands there is their average and the curve moves. Averaging by incidence count does the
   * right thing at both ends of the chain and at a split joint without a special case, since
   * a split entry is only written by the one end that owns it.
   */
  seedFrom(prev: ChainSystem) {
    const { verts, edges } = this.chain;
    const carry = Math.min(prev.block, this.block);
    const counts = new Float64Array(this.n);

    for (let i = 0; i < verts.length; i++) {
      for (let n = 0; n < carry; n++) {
        for (const arriving of [true, false]) {
          this.z[this.slot(i, n, arriving)] = prev.z[prev.slot(i, n, arriving)];
        }
      }
    }

    if (this.block > prev.block) {
      const scratch = new Float64Array(2 * prev.block);

      for (let e = 0; e < edges.length; e++) {
        for (let k = 0; k < scratch.length; k++) {
          scratch[k] = prev.z[prev.at(e, k)];
        }

        for (const earlier of [true, false]) {
          const at = this.slot(earlier ? e : e + 1, prev.block, !earlier);
          const value = continuationEntry(prev.p, prev.frames[e], scratch, earlier);

          this.z[at] += value;
          counts[at]++;
        }
      }

      for (let i = 0; i < this.n; i++) {
        if (counts[i] > 1) {
          this.z[i] /= counts[i];
        }
      }
    }

    for (let e = 0; e < edges.length; e++) {
      this.frames[e].arclength = prev.frames[e].arclength;
    }
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

      this.residuals[i] = wrapAngle(this.endAngle(i, false) - this.endAngle(i - 1, true));
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
    const row = this.g1;

    row.clear();

    const add = (column: number, value: number) => {
      row.set(column, (row.get(column) ?? 0.0) + value);
    };

    // The arriving edge enters negatively, the leaving edge positively — as in `measure`.
    this.endRow(i - 1, true, -1.0, add);
    this.endRow(i, false, 1.0, add);

    return row;
  }

  /**
   * `sign · ∂(endAngle(e, later))/∂x`, emitted one local column at a time.
   *
   * A column at a time rather than into a row object, because a *node* row spans edge-ends in
   * different chains and there is no one local index space to write it into: the caller routes
   * each column to wherever that unknown lives. {@link g1Row} is the two-call chain case.
   *
   * `∂KTH/∂a` is what the frozen version throws away, and it is the larger of the two terms.
   * `L_e` still lags: it enters only through `M_e`, so the transform used here is the one the
   * coefficients were built with, refreshed between passes. See §14.
   */
  endRow(e: number, later: boolean, sign: number, emit: (column: number, value: number) => void) {
    const { p, block, frames } = this;
    const dof = this.options.dof;
    const frame = frames[e];
    const far = seesTurning(frame, later);

    let entries: Float64Array;

    if (this.options.jacobian === "frozen") {
      if (!far) {
        return;
      }

      entries = dof.turningRow(p, frame);
    } else {
      const w = this.weights;

      w.set(this.jacobians[e].dKth);

      if (far) {
        const turning = dof.turningWeights(p);

        for (let c = 0; c < w.length; c++) {
          w[c] += turning[c];
        }
      }

      entries = dof.pullback(p, frame, w);
    }

    for (let k = 0; k < 2 * block; k++) {
      emit(this.at(e, k), sign * entries[k]);
    }
  }

  /** Solve this chain on its own. Shorthand for a {@link ComponentSystem} with no interface. */
  run(): ComponentRun {
    return new ComponentSystem([this], []).run();
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
  faults(run: ComponentRun, broken?: Map<SolvableVertex, number>): Fault[] {
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
  diagnose(chain: number, run: ComponentRun, into: Diagnostic[]) {
    const limits = diagnosticThresholds;
    const opts = this.options;

    this.chordDiagnostics(chain, into);

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
   * The half of {@link diagnose} that is about this chain rather than about the run.
   *
   * A component's chains share one Newton loop, so the run-level records — the factorization,
   * the rate, the line search — are one event and are reported once, against the component's
   * first chain. These are per-chain measurements and are reported for every chain.
   */
  localDiagnostics(chain: number, into: Diagnostic[]) {
    this.chordDiagnostics(chain, into);
    this.multipliers(chain, into);

    return into;
  }

  /** Edges whose canonical chord has collapsed, or which the quadrature has lost outright. */
  chordDiagnostics(chain: number, into: Diagnostic[]) {
    const limits = diagnosticThresholds;

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
   * `zᵀHz`, using the element Hessians as of the last {@link assemble}.
   *
   * Reading `H` from the cached elements rather than from the band keeps the multiplier slots
   * out of it — `kkt.apply` would fold `Jᵀλ` into the same numbers.
   */
  energy() {
    const d = 2 * this.block;
    let total = 0.0;

    for (let e = 0; e < this.elements.length; e++) {
      const h = this.elements[e];

      for (let r = 0; r < d; r++) {
        const zr = this.z[this.at(e, r)];

        for (let c = 0; c < d; c++) {
          total += h[r * d + c] * zr * this.z[this.at(e, c)];
        }
      }
    }

    return total;
  }

  /** `Σ|c_i|` over this chain's interior joints. */
  infeasibility() {
    let total = 0.0;

    for (let i = 1; i < this.frames.length; i++) {
      total += Math.abs(this.residuals[i]);
    }

    return total;
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
 * Several chains solved as one system, coupled only through the node DOF they share — §5's
 * substructuring, and the whole of Phase 8.
 *
 * ## The split
 *
 * Order the unknowns as every chain's private DOF followed by the interface `γ`, and the
 * global KKT matrix is
 *
 * ```
 *   ⎡ A₁        B₁ ⎤        A_i  chain `i`'s band, on its interior unknowns only
 *   ⎢    A₂     B₂ ⎥        B_i  what chain `i` contributes to the interface
 *   ⎢       ⋱   ⋮  ⎥        D    the interface's own block, summed over chains
 *   ⎣ B₁ᵀ B₂ᵀ ⋯ D  ⎦             plus the node G1 rows
 * ```
 *
 * with `S = D − Σ B_iᵀ A_i⁻¹ B_i` the Schur complement, solved densely. `A_i` is still banded
 * and still quasi-definite: a principal submatrix taken in increasing index order has
 * `|pos(a) − pos(b)| ≤ |a − b|`, so it inherits the chain's bandwidth, and a principal
 * submatrix of an SPD (1,1) block is SPD while every chain-interior multiplier row stays a
 * multiplier row. So the expensive half is unchanged and `math/banded.ts` still applies to it;
 * only `S` needs a general method, and `S` is `O(#shared groups × (p+1))` on a side.
 *
 * ## What lands in `γ`
 *
 * Both of §5's coupling caveats, which `junctions.ts` derives: a block entry two or more
 * edge-ends read out of the same DOF, *and* a node G1 row's multiplier. The second is easy to
 * forget because it shares no block — a G1-only junction couples its chains through a single
 * row and nothing else — and forgetting it is silent.
 *
 * A node row's multiplier has to be an interface unknown even though a chain-interior one need
 * not be: it touches the far block of each incident edge, which is interior to a *different*
 * chain, so there is no chain whose band could hold it.
 *
 * ## What does not change
 *
 * The chains keep the layout, band, assembly and DOF numbering they had in Phase 2. An
 * interface unknown is an ordinary local slot in every chain that reads it, mirrored between
 * them by {@link sync}; {@link ChainSystem.localToGlobal} is the only new state, and it is
 * read only here. With an empty interface — one chain, or a mesh whose junctions are all
 * fully split, which §3's default at valence ≥ 3 makes the common case — the Schur complement
 * is `0 × 0` and the step reduces to exactly the single banded solve it always was.
 */
export class ComponentSystem {
  options: SPowerSolverOptions;

  /** Global index of each local unknown that is private, `-1` for the interface ones. */
  interior: Int32Array[] = [];

  /** Where each chain's private unknowns start in the global vector, and how many there are. */
  offsets: number[] = [];
  counts: number[] = [];

  /** First global index of the interface block, and the total unknown count. */
  gammaAt = 0;
  n = 0;

  /** Current global state — chain-private unknowns then `γ`. Mirrored into the chains by {@link take}. */
  z: Float64Array;

  /** Interface values, and the `γ` part of {@link z}. */
  zGamma: Float64Array;

  /** Global index of every multiplier: each chain's, renumbered, then the node rows'. */
  multiplierRows: number[] = [];

  /** Node G1 gaps as of the last {@link measure}, indexed by `γ` slot. `0` at an entry slot. */
  nodeResiduals: Float64Array;

  /** `A_i`, `B_i` row-major `counts[i] × γ`, and the dense `D`. */
  blocks: BandedSymmetric[] = [];
  coupling: Float64Array[] = [];
  d: Float64Array;

  /** The Schur complement, its factorization, and `A_i⁻¹B_i` per chain. */
  schur: Float64Array;
  lu?: DenseLU;
  ys: Float64Array[] = [];

  private factors: KKTFactorization[] = [];
  private rowsOf: number[][] = [];
  private solution: Float64Array[] = [];
  private column: Float64Array[] = [];
  private reduced: Float64Array;
  private work: Float64Array;
  private correction: Float64Array;

  constructor(
    public systems: ChainSystem[],
    public gamma: GammaSlot[]
  ) {
    this.options = systems[0].options;

    for (const s of systems) {
      s.localToGlobal.fill(-1);
      s.localSign.fill(1.0);
    }

    for (let g = 0; g < gamma.length; g++) {
      const slot = gamma[g];

      if (slot.kind !== "entry") {
        continue;
      }

      for (const end of slot.ends) {
        const s = systems[end.chain];
        const at = s.endSlot(end.end, slot.order);

        s.localToGlobal[at] = g;
        s.localSign[at] = ComponentSystem.orientation(end, slot.order);
      }
    }

    const m = gamma.length;
    let at = 0;

    for (const s of systems) {
      const map = new Int32Array(s.n).fill(-1);
      let count = 0;

      for (let i = 0; i < s.n; i++) {
        if (s.localToGlobal[i] < 0) {
          map[i] = count++;
        }
      }

      this.interior.push(map);
      this.offsets.push(at);
      this.counts.push(count);
      this.rowsOf.push(s.rows.map((row) => map[row]));

      for (const row of s.rows) {
        this.multiplierRows.push(at + map[row]);
      }

      this.blocks.push(new BandedSymmetric(count, Math.min(s.kkt.bandwidth, Math.max(0, count - 1))));
      this.coupling.push(new Float64Array(count * m));
      this.ys.push(new Float64Array(count * m));
      this.solution.push(new Float64Array(count));
      this.column.push(new Float64Array(count));

      at += count;
    }

    this.gammaAt = at;
    this.n = at + m;

    for (let g = 0; g < m; g++) {
      if (gamma[g].kind === "row") {
        this.multiplierRows.push(at + g);
      }
    }

    this.z = new Float64Array(this.n);
    this.zGamma = new Float64Array(m);
    this.nodeResiduals = new Float64Array(m);
    this.d = new Float64Array(m * m);
    this.schur = new Float64Array(m * m);
    this.reduced = new Float64Array(m);
    this.work = new Float64Array(this.n);
    this.correction = new Float64Array(this.n);
  }

  /**
   * The sign relating entry `order` at `end` to the interface value it mirrors — §5's
   * orientation, the caveat `transformEntry` handles inside a chain.
   *
   * A vertex block holds `Rᵥⁿ/n! · dⁿκ/dsⁿ` measured along the direction the chain travels
   * *through* the vertex. Reverse that direction and `κ` flips sign while `d/ds` flips too, so
   * entry `n` picks up `(−1)^{n+1}` and the odd entries come through unchanged. Which ends are
   * reversed is {@link ChainEnd.reversed}, decided per group in `junctions.ts`.
   *
   * Getting this wrong is quiet: the tangent row still closes, because it carries its own `π`,
   * and the solve still converges — onto `κ_a = −κ_b`, a curvature *jump* twice the size of the
   * one the pairing asked to remove.
   */
  static orientation(end: ChainEnd, order: number) {
    return end.reversed && (order & 1) === 0 ? -1.0 : 1.0;
  }

  /** The chain system, edge and end-orientation a {@link ChainEnd} names. */
  endOf(end: ChainEnd) {
    const s = this.systems[end.chain];

    return { system: s, edge: s.endEdge(end.end), later: end.end === 1 };
  }

  /** Push every interface value out to the chain-local slots that mirror it. */
  sync() {
    for (const s of this.systems) {
      for (let i = 0; i < s.n; i++) {
        const g = s.localToGlobal[i];

        if (g >= 0) {
          s.z[i] = s.localSign[i] * this.zGamma[g];
        }
      }
    }
  }

  /**
   * Take the interface values *from* the chains, averaging where they disagree, then mirror
   * the result back.
   *
   * Disagreement is not an error state: §9's continuation seeds each chain separately and the
   * two ends at a junction land on different numbers for the same shared entry, exactly as the
   * two edges at a joint do within one chain. Averaging is what {@link ChainSystem.seedFrom}
   * already does about it, for the same reason and by the same argument.
   */
  gather() {
    const counts = new Float64Array(this.gamma.length);

    this.zGamma.fill(0.0);

    for (const s of this.systems) {
      for (let i = 0; i < s.n; i++) {
        const g = s.localToGlobal[i];

        if (g >= 0) {
          this.zGamma[g] += s.localSign[i] * s.z[i];
          counts[g]++;
        }
      }
    }

    for (let g = 0; g < counts.length; g++) {
      if (counts[g] > 1) {
        this.zGamma[g] /= counts[g];
      }
    }

    this.sync();
    this.read(this.z);
  }

  /** Collect the current state into a global vector. */
  read(out: Float64Array) {
    for (let c = 0; c < this.systems.length; c++) {
      const s = this.systems[c];
      const map = this.interior[c];
      const off = this.offsets[c];

      for (let i = 0; i < s.n; i++) {
        if (map[i] >= 0) {
          out[off + map[i]] = s.z[i];
        }
      }
    }

    out.set(this.zGamma, this.gammaAt);

    return out;
  }

  /** Distribute a global vector back into the chains and the interface. */
  take(x: Float64Array) {
    for (let c = 0; c < this.systems.length; c++) {
      const s = this.systems[c];
      const map = this.interior[c];
      const off = this.offsets[c];

      for (let i = 0; i < s.n; i++) {
        if (map[i] >= 0) {
          s.z[i] = x[off + map[i]];
        }
      }
    }

    for (let g = 0; g < this.zGamma.length; g++) {
      this.zGamma[g] = x[this.gammaAt + g];
    }

    this.sync();
  }

  write() {
    for (const s of this.systems) {
      s.write();
    }
  }

  /**
   * Measure every chain's joints and every node's, and return the worst gap.
   *
   * A node gap is `wrap(outgoing(a) − outgoing(b) − π)`, which is the same quantity a chain
   * joint measures — see {@link ChainSystem.outgoing} — written without reference to which
   * edge arrives and which leaves, since at a node neither does.
   */
  measure() {
    let worst = 0.0;

    for (const s of this.systems) {
      worst = Math.max(worst, s.measure());
    }

    for (let g = 0; g < this.gamma.length; g++) {
      const slot = this.gamma[g];

      if (slot.kind !== "row") {
        this.nodeResiduals[g] = 0.0;
        continue;
      }

      const a = this.endOf(slot.a);
      const b = this.endOf(slot.b);

      const gap = wrapAngle(a.system.outgoing(a.edge, a.later) - b.system.outgoing(b.edge, b.later) - Math.PI);

      this.nodeResiduals[g] = gap;
      worst = Math.max(worst, Math.abs(gap));
    }

    return worst;
  }

  /** `½zᵀHz + μ·Σ|c_i|` over the whole component — §5's `ℓ1` merit, summed. */
  merit(mu: number) {
    let energy = 0.0;
    let infeasible = 0.0;

    for (const s of this.systems) {
      energy += s.energy();
      infeasible += s.infeasibility();
    }

    for (const gap of this.nodeResiduals) {
      infeasible += Math.abs(gap);
    }

    return 0.5 * energy + mu * infeasible;
  }

  /** `max|λ|` over every multiplier in the component, chain-interior and node alike. */
  largestMultiplier() {
    let worst = 0.0;

    for (const row of this.multiplierRows) {
      worst = Math.max(worst, Math.abs(this.z[row]));
    }

    return worst;
  }

  /**
   * Assemble each chain, then split its band into `A_i`, `B_i` and `D`, and add the node rows.
   *
   * The split is by lookup rather than by re-derivation: a chain assembles into its own local
   * numbering exactly as before and every entry is then routed by {@link
   * ChainSystem.localToGlobal}. Two ends reading the same interface DOF therefore *sum* into
   * the same row and column of `D`, which is what identifying two variables means.
   */
  build() {
    const m = this.gamma.length;

    this.d.fill(0.0);

    for (let c = 0; c < this.systems.length; c++) {
      const s = this.systems[c];

      s.assemble();

      const a = this.blocks[c].zero();
      const b = this.coupling[c].fill(0.0);
      const map = this.interior[c];
      const bw = s.kkt.bandwidth;

      for (let j = 0; j < s.n; j++) {
        const hi = Math.min(s.n - 1, j + bw);

        for (let i = j; i <= hi; i++) {
          const v = s.kkt.get(i, j);

          if (v === 0.0) {
            continue;
          }

          const gi = s.localToGlobal[i];
          const gj = s.localToGlobal[j];

          if (gi < 0 && gj < 0) {
            a.add(map[i], map[j], v);
          } else if (gi < 0) {
            b[map[i] * m + gj] += s.localSign[j] * v;
          } else if (gj < 0) {
            b[map[j] * m + gi] += s.localSign[i] * v;
          } else {
            const both = s.localSign[i] * s.localSign[j] * v;

            this.d[gi * m + gj] += both;

            if (i !== j) {
              this.d[gj * m + gi] += both;
            }
          }
        }
      }
    }

    for (let g = 0; g < m; g++) {
      const slot = this.gamma[g];

      if (slot.kind === "row") {
        this.nodeRow(g, slot.a, 1.0);
        this.nodeRow(g, slot.b, -1.0);
      }
    }
  }

  /** Scatter one end's contribution to node row `g`, by where each column's unknown lives. */
  nodeRow(g: number, at: ChainEnd, sign: number) {
    const m = this.gamma.length;
    const { system, edge, later } = this.endOf(at);
    const map = this.interior[at.chain];
    const b = this.coupling[at.chain];

    system.endRow(edge, later, sign, (column, value) => {
      const local = system.localToGlobal[column];

      if (local < 0) {
        b[map[column] * m + g] += value;
      } else {
        const mirrored = system.localSign[column] * value;

        this.d[local * m + g] += mirrored;
        this.d[g * m + local] += mirrored;
      }
    });
  }

  /** `out = K x` for the assembled global matrix — what the global refinement measures against. */
  applyGlobal(x: Float64Array, out: Float64Array) {
    const m = this.gamma.length;

    out.fill(0.0);

    for (let c = 0; c < this.systems.length; c++) {
      const off = this.offsets[c];
      const k = this.counts[c];
      const b = this.coupling[c];
      const xi = x.subarray(off, off + k);
      const ai = this.solution[c];

      this.blocks[c].apply(xi, ai);

      for (let i = 0; i < k; i++) {
        out[off + i] += ai[i];

        const row = i * m;

        for (let g = 0; g < m; g++) {
          const v = b[row + g];

          if (v !== 0.0) {
            out[off + i] += v * x[this.gammaAt + g];
            out[this.gammaAt + g] += v * xi[i];
          }
        }
      }
    }

    for (let r = 0; r < m; r++) {
      let sum = 0.0;

      for (let c = 0; c < m; c++) {
        sum += this.d[r * m + c] * x[this.gammaAt + c];
      }

      out[this.gammaAt + r] += sum;
    }

    return out;
  }

  /** Factor every `A_i`, then form and factor `S`. False if any factorization underflowed. */
  factor() {
    const m = this.gamma.length;

    this.factors.length = 0;

    for (let c = 0; c < this.systems.length; c++) {
      const f = new KKTFactorization(this.blocks[c], this.rowsOf[c], this.options.kkt);

      if (!f.ok) {
        return false;
      }

      this.factors.push(f);
    }

    if (m === 0) {
      return true;
    }

    this.schur.set(this.d);

    for (let c = 0; c < this.systems.length; c++) {
      const k = this.counts[c];
      const b = this.coupling[c];
      const y = this.ys[c];
      const col = this.column[c];
      const sol = this.solution[c];

      for (let g = 0; g < m; g++) {
        for (let i = 0; i < k; i++) {
          col[i] = b[i * m + g];
        }

        this.factors[c].solve(col, sol);

        for (let i = 0; i < k; i++) {
          y[i * m + g] = sol[i];
        }
      }

      for (let r = 0; r < m; r++) {
        for (let cc = 0; cc < m; cc++) {
          let sum = 0.0;

          for (let i = 0; i < k; i++) {
            sum += b[i * m + r] * y[i * m + cc];
          }

          this.schur[r * m + cc] -= sum;
        }
      }
    }

    this.lu = new DenseLU(this.schur, m);

    return this.lu.ok;
  }

  /** One pass of the substructured solve, against the factorization {@link factor} left behind. */
  solveOnce(b: Float64Array, out: Float64Array) {
    const m = this.gamma.length;

    for (let c = 0; c < this.systems.length; c++) {
      const off = this.offsets[c];
      const k = this.counts[c];

      this.factors[c].solve(b.subarray(off, off + k), this.solution[c]);
      out.set(this.solution[c], off);
    }

    if (m === 0) {
      return;
    }

    const rs = this.reduced;

    for (let g = 0; g < m; g++) {
      rs[g] = b[this.gammaAt + g];
    }

    for (let c = 0; c < this.systems.length; c++) {
      const off = this.offsets[c];
      const k = this.counts[c];
      const bc = this.coupling[c];

      for (let i = 0; i < k; i++) {
        const xi = out[off + i];

        for (let g = 0; g < m; g++) {
          rs[g] -= bc[i * m + g] * xi;
        }
      }
    }

    this.lu?.solveInPlace(rs);

    for (let g = 0; g < m; g++) {
      out[this.gammaAt + g] = rs[g];
    }

    for (let c = 0; c < this.systems.length; c++) {
      const off = this.offsets[c];
      const k = this.counts[c];
      const y = this.ys[c];

      for (let i = 0; i < k; i++) {
        let sum = 0.0;

        for (let g = 0; g < m; g++) {
          sum += y[i * m + g] * rs[g];
        }

        out[off + i] -= sum;
      }
    }
  }

  /**
   * Solve the global step, refining against the assembled matrix.
   *
   * With no interface this is `solveKKT` on the one chain, down to the object doing the
   * factoring, so a single chain takes the Phase 2 path and reports the Phase 2 refinement
   * residual — in equilibrated units, as {@link KKTSolveResult.residual} always is. The
   * substructured path cannot report those units, having several equilibrations and a dense
   * complement between it and the answer, so its residual is the plain `‖Kx − b‖∞`.
   */
  solve(b: Float64Array, out: Float64Array) {
    if (!this.factor()) {
      out.fill(0.0);

      return { residual: Infinity, ok: false };
    }

    if (this.gamma.length === 0) {
      return { residual: this.factors[0].solve(b, out).residual, ok: true };
    }

    const opts = { ...defaultKKTSolveOptions, ...this.options.kkt };

    this.solveOnce(b, out);

    let residual = Infinity;

    for (let pass = 0; pass <= opts.refinement; pass++) {
      this.applyGlobal(out, this.work);

      residual = 0.0;

      for (let i = 0; i < this.n; i++) {
        this.work[i] = b[i] - this.work[i];
        residual = Math.max(residual, Math.abs(this.work[i]));
      }

      if (pass === opts.refinement) {
        break;
      }

      this.solveOnce(this.work, this.correction);

      for (let i = 0; i < this.n; i++) {
        out[i] += this.correction[i];
      }
    }

    return { residual, ok: true };
  }

  /** `−(Kz)` with the multiplier rows replaced by `−c` — the Newton right-hand side. */
  gradient(rhs: Float64Array) {
    this.applyGlobal(this.z, rhs);

    for (let i = 0; i < this.n; i++) {
      rhs[i] = -rhs[i];
    }

    for (let c = 0; c < this.systems.length; c++) {
      const s = this.systems[c];
      const map = this.interior[c];
      const off = this.offsets[c];

      for (let i = 1; i < s.frames.length; i++) {
        if (s.lambdaAt[i] >= 0) {
          rhs[off + map[s.lambdaAt[i]]] = -s.residuals[i];
        }
      }
    }

    for (let g = 0; g < this.gamma.length; g++) {
      if (this.gamma[g].kind === "row") {
        rhs[this.gammaAt + g] = -this.nodeResiduals[g];
      }
    }

    return rhs;
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
   *
   * The search is over the *component*, not over one chain at a time. A per-chain search would
   * be a different algorithm: a coupled interface makes one direction, and cutting one chain's
   * share of it while leaving another's would leave the shared DOF holding two values.
   */
  run(): ComponentRun {
    const opts = this.options;
    const rhs = new Float64Array(this.n);
    const step = new Float64Array(this.n);
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

    this.gather();
    this.write();

    let residual = this.measure();

    for (let iter = 0; ; iter++) {
      if (iter >= opts.iterations || residual < opts.tolerance || this.multiplierRows.length === 0) {
        break;
      }

      this.build();
      this.gradient(rhs);

      const solved = this.solve(rhs, step);

      refinement = Math.max(refinement, solved.residual);

      if (!solved.ok) {
        factored = false;
        break;
      }

      for (const row of this.multiplierRows) {
        mu = Math.max(mu, 2.0 * Math.abs(this.z[row] + step[row]));
      }

      base.set(this.z);

      const start = this.merit(mu);
      let t = opts.relaxation;
      let cuts = 0;

      for (;;) {
        for (let i = 0; i < this.n; i++) {
          this.z[i] = base[i] + t * step[i];
        }

        this.take(this.z);
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
      for (const s of this.systems) {
        for (let i = 0; i < s.frames.length; i++) {
          s.frames[i].arclength = s.arclength[i];
        }
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
        stepNorm     : t * infinityNorm(step),
        maxMultiplier: largest,
        relaxation   : t,
        refinement   : solved.residual,
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
}

/**
 * Fits every {@link SPowerClothoid} in a mesh so the segments meet tangentially.
 *
 * The mesh is cut into chains, the chains grouped into the components of §5 by what they
 * share at their ends, and each component solved as one system. Two chains are in different
 * components exactly when no vertex block and no G1 row spans them, which is what §5 requires
 * before independent solves are licensed — and which §3's default at valence ≥ 3 makes true
 * of every unauthored junction, so the ordinary mesh still comes apart into one chain per
 * component. `Rᵥ` couples every incident edge at every node and is not a component boundary,
 * benignly: it is built from chord lengths and so does not move during the solve.
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

    const all = chains(this.mesh).map(cutOpen);

    report.chains = all.length;

    for (const chain of all) {
      if (chain.closed) {
        report.cutClosed++;
      }
    }

    for (const component of components(all, this.options.order)) {
      this.solveComponent(component, refs, report);
    }

    this.report = report;

    return report;
  }

  /**
   * Solve one component, lowering a joint and retrying for as long as §8 says the answer
   * cannot be trusted.
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
   *
   * Every joint reachable this way is chain-interior. §8's localization does not yet run at a
   * junction, which is a gap and not a decision the machinery forces: the criteria are all
   * measurable there. It costs nothing today because a junction §3 has not been told about is
   * fully split already, and so has no level left to lose.
   */
  solveComponent(component: Component, refs: Map<SolvableVertex, number>, report: SPowerSolverReport) {
    const opts = this.options;
    const caps = new Map<SolvableVertex, number>();
    const actions: Diagnostic[] = [];

    let { systems, run } = this.ladder(component, refs, caps, report);
    let refused = false;

    for (let attempt = 0; attempt < opts.breaks; attempt++) {
      const faults = this.rank(component, systems, run);

      if (faults.length === 0) {
        break;
      }

      if (opts.mode === "engineering") {
        refused = true;

        for (const { system, index, fault } of faults) {
          actions.push(this.record(fault, "refused", index, system, system.delivered[fault.at]));
        }

        break;
      }

      // A joint already at 0 has nothing left to give; look past it rather than giving up,
      // since the fault it still reports may be the second-worst one and fixable.
      const found = faults.find(({ system, fault }) => system.delivered[fault.at] > 0);

      if (!found) {
        break;
      }

      const to = found.system.delivered[found.fault.at] - 1;

      actions.push(this.record(found.fault, "degraded", found.index, found.system, to));
      caps.set(found.fault.vertex, to);

      ({ systems, run } = this.ladder(component, refs, caps, report));
    }

    for (const { verts, edges } of component.chains) {
      for (let i = 1; i < edges.length; i++) {
        const v = verts[i];
        const cap = caps.get(v);

        if (cap === undefined) {
          this.broken.delete(v);
        } else {
          this.broken.set(v, cap);
          report.degraded++;
        }
      }
    }

    report.steps += run.steps;
    report.maxResidual = Math.max(report.maxResidual, run.residual);
    report.ok &&= run.ok && !refused;

    for (let k = 0; k < systems.length; k++) {
      const index = component.indices[k];

      report.unenforced += systems[k].unenforced();

      if (k === 0) {
        systems[k].diagnose(index, run, report.diagnostics);
      } else {
        systems[k].localDiagnostics(index, report.diagnostics);
      }
    }

    report.diagnostics.push(...actions);

    if (run.trace) {
      report.traces?.push({ chain: component.indices[0], steps: run.trace });
    }
  }

  /**
   * Every chain's faults, merged into one list in §8's order: rank, then chord, then
   * divergence — cause before symptom.
   *
   * The sort is across chains and not within one, because concatenation would let the first
   * chain's *symptom* outrank a later chain's *cause* and the whole point of the ordering is
   * that it does not. A single-chain component sorts to the order {@link ChainSystem.faults}
   * already emitted, the sort being stable.
   */
  rank(component: Component, systems: ChainSystem[], run: ComponentRun) {
    const found = systems.flatMap((system, k) =>
      system.faults(run, this.broken).map((fault) => ({ system, index: component.indices[k], fault }))
    );

    return found.sort((a, b) => FAULT_ORDER.indexOf(a.fault.condition) - FAULT_ORDER.indexOf(b.fault.condition));
  }

  /**
   * Build and solve one component at the requested order, climbing to it if §9's continuation
   * is on. Returns the top rung, which is the only one whose answer counts.
   *
   * A rung that failed is not seeded from: a warm start built out of a broken solve is worse
   * than no warm start, and the failure is the top rung's to rediscover and report. The ladder
   * is rebuilt from `p = 0` after every break rather than resumed, because a lowered level
   * changes the layout at every order and the rungs below have to answer the same question the
   * top one is being asked.
   *
   * The interface is rebuilt per rung and the *grouping* is not — see {@link components} for
   * why lowering `p` cannot split a component, only shrink its interface.
   */
  ladder(
    component: Component,
    refs: Map<SolvableVertex, number>,
    caps: Map<SolvableVertex, number>,
    report: SPowerSolverReport
  ) {
    const opts = this.options;
    const first = opts.continuation ? 0 : opts.order;

    let systems!: ChainSystem[];
    let run!: ComponentRun;
    let seed: ChainSystem[] | undefined;

    for (let p = first; p <= opts.order; p++) {
      const rung = p === opts.order ? opts : { ...opts, order: p };

      systems = component.chains.map((chain) => new ChainSystem(chain, refs, rung, caps));

      if (seed) {
        for (let k = 0; k < systems.length; k++) {
          systems[k].seedFrom(seed[k]);
        }
      }

      const slots = p === opts.order ? component.gamma : interfaceSlots(component.chains, p);

      run = new ComponentSystem(systems, slots).run();

      if (p < opts.order) {
        report.seedSteps += run.steps;

        // §8's criteria are exactly "this answer cannot be trusted", and a warm start is a way
        // of trusting one. A rung that trips any of them, or that stopped short of its own
        // tolerance, is not seeded from.
        const converged = run.ok && run.residual < opts.tolerance;

        seed = converged && systems.every((s) => s.faults(run).length === 0) ? systems : seed;
      }
    }

    report.attempts++;

    return { systems, run };
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
