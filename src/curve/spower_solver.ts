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
import { BandedSymmetric, type KKTSolveOptions, solveKKT } from "../math/index.js";
import { type EdgeFrame, type ProfileDOF, referenceLength, sPowerDOF } from "./blocks.js";
import { type CurveSolver, type SolvableEdge, type SolvableMesh, type SolvableVertex } from "./mesh_types.js";
import { SPOWER_ORDER, type SPowerClothoid } from "./spower_clothoid.js";
import { type TransformJacobian, transformJacobian } from "./quadrature.js";
import { type Chain, chains } from "./topology.js";

const TAU = Math.PI * 2.0;

/** Backtracking line search on the merit — see {@link ChainSystem.run}. */
const SHRINK = 0.5;
const MIN_RELAXATION = 1e-4;

/** Signed angle in `(−π, π]`. §4: the residual has to carry a sign and stay differentiable. */
function wrapAngle(a: number) {
  return a - TAU * Math.round(a / TAU);
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
  kkt       : {},
};

export interface SPowerSolverReport {
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

  /** False if any KKT factorization underflowed — `kkt.delta` is too small for the system. */
  ok: boolean;
}

function emptyReport(): SPowerSolverReport {
  return { chains: 0, skippedClosed: 0, unenforced: 0, maxResidual: 0.0, steps: 0, ok: true };
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

  /** Index of `block(verts[i])` in the solution vector. */
  blockAt: Int32Array;

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

  /** Element Hessians as of the last {@link assemble}, kept so the merit can reuse them. */
  elements: Float64Array[] = [];

  residuals: Float64Array;

  constructor(
    public chain: Chain,
    refs: Map<SolvableVertex, number>,
    public options: SPowerSolverOptions
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
    this.lambdaAt = new Int32Array(verts.length).fill(-1);

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

    return this.options.jacobian === "exact" || seesTurning(f[i - 1], true) || seesTurning(f[i], false);
  }

  /**
   * Interleave multipliers with blocks: `[b₀] [λ₁] [b₁] … [λ_{m−1}] [b_m]`.
   *
   * Placing `λ_i` immediately before `block(v_i)` is what keeps the system banded — its row
   * reaches the two adjacent blocks in either direction, and no further.
   */
  layout() {
    const m = this.chain.edges.length;
    let at = 0;

    for (let i = 0; i < this.chain.verts.length; i++) {
      if (i > 0 && i < m && this.enforced(i)) {
        this.lambdaAt[i] = at++;
        this.rows.push(this.lambdaAt[i]);
      }

      this.blockAt[i] = at;
      at += this.block;
    }

    this.n = at;
  }

  /**
   * Widest offset any assembled entry spans.
   *
   * Sized for a row touching *three* blocks even though the frozen row touches two, so that
   * Phase 3's exact Jacobian needs no change here.
   */
  bandwidth() {
    const m = this.chain.edges.length;
    const b = this.block;
    let w = 0;

    for (let i = 0; i < m; i++) {
      w = Math.max(w, this.blockAt[i + 1] + b - 1 - this.blockAt[i]);
    }

    for (let i = 1; i < m; i++) {
      const at = this.lambdaAt[i];

      if (at >= 0) {
        w = Math.max(w, at - this.blockAt[i - 1], this.blockAt[i + 1] + b - 1 - at);
      }
    }

    return Math.min(w, Math.max(0, this.n - 1));
  }

  /** Global index of local DOF `k` of edge `i`, whose DOF are its two blocks in chain order. */
  at(i: number, k: number) {
    return k < this.block ? this.blockAt[i] + k : this.blockAt[i + 1] + k - this.block;
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
    const { p, block, frames, kkt } = this;
    const dof = this.options.dof;
    const at = this.lambdaAt[i];
    const frozen = this.options.jacobian === "frozen";

    // The arriving edge enters negatively, the leaving edge positively — as in `measure`.
    for (const arriving of [true, false]) {
      const edge = arriving ? i - 1 : i;
      const frame = frames[edge];
      const sign = arriving ? -1.0 : 1.0;
      const far = seesTurning(frame, arriving);

      let row: Float64Array;

      if (frozen) {
        if (!far) {
          continue;
        }

        row = dof.turningRow(p, frame);
      } else {
        const w = this.weights;

        w.set(this.jacobians[edge].dKth);

        if (far) {
          const turning = dof.turningWeights(p);

          for (let c = 0; c < w.length; c++) {
            w[c] += turning[c];
          }
        }

        row = dof.pullback(p, frame, w);
      }

      for (let k = 0; k < 2 * block; k++) {
        kkt.add(at, this.at(edge, k), sign * row[k]);
      }
    }
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
  run() {
    const opts = this.options;
    const grad = new Float64Array(this.n);
    const rhs = new Float64Array(this.n);
    const base = new Float64Array(this.n);

    let steps = 0;
    let ok = true;
    let mu = 1.0;

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

      if (!step.ok) {
        ok = false;
        break;
      }

      for (const row of this.rows) {
        mu = Math.max(mu, 2.0 * Math.abs(this.z[row] + step.x[row]));
      }

      base.set(this.z);

      const start = this.merit(mu);
      let t = opts.relaxation;

      for (;;) {
        for (let i = 0; i < this.n; i++) {
          this.z[i] = base[i] + t * step.x[i];
        }

        this.write();
        this.measure();

        if (this.merit(mu) < start || t < MIN_RELAXATION) {
          break;
        }

        t *= SHRINK;
      }

      // `L_e` only moves here, between passes — the coefficients are rebuilt through the
      // updated `M_e` so the next linearization point is a self-consistent one.
      for (let i = 0; i < this.frames.length; i++) {
        this.frames[i].arclength = this.arclength[i];
      }

      this.write();
      residual = this.measure();

      steps++;
    }

    return { steps, residual, ok };
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

  /** Interior joints this system chose not to constrain. See the module doc comment. */
  unenforced() {
    let count = 0;

    for (let i = 1; i < this.frames.length; i++) {
      if (this.lambdaAt[i] < 0) {
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

    for (const chain of chains(this.mesh)) {
      report.chains++;

      if (chain.closed) {
        report.skippedClosed++;
        continue;
      }

      const system = new ChainSystem(chain, refs, this.options);
      const result = system.run();

      report.steps += result.steps;
      report.unenforced += system.unenforced();
      report.maxResidual = Math.max(report.maxResidual, result.residual);
      report.ok &&= result.ok;
    }

    this.report = report;

    return report;
  }
}
