/**
 * The solver's diagnostic channel — `docs/plans/spower-solver.md` §8.
 *
 * A diagnostic is a **located, typed record, not a boolean**: where it was measured, what
 * condition, the measured value against the threshold it was compared with, and what the
 * solver did about it. Engineering mode (§9) is what makes the shape load-bearing — if the
 * answer to a degenerate configuration is "refuse", the user has to learn *what* and
 * *where* well enough to fix it — but the stronger reason to build it is that it is the
 * debugging surface for development, which until now was the demo and a pair of eyes.
 *
 * Three rules from §8 govern what belongs here.
 *
 * 1. **Measurement is not gated on failure.** A chain that converges slowly still emits a
 *    record, because a sluggish stroke is exactly the thing you want to see and it is not a
 *    fault. What emission *is* gated on is the threshold being crossed at all — a healthy
 *    edge produces nothing, or a thousand-edge stroke would allocate a thousand records per
 *    solve. {@link DiagnosticSeverity} then grades how bad, and gates only what a client
 *    chooses to surface.
 * 2. **It is nearly free**, which is what makes (1) affordable. Every condition emitted
 *    today is a byproduct of a computation the solve already performs: the canonical chord
 *    from `_update`, the backtrack count from the line search, the refinement residual from
 *    the refinement, the residual history from the convergence test.
 * 3. **It is deterministic.** Identical input, identical diagnostics — the same requirement
 *    §8 puts on hysteresis, and for the same reason: a report that varies run to run is
 *    worse than none when bisecting.
 *
 * Per-iteration history is the thing that actually diagnoses non-convergence, and it is the
 * one part that allocates per step, so it lives behind a trace flag as {@link SolveTrace}.
 */
import { type CurveEdge, type CurveVertex } from "./curve.js";

/**
 * What was measured. One member per row of §8's table.
 *
 * Three are declared but not yet emitted, because they are the two measurements §8 admits
 * are not free plus the outcome of acting on them, and acting is Phase 6's job:
 * `g1-rank-deficiency` (a pivoted QR of the local `J` block), `ill-conditioned-hessian` (a
 * condition estimate of the `(1,1)` block), and `level-lowered`. `large-multiplier` is
 * likewise Phase 6 — §8 is explicit that `|λ_v|` is corner evidence for the client and not
 * a fault the solver may act on.
 */
export type DiagnosticCondition =
  | "chord-degeneracy"
  | "newton-not-converging"
  | "line-search-failing"
  | "refinement-not-converging"
  | "factorization-underflow"
  | "g1-rank-deficiency"
  | "ill-conditioned-hessian"
  | "level-lowered"
  | "large-multiplier";

/**
 * How bad, which per §8 gates what a client *surfaces* and never what the solver *records*.
 *
 * `info` is a measurement worth having, `warning` is a degradation the answer survived, and
 * `error` means this part of the solve did not do what was asked of it.
 */
export type DiagnosticSeverity = "info" | "warning" | "error";

/** What the solver did about it: nothing, degraded from level `k` to `j`, or refused. */
export type DiagnosticAction = "none" | "degraded" | "refused";

export interface Diagnostic {
  condition: DiagnosticCondition;
  severity: DiagnosticSeverity;
  action: DiagnosticAction;

  /**
   * Index of the chain, in the order the solve walked them, counting skipped ones.
   *
   * `-1` for a record that is not chain-local.
   */
  chain: number;

  /** Vertex or edge index within that chain, or `-1` when the record is about the chain. */
  at: number;

  vertex?: CurveVertex;
  edge?: CurveEdge;

  /** The measured quantity. Its units depend on {@link condition} — see {@link diagnosticThresholds}. */
  measured: number;

  /** What {@link measured} was compared with. */
  threshold: number;

  /** Pairing level asked for and pairing level delivered, when {@link action} is `"degraded"`. */
  requested?: number;
  delivered?: number;
}

/**
 * Where each condition starts being worth a record.
 *
 * These are *reporting* thresholds, not failure thresholds — §8's first rule. Each is set
 * where the measurement becomes interesting rather than where the solve breaks, and the
 * severity assigned alongside is what says which of the two happened.
 */
export const diagnosticThresholds = {
  /**
   * `|C(1)|` of the canonical curve, whose reciprocal is `L_e / C_e`. At `0.05` the segment
   * is twenty times its own chord and the similarity transform is already stiff.
   *
   * The measurement is bounded above by `1` as well — the canonical curve has unit arclength
   * — so a value over it is not a tight segment but a lost one, and is reported through the
   * same condition at `error` rather than through a threshold of its own.
   */
  chordDegeneracy: 0.05,

  /** Fitted per-step residual factor. Above this the iteration is crawling rather than converging. */
  newtonRate: 0.9,

  /** Backtracks within one line search. Four is a step cut to a sixteenth. */
  lineSearch: 4,

  /** `‖Ax − b‖∞` after iterative refinement, in equilibrated units — §5's `δ` left too much behind. */
  refinement: 1e-8,
};

/** Residuals kept for {@link geometricRate}. Six, so the fit sees a trend and not the whole run. */
export const RATE_WINDOW = 6;

/**
 * The `ρ` in `r_k ≈ A·ρ^k`, by least squares on `log r_k`.
 *
 * §8 asks for a geometric fit over at least three iterations rather than a one-step ratio,
 * because the Gauss-Newton iteration is highly non-normal: a single step can grow while the
 * sequence converges, and shrink while it diverges. So `ρ < 1` is convergence, `ρ ≈ 1` is a
 * stall, and `ρ > 1` is divergence — read off a trend rather than off one pair of numbers.
 *
 * Returns `NaN` when fewer than three positive samples are available to fit, which is not
 * the same as a rate of zero and should not be reported as one.
 */
export function geometricRate(history: ArrayLike<number>) {
  let n = 0;
  let sx = 0.0;
  let sy = 0.0;
  let sxx = 0.0;
  let sxy = 0.0;

  for (let i = 0; i < history.length; i++) {
    const r = history[i];

    if (!(r > 0.0) || !Number.isFinite(r)) {
      continue;
    }

    const y = Math.log(r);

    n++;
    sx += i;
    sy += y;
    sxx += i * i;
    sxy += i * y;
  }

  const denom = n * sxx - sx * sx;

  if (n < 3 || denom === 0.0) {
    return NaN;
  }

  return Math.exp((n * sxy - sx * sy) / denom);
}

/** One Gauss-Newton step, as it looked from inside. See {@link SolveTrace}. */
export interface TraceStep {
  /** §5's `ℓ1` merit at the accepted point. */
  merit: number;

  /** Largest `|wrapped angle gap|` over the chain's enforced joints, after the step. */
  residual: number;

  /** `‖t·Δz‖∞` — the accepted step, not the Newton step it was cut from. */
  stepNorm: number;

  /** Largest `|λ_v|` over the chain, after the step. */
  maxMultiplier: number;

  /** The `t` the line search settled on, in `(0, relaxation]`. */
  relaxation: number;

  /** `‖Ax − b‖∞` after refinement, in equilibrated units. */
  refinement: number;
}

/**
 * Per-iteration history for one chain. Opt-in, because it is the only part of the
 * diagnostics that allocates per step.
 */
export interface SolveTrace {
  /** Matches {@link Diagnostic.chain}. */
  chain: number;

  steps: TraceStep[];
}

/** What every {@link CurveSolver} hands back. */
export interface SolveReport {
  /**
   * False if the solve produced nothing usable — a factorization that failed, or geometry
   * that went non-finite.
   *
   * Not a convergence test. A chain that stalls at a visible angle gap has still produced an
   * answer, and reports that through an `error` diagnostic; `ok` stays true. Anything that
   * clears `ok` also leaves at least one located record saying where.
   */
  ok: boolean;

  diagnostics: Diagnostic[];

  /** Present only when the solve was run with tracing on. */
  traces?: SolveTrace[];
}

/** A clean report, for solvers that have nothing to measure yet. */
export function emptySolveReport(): SolveReport {
  return { ok: true, diagnostics: [] };
}
