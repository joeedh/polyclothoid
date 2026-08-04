/**
 * The piecewise-linear profile as a {@link ProfileDOF} — §12's control experiment.
 *
 * §1 claims the s-power basis is the win. The new solver changed four things at once
 * relative to `ClothoidSolver`, though: the basis, where the unknowns live, the energy that
 * regularizes them, and the banded KKT step. Swapping only the first of those is what tells
 * "the basis was the win" apart from "the KKT and the energy were the win," and §12 is
 * explicit that §1 is not settled without it.
 *
 * ## What is held fixed
 *
 * Everything except the representation inside an edge. The unknowns are still the vertex
 * blocks of §3 — `p + 1` scaled derivatives of world curvature per vertex — so the DOF count,
 * the constraint, the sparsity pattern and the bandwidth are identical. What changes is only
 * that an edge carries `SAMPLE_COUNT` point values of `q(u)` instead of `2p + 2` s-power
 * coefficients, and reads them back with {@link piecewiseLinear}.
 *
 * The transform is therefore the s-power one with a sampling matrix in front,
 *
 * ```
 *   T_e = S · M_e ,      S_{ij} = ψ_j(i/(m−1))
 * ```
 *
 * with `ψ_j` the `j`-th s-power basis function and `m` the sample count. Going through `M_e`
 * rather than deriving a Hermite sampling directly is deliberate: it guarantees the two runs
 * differ *only* by `S`, so any difference in the measurements is attributable to the sampling
 * and to nothing else. `S` has full column rank for `m ≥ 2p + 2`, so no DOF is lost.
 *
 * ## What that costs, predicted
 *
 * Two things, both of which §1 asserts and this module exists to measure.
 *
 * 1. **Quadrature order.** `piecewiseLinear` has no `d2Curvature`, and cannot: it is zero on
 *    every open interval and a delta at every knot. So the fourth-order terms of §7 are
 *    unavailable, and `docs/research/clothoids.md` §4 already measures the knots as the cap
 *    on the observed order of the integrator. This is not a tuning difference — the option
 *    is not reachable.
 * 2. **Interpolation error.** `q` is a polynomial of degree `2p + 1`; sampling it and
 *    rejoining the samples with straight lines loses `O(h²·q″)`. That is a *bias*, not a
 *    solver failure: the solve converges to a curve whose profile is the polyline, and the
 *    G1 residual goes to zero all the same.
 */
import {
  type EdgeFrame,
  type ProfileDOF,
  blockLength,
  congruence,
  edgeTransform,
  matrixTimesVector,
  multiply,
  rowTimesMatrix,
} from "./blocks.js";
import { type QuadratureOptions, defaultQuadratureFor } from "./quadrature.js";
import { piecewiseLinear } from "./profile.js";
import { evalSPower, sPowerLength } from "./spower.js";

/**
 * Samples an edge carries, matching `KORDER` in `clothoid.ts`.
 *
 * The same twelve the existing solver uses, so the control is the existing profile at the
 * existing resolution rather than a straw man. Note that it is not a DOF count here — it is
 * a *rendering* of `2p + 2` DOF, and twelve of them describe four numbers.
 */
export const SAMPLE_COUNT = 12;

const sampleCache = new Map<string, Float64Array>();

/**
 * The `m × n` matrix evaluating the s-power basis at `m` uniform nodes over `[0, 1]`.
 *
 * Column `j` is the `j`-th basis function sampled, which is what {@link evalSPower} returns
 * on the `j`-th unit vector.
 */
export function sampleMatrix(p: number, samples = SAMPLE_COUNT) {
  const key = `${p}:${samples}`;
  const hit = sampleCache.get(key);

  if (hit !== undefined) {
    return hit;
  }

  const n = sPowerLength(p);
  const s = new Float64Array(samples * n);
  const unit = new Float64Array(n);

  for (let j = 0; j < n; j++) {
    unit.fill(0.0);
    unit[j] = 1.0;

    for (let i = 0; i < samples; i++) {
      s[i * n + j] = evalSPower(unit, n, i / (samples - 1));
    }
  }

  sampleCache.set(key, s);

  return s;
}

/** `T_e = S · M_e`, the `m × (2p+2)` map from vertex DOF to sampled canonical curvature. */
export function sampleTransform(p: number, frame: EdgeFrame, samples = SAMPLE_COUNT, out?: Float64Array) {
  const n = sPowerLength(p);

  return multiply(sampleMatrix(p, samples), edgeTransform(p, frame), samples, n, n, out);
}

/**
 * P1 finite-element `(K + α²·M) / C_e³` over `m` uniform nodes, row-major `m²`.
 *
 * The same energy as {@link edgeEnergy}, assembled against hat functions instead of the
 * s-power basis: the standard element pair `(1/h)[[1,−1],[−1,1]]` for `∫q′²` and
 * `(h/6)[[2,1],[1,2]]` for `∫q²`. Both are exact for the piecewise-linear interpolant, so
 * this is the true energy *of the represented profile* rather than an approximation of the
 * s-power one — which is the honest comparison, since the represented profile is what gets
 * integrated.
 */
export function sampleEnergy(chord: number, alpha: number, samples = SAMPLE_COUNT, out?: Float64Array) {
  const e = out ?? new Float64Array(samples * samples);

  const h = 1.0 / (samples - 1);
  const w = 1.0 / (chord * chord * chord);
  const a2 = alpha * alpha;

  e.fill(0.0);

  for (let i = 0; i < samples - 1; i++) {
    const k = w / h;
    const m = (w * a2 * h) / 6.0;

    e[i * samples + i] += k + 2.0 * m;
    e[(i + 1) * samples + i + 1] += k + 2.0 * m;
    e[i * samples + i + 1] += -k + m;
    e[(i + 1) * samples + i] += -k + m;
  }

  return e;
}

/**
 * Trapezoid weights, which integrate the piecewise-linear interpolant *exactly*.
 *
 * Not an approximation and not a choice: `∫₀¹` of a polyline through the samples is the
 * trapezoid sum, and this has to agree with `piecewiseLinear.integral` to the last bit or
 * the constraint row and the measured residual would be describing different curves.
 */
export function trapezoidWeights(samples = SAMPLE_COUNT) {
  const h = 1.0 / (samples - 1);
  const w = new Float64Array(samples).fill(h);

  w[0] = h * 0.5;
  w[samples - 1] = h * 0.5;

  return w;
}

/**
 * The piecewise-linear profile as a {@link ProfileDOF}, over `samples` nodes per edge.
 *
 * A factory rather than a constant so the sample count can be swept — §14's control table
 * reports what happens as it rises, which is the quantitative form of "the polyline is the
 * only thing that changed".
 */
export function sampledDOF(samples = SAMPLE_COUNT): ProfileDOF {
  const weights = trapezoidWeights(samples);
  const quadrature: QuadratureOptions = defaultQuadratureFor(piecewiseLinear);

  return {
    name   : `piecewise-linear(${samples})`,
    profile: piecewiseLinear,
    quadrature,

    blockLength,

    coefficientLength() {
      return samples;
    },

    hessian(p, frame, alpha, out) {
      const t = sampleTransform(p, frame, samples);

      return congruence(sampleEnergy(frame.chord, alpha, samples), t, samples, sPowerLength(p), out);
    },

    turningRow(p, frame, out) {
      return this.pullback(p, frame, weights, out);
    },

    turningWeights() {
      return weights;
    },

    pullback(p, frame, w, out) {
      return rowTimesMatrix(w, sampleTransform(p, frame, samples), samples, sPowerLength(p), out);
    },

    coefficients(p, frame, dof, out) {
      return matrixTimesVector(sampleTransform(p, frame, samples), dof, samples, sPowerLength(p), out);
    },
  };
}

/** The control at the default resolution. Swap it into `SPowerSolverOptions.dof`. */
export const sampleDOF = sampledDOF();
