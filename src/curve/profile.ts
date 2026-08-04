/**
 * Curvature profiles: the sampled `k(s)` a clothoid segment integrates to recover its
 * geometry.
 *
 * Split out of `clothoid.ts` so that anything integrating a profile — the segment itself,
 * the standalone quadrature harness, and eventually the s-power solver — can reach them
 * without importing the curve. `docs/plans/spower-solver.md` §10 is the reason: the
 * solver wants to be generic over *scalar profiles*, with the piecewise-linear one as
 * just another member.
 */
import { binomial, clamp, fract } from "../math/index.js";
import { evalSPower, sPowerCurvature2, sPowerDerivative, sPowerIntegral } from "./spower.js";

/**
 * How a curvature profile is sampled between its control values.
 *
 * The members must be mutually consistent: `curvature` is the profile, `dCurvature` its
 * derivative, `d2Curvature` its second derivative, and `integral` its exact integral from
 * 0 to s. Getting `integral` wrong bends the curve without any visible error in the
 * curvature plot, so prefer exact integrals over numeric ones.
 *
 * `d2Curvature` is optional because a piecewise profile does not have one in any useful
 * sense — it is a sum of deltas at the knots. It is what the fourth-order quadrature
 * terms need, so a profile without it cannot drive that scheme; see
 * {@link integrateProfile}.
 *
 * `klen` is the length of the meaningful prefix of `ks`, not `ks.length`: for the sampled
 * profiles it is the sample count, and for {@link sPowerProfile} it is the s-power
 * coefficient count.
 */
export interface CurvatureProfile {
  readonly name: string;
  dCurvature(ks: Float64Array, klen: number, s: number): number;
  curvature(ks: Float64Array, klen: number, s: number): number;
  d2Curvature?(ks: Float64Array, klen: number, s: number): number;
  integral(ks: Float64Array, klen: number, s: number): number;
}

/*
  `t` is derived from the same `i1` used to index, so the piece boundary is continuous and
  needs no snapping epsilon. Deriving it independently (`fract` of the unrounded index)
  pairs a t near 1 with an i1 already bumped past the knot.
*/
function linearCurvature(ks: Float64Array, klen: number, s: number) {
  const si = s * (klen - 1);
  const i1 = ~~si;
  const i2 = i1 + 1;

  if (i2 <= klen - 1) {
    return ks[i1] + (ks[i2] - ks[i1]) * clamp(si - i1, 0.0, 1.0);
  }

  return ks[i1];
}

/** Exact integral of a linear ramp from `a` to `b` over `[0, s]`, with s in `[0, 1]`. */
function rampIntegral(a: number, b: number, s: number) {
  return -((s - 2.0) * a - b * s) * s * 0.5;
}

/**
 * Piecewise-linear curvature: interpolate linearly between samples.
 *
 * This is the profile the existing solver runs on. The integral is computed exactly by
 * trapezoid accumulation rather than by quadrature, so it carries no step error.
 *
 * No `d2Curvature`: it is zero on every open interval and a delta at every knot, and the
 * ten interior knots are exactly what `docs/research/clothoids.md` §4 measures as the
 * limit on the integrator's observed order.
 */
export const piecewiseLinear: CurvatureProfile = {
  name: "piecewise-linear",

  curvature: linearCurvature,

  dCurvature(ks, klen, s) {
    const df = 0.00001;

    return (linearCurvature(ks, klen, s + df) - linearCurvature(ks, klen, s)) / df;
  },

  integral(ks, klen, s) {
    const klen2 = klen - 1;

    const si = s * klen2;
    const i1 = ~~si;
    const i2 = clamp(i1 + 1, 0, klen2);

    let sum = 0.0;

    for (let i = 0; i < i1; i++) {
      sum += rampIntegral(ks[i], ks[i + 1], 1.0) / klen2;
    }

    if (i2 !== i1) {
      sum += rampIntegral(ks[i1], ks[i2], clamp(si - i1, 0.0, 1.0)) / klen2;
    }

    return sum;
  },
};

/**
 * Piecewise-constant curvature — a chain of circular arcs.
 *
 * **Inconsistent, and slated for removal.** `curvature` indexes intervals of width
 * `1/(klen-1)` while `integral` accumulates with `ds = 1/klen`, an 8.3% error in theta at
 * `klen = 12`. `docs/research/clothoids.md` §8 records it as a won't-fix on an unused path.
 * Do not use it as a reference shape for quadrature measurements.
 */
export const circleArc: CurvatureProfile = {
  name: "circle-arc",

  dCurvature() {
    return 0.0;
  },

  curvature(ks, klen, s) {
    return ks[~~(s * (klen - 1))];
  },

  integral(ks, klen, s) {
    let si = s * (klen - 1);
    const t = fract(si);

    si = ~~si;

    const ds = 1.0 / klen;
    let sum = 0.0;

    for (let i = 0; i < si; i++) {
      sum += ks[i] * ds;
    }

    return sum + t * ks[si] * ds;
  },
};

/**
 * s-power curvature: a smooth polynomial of degree `klen - 1`, with `klen` the s-power
 * coefficient count `2p + 2`.
 *
 * All four members are exact and closed-form, which is the representational half of
 * `docs/plans/spower-solver.md` §1: no knots, so nothing caps the integrator's order, and
 * an exact `d2Curvature` so the fourth-order quadrature terms are available both as extra
 * accuracy and as an error estimator (§7).
 */
export const sPowerProfile: CurvatureProfile = {
  name: "s-power",

  curvature  : evalSPower,
  dCurvature : sPowerDerivative,
  d2Curvature: sPowerCurvature2,
  integral   : sPowerIntegral,
};

/**
 * Bernstein-basis curvature profile.
 *
 * Incomplete: the integral is not derived, so this cannot drive `evaluate`. path.ux's
 * original also failed to weight the basis functions by `ks[i]` at all, which made it a
 * constant-1 profile regardless of input. Weighted correctly here, but still parked.
 */
export function bernsteinCurvature(ks: Float64Array, klen: number, s: number) {
  let sum = 0.0;

  for (let i = 0; i < klen; i++) {
    sum += ks[i] * binomial(klen - 1, i) * Math.pow(s, i) * Math.pow(1.0 - s, klen - 1 - i);
  }

  return sum;
}

/** Which profile the clothoid integrator uses. Swap to compare formulations. */
export let activeProfile: CurvatureProfile = piecewiseLinear;

export function setCurvatureProfile(profile: CurvatureProfile) {
  activeProfile = profile;
}
