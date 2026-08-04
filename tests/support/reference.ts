/**
 * An independent reference for the position integrals, for measuring the Taylor
 * integrator against.
 *
 * Deliberately *not* the same scheme at a high step count: refining the scheme under test
 * gives correlated errors, and the order it reports is then partly its own. This
 * integrates `cos theta` and `sin theta` directly by composite Gauss-Legendre, taking
 * `theta` from the profile's exact closed-form integral. The only shared assumption is
 * that `profile.integral` is right, which the profile's own tests cover.
 */
import { type CurvatureProfile } from "../../src/curve/profile.js";

/** Gauss-Legendre nodes and weights on `[-1, 1]`, by Newton on the Legendre polynomial. */
export function gaussLegendre(n: number) {
  const x = new Float64Array(n);
  const w = new Float64Array(n);

  for (let i = 0; i < n; i++) {
    let z = Math.cos((Math.PI * (i + 0.75)) / (n + 0.5));
    let dp = 0.0;

    for (let iter = 0; iter < 100; iter++) {
      let p0 = 1.0;
      let p1 = 0.0;

      for (let j = 0; j < n; j++) {
        const p2 = p1;
        p1 = p0;
        p0 = ((2 * j + 1) * z * p1 - j * p2) / (j + 1);
      }

      dp = (n * (z * p0 - p1)) / (z * z - 1.0);

      const dz = p0 / dp;
      z -= dz;

      if (Math.abs(dz) < 1e-16) {
        break;
      }
    }

    x[i] = z;
    w[i] = 2.0 / ((1.0 - z * z) * dp * dp);
  }

  return { x, w };
}

const gl = gaussLegendre(12);

/**
 * `(integral of cos theta, integral of sin theta)` over `[0, 1]`, to machine precision.
 *
 * `panels` should be a multiple of any knot count in the profile so that knots land on
 * panel boundaries — Gauss-Legendre converges spectrally on an analytic integrand and
 * merely algebraically across a kink.
 */
export function referenceEndpoint(profile: CurvatureProfile, ks: Float64Array, klen: number, panels: number) {
  const h = 1.0 / panels;

  let x = 0.0;
  let y = 0.0;

  for (let panel = 0; panel < panels; panel++) {
    const mid = (panel + 0.5) * h;

    for (let i = 0; i < gl.x.length; i++) {
      const u = mid + (gl.x[i] * h) / 2.0;
      const th = profile.integral(ks, klen, u);

      x += gl.w[i] * Math.cos(th);
      y += gl.w[i] * Math.sin(th);
    }
  }

  return [(x * h) / 2.0, (y * h) / 2.0] as const;
}
