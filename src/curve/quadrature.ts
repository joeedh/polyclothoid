/**
 * The Taylor position integrator, with the step count as a parameter.
 *
 * `Clothoid` carries its own copy with `QUADRATURE_STEPS` fixed at 19 and the integration
 * interval centred on the origin. This one is the general form: an explicit step count, an
 * optional fourth-order pair of terms, and an arbitrary sub-interval. It exists because
 * `docs/plans/spower-solver.md` §12's Phase 0a has to vary the step count, which the
 * private constant in `clothoid.ts` does not allow, and because Phase 3 needs a quadrature
 * loop it can hang derivative accumulators off.
 *
 * The parameter here is canonical arclength: the canonical curve has unit total arclength,
 * so `u` runs over `[0, 1]` and `du` *is* `ds`.
 */
import { Vector2 } from "../math/index.js";
import { type CurvatureProfile } from "./profile.js";

/**
 * Floor on `|C(1)|` before anything divides by it.
 *
 * A canonical curve that closes on itself has no similarity image, and the transform is
 * genuinely undefined there rather than merely large. Clamping keeps a solver step finite
 * long enough for a line search to reject it.
 */
export const MIN_CANONICAL_CHORD = 1e-9;

export interface QuadratureOptions {
  /** Number of Taylor steps over the whole `[0, 1]` interval. */
  steps: number;

  /**
   * Include the `ds^3` terms, making the scheme `O(ds^4)`.
   *
   * Requires the profile to supply `d2Curvature`; a profile without one throws rather
   * than silently dropping to third order, because a piecewise profile's second
   * derivative is a train of deltas and "zero away from the knots" is the wrong answer.
   */
  fourthOrder: boolean;
}

export const defaultQuadratureOptions: QuadratureOptions = {
  steps      : 19,
  fourthOrder: false,
};

/**
 * Total turning over `[0, 1]`, i.e. `integral of kappa` — the natural step-schedule
 * currency (`spower-solver.md` §7) and the quantity to hold fixed when comparing profiles.
 */
export function totalTurning(profile: CurvatureProfile, ks: Float64Array, klen: number) {
  return profile.integral(ks, klen, 1.0);
}

/**
 * Integrate `(cos theta, sin theta)` from `u0` to `u1`, expanding about the left end of
 * each step.
 *
 * Per step, with `t` the offset into the step:
 *
 * ```
 * theta(s + t) = theta + k t + dk t^2 / 2 + d2k t^3 / 6
 * ```
 *
 * and expanding the cosine and sine of that, then integrating `t` over `[0, ds]`:
 *
 * ```
 * dx = cos - k sin (ds/2) - (cos k^2 + dk sin)(ds^2/6)
 *    + ((k^3 - d2k) sin - 3 k dk cos)(ds^3/24)
 * dy = sin + k cos (ds/2) + (cos dk - k^2 sin)(ds^2/6)
 *    + ((d2k - k^3) cos - 3 k dk sin)(ds^3/24)
 * ```
 *
 * with the whole accumulation scaled by `ds` at the end — which is why each term carries
 * one power of `ds` less than its order. The first two lines are what `Clothoid`
 * integrates; the `ds^3` pair is `spower-solver.md` §7's fourth-order extension, and its
 * magnitude doubles as a local truncation error estimate for the scheme without it.
 *
 * `steps` counts over the full unit interval, so a partial interval takes proportionally
 * fewer — that is what keeps `ds` constant regardless of where the query lands, and hence
 * keeps `evaluate` consistent along a segment.
 */
export function integrateProfile(
  profile: CurvatureProfile,
  ks: Float64Array,
  klen: number,
  u0: number,
  u1: number,
  opts: QuadratureOptions,
  out = new Vector2()
) {
  const d2 = profile.d2Curvature;

  if (opts.fourthOrder && d2 === undefined) {
    throw new Error(`profile "${profile.name}" has no d2Curvature, so it cannot drive the fourth-order scheme`);
  }

  const steps = Math.max(1, Math.round(opts.steps * Math.abs(u1 - u0)));
  const ds = (u1 - u0) / steps;

  const ds2 = ds * ds;
  const ds3 = ds2 * ds;

  let x = 0.0;
  let y = 0.0;

  for (let i = 0; i < steps; i++) {
    const u = u0 + i * ds;

    const k = profile.curvature(ks, klen, u);
    const dk = profile.dCurvature(ks, klen, u);
    const th = profile.integral(ks, klen, u);

    const cos = Math.cos(th);
    const sin = Math.sin(th);

    x += cos - k * sin * ds * 0.5 - (cos * k * k + dk * sin) * ds2 * (1.0 / 6.0);
    y += sin + k * cos * ds * 0.5 + (cos * dk - k * k * sin) * ds2 * (1.0 / 6.0);

    if (d2 !== undefined && opts.fourthOrder) {
      const cube = k * k * k - d2(ks, klen, u);
      const cross = 3.0 * k * dk;

      x += (cube * sin - cross * cos) * ds3 * (1.0 / 24.0);
      y += (-cube * cos - cross * sin) * ds3 * (1.0 / 24.0);
    }
  }

  out[0] = x * ds;
  out[1] = y * ds;

  return out;
}

const unitCache = new Map<number, Float64Array[]>();

/** The coefficient-space basis, cached: `units[j]` selects `a_j` and nothing else. */
function unitVectors(klen: number) {
  let units = unitCache.get(klen);

  if (units === undefined) {
    units = [];

    for (let j = 0; j < klen; j++) {
      const e = new Float64Array(klen);

      e[j] = 1.0;
      units.push(e);
    }

    unitCache.set(klen, units);
  }

  return units;
}

/**
 * The canonical endpoint `C(1)` together with `∂C(1)/∂a_j` for every coefficient — §5's
 * exact Jacobian, accumulated in the same loop that produces the endpoint.
 *
 * Continuously the derivative is
 *
 * ```
 *   ∂C(1)/∂a_j = ∫₀¹ i · φ_j(v) · e^{iθ(v)} dv ,      φ_j = ∂θ/∂a_j = ∫₀^v ∂κ/∂a_j
 * ```
 *
 * but what is differentiated here is the *discrete* sum of {@link integrateProfile}, term by
 * term, not that integral requadratured. The distinction is not pedantic: the residual
 * `SPowerSolver` drives to zero is measured by the discrete scheme, so only the discrete
 * derivative is the Jacobian of the function actually being solved. Quadraturing the
 * continuous integrand instead leaves an `O(ds³)` mismatch — around `3e-4` relative at the
 * default nineteen steps — which is invisible in a finite-difference check only if the
 * tolerance is loose enough to hide it. This version agrees with central differences to
 * roundoff at any step count.
 *
 * `κ` is linear in the coefficients for every profile in this codebase, so `∂κ/∂a_j` is the
 * profile evaluated on the `j`-th unit vector, and likewise for its two derivatives and its
 * integral: writing `φ = ∂θ/∂a_j`, `ψ = ∂κ/∂a_j`, `ψ′ = ∂κ′/∂a_j`, `ψ″ = ∂κ″/∂a_j`, every
 * step's contribution differentiates by the product rule with `∂cos θ = −sin θ · φ`. That
 * keeps this generic over {@link CurvatureProfile} at the cost of `klen` extra profile
 * evaluations per step.
 */
export function integrateProfileJacobian(
  profile: CurvatureProfile,
  ks: Float64Array,
  klen: number,
  opts: QuadratureOptions,
  end = new Vector2(),
  out?: Float64Array
) {
  const d2 = profile.d2Curvature;

  if (opts.fourthOrder && d2 === undefined) {
    throw new Error(`profile "${profile.name}" has no d2Curvature, so it cannot drive the fourth-order scheme`);
  }

  const fourth = opts.fourthOrder && d2 !== undefined;

  const steps = Math.max(1, Math.round(opts.steps));
  const ds = 1.0 / steps;

  const ds2 = ds * ds;
  const ds3 = ds2 * ds;

  const units = unitVectors(klen);
  const d = out ?? new Float64Array(2 * klen);

  d.fill(0.0);

  let x = 0.0;
  let y = 0.0;

  for (let i = 0; i < steps; i++) {
    const u = i * ds;

    const k = profile.curvature(ks, klen, u);
    const dk = profile.dCurvature(ks, klen, u);
    const th = profile.integral(ks, klen, u);

    const cos = Math.cos(th);
    const sin = Math.sin(th);

    x += cos - k * sin * ds * 0.5 - (cos * k * k + dk * sin) * ds2 * (1.0 / 6.0);
    y += sin + k * cos * ds * 0.5 + (cos * dk - k * k * sin) * ds2 * (1.0 / 6.0);

    let cube = 0.0;
    let cross = 0.0;

    if (fourth) {
      cube = k * k * k - d2(ks, klen, u);
      cross = 3.0 * k * dk;

      x += (cube * sin - cross * cos) * ds3 * (1.0 / 24.0);
      y += (-cube * cos - cross * sin) * ds3 * (1.0 / 24.0);
    }

    for (let j = 0; j < klen; j++) {
      const unit = units[j];

      const phi = profile.integral(unit, klen, u);
      const psi = profile.curvature(unit, klen, u);
      const dpsi = profile.dCurvature(unit, klen, u);

      const dcos = -sin * phi;
      const dsin = cos * phi;

      let dx = dcos - (psi * sin + k * dsin) * ds * 0.5;
      let dy = dsin + (psi * cos + k * dcos) * ds * 0.5;

      dx -= (dcos * k * k + cos * 2.0 * k * psi + dpsi * sin + dk * dsin) * ds2 * (1.0 / 6.0);
      dy += (dcos * dk + cos * dpsi - 2.0 * k * psi * sin - k * k * dsin) * ds2 * (1.0 / 6.0);

      if (fourth) {
        const dcube = 3.0 * k * k * psi - d2(unit, klen, u);
        const dcross = 3.0 * (psi * dk + k * dpsi);

        dx += (dcube * sin + cube * dsin - dcross * cos - cross * dcos) * ds3 * (1.0 / 24.0);
        dy += (-dcube * cos - cube * dcos - dcross * sin - cross * dsin) * ds3 * (1.0 / 24.0);
      }

      d[j * 2] += dx;
      d[j * 2 + 1] += dy;
    }
  }

  end[0] = x * ds;
  end[1] = y * ds;

  for (let i = 0; i < 2 * klen; i++) {
    d[i] *= ds;
  }

  return { end, d };
}

/** The similarity transform an edge's coefficients imply, and how it moves with them. */
export interface TransformJacobian {
  /** `KSCALE_e`, which is the arclength `L_e`. */
  scale: number;

  /** `arg C(1)`, the canonical chord's direction. `KTH_e` is the edge's chord angle minus this. */
  angle: number;

  /** `∂L_e/∂a_j`, length `klen`. */
  dScale: Float64Array;

  /** `∂KTH_e/∂a_j`, length `klen` — already negated, so it adds into a tangent-angle row directly. */
  dKth: Float64Array;
}

/**
 * `KSCALE_e` and `KTH_e` and their coefficient derivatives, from one pass of
 * {@link integrateProfileJacobian}.
 *
 * §5: "`KTH` and `KSCALE` are then `−arg` and `1/|·|` of `C(1)` scaled by the chord, so their
 * derivatives follow from `∂C(1)` by two lines of complex arithmetic."
 *
 * ```
 *   L = chord/|C| ,      ∂L/∂a_j = −chord·(C·D_j)/|C|³
 *   KTH = θ_chord − arg C ,    ∂KTH/∂a_j = −(C.x·D_j.y − C.y·D_j.x)/|C|²
 * ```
 *
 * The chord angle itself is a constant of the input polyline and so contributes nothing.
 */
export function transformJacobian(
  profile: CurvatureProfile,
  ks: Float64Array,
  klen: number,
  chord: number,
  opts: QuadratureOptions,
  out?: TransformJacobian
): TransformJacobian {
  const { end, d } = integrateProfileJacobian(profile, ks, klen, opts);

  const x = end[0];
  const y = end[1];

  const r2 = Math.max(x * x + y * y, MIN_CANONICAL_CHORD * MIN_CANONICAL_CHORD);
  const r = Math.sqrt(r2);

  const j = out ?? { scale: 0.0, angle: 0.0, dScale: new Float64Array(klen), dKth: new Float64Array(klen) };

  j.scale = chord / r;
  j.angle = Math.atan2(y, x);

  for (let i = 0; i < klen; i++) {
    const dx = d[i * 2];
    const dy = d[i * 2 + 1];

    j.dScale[i] = (-chord * (x * dx + y * dy)) / (r2 * r);
    j.dKth[i] = -(x * dy - y * dx) / r2;
  }

  return j;
}
