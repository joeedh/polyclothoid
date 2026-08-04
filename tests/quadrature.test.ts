/**
 * Phase 0a — the quadrature gate of `docs/plans/spower-solver.md` §12.
 *
 * The question is representational and has nothing to do with the solver: does removing
 * the curvature profile's knots let the Taylor integrator reach its asymptotic `O(ds^3)`,
 * and is the error constant at low step counts competitive at realistic turning?
 * `docs/research/clothoids.md` §4 records the piecewise-linear answer — order oscillating
 * between ~1 and ~4, error not even monotonic below `N ~ 11` — and blames the ten interior
 * knots. A polynomial profile has none.
 *
 * Both representations carry the *same* underlying curvature function, rescaled so their
 * total turning matches exactly, so the comparison is of representations rather than of
 * shapes. Each is measured against its own Gauss-Legendre reference (`support/reference.ts`),
 * not against a refinement of the scheme under test.
 *
 * `circleArc` is deliberately not involved: its `curvature` and `integral` disagree by 8.3%
 * (`clothoids.md` §8), which would bias exactly this measurement.
 */
import assert from "node:assert/strict";
import { describe, it, type TestContext } from "node:test";

import { KORDER } from "../src/curve/clothoid.js";
import { integrateProfile } from "../src/curve/quadrature.js";
import { piecewiseLinear, sPowerProfile, type CurvatureProfile } from "../src/curve/profile.js";
import { sPowerLength, taylorToPairs } from "../src/curve/spower.js";
import { referenceEndpoint } from "./support/reference.js";

/**
 * A degree-7 curvature polynomial in `u`, monomial coefficients.
 *
 * Chosen to have an interior sign change and a non-constant `kappa'` — a profile that a
 * low-order scheme cannot get right by accident. Degree 7 is exactly order `p = 3` in the
 * s-power basis, so the s-power representation of it is exact rather than an approximation,
 * and the measurement is of the integrator alone.
 */
const KAPPA = [1.0, 2.4, -6.0, 3.2, 1.5, -2.0, 0.8, -0.5];

function polyEval(c: readonly number[], u: number) {
  let acc = 0.0;

  for (let i = c.length - 1; i >= 0; i--) {
    acc = acc * u + c[i];
  }

  return acc;
}

function polyDeriv(c: readonly number[]) {
  return c.slice(1).map((v, i) => v * (i + 1));
}

/** Exact `integral of KAPPA over [0, 1]`. */
function polyTurning(c: readonly number[]) {
  let sum = 0.0;

  for (let i = 0; i < c.length; i++) {
    sum += c[i] / (i + 1);
  }

  return sum;
}

function endpointTaylor(c: readonly number[], count: number) {
  const f = new Float64Array(count);
  const g = new Float64Array(count);

  let d = c;
  let fac = 1.0;

  for (let i = 0; i < count; i++) {
    f[i] = polyEval(d, 0.0) / fac;
    g[i] = ((i & 1 ? -1 : 1) * polyEval(d, 1.0)) / fac;

    d = polyDeriv(d);
    fac *= i + 1;
  }

  return { f, g };
}

interface Case {
  profile: CurvatureProfile;
  ks: Float64Array;
  klen: number;
  panels: number;
}

/** The polynomial, scaled to total turning `phi`, exactly in the s-power basis. */
function sPowerCase(phi: number): Case {
  const scale = phi / polyTurning(KAPPA);
  const c = KAPPA.map((v) => v * scale);

  const p = 3;
  const klen = sPowerLength(p);
  const { f, g } = endpointTaylor(c, p + 1);

  return { profile: sPowerProfile, ks: taylorToPairs(f, g, p, new Float64Array(klen)), klen, panels: 512 };
}

/** The same polynomial sampled at `KORDER` points, rescaled to the same total turning. */
function piecewiseLinearCase(phi: number): Case {
  const ks = new Float64Array(KORDER);

  for (let i = 0; i < KORDER; i++) {
    ks[i] = polyEval(KAPPA, i / (KORDER - 1));
  }

  const scale = phi / piecewiseLinear.integral(ks, KORDER, 1.0);

  for (let i = 0; i < KORDER; i++) {
    ks[i] *= scale;
  }

  // Panels are a multiple of the 11 knot intervals so no panel straddles a kink.
  return { profile: piecewiseLinear, ks, klen: KORDER, panels: 11 * 64 };
}

function endpointError(c: Case, steps: number, fourthOrder: boolean) {
  const ref = referenceEndpoint(c.profile, c.ks, c.klen, c.panels);
  const got = integrateProfile(c.profile, c.ks, c.klen, 0.0, 1.0, { steps, fourthOrder });

  return Math.hypot(got[0] - ref[0], got[1] - ref[1]);
}

const STEP_COUNTS = [2, 3, 4, 6, 8, 11, 16, 19, 32, 64, 128];

function errorCurve(c: Case, fourthOrder = false) {
  return STEP_COUNTS.map((n) => ({ n, err: endpointError(c, n, fourthOrder) }));
}

/** Observed order between successive step counts. */
function observedOrders(curve: { n: number; err: number }[]) {
  const out: { from: number; to: number; order: number }[] = [];

  for (let i = 1; i < curve.length; i++) {
    const a = curve[i - 1];
    const b = curve[i];

    out.push({ from: a.n, to: b.n, order: Math.log(a.err / b.err) / Math.log(b.n / a.n) });
  }

  return out;
}

function report(t: TestContext, label: string, curve: { n: number; err: number }[]) {
  const orders = observedOrders(curve);

  t.diagnostic(label);

  for (let i = 0; i < curve.length; i++) {
    const order = i === 0 ? "" : `  order ${orders[i - 1].order.toFixed(2)}`;

    t.diagnostic(`  N = ${String(curve[i].n).padStart(3)}   err ${curve[i].err.toExponential(2)}${order}`);
  }
}

const TURNINGS = [
  { name: "phi = pi/4", phi: Math.PI / 4 },
  { name: "phi = pi/2", phi: Math.PI / 2 },
  { name: "phi = pi", phi: Math.PI },
];

describe("Phase 0a: quadrature order", () => {
  it("reproduces the piecewise-linear behaviour clothoids.md §4 records", (t) => {
    const curve = errorCurve(piecewiseLinearCase(Math.PI / 2));

    report(t, "piecewise-linear, phi = pi/2", curve);

    const orders = observedOrders(curve).map((o) => o.order);
    const low = orders.slice(0, 5);

    assert.ok(
      Math.max(...low) - Math.min(...low) > 1.0,
      `expected the observed order to oscillate at low N, got ${low.map((o) => o.toFixed(2)).join(", ")}`
    );
  });

  for (const { name, phi } of TURNINGS) {
    it(`holds a clean third order on the s-power profile, ${name}`, (t) => {
      const curve = errorCurve(sPowerCase(phi));

      report(t, `s-power, ${name}`, curve);

      for (let i = 1; i < curve.length; i++) {
        assert.ok(
          curve[i].err < curve[i - 1].err,
          `error must decrease monotonically: N=${curve[i - 1].n} gave ${curve[i - 1].err}, N=${curve[i].n} gave ${curve[i].err}`
        );
      }

      // Read the order away from the very coarse end, where the asymptotic regime has not
      // started, and away from the fine end, where the reference's own error shows up.
      for (const o of observedOrders(curve).filter((o) => o.from >= 8 && o.to <= 64)) {
        assert.ok(o.order > 2.8 && o.order < 3.2, `order ${o.from}->${o.to} should be ~3, got ${o.order.toFixed(2)}`);
      }
    });
  }

  it("reaches fourth order with the ds^3 terms", (t) => {
    const curve = errorCurve(sPowerCase(Math.PI / 2), true);

    report(t, "s-power with fourth-order terms, phi = pi/2", curve);

    for (const o of observedOrders(curve).filter((o) => o.from >= 8 && o.to <= 32)) {
      assert.ok(o.order > 3.7 && o.order < 4.3, `order ${o.from}->${o.to} should be ~4, got ${o.order.toFixed(2)}`);
    }
  });

  /**
   * The plan's gate asks for the low-`N` constant to be *competitive*, not dominant, and
   * that distinction turns out to matter: at third order the piecewise-linear profile is
   * consistently a little better. That is not an accident. Its second derivative vanishes
   * inside every piece, so the `ds^2` term of the scheme is exact almost everywhere and
   * the constant it pays is only the kink contribution — which is also precisely why its
   * observed *order* refuses to settle. The s-power profile pays a genuine `kappa''`
   * everywhere and gets a clean order in exchange.
   *
   * Turning the fourth-order terms on is where the representation earns its keep: they
   * need `kappa''`, so they are available to s-power and not to a piecewise-linear profile
   * at all, and with them the same step count is several times more accurate.
   */
  it("is competitive at low N and wins outright once the ds^3 terms are on", (t) => {
    const rows: string[] = [];
    let worstRatio = 0.0;
    let worstGain = Infinity;
    let worstGainAt8 = Infinity;

    for (const { name, phi } of TURNINGS) {
      const sp = sPowerCase(phi);
      const pl = piecewiseLinearCase(phi);

      for (const n of [3, 4, 6, 8]) {
        const third = endpointError(sp, n, false);
        const fourth = endpointError(sp, n, true);
        const linear = endpointError(pl, n, false);

        worstRatio = Math.max(worstRatio, third / linear);
        worstGain = Math.min(worstGain, linear / fourth);

        if (n === 8) {
          worstGainAt8 = Math.min(worstGainAt8, linear / fourth);
        }

        rows.push(
          `  ${name.padEnd(11)} N = ${n}   s-power ${third.toExponential(2)}` +
            `   +ds^3 ${fourth.toExponential(2)}   pw-linear ${linear.toExponential(2)}`
        );
      }
    }

    t.diagnostic("low-N error at matched total turning");
    rows.forEach((r) => t.diagnostic(r));
    t.diagnostic(
      `  worst s-power/pw-linear ratio ${worstRatio.toFixed(2)}, ` +
        `worst ds^3 gain ${worstGain.toFixed(2)}x (${worstGainAt8.toFixed(2)}x at N = 8)`
    );

    assert.ok(
      worstRatio < 2.0,
      `third-order s-power error should stay within 2x of piecewise-linear, got ${worstRatio}`
    );

    // At N = 3 and half a turn per segment the fourth-order terms are only just ahead;
    // the separation opens up as soon as the step is small enough to be asymptotic.
    assert.ok(worstGain > 1.25, `the fourth-order scheme should never lose to piecewise-linear, got ${worstGain}x`);
    assert.ok(
      worstGainAt8 > 3.0,
      `the fourth-order scheme should be several times better by N = 8, got ${worstGainAt8}x`
    );
  });
});
