/**
 * Phase 0b acceptance tests for the s-power kernel.
 *
 * The two properties `docs/research/spower.md` §5 checks against `s-power.reduce` —
 * reproduction and Hermite contact — plus the calculus rules, the reversal convention,
 * and the Gram matrices against independent numerical quadrature.
 */
import assert from "node:assert/strict";
import { describe, it } from "node:test";

import {
  differentiateSPower,
  differentiationMatrix,
  evalSPower,
  hik,
  integralWeights,
  integrateSPower,
  massMatrix,
  pairsToTaylor,
  reverseCurvature,
  reverseSPower,
  sPowerLength,
  stiffnessMatrix,
  taylorToPairs,
} from "../src/curve/spower.js";

/** Monomial-basis polynomial, `c[i]` weighting `u^i`. */
function polyEval(c: number[], u: number) {
  let acc = 0.0;

  for (let i = c.length - 1; i >= 0; i--) {
    acc = acc * u + c[i];
  }

  return acc;
}

function polyDeriv(c: number[]) {
  return c.slice(1).map((v, i) => v * (i + 1));
}

/**
 * Endpoint Taylor data in the convention {@link taylorToPairs} wants:
 * `f_i = f^(i)(0)/i!`, `g_i = (-1)^i f^(i)(1)/i!`.
 */
function endpointTaylor(c: number[], count: number) {
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

/** Sample points that avoid the endpoints and the symmetry axis. */
const SAMPLES = [0.0, 0.03, 0.17, 0.31, 0.5, 0.62, 0.79, 0.94, 1.0];

function closeTo(a: number, b: number, eps: number, msg: string) {
  assert.ok(Math.abs(a - b) <= eps, `${msg}: ${a} vs ${b} (delta ${Math.abs(a - b)})`);
}

describe("s-power basis", () => {
  it("evaluates the basis functions as u^k (1-u)^{k+1} and u^{k+1} (1-u)^k", () => {
    const p = 3;
    const n = sPowerLength(p);
    const a = new Float64Array(n);

    for (let j = 0; j < n; j++) {
      a.fill(0.0);
      a[j] = 1.0;

      const k = j >> 1;
      const e = j & 1;

      for (const u of SAMPLES) {
        const want = Math.pow(u, k + e) * Math.pow(1.0 - u, k + 1 - e);

        closeTo(evalSPower(a, n, u), want, 1e-14, `basis ${j} at u=${u}`);
      }
    }
  });

  it("reproduces a degree-7 polynomial at order 3", () => {
    const c = [0.7, -1.3, 2.1, 0.4, -3.2, 1.9, 0.6, -2.5];
    const p = 3;
    const n = sPowerLength(p);

    const { f, g } = endpointTaylor(c, p + 1);
    const a = taylorToPairs(f, g, p, new Float64Array(n));

    for (const u of SAMPLES) {
      closeTo(evalSPower(a, n, u), polyEval(c, u), 1e-12, `reproduction at u=${u}`);
    }
  });

  it("round-trips pairs through endpoint Taylor data", () => {
    const p = 3;
    const n = sPowerLength(p);
    const a = new Float64Array([0.3, -1.1, 2.0, 0.5, -0.8, 1.4, 0.25, -0.6]);

    const f = new Float64Array(p + 1);
    const g = new Float64Array(p + 1);

    pairsToTaylor(a, p, f, g);

    const back = taylorToPairs(f, g, p, new Float64Array(n));

    for (let i = 0; i < n; i++) {
      closeTo(back[i], a[i], 1e-12, `pair ${i}`);
    }
  });

  it("matches the h(i,k) values spower.md quotes", () => {
    assert.deepEqual(hik(0, 0), [1.0, 0.0]);
    assert.deepEqual(hik(0, 1), [1.0, -1.0]);
    assert.deepEqual(hik(1, 1), [1.0, 0.0]);
    assert.deepEqual(hik(0, 2), [binom(3, 2), -binom(3, 1)]);
    assert.deepEqual(hik(1, 2), [binom(2, 1), -binom(2, 0)]);
  });
});

function binom(n: number, i: number) {
  let r = 1.0;

  for (let j = 0; j < i; j++) {
    r = (r * (n - j)) / (j + 1);
  }

  return r;
}

describe("s-power truncation is two-point Hermite interpolation", () => {
  const c = [1.2, -0.7, 3.1, 0.9, -2.2, 1.5, 0.3, -1.8, 2.4, -0.5, 1.1];
  const full = 5;

  /** `d^n a/du^n` at `u`, by applying the exact differentiation rule `n` times. */
  function nthDerivative(a: Float64Array, len: number, n: number, u: number) {
    const work = Float64Array.from(a.subarray(0, len));

    for (let i = 0; i < n; i++) {
      differentiateSPower(work, len, work);
    }

    return evalSPower(work, len, u);
  }

  it("matches derivatives 0..k at both ends and diverges at k+1", () => {
    const { f, g } = endpointTaylor(c, full + 1);
    const nFull = sPowerLength(full);
    const exact = taylorToPairs(f, g, full, new Float64Array(nFull));

    for (let k = 0; k < full; k++) {
      const nk = sPowerLength(k);
      const trunc = new Float64Array(nFull);
      trunc.set(exact.subarray(0, nk));

      for (let d = 0; d <= k; d++) {
        for (const u of [0.0, 1.0]) {
          const want = nthDerivative(exact, nFull, d, u);
          const got = nthDerivative(trunc, nFull, d, u);

          closeTo(got, want, 1e-9, `truncation ${k}, derivative ${d} at u=${u}`);
        }
      }

      const wantNext = nthDerivative(exact, nFull, k + 1, 0.0);
      const gotNext = nthDerivative(trunc, nFull, k + 1, 0.0);

      assert.ok(
        Math.abs(gotNext - wantNext) > 1e-6,
        `truncation ${k} should not match derivative ${k + 1}, got ${gotNext} vs ${wantNext}`
      );
    }
  });
});

describe("s-power calculus", () => {
  const p = 3;
  const n = sPowerLength(p);
  const a = new Float64Array([0.4, -1.7, 2.3, 0.8, -1.1, 0.9, 1.6, -0.3]);

  it("differentiates against a central finite difference", () => {
    const d = differentiateSPower(a, n, new Float64Array(n));
    const h = 1e-5;

    for (const u of SAMPLES.slice(1, -1)) {
      const fd = (evalSPower(a, n, u + h) - evalSPower(a, n, u - h)) / (2 * h);

      closeTo(evalSPower(d, n, u), fd, 1e-6, `derivative at u=${u}`);
    }
  });

  it("leaves the top pair symmetric, i.e. drops one degree", () => {
    const d = differentiateSPower(a, n, new Float64Array(n));

    closeTo(d[n - 1], d[n - 2], 1e-14, "top pair of the derivative");
  });

  it("integrates as the exact inverse of differentiation", () => {
    const raised = new Float64Array(n + 2);
    raised.set(a);

    const d = differentiateSPower(raised, n + 2, new Float64Array(n + 2));
    const back = integrateSPower(d, n + 2, new Float64Array(n + 4), evalSPower(a, n, 0.0));

    for (const u of SAMPLES) {
      closeTo(evalSPower(back, n + 4, u), evalSPower(a, n, u), 1e-12, `integrate(differentiate) at u=${u}`);
    }
  });

  it("integrates to the antiderivative pinned at u=0", () => {
    const anti = integrateSPower(a, n, new Float64Array(n + 2), 0.0);

    closeTo(evalSPower(anti, n + 2, 0.0), 0.0, 1e-14, "antiderivative at 0");

    for (const u of SAMPLES.slice(1)) {
      let ref = 0.0;
      const steps = 40000;

      for (let i = 0; i < steps; i++) {
        ref += evalSPower(a, n, (u * (i + 0.5)) / steps);
      }

      ref *= u / steps;

      closeTo(evalSPower(anti, n + 2, u), ref, 1e-8, `antiderivative at u=${u}`);
    }
  });

  it("gives an antiderivative whose top pair is skew-free", () => {
    const anti = integrateSPower(a, n, new Float64Array(n + 2), 0.0);

    closeTo(anti[n + 1], anti[n], 1e-14, "top pair of the antiderivative");
  });
});

describe("s-power reversal", () => {
  const p = 2;
  const n = sPowerLength(p);
  const a = new Float64Array([0.4, -1.7, 2.3, 0.8, -1.1, 0.9]);

  it("swaps pairs, giving a(1 - u)", () => {
    const r = reverseSPower(a, n, new Float64Array(n));

    for (const u of SAMPLES) {
      closeTo(evalSPower(r, n, u), evalSPower(a, n, 1.0 - u), 1e-14, `reversal at u=${u}`);
    }
  });

  it("round-trips", () => {
    const r = reverseSPower(reverseSPower(a, n, new Float64Array(n)), n, new Float64Array(n));

    for (let i = 0; i < n; i++) {
      closeTo(r[i], a[i], 0.0, `pair ${i}`);
    }
  });

  it("negates as well for signed curvature, and round-trips", () => {
    const r = reverseCurvature(a, n, new Float64Array(n));

    for (const u of SAMPLES) {
      closeTo(evalSPower(r, n, u), -evalSPower(a, n, 1.0 - u), 1e-14, `curvature reversal at u=${u}`);
    }

    const back = reverseCurvature(r, n, new Float64Array(n));

    for (let i = 0; i < n; i++) {
      closeTo(back[i], a[i], 1e-15, `pair ${i}`);
    }
  });

  it("reverses the derivative with a sign, and leaves the second derivative alone", () => {
    const rev = reverseCurvature(a, n, new Float64Array(n));

    const d = differentiateSPower(a, n, new Float64Array(n));
    const dRev = differentiateSPower(rev, n, new Float64Array(n));

    // d/d(1-u) picks up a minus from the chain rule, cancelling the curvature flip.
    for (const u of SAMPLES) {
      closeTo(evalSPower(dRev, n, u), evalSPower(d, n, 1.0 - u), 1e-13, `derivative reversal at u=${u}`);
    }
  });
});

describe("s-power Gram matrices", () => {
  /** Gauss-Legendre would be tidier; a fine midpoint rule is enough to catch a wrong entry. */
  function numericInner(
    p: number,
    j: number,
    l: number,
    transform: (a: Float64Array, len: number) => Float64Array | ArrayLike<number>
  ) {
    const n = sPowerLength(p);

    const ej = new Float64Array(n);
    const el = new Float64Array(n);
    ej[j] = 1.0;
    el[l] = 1.0;

    const tj = transform(ej, n);
    const tl = transform(el, n);

    const steps = 200000;
    let sum = 0.0;

    for (let i = 0; i < steps; i++) {
      const u = (i + 0.5) / steps;

      sum += evalSPower(tj, n, u) * evalSPower(tl, n, u);
    }

    return sum / steps;
  }

  const identity = (a: Float64Array, _len: number) => a;
  const derivative = (a: Float64Array, len: number) => differentiateSPower(a, len, new Float64Array(len));

  for (const p of [0, 1, 2]) {
    it(`mass matrix at order ${p} matches quadrature and is symmetric`, () => {
      const n = sPowerLength(p);
      const m = massMatrix(p);

      for (let j = 0; j < n; j++) {
        for (let l = 0; l < n; l++) {
          closeTo(m[j * n + l], m[l * n + j], 1e-15, `mass symmetry ${j},${l}`);
          closeTo(m[j * n + l], numericInner(p, j, l, identity), 1e-9, `mass ${j},${l}`);
        }
      }
    });

    it(`stiffness matrix at order ${p} matches quadrature and is symmetric`, () => {
      const n = sPowerLength(p);
      const k = stiffnessMatrix(p);

      for (let j = 0; j < n; j++) {
        for (let l = 0; l < n; l++) {
          closeTo(k[j * n + l], k[l * n + j], 1e-12, `stiffness symmetry ${j},${l}`);
          closeTo(k[j * n + l], numericInner(p, j, l, derivative), 1e-8, `stiffness ${j},${l}`);
        }
      }
    });
  }

  it("has constant curvature in the null space of the stiffness matrix", () => {
    for (const p of [0, 1, 2]) {
      const n = sPowerLength(p);
      const k = stiffnessMatrix(p);

      // The constant function 1 is (1, 1) in pair 0 and zero above.
      const c = new Float64Array(n);
      c[0] = 1.0;
      c[1] = 1.0;

      for (let i = 0; i < n; i++) {
        let sum = 0.0;

        for (let j = 0; j < n; j++) {
          sum += k[i * n + j] * c[j];
        }

        closeTo(sum, 0.0, 1e-13, `stiffness null space, order ${p}, row ${i}`);
      }
    }
  });

  it("agrees with the differentiation matrix", () => {
    const p = 2;
    const n = sPowerLength(p);
    const d = differentiationMatrix(p);
    const a = new Float64Array([0.4, -1.7, 2.3, 0.8, -1.1, 0.9]);

    const want = differentiateSPower(a, n, new Float64Array(n));

    for (let i = 0; i < n; i++) {
      let sum = 0.0;

      for (let j = 0; j < n; j++) {
        sum += d[i * n + j] * a[j];
      }

      closeTo(sum, want[i], 1e-14, `differentiation matrix row ${i}`);
    }
  });

  it("gives integral weights that reproduce the antiderivative at u=1", () => {
    for (const p of [0, 1, 2, 3]) {
      const n = sPowerLength(p);
      const w = integralWeights(p);
      const a = new Float64Array(n);

      for (let i = 0; i < n; i++) {
        a[i] = Math.sin(1.7 * i + 0.3);
      }

      const anti = integrateSPower(a, n, new Float64Array(n + 2), 0.0);

      let dot = 0.0;
      for (let i = 0; i < n; i++) {
        dot += w[i] * a[i];
      }

      closeTo(dot, evalSPower(anti, n + 2, 1.0), 1e-13, `integral weights at order ${p}`);
    }
  });
});
