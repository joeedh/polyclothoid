/**
 * Phase 2 of `docs/plans/spower-solver.md` §12, layer 2 — the vertex-block transform.
 *
 * The headline claim of §3 is that `Gᵖ` continuity is *structural*: two edges that read the
 * same vertex block reconstruct the same curvature derivatives at that vertex, whatever
 * their lengths and whichever way they happen to be parameterized. The tests below check
 * that against ground truth rather than against the transform's own bookkeeping.
 *
 * The trick is to start from a single global curvature polynomial `κ(s)` in world
 * arclength, of degree at most `2p+1`, and lay two edges end to end along it. Since a
 * degree-`2p+1` polynomial is its own two-point Hermite interpolant, each edge must
 * reproduce its slice of `κ` *exactly*, and continuity at the joint follows from both
 * matching the same function rather than from any assertion about signs.
 *
 * Two things are deliberately adversarial. The second edge is laid down backwards, so the
 * `(−1)ⁿ⁺¹` orientation rule is on the critical path. And the two edges are given
 * *different* arclength-to-chord ratios — 1.46 against 1.09 — which is what discriminates
 * the arclength-based rescale from the chord-based one §3 specifies. With chords the joint
 * gaps here run to several percent instead of `1e-16`.
 */
import assert from "node:assert/strict";
import { describe, it } from "node:test";

import {
  congruence,
  edgeEnergy,
  edgeHessian,
  edgeTransform,
  edgeTurningRow,
  referenceLength,
  type EdgeFrame,
} from "../src/curve/blocks.js";
import { evalSPower, reverseCurvature, sPowerDerivative, sPowerIntegral, sPowerLength } from "../src/curve/spower.js";

/** Deterministic LCG — tests must not depend on `Math.random()`. */
function lcg(seed: number) {
  let state = seed >>> 0;

  return () => {
    state = (Math.imul(state, 1664525) + 1013904223) >>> 0;

    return state / 4294967296;
  };
}

const KAPPA = [0.7, -1.3, 2.1, 0.9, -1.7, 0.4, 1.1, -0.6];

/** Chord, arclength: two edges with visibly different `L/C`, which is the point. */
const C1 = 1.3;
const L1 = 1.9;
const C2 = 0.55;
const L2 = 0.6;

function polyEval(c: ArrayLike<number>, x: number) {
  let sum = 0.0;

  for (let i = c.length - 1; i >= 0; i--) {
    sum = sum * x + c[i];
  }

  return sum;
}

function polyDerivative(c: ArrayLike<number>) {
  const d = new Float64Array(Math.max(1, c.length - 1));

  for (let i = 1; i < c.length; i++) {
    d[i - 1] = i * c[i];
  }

  return d;
}

/** `dⁿκ/dsⁿ` at `s`, for `n = 0 … p`. */
function derivatives(c: ArrayLike<number>, s: number, p: number) {
  const out: number[] = [];
  let cur: ArrayLike<number> = c;

  for (let n = 0; n <= p; n++) {
    out.push(polyEval(cur, s));
    cur = polyDerivative(cur);
  }

  return out;
}

/** `block(v)_n = (Rᵥⁿ / n!) · dⁿκ/dsⁿ|_v`, straight from the definition in §3. */
function makeBlock(c: ArrayLike<number>, s: number, r: number, p: number) {
  const d = derivatives(c, s, p);
  const out = new Float64Array(p + 1);

  let fac = 1.0;

  for (let n = 0; n <= p; n++) {
    if (n > 0) {
      fac *= n;
    }

    out[n] = (Math.pow(r, n) / fac) * d[n];
  }

  return out;
}

function applyMatrix(m: ArrayLike<number>, x: ArrayLike<number>, n: number) {
  const out = new Float64Array(n);

  for (let i = 0; i < n; i++) {
    let sum = 0.0;

    for (let j = 0; j < n; j++) {
      sum += m[i * n + j] * x[j];
    }

    out[i] = sum;
  }

  return out;
}

function quadratic(m: ArrayLike<number>, x: ArrayLike<number>, n: number) {
  let sum = 0.0;

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      sum += x[i] * m[i * n + j] * x[j];
    }
  }

  return sum;
}

/** Composite Simpson on `[0, 1]`. The integrands here are polynomials, so this converges hard. */
function simpson(f: (u: number) => number, panels = 4096) {
  const h = 1.0 / panels;
  let sum = f(0.0) + f(1.0);

  for (let i = 1; i < panels; i++) {
    sum += (i & 1 ? 4.0 : 2.0) * f(i * h);
  }

  return (sum * h) / 3.0;
}

/**
 * The chain `A -> V -> B` laid along `KAPPA`, with edge 2 parameterized against the chain.
 *
 * Returns each edge's canonical coefficients. `q(u) = L·κ(L·u)` is what they should hold,
 * negated and flipped on the reversed edge.
 */
function layChain(p: number) {
  const n = sPowerLength(p);
  const kappa = KAPPA.slice(0, 2 * p + 2);

  const rA = referenceLength([C1]);
  const rV = referenceLength([C1, C2]);
  const rB = referenceLength([C2]);

  const frame1: EdgeFrame = { chord: C1, arclength: L1, rEarlier: rA, rLater: rV, forward: true };
  const frame2: EdgeFrame = { chord: C2, arclength: L2, rEarlier: rV, rLater: rB, forward: false };

  const dof1 = new Float64Array([...makeBlock(kappa, 0.0, rA, p), ...makeBlock(kappa, L1, rV, p)]);
  const dof2 = new Float64Array([...makeBlock(kappa, L1, rV, p), ...makeBlock(kappa, L1 + L2, rB, p)]);

  return {
    n,
    kappa,
    a1: applyMatrix(edgeTransform(p, frame1), dof1, n),
    a2: applyMatrix(edgeTransform(p, frame2), dof2, n),
  };
}

describe("vertex blocks: the DOF transform", () => {
  it("reproduces a degree-2p+1 curvature polynomial exactly on both edges", (t) => {
    for (let p = 0; p <= 3; p++) {
      const { n, kappa, a1, a2 } = layChain(p);

      let worst = 0.0;

      for (let i = 0; i <= 20; i++) {
        const u = i / 20;

        // Edge 1 runs with the chain; edge 2 runs against it, so it is negated and flipped.
        worst = Math.max(worst, Math.abs(evalSPower(a1, n, u) - L1 * polyEval(kappa, u * L1)));
        worst = Math.max(worst, Math.abs(evalSPower(a2, n, u) + L2 * polyEval(kappa, L1 + L2 * (1.0 - u))));
      }

      t.diagnostic(`  p=${p}  worst reconstruction error ${worst.toExponential(2)}`);

      assert.ok(worst < 1e-12, `p=${p}: reconstruction should be exact, worst ${worst}`);
    }
  });

  it("meets to order p at the shared vertex whatever the lengths and directions", (t) => {
    for (let p = 0; p <= 3; p++) {
      const { n, a1, a2 } = layChain(p);

      // World curvature at the joint, read from each edge and put back in chain orientation.
      // Both edges *end* there — edge 2 because it runs backwards — so both are read at u=1,
      // and edge 2 carries the reversal sign. `q = L·κ`, so undoing it divides by L.
      const k1 = evalSPower(a1, n, 1.0) / L1;
      const k2 = -evalSPower(a2, n, 1.0) / L2;

      // For dκ/ds the parameter reverses too, and the two flips cancel: du = ds/L both times.
      const d1 = sPowerDerivative(a1, n, 1.0) / (L1 * L1);
      const d2 = sPowerDerivative(a2, n, 1.0) / (L2 * L2);

      t.diagnostic(
        `  p=${p}  G0 gap ${Math.abs(k1 - k2).toExponential(2)}  G1 gap ${Math.abs(d1 - d2).toExponential(2)}`
      );

      assert.ok(Math.abs(k1 - k2) < 1e-12, `p=${p}: curvature should be continuous, gap ${k1 - k2}`);

      if (p >= 1) {
        assert.ok(Math.abs(d1 - d2) < 1e-12, `p=${p}: dκ/ds should be continuous, gap ${d1 - d2}`);
      }
    }
  });

  /**
   * The guard on the §3 deviation, and the reason the test above is worth anything.
   *
   * §3 specifies "the rescale from `Rᵥ` to `C_e`". Feeding the chord in where the arclength
   * belongs is what that would amount to, and it puts a visible gap in the joint — so the
   * clean result above is a property of the arclength rescale, not of the setup. If someone
   * restores the letter of §3 this fails rather than quietly degrading continuity.
   *
   * The gap only appears when the coefficients are read back through the *realized*
   * geometry, which is the subtlety that makes the wrong version look fine. Substituting a
   * constant `S_e` for `L_e` and then reading the derivatives back with the same `S_e`
   * cancels exactly — it is just a change of normalization. What does not cancel is the
   * curve: coefficients drive `integrateProfile`, whose output is scaled to arclength
   * `L_e`, so the derivative the geometry actually carries at the vertex is the intended
   * one times `(S_e/L_e)ⁿ⁺¹`. Continuity then needs `C₁/L₁ = C₂/L₂`, which no polyline
   * owes anybody.
   */
  it("loses continuity at the joint if the rescale uses the chord", (t) => {
    for (let p = 1; p <= 3; p++) {
      const n = sPowerLength(p);
      const kappa = KAPPA.slice(0, 2 * p + 2);

      const rV = referenceLength([C1, C2]);
      const rA = referenceLength([C1]);
      const rB = referenceLength([C2]);

      // The only change from `layChain`: arclength := chord.
      const f1: EdgeFrame = { chord: C1, arclength: C1, rEarlier: rA, rLater: rV, forward: true };
      const f2: EdgeFrame = { chord: C2, arclength: C2, rEarlier: rV, rLater: rB, forward: false };

      const a1 = applyMatrix(
        edgeTransform(p, f1),
        new Float64Array([...makeBlock(kappa, 0.0, rA, p), ...makeBlock(kappa, L1, rV, p)]),
        n
      );
      const a2 = applyMatrix(
        edgeTransform(p, f2),
        new Float64Array([...makeBlock(kappa, L1, rV, p), ...makeBlock(kappa, L1 + L2, rB, p)]),
        n
      );

      // Read back through the realized arclengths, which is what the geometry carries.
      const k1 = evalSPower(a1, n, 1.0) / L1;
      const k2 = -evalSPower(a2, n, 1.0) / L2;

      const relative = Math.abs(k1 - k2) / Math.max(Math.abs(k1), Math.abs(k2));

      t.diagnostic(`  p=${p}  chord-based rescale leaves a ${(relative * 100).toFixed(1)}% curvature jump`);

      assert.ok(relative > 0.05, `p=${p}: expected the chord rescale to break continuity, got ${relative}`);
    }
  });

  it("reverses to exactly the negated, flipped profile", () => {
    const rand = lcg(0x5eed);

    for (let p = 0; p <= 3; p++) {
      const n = sPowerLength(p);
      const base = { chord: 0.8, arclength: 0.83, rEarlier: 1.1, rLater: 0.6 };

      const dof = new Float64Array(n);

      for (let i = 0; i < n; i++) {
        dof[i] = rand() * 2.0 - 1.0;
      }

      const fwd = applyMatrix(edgeTransform(p, { ...base, forward: true }), dof, n);
      const rev = applyMatrix(edgeTransform(p, { ...base, forward: false }), dof, n);

      const expect = reverseCurvature(fwd, n, new Float64Array(n));

      for (let i = 0; i < n; i++) {
        assert.ok(
          Math.abs(rev[i] - expect[i]) < 1e-12,
          `p=${p}, coefficient ${i}: ${rev[i]} should be the reversal ${expect[i]}`
        );
      }
    }
  });

  it("scales the DOF but not the coefficients under a uniform zoom", () => {
    const rand = lcg(0xc0ffee);
    const lambda = 7.25;

    for (let p = 0; p <= 3; p++) {
      const n = sPowerLength(p);
      const frame: EdgeFrame = { chord: 0.9, arclength: 1.05, rEarlier: 1.4, rLater: 0.7, forward: true };

      const dof = new Float64Array(n);
      const zoomed = new Float64Array(n);

      for (let i = 0; i < n; i++) {
        dof[i] = rand() * 2.0 - 1.0;
        zoomed[i] = dof[i] / lambda;
      }

      const a = applyMatrix(edgeTransform(p, frame), dof, n);
      const b = applyMatrix(
        edgeTransform(p, {
          chord    : frame.chord * lambda,
          arclength: frame.arclength * lambda,
          rEarlier : frame.rEarlier * lambda,
          rLater   : frame.rLater * lambda,
          forward  : true,
        }),
        zoomed,
        n
      );

      for (let i = 0; i < n; i++) {
        assert.ok(Math.abs(a[i] - b[i]) < 1e-12, `p=${p}, coefficient ${i}: ${a[i]} vs ${b[i]} under zoom`);
      }
    }
  });

  it("uses the geometric mean for the reference length", () => {
    assert.ok(Math.abs(referenceLength([4.0]) - 4.0) < 1e-15);
    assert.ok(Math.abs(referenceLength([2.0, 8.0]) - 4.0) < 1e-14);
    assert.ok(Math.abs(referenceLength([1.0, 2.0, 4.0]) - 2.0) < 1e-14);

    // A degenerate incident edge must not take the mean to zero or NaN.
    assert.ok(Math.abs(referenceLength([0.0, 3.0]) - 3.0) < 1e-14);
    assert.equal(referenceLength([]), 1.0);
  });
});

describe("vertex blocks: energy and turning", () => {
  it("matches the quadrature of (1/C³)(∫q′² + α²∫q²)", (t) => {
    const rand = lcg(0xbeef);

    for (let p = 1; p <= 3; p++) {
      const n = sPowerLength(p);
      const frame: EdgeFrame = { chord: 1.7, arclength: 1.85, rEarlier: 1.2, rLater: 0.9, forward: true };
      const alpha = 0.75;

      const dof = new Float64Array(n);

      for (let i = 0; i < n; i++) {
        dof[i] = rand() * 2.0 - 1.0;
      }

      const a = applyMatrix(edgeTransform(p, frame), dof, n);

      const bend = simpson((u) => sPowerDerivative(a, n, u) ** 2);
      const foundation = simpson((u) => evalSPower(a, n, u) ** 2);
      const expect = (bend + alpha * alpha * foundation) / frame.chord ** 3;

      const viaCoefficients = quadratic(edgeEnergy(p, frame.chord, alpha), a, n);
      const viaDOF = quadratic(edgeHessian(p, frame, alpha), dof, n);

      t.diagnostic(`  p=${p}  quadrature ${expect.toExponential(6)}  Gram ${viaCoefficients.toExponential(6)}`);

      assert.ok(
        Math.abs(viaCoefficients - expect) < 1e-9 * Math.abs(expect),
        `p=${p}: coefficient-space energy ${viaCoefficients} vs quadrature ${expect}`
      );
      assert.ok(
        Math.abs(viaDOF - expect) < 1e-9 * Math.abs(expect),
        `p=${p}: DOF-space energy ${viaDOF} vs quadrature ${expect}`
      );
    }
  });

  it("gives a symmetric Hessian that is positive semi-definite with a constant null space", () => {
    const rand = lcg(0x1234);

    for (let p = 0; p <= 3; p++) {
      const n = sPowerLength(p);
      const frame: EdgeFrame = { chord: 1.1, arclength: 1.3, rEarlier: 0.8, rLater: 1.6, forward: false };
      const h = edgeHessian(p, frame, 0.0);

      for (let i = 0; i < n; i++) {
        for (let j = 0; j < i; j++) {
          const d = Math.abs(h[i * n + j] - h[j * n + i]);

          assert.ok(d < 1e-12 * (1.0 + Math.abs(h[i * n + j])), `p=${p}: H is not symmetric at (${i}, ${j})`);
        }
      }

      for (let trial = 0; trial < 8; trial++) {
        const x = new Float64Array(n);

        for (let i = 0; i < n; i++) {
          x[i] = rand() * 2.0 - 1.0;
        }

        assert.ok(quadratic(h, x, n) > -1e-12, `p=${p}: H should be positive semi-definite`);
      }

      // With α = 0 the energy is pure bending, so a constant profile costs nothing. Constant
      // q means f₀ = g₀ and no higher Taylor data, i.e. both blocks are (v/L, 0, …).
      const constant = new Float64Array(n);
      constant[0] = 1.0 / frame.arclength;
      constant[p + 1] = 1.0 / frame.arclength;

      assert.ok(
        Math.abs(quadratic(h, constant, n)) < 1e-12,
        `p=${p}: constant curvature should have zero bending energy, got ${quadratic(h, constant, n)}`
      );
    }
  });

  it("reads total turning off the DOF directly", () => {
    const rand = lcg(0xfeed);

    for (let p = 0; p <= 3; p++) {
      const n = sPowerLength(p);

      for (const forward of [true, false]) {
        const frame: EdgeFrame = { chord: 0.62, arclength: 0.71, rEarlier: 1.05, rLater: 0.44, forward };

        const dof = new Float64Array(n);

        for (let i = 0; i < n; i++) {
          dof[i] = rand() * 2.0 - 1.0;
        }

        const a = applyMatrix(edgeTransform(p, frame), dof, n);
        const row = edgeTurningRow(p, frame);

        let dot = 0.0;

        for (let i = 0; i < n; i++) {
          dot += row[i] * dof[i];
        }

        const expect = sPowerIntegral(a, n, 1.0);

        assert.ok(
          Math.abs(dot - expect) < 1e-12 * (1.0 + Math.abs(expect)),
          `p=${p}, forward=${forward}: turning row gives ${dot}, integral gives ${expect}`
        );
      }
    }
  });

  it("computes BᵀAB", () => {
    // A = [[1,2],[3,4]], B = [[5,6],[7,8]]: AB = [[19,22],[43,50]], BᵀAB = [[396,460],[458,532]].
    const r = congruence([1, 2, 3, 4], [5, 6, 7, 8], 2);

    assert.deepEqual(Array.from(r), [396, 460, 458, 532]);
  });
});
