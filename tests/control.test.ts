/**
 * Phase 4 of `docs/plans/spower-solver.md` §12 — the control experiment.
 *
 * The new solver changed four things at once: the basis, where the unknowns live, the energy,
 * and the banded KKT step. §12 asks for the old piecewise-linear profile run through the *new*
 * solver, so that "the basis was the win" can be told apart from "the KKT and the energy were
 * the win." `src/curve/samples.ts` is that second {@link ProfileDOF}; this file is what says
 * it is a fair one.
 *
 * Three things have to hold before any measurement means anything:
 *
 * 1. **The sampling matrix samples.** `S·a` has to be the s-power polynomial evaluated at the
 *    nodes, or the control is running on a different curve rather than a different
 *    representation of the same one.
 * 2. **The turning row is exact.** Trapezoid weights integrate a polyline exactly, so the
 *    constraint the solver linearizes and the residual it measures have to agree to the last
 *    bit — the same requirement §4 puts on the s-power weights.
 * 3. **The Jacobian is still right.** §10's interface is generic, but the finite-difference
 *    check of Phase 3 is not something to take on faith across a substitution.
 *
 * Only then are the convergence and shape comparisons at the bottom worth reading.
 */
import assert from "node:assert/strict";
import { describe, it } from "node:test";

import { Vector2 } from "../src/math/index.js";
import {
  SAMPLE_COUNT,
  SPowerClothoid,
  SPowerSolver,
  chains,
  evalSPower,
  integrateProfile,
  piecewiseLinear,
  sPowerDOF,
  sPowerLength,
  sampleDOF,
  sampleEnergy,
  sampleMatrix,
  sampleTransform,
  sampledDOF,
  trapezoidWeights,
} from "../src/curve/index.js";
import { ChainSystem, defaultSPowerSolverOptions, referenceLengths } from "../src/curve/spower_solver.js";
import { Mesh, type Edge, type Vertex } from "../src/mesh/index.js";

const TAU = Math.PI * 2.0;
const P = 1;
const N = sPowerLength(P);

function wrap(a: number) {
  return a - TAU * Math.round(a / TAU);
}

const ZIGZAG = [
  [0, 0],
  [1, 0.35],
  [2, 0.2],
  [3, 0.9],
  [4, 0.7],
];

const ELBOW = ZIGZAG.slice(0, 3);

/** Nine points off a widening spiral — the strongest turning of the three. */
const SPIRAL = [0, 1, 2, 3, 4, 5, 6, 7, 8].map((i) => {
  const r = 1.0 + 0.45 * i;
  const th = 0.9 * i;

  return [r * Math.cos(th), r * Math.sin(th)];
});

const SHAPES: [string, number[][]][] = [
  ["elbow", ELBOW],
  ["zigzag", ZIGZAG],
  ["spiral", SPIRAL],
];

function polyline(points: number[][], flip: boolean[] = []) {
  const mesh = new Mesh();

  mesh.CurveCls = SPowerClothoid;
  mesh.SolverCls = SPowerSolver;

  const verts = points.map((p) => mesh.makeVertex(p));
  const edges: Edge[] = [];

  for (let i = 0; i + 1 < verts.length; i++) {
    edges.push(flip[i] ? mesh.makeEdge(verts[i + 1], verts[i]) : mesh.makeEdge(verts[i], verts[i + 1]));
  }

  return { mesh, verts, edges };
}

function travelAngle(e: Edge, v: Vertex, arriving: boolean) {
  const c = e.curve as SPowerClothoid;
  const atV2 = v === e.v2;

  const t = new Vector2(atV2 ? c.derivative(c.length) : c.derivative(0.0));

  if (atV2 !== arriving) {
    t.negate();
  }

  return Math.atan2(t[1], t[0]);
}

function worstGap(edges: Edge[], verts: Vertex[]) {
  let worst = 0.0;

  for (let i = 1; i < edges.length; i++) {
    const gap = wrap(travelAngle(edges[i], verts[i], false) - travelAngle(edges[i - 1], verts[i], true));

    worst = Math.max(worst, Math.abs(gap));
  }

  return worst;
}

describe("Phase 4: the sampling matrix", () => {
  it("evaluates the s-power basis at the nodes", () => {
    for (const p of [0, 1, 2]) {
      const n = sPowerLength(p);
      const s = sampleMatrix(p, SAMPLE_COUNT);
      const a = new Float64Array(n);

      for (let j = 0; j < n; j++) {
        a[j] = 1.0 + 0.7 * j - 0.3 * j * j;
      }

      for (let i = 0; i < SAMPLE_COUNT; i++) {
        let got = 0.0;

        for (let j = 0; j < n; j++) {
          got += s[i * n + j] * a[j];
        }

        const want = evalSPower(a, n, i / (SAMPLE_COUNT - 1));

        assert.ok(Math.abs(got - want) < 1e-12, `p=${p}, node ${i}: ${got} vs ${want}`);
      }
    }
  });

  /*
    The whole point of routing through `edgeTransform` rather than deriving a Hermite
    sampling directly: the two DOF interfaces differ by `S` and by nothing else.
  */
  it("is the s-power transform with the sampling in front", () => {
    const frame = { chord: 1.4, arclength: 1.55, rEarlier: 1.1, rLater: 0.8, forward: false };
    const dof = new Float64Array([0.4, -0.9, 1.3, 0.25]);

    const a = sPowerDOF.coefficients(P, frame, dof);
    const sampled = sampleDOF.coefficients(P, frame, dof);

    for (let i = 0; i < SAMPLE_COUNT; i++) {
      const want = evalSPower(a, N, i / (SAMPLE_COUNT - 1));

      assert.ok(Math.abs(sampled[i] - want) < 1e-12, `node ${i}: ${sampled[i]} vs ${want}`);
    }
  });

  it("has full column rank, so no degree of freedom is lost", () => {
    const frame = { chord: 1.0, arclength: 1.0, rEarlier: 1.0, rLater: 1.0, forward: true };
    const t = sampleTransform(P, frame, SAMPLE_COUNT);

    // Gram determinant by elimination on TᵀT; zero would mean two DOF are indistinguishable.
    const g: number[][] = [];

    for (let i = 0; i < N; i++) {
      g.push([]);

      for (let j = 0; j < N; j++) {
        let sum = 0.0;

        for (let k = 0; k < SAMPLE_COUNT; k++) {
          sum += t[k * N + i] * t[k * N + j];
        }

        g[i].push(sum);
      }
    }

    let det = 1.0;

    for (let i = 0; i < N; i++) {
      det *= g[i][i];

      for (let r = i + 1; r < N; r++) {
        const f = g[r][i] / g[i][i];

        for (let c = i; c < N; c++) {
          g[r][c] -= f * g[i][c];
        }
      }
    }

    assert.ok(det > 1e-6, `TᵀT is near singular, det ${det}`);
  });
});

describe("Phase 4: the sampled energy and turning", () => {
  it("integrates a polyline's turning exactly", () => {
    const w = trapezoidWeights(SAMPLE_COUNT);
    const ks = new Float64Array(SAMPLE_COUNT);

    for (let i = 0; i < SAMPLE_COUNT; i++) {
      ks[i] = Math.sin(2.3 * i) + 0.4 * i;
    }

    let dot = 0.0;

    for (let i = 0; i < SAMPLE_COUNT; i++) {
      dot += w[i] * ks[i];
    }

    const want = piecewiseLinear.integral(ks, SAMPLE_COUNT, 1.0);

    assert.ok(Math.abs(dot - want) < 1e-14, `${dot} vs ${want}`);
  });

  /*
    A ramp `q(u) = u` has `∫q′² = 1` and `∫q² = 1/3`, and P1 elements are exact on it. Both
    terms are checked because the mass matrix is what removes the constant-curvature null
    space, and a wrong one leaves `H` only semi-definite.
  */
  it("assembles the P1 stiffness and mass of a ramp", () => {
    const ks = new Float64Array(SAMPLE_COUNT);

    for (let i = 0; i < SAMPLE_COUNT; i++) {
      ks[i] = i / (SAMPLE_COUNT - 1);
    }

    const quad = (e: Float64Array) => {
      let sum = 0.0;

      for (let i = 0; i < SAMPLE_COUNT; i++) {
        for (let j = 0; j < SAMPLE_COUNT; j++) {
          sum += e[i * SAMPLE_COUNT + j] * ks[i] * ks[j];
        }
      }

      return sum;
    };

    const stiffness = quad(sampleEnergy(1.0, 0.0, SAMPLE_COUNT));
    const mass = quad(sampleEnergy(1.0, 1.0, SAMPLE_COUNT)) - stiffness;

    assert.ok(Math.abs(stiffness - 1.0) < 1e-12, `∫q′² = ${stiffness}`);
    assert.ok(Math.abs(mass - 1.0 / 3.0) < 1e-12, `∫q² = ${mass}`);
  });

  it("weights the energy by the chord cubed, as the s-power one does", () => {
    const a = sampleEnergy(1.0, 0.3, SAMPLE_COUNT);
    const b = sampleEnergy(2.0, 0.3, SAMPLE_COUNT);

    for (let i = 0; i < a.length; i++) {
      assert.ok(Math.abs(b[i] * 8.0 - a[i]) < 1e-12, `entry ${i} does not scale as 1/C³`);
    }
  });

  /*
    §1's representational claim, as an assertion rather than as prose. The fourth-order terms
    are not a tuning choice the control declined to make — they are unreachable, because a
    polyline's second derivative is a train of deltas and zero is the wrong answer for it.
  */
  it("cannot reach the fourth-order quadrature at all", () => {
    assert.equal(sampleDOF.quadrature.fourthOrder, false);
    assert.equal(sPowerDOF.quadrature.fourthOrder, true);

    assert.throws(
      () =>
        integrateProfile(piecewiseLinear, new Float64Array(SAMPLE_COUNT), SAMPLE_COUNT, 0.0, 1.0, {
          steps      : 19,
          fourthOrder: true,
        }),
      /d2Curvature/
    );
  });
});

function system(points: number[][], flip: boolean[] = [], seed = 4093) {
  const { mesh, edges } = polyline(points, flip);

  const chain = chains(mesh)[0];
  const out = new ChainSystem(chain, referenceLengths(mesh), { ...defaultSPowerSolverOptions, dof: sampleDOF });

  let state = seed >>> 0;

  for (let i = 0; i < chain.verts.length; i++) {
    for (let k = 0; k < out.block; k++) {
      state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
      out.z[out.blockAt[i] + k] = (state / 0x100000000) * 1.2 - 0.6;
    }
  }

  out.write();
  out.measure();

  return { mesh, edges, system: out };
}

describe("Phase 4: the control's constraint row", () => {
  /*
    Phase 3's check, rerun through the substituted `ProfileDOF`. The row is assembled by the
    same code, but it crosses a rectangular pullback and reads derivatives of a profile whose
    `dCurvature` is itself a finite difference, so "generic" is a claim to verify rather than
    assume.
  */
  for (const [name, flip] of [
    ["all forward", []],
    ["mixed", [false, true, true, false]],
  ] as [string, boolean[]][]) {
    it(`matches a differenced residual, ${name}`, () => {
      const { system: s } = system(ZIGZAG, flip);
      const h = 1e-6;

      s.assemble();

      let checked = 0;
      let peak = 0.0;

      for (let i = 1; i < s.frames.length; i++) {
        const at = s.lambdaAt[i];

        assert.ok(at >= 0, `joint ${i} has no row`);

        for (let c = 0; c < s.n; c++) {
          if (s.rows.includes(c)) {
            continue;
          }

          const keep = s.z[c];

          s.z[c] = keep + h;
          s.write();
          s.measure();
          const plus = s.residuals[i];

          s.z[c] = keep - h;
          s.write();
          s.measure();
          const minus = s.residuals[i];

          s.z[c] = keep;

          const fd = (plus - minus) / (2.0 * h);
          const row = s.kkt.get(at, c);

          peak = Math.max(peak, Math.abs(fd));

          assert.ok(Math.abs(row - fd) < 1e-7, `joint ${i}, dof ${c}: ${row} vs ${fd}`);
          checked++;
        }
      }

      s.write();
      s.measure();

      assert.ok(checked >= 30, `only ${checked} entries compared`);
      assert.ok(peak > 0.1, `the whole Jacobian is near zero, peak ${peak}`);
    });
  }
});

describe("Phase 4: the control experiment", () => {
  /*
    The headline. If the piecewise-linear profile converges in comparable step counts, then
    what §14 measured is the solver rather than the basis, and §1's claim has to rest on the
    representational half — the reachable quadrature order and the smoothness — instead.
  */
  for (const [name, points] of SHAPES) {
    it(`reaches tangent continuity on the ${name} with the sampled profile`, () => {
      const { mesh, verts, edges } = polyline(points);

      const report = new SPowerSolver(mesh, { dof: sampleDOF }).solve();

      assert.equal(report.ok, true);
      assert.equal(report.unenforced, 0);
      assert.ok(report.maxResidual < 1e-9, `${name}: solver residual ${report.maxResidual}`);
      assert.ok(worstGap(edges, verts) < 1e-8, `${name}: world G1 gap ${worstGap(edges, verts)}`);
    });
  }

  it("takes the same steps whichever way the edges were built", () => {
    const counts = new Set<number>();

    for (const flip of [[], [true, true, true, true], [false, true, true, false]]) {
      const { mesh } = polyline(ZIGZAG, flip);

      counts.add(new SPowerSolver(mesh, { dof: sampleDOF }).solve().steps);
    }

    assert.equal(counts.size, 1, `step counts differ by orientation: ${[...counts].join(", ")}`);
  });

  /*
    The bias §1 predicts, bounded rather than merely asserted. Sampling a degree-`2p+1`
    polynomial and rejoining the samples with straight lines loses `O(h²q″)`, so the two
    solves land on genuinely different curves — both G1, one of them not the profile it was
    written in. The gap has to be small enough that the control is a control, and nonzero, or
    the substitution did nothing.
  */
  it("lands near the s-power solve, but not on it", () => {
    const exact = polyline(ZIGZAG);
    const control = polyline(ZIGZAG);

    new SPowerSolver(exact.mesh, { dof: sPowerDOF }).solve();
    new SPowerSolver(control.mesh, { dof: sampleDOF }).solve();

    let worst = 0.0;

    for (let i = 0; i < exact.edges.length; i++) {
      const a = exact.edges[i].curve as SPowerClothoid;
      const b = control.edges[i].curve as SPowerClothoid;

      for (let k = 0; k <= 8; k++) {
        const t = k / 8;

        worst = Math.max(worst, new Vector2(a.evaluate(t * a.length)).vectorDistance(b.evaluate(t * b.length)));
      }
    }

    assert.ok(worst < 0.02, `the two solves are ${worst} apart, which is not a control`);
    assert.ok(worst > 1e-9, `the sampled profile reproduced the s-power one exactly, at ${worst}`);
  });

  /*
    That bias is interpolation error and nothing else, so it has to fall like `h²` as the
    sample count rises. If it did not, the difference would be coming from the energy or the
    constraint instead, and the experiment would not be isolating what it claims to.
  */
  it("converges on the s-power solve as the sampling refines", () => {
    const exact = polyline(ZIGZAG);

    new SPowerSolver(exact.mesh, { dof: sPowerDOF }).solve();

    const deviation = (samples: number) => {
      const control = polyline(ZIGZAG);

      new SPowerSolver(control.mesh, { dof: sampledDOF(samples) }).solve();

      let worst = 0.0;

      for (let i = 0; i < exact.edges.length; i++) {
        const a = exact.edges[i].curve as SPowerClothoid;
        const b = control.edges[i].curve as SPowerClothoid;

        for (let k = 0; k <= 8; k++) {
          const t = k / 8;

          worst = Math.max(worst, new Vector2(a.evaluate(t * a.length)).vectorDistance(b.evaluate(t * b.length)));
        }
      }

      return worst;
    };

    const coarse = deviation(6);
    const fine = deviation(24);

    assert.ok(fine < coarse * 0.1, `refining 4× cut the deviation only from ${coarse} to ${fine}`);
  });
});
