/**
 * Phase 7 of `docs/plans/spower-solver.md` §12 — degree continuation, `p = 0 → 2`.
 *
 * §9 makes a precise claim and a vague one, and they need different tests.
 *
 * The precise one is about *one edge*: the order-`p` profile is a polynomial of degree
 * `2p + 1`, so the order-`p+1` blocks that reproduce it exist and are computable. That has an
 * exact answer, and the first suite pins it to roundoff — a seed that is merely close would
 * mean a sign or a scale is wrong in a way the solver would quietly absorb.
 *
 * The vague one is about *a chain*: the two edges at a joint disagree about the new
 * derivative, the ceiling makes them share one, so the rung below is a starting point and not
 * a solution. §12 says the phase is worth keeping only if the step count improves, so the
 * second suite measures the step count. It also fixes the two things that must hold whatever
 * the counts say: continuation may not change the answer, and it may not build a warm start
 * out of a rung §8 says cannot be trusted.
 */
import assert from "node:assert/strict";
import { describe, it } from "node:test";

import {
  type EdgeFrame,
  SPowerClothoid,
  SPowerSolver,
  continuationEntry,
  edgeCoefficients,
  endpointTaylor,
  evalSPower,
  pairsToTaylor,
  sPowerLength,
} from "../src/curve/index.js";
import { ChainSystem, defaultSPowerSolverOptions, referenceLengths } from "../src/curve/spower_solver.js";
import { chains } from "../src/curve/topology.js";
import { Mesh, type Edge } from "../src/mesh/index.js";

const ZIGZAG = [
  [0, 0],
  [1, 0.4],
  [2, -0.4],
  [3, 0.4],
  [4, 0],
];

const ARC = [
  [0, 0],
  [1, 0.5],
  [2, 1.4],
  [2.8, 2.6],
  [3.2, 4.0],
];

/** A hairpin pair, which §14's Phase 6 table has degrading rather than converging cleanly. */
const TIGHT = [
  [0, 0],
  [1, 1.2],
  [2, 0],
  [3, 1.2],
  [4, 0],
  [5, 1.2],
];

function polyline(points: number[][]) {
  const mesh = new Mesh();

  mesh.CurveCls = SPowerClothoid;
  mesh.SolverCls = SPowerSolver;

  const verts = points.map((p) => mesh.makeVertex(p));
  const edges: Edge[] = [];

  for (let i = 0; i + 1 < verts.length; i++) {
    edges.push(mesh.makeEdge(verts[i], verts[i + 1]));
  }

  return { mesh, verts, edges };
}

function solve(points: number[][], options = {}) {
  const { mesh, edges } = polyline(points);
  const report = new SPowerSolver(mesh, options).solve();

  return { report, edges };
}

/** Twenty points along each curve of a solved chain, as one flat list. */
function trace(edges: Edge[]) {
  const out: number[] = [];

  for (const e of edges) {
    const curve = e.curve as SPowerClothoid;

    for (let i = 0; i <= 20; i++) {
      const p = curve.evaluate((i / 20) * curve.length);

      out.push(p[0], p[1]);
    }
  }

  return out;
}

function farthest(a: number[], b: number[]) {
  let worst = 0.0;

  for (let i = 0; i < a.length; i += 2) {
    worst = Math.max(worst, Math.hypot(a[i] - b[i], a[i + 1] - b[i + 1]));
  }

  return worst;
}

/** DOF that are nothing in particular, so no accidental symmetry carries the test. */
function arbitraryDOF(p: number) {
  const dof = new Float64Array(2 * (p + 1));

  for (let i = 0; i < dof.length; i++) {
    dof[i] = Math.cos(i * 2.3 + p) * (1.0 + i);
  }

  return dof;
}

/** The order-`p+1` DOF of the same edge: the old entries, plus one new one per end. */
function raise(p: number, frame: EdgeFrame, dof: Float64Array) {
  const next = new Float64Array(2 * (p + 2));

  for (let end = 0; end < 2; end++) {
    for (let n = 0; n <= p; n++) {
      next[end * (p + 2) + n] = dof[end * (p + 1) + n];
    }

    next[end * (p + 2) + p + 1] = continuationEntry(p, frame, dof, end === 0);
  }

  return next;
}

const FRAMES: [string, EdgeFrame][] = [
  ["forward", { chord: 1.3, arclength: 1.45, rEarlier: 1.1, rLater: 0.9, forward: true }],
  ["reversed", { chord: 1.3, arclength: 1.45, rEarlier: 1.1, rLater: 0.9, forward: false }],
  ["unit", { chord: 1.0, arclength: 1.0, rEarlier: 1.0, rLater: 1.0, forward: true }],
];

describe("Phase 7: the seed, on one edge", () => {
  it("reads the endpoint derivatives the basis already carries", () => {
    const p = 2;
    const a = edgeCoefficients(p, FRAMES[0][1], arbitraryDOF(p));

    const f = new Float64Array(p + 1);
    const g = new Float64Array(p + 1);

    pairsToTaylor(a, p, f, g);

    for (let n = 0; n <= p; n++) {
      const [ef, eg] = endpointTaylor(a, sPowerLength(p), n);

      assert.ok(Math.abs(ef - f[n]) < 1e-12, `f[${n}]: ${ef} vs ${f[n]}`);
      assert.ok(Math.abs(eg - g[n]) < 1e-12, `g[${n}]: ${eg} vs ${g[n]}`);
    }
  });

  it("returns zero past the degree, rather than refusing", () => {
    const p = 1;
    const a = edgeCoefficients(p, FRAMES[0][1], arbitraryDOF(p));

    // Degree is 2p+1 = 3, so the fourth derivative is genuinely zero and not an error case.
    for (const value of endpointTaylor(a, sPowerLength(p), 2 * p + 2)) {
      assert.ok(Math.abs(value) < 1e-10, `expected 0, got ${value}`);
    }
  });

  for (const [name, frame] of FRAMES) {
    it(`reproduces the order-p curve exactly at order p+1 (${name})`, () => {
      for (const p of [0, 1, 2]) {
        const dof = arbitraryDOF(p);

        const a = edgeCoefficients(p, frame, dof);
        const b = edgeCoefficients(p + 1, frame, raise(p, frame, dof));

        let worst = 0.0;

        for (let i = 0; i <= 20; i++) {
          const u = i / 20;

          worst = Math.max(worst, Math.abs(evalSPower(a, sPowerLength(p), u) - evalSPower(b, sPowerLength(p + 1), u)));
        }

        // Exact, not close: a degree-2p+1 polynomial *is* its own Hermite interpolant at
        // order p+1, so anything above roundoff is a wrong sign or a wrong scale.
        assert.ok(worst < 1e-12, `p=${p}: worst |Δq| = ${worst}`);
      }
    });
  }
});

describe("Phase 7: the ladder", () => {
  it("solves at every rung of §9's range", () => {
    for (const order of [0, 1, 2]) {
      const { report } = solve(ZIGZAG, { order });

      assert.ok(report.ok, `p=${order} did not solve`);
      assert.ok(report.maxResidual < 1e-9, `p=${order} residual ${report.maxResidual}`);
    }
  });

  it("costs nothing and reports nothing when it is off", () => {
    const { report } = solve(ZIGZAG, { order: 2 });

    assert.equal(report.seedSteps, 0);
  });

  it("charges the rungs below to seedSteps and not to steps", () => {
    const { report } = solve(ZIGZAG, { order: 2, continuation: true });

    assert.ok(report.seedSteps > 0, "the ladder ran no lower rungs");
    assert.equal(report.attempts, 1, "the rungs are one attempt, not one each");
  });

  it("reaches the same curve as a cold start", () => {
    const cold = solve(ARC, { order: 2 });
    const warm = solve(ARC, { order: 2, continuation: true });

    assert.ok(cold.report.ok && warm.report.ok);

    // §9: a warm start moves where the solve *begins*, never where it ends. Both runs drive
    // the same residual below the same tolerance, so the curves agree far tighter than the
    // eye — the bound here is the tolerance's geometric footprint, not a fudge factor.
    assert.ok(farthest(trace(cold.edges), trace(warm.edges)) < 1e-6);
  });

  it("takes fewer top-rung steps, which is the whole of §12's test", () => {
    for (const points of [ZIGZAG, ARC]) {
      const cold = solve(points, { order: 2 }).report;
      const warm = solve(points, { order: 2, continuation: true }).report;

      assert.ok(warm.steps <= cold.steps, `${warm.steps} steps warm against ${cold.steps} cold`);
      assert.ok(warm.seedSteps > 0);
    }

    // At least one of them has to actually improve, or the phase measured nothing.
    const cold = solve(ARC, { order: 2 }).report;
    const warm = solve(ARC, { order: 2, continuation: true }).report;

    assert.ok(warm.steps < cold.steps, `no improvement: ${warm.steps} against ${cold.steps}`);
  });

  it("will not seed from a rung §8 says cannot be trusted", () => {
    const { mesh } = polyline(TIGHT);
    const refs = referenceLengths(mesh);
    const [chain] = [...chains(mesh)];

    // `ladder`'s seeding predicate, evaluated on the rung it would seed *from*, at the caps it
    // has on the first attempt — none. Asserted here rather than inferred from two whole-solve
    // reports agreeing, which they can do for reasons that have nothing to do with seeding.
    const rung = new ChainSystem(chain, refs, { ...defaultSPowerSolverOptions, order: 0 });
    const run = rung.run();

    assert.equal(run.ok && run.residual < defaultSPowerSolverOptions.tolerance, false);
    assert.ok(rung.faults(run).length > 0, "the fixture stopped tripping any criterion");

    // What a discarded rung costs is steps, and only steps. Both solves reach the same curve
    // and break the same joints; the warm one has paid `seedSteps` for a rung it threw away.
    const cold = solve(TIGHT, { order: 1 });
    const warm = solve(TIGHT, { order: 1, continuation: true });

    assert.ok(cold.report.ok && warm.report.ok);
    assert.ok(warm.report.seedSteps > 0, "the ladder ran no lower rungs");
    assert.equal(warm.report.degraded, cold.report.degraded);
    assert.ok(farthest(trace(cold.edges), trace(warm.edges)) < 1e-6);
  });
});
