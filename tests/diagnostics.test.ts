/**
 * Phase 5 of `docs/plans/spower-solver.md` §12 — the §8 diagnostic record.
 *
 * Two things are under test and they fail differently.
 *
 * The first is {@link geometricRate}, which is arithmetic with an exact answer: feed it a
 * geometric sequence and it must return the ratio. §8 asks for a fit rather than a one-step
 * ratio precisely because the Gauss-Newton iteration is non-normal, so the interesting case
 * is the sequence that grows before it converges.
 *
 * The second is the emission rule, which has no exact answer — a threshold is a judgement —
 * so what is pinned here is the *shape* of the record rather than its numbers: that a
 * healthy solve stays silent, that a solve which quietly produces garbage does not, that a
 * record is not the same thing as a failure, and that the same input twice gives the same
 * records. §8's third rule is the one worth a test of its own; a report that varies run to
 * run is worse than none when bisecting.
 *
 * The two conditions §8 admits are not free — the local pivoted QR for G1 rank deficiency
 * and the condition estimate of the `(1,1)` block — are Phase 6's, and are not exercised.
 */
import assert from "node:assert/strict";
import { describe, it } from "node:test";

import {
  BSplineSolver,
  BezierSolver,
  ClothoidSolver,
  Clothoid,
  CubicBezier,
  BSpline,
  type ChainRun,
  type Diagnostic,
  SPowerClothoid,
  SPowerSolver,
  diagnosticThresholds,
  geometricRate,
} from "../src/curve/index.js";
import { ChainSystem, referenceLengths, defaultSPowerSolverOptions } from "../src/curve/spower_solver.js";
import { chains } from "../src/curve/topology.js";
import { Mesh, type Edge } from "../src/mesh/index.js";

const ZIGZAG = [
  [0, 0],
  [1, 0.35],
  [2, 0.2],
  [3, 0.9],
  [4, 0.7],
];

/**
 * A near-doubling-back chain, which the solver drives to coefficients around `1e+80`.
 *
 * It then reports a zero angle gap, because every segment has collapsed to a point and two
 * points agree tangentially for free. This is the "looks right and is not" case §8 exists
 * for, and it is why the canonical chord is checked against `1` as well as against zero.
 */
const HAIRPIN = [
  [0, 0],
  [1, 0],
  [0.02, 0.1],
  [1, 0.2],
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

function only(diagnostics: Diagnostic[], condition: string) {
  const found = diagnostics.filter((d) => d.condition === condition);

  assert.equal(
    found.length,
    1,
    `expected one ${condition}, got ${JSON.stringify(diagnostics.map((d) => d.condition))}`
  );

  return found[0];
}

/** Everything about a record except the element references, which are identity-compared. */
function shape(diagnostics: Diagnostic[]) {
  return diagnostics.map((d) => [d.condition, d.severity, d.action, d.chain, d.at, d.measured, d.threshold]);
}

describe("Phase 5: geometricRate", () => {
  it("recovers the ratio of an exact geometric sequence", () => {
    const history = [1.0, 0.3, 0.09, 0.027, 0.0081];

    assert.ok(Math.abs(geometricRate(history) - 0.3) < 1e-12, `${geometricRate(history)}`);
  });

  it("fits across holes, since a residual of zero carries no logarithm", () => {
    // 0.5^k at k = 0, 2, 4, with the odd entries unusable.
    assert.ok(Math.abs(geometricRate([1.0, 0.0, 0.25, NaN, 0.0625]) - 0.5) < 1e-12);
  });

  it("returns NaN rather than a rate when fewer than three samples survive", () => {
    assert.ok(Number.isNaN(geometricRate([1.0, 0.5])));
    assert.ok(Number.isNaN(geometricRate([1.0, 0.0, 0.5, -1.0])));
    assert.ok(Number.isNaN(geometricRate([])));
  });

  /*
    §8's reason for asking for a fit at all: one step grew by four and the sequence still
    converges, so any single-step ratio test would have called this divergence.
  */
  it("calls a sequence that grows before it converges convergent", () => {
    const history = [1.0, 4.0, 0.5, 0.06, 0.008, 0.001];

    assert.ok(history[1] > history[0]);
    assert.ok(geometricRate(history) < 1.0, `${geometricRate(history)}`);
  });

  it("calls a growing sequence divergent", () => {
    assert.ok(geometricRate([1e-3, 3e-3, 9e-3, 2.7e-2]) > 1.0);
  });
});

describe("Phase 5: the report", () => {
  it("says nothing about a solve that went well", () => {
    const report = new SPowerSolver(polyline(ZIGZAG).mesh).solve();

    assert.equal(report.ok, true);
    assert.deepEqual(report.diagnostics, []);
    assert.equal(report.traces, undefined);
  });

  it("comes back clean from the solvers that have nothing to measure", () => {
    for (const [CurveCls, SolverCls] of [
      [CubicBezier, BezierSolver],
      [BSpline, BSplineSolver],
      [Clothoid, ClothoidSolver],
    ] as const) {
      const { mesh } = polyline(ZIGZAG);

      mesh.switchSplineType(CurveCls, SolverCls);

      const report = new SolverCls(mesh).solve();

      assert.equal(report.ok, true, CurveCls.name);
      assert.deepEqual(report.diagnostics, [], CurveCls.name);
    }
  });

  it("reaches the client through the mesh, which also keeps it", () => {
    const { mesh } = polyline(ZIGZAG);

    const report = mesh.solve();

    assert.equal(report.ok, true);
    assert.equal(mesh.report, report);
  });

  it("gives identical records for identical input", () => {
    const first = new SPowerSolver(polyline(HAIRPIN).mesh).solve();
    const second = new SPowerSolver(polyline(HAIRPIN).mesh).solve();

    assert.ok(first.diagnostics.length > 0, "fixture stopped being diagnostic");
    assert.deepEqual(shape(first.diagnostics), shape(second.diagnostics));
  });
});

describe("Phase 5: what gets recorded", () => {
  /*
    Phase 2's frozen Jacobian on a four-edge chain — §14's stall. It is the one fixture in
    the suite that genuinely fails to converge, so it is what pins both Newton criteria.
  */
  it("reports the frozen Jacobian's stall as not converging", () => {
    const report = new SPowerSolver(polyline(ZIGZAG).mesh, { jacobian: "frozen" }).solve();

    assert.ok(report.maxResidual > 1e-3, `unexpectedly converged to ${report.maxResidual}`);

    const newton = only(report.diagnostics, "newton-not-converging");

    assert.equal(newton.severity, "error");
    assert.equal(newton.action, "none");
    assert.equal(newton.at, -1);
    assert.ok(newton.measured >= diagnosticThresholds.newtonRate, `rate ${newton.measured}`);

    // The line search runs out of room every pass, which is §8's third Newton signal.
    const search = only(report.diagnostics, "line-search-failing");

    assert.equal(search.severity, "error");
    assert.ok(search.measured >= diagnosticThresholds.lineSearch);
  });

  it("reports a solve cut off before it could converge", () => {
    const report = new SPowerSolver(polyline(ZIGZAG).mesh, { iterations: 2 }).solve();

    const newton = only(report.diagnostics, "newton-not-converging");

    assert.equal(newton.severity, "error");
    assert.ok(Number.isNaN(newton.measured), "two steps is not enough to fit a rate");
  });

  /*
    A record is not a failure. Refinement off and `δ` two hundred thousand times its default
    leaves a visible linear-solve residual, and the outer iteration absorbs it and converges
    anyway — which is exactly the measurement §8's first rule wants recorded regardless.
  */
  it("reports a linear solve it nonetheless recovered from", () => {
    const report = new SPowerSolver(polyline(ZIGZAG).mesh, { kkt: { delta: 1e-2, refinement: 0 } }).solve();

    assert.equal(report.ok, true);
    assert.ok(report.maxResidual < 1e-10, `residual ${report.maxResidual}`);

    const refinement = only(report.diagnostics, "refinement-not-converging");

    assert.equal(refinement.severity, "warning");
    assert.ok(refinement.measured > diagnosticThresholds.refinement);
  });

  it("locates a canonical chord that has closed on itself", () => {
    const { mesh, edges } = polyline(ZIGZAG);
    const refs = referenceLengths(mesh);
    const [chain] = [...chains(mesh)];

    const system = new ChainSystem(chain, refs, defaultSPowerSolverOptions);
    const run = system.run();

    // A constant canonical profile is a circular arc; at five radians it has nearly closed.
    const curve = edges[1].curve as SPowerClothoid;

    curve.ks.fill(5.0);
    curve.recalc = 1;

    assert.ok(curve.length > 0.0, "reading the length is what forces the update");
    assert.ok(curve.canonical < diagnosticThresholds.chordDegeneracy, `canonical ${curve.canonical}`);

    const found: Diagnostic[] = [];

    system.diagnose(0, run, found);

    const chord = only(found, "chord-degeneracy");

    assert.equal(chord.severity, "warning");
    assert.equal(chord.at, 1);
    assert.equal(chord.edge, edges[1]);
    assert.equal(chord.measured, curve.canonical);
    assert.equal(chord.threshold, diagnosticThresholds.chordDegeneracy);
  });

  /*
    The reason the canonical chord is checked from both ends. Here the coefficients run away,
    every segment collapses, and the angle gap goes to zero because two points meet
    tangentially for free — a solve that reports success and has produced nothing.
  */
  it("catches a runaway from the other side, where the residual claims success", () => {
    const { mesh, edges } = polyline(HAIRPIN);

    const report = new SPowerSolver(mesh).solve();

    assert.ok(report.maxResidual < 1e-12, `fixture stopped claiming success: ${report.maxResidual}`);

    const chord = report.diagnostics.filter((d) => d.condition === "chord-degeneracy");

    assert.equal(chord.length, edges.length);

    for (const d of chord) {
      assert.equal(d.severity, "error");
      assert.ok(d.measured > 1.0, `canonical chord ${d.measured} is not out of range`);
    }
  });

  it("reports a failed factorization once, and not as a refinement fault", () => {
    const { mesh } = polyline(ZIGZAG);
    const refs = referenceLengths(mesh);
    const [chain] = [...chains(mesh)];

    const system = new ChainSystem(chain, refs, defaultSPowerSolverOptions);

    const failed: ChainRun = {
      steps     : 1,
      residual  : 0.5,
      factored  : false,
      ok        : false,
      history   : [0.5],
      backtracks: 0,
      starved   : false,
      refinement: Infinity,
      multiplier: 0.0,
    };

    const found: Diagnostic[] = [];

    system.diagnose(3, failed, found);

    const underflow = only(found, "factorization-underflow");

    assert.equal(underflow.severity, "error");
    assert.equal(underflow.chain, 3);
    assert.equal(underflow.at, -1);
    assert.equal(found.filter((d) => d.condition === "refinement-not-converging").length, 0);
  });
});

describe("Phase 5: the trace", () => {
  it("records one step per step, and only when asked", () => {
    const { mesh } = polyline(ZIGZAG);

    const report = new SPowerSolver(mesh, { trace: true }).solve();

    assert.equal(report.traces?.length, 1);

    const trace = report.traces![0];

    assert.equal(trace.chain, 0);
    assert.equal(trace.steps.length, report.steps);

    for (const step of trace.steps) {
      assert.ok(Number.isFinite(step.merit), "merit");
      assert.ok(Number.isFinite(step.stepNorm), "stepNorm");
      assert.ok(step.maxMultiplier > 0.0, "a joint held by no force is not being held");
      assert.ok(step.relaxation > 0.0 && step.relaxation <= 1.0, `relaxation ${step.relaxation}`);
      assert.ok(step.refinement >= 0.0, "refinement");
    }

    // The step count is the convergence claim, so the trace has to end where the report does.
    assert.equal(trace.steps[trace.steps.length - 1].residual, report.maxResidual);
  });

  it("carries the history that shows what a stall looks like", () => {
    const { mesh } = polyline(ZIGZAG);

    const report = new SPowerSolver(mesh, { jacobian: "frozen", trace: true, iterations: 12 }).solve();
    const steps = report.traces![0].steps;

    assert.equal(steps.length, 12);

    // Flat, not falling: the last third of the run is doing nothing the first third did not.
    assert.ok(steps[11].residual > 0.5 * steps[7].residual, "the fixture started converging");
  });
});
