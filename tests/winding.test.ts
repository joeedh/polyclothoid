import assert from "node:assert/strict";
import { describe, test } from "node:test";
import { Clothoid, ClothoidSolver, Mesh, SPowerClothoid, SPowerSolver, Stroker, type Vertex } from "../src/index.js";
import { defaultSPowerSolverOptions } from "../src/curve/spower_solver.js";

/**
 * Total turning of an edge, differenced out of `derivative()`.
 *
 * Deliberately not `SPowerClothoid.turning`: the point of these tests is the *geometry* the
 * solve produced, and reading the profile back would be asking the same object that produced
 * the shape whether the shape is right.
 */
function turned(edge: { curve: { derivative(s: number): ArrayLike<number>; length: number } }) {
  const steps = 256;
  const angle = (s: number) => {
    const d = edge.curve.derivative(s);
    return Math.atan2(d[1], d[0]);
  };

  let prev = angle(0.0);
  let total = 0.0;

  for (let i = 1; i <= steps; i++) {
    const th = angle((i / steps) * edge.curve.length);
    let delta = th - prev;

    while (delta > Math.PI) delta -= Math.PI * 2.0;
    while (delta < -Math.PI) delta += Math.PI * 2.0;

    total += delta;
    prev = th;
  }

  return total;
}

function solverWith(options: Partial<ConstructorParameters<typeof SPowerSolver>[1]>) {
  return class extends SPowerSolver {
    constructor(mesh: ConstructorParameters<typeof SPowerSolver>[0]) {
      super(mesh, options);
    }
  };
}

/** The three-vertex chain of the bug report, ready to be dragged. */
function folding(branchLimit?: number) {
  const mesh = new Mesh();

  mesh.switchSplineType(SPowerClothoid, branchLimit === undefined ? SPowerSolver : solverWith({ branchLimit }));

  const verts: Vertex[] = [mesh.makeVertex([0, 0, 0]), mesh.makeVertex([100, 0, 0]), mesh.makeVertex([200, 10, 0])];

  mesh.makeEdge(verts[0], verts[1]);
  mesh.makeEdge(verts[1], verts[2]);
  mesh.solve();

  /**
   * Walk the far vertex back past the near one, re-solving every frame.
   *
   * Warm-started throughout, because that is the whole point — a cold solve of any single
   * frame here is well behaved, and it is the drag that walks the iteration off its branch.
   */
  return (report: (m: Mesh) => void) => {
    for (let i = 0; i <= 250; i++) {
      verts[2][0] = 200 - i;
      mesh.regenSolve().ensureSolve();
      report(mesh);
    }

    return mesh;
  };
}

/** One pass of the drag, since it is slow enough that the tests should share it. */
function measured(branchLimit?: number) {
  let peak = 0.0;
  let cuts = 0;
  let capped = 0;

  folding(branchLimit)((m) => {
    const report = m.report as unknown as { branchCuts: number; steps: number };

    cuts += report.branchCuts;

    if (report.steps >= 400) {
      capped++;
    }

    for (const e of m.edges) {
      peak = Math.max(peak, Math.abs(turned(e)));
    }
  });

  return { peak, cuts, capped };
}

const guarded = measured();
const free = measured(Infinity);

describe("the winding branch", () => {
  test("a fold through collinearity does not wind the spiral", () => {
    // Under a full turn is the claim; in practice this fixture stays near 0.9.
    assert.ok(guarded.peak < Math.PI * 2.0, `wound to ${(guarded.peak / (Math.PI * 2)).toFixed(2)} turns`);
  });

  test("and without the guard it winds nine times, which is what the bound is for", () => {
    // Pinned deliberately: if this ever stops reproducing, the guard is no longer being
    // exercised by the test above, and the fixture needs rebuilding rather than deleting.
    assert.ok(
      free.peak > Math.PI * 2.0 * 8.0,
      `expected a multi-turn spiral, got ${(free.peak / (Math.PI * 2)).toFixed(2)} turns`
    );
  });

  test("the guard reports the steps it cut", () => {
    assert.ok(guarded.cuts > 0, "the fold should have aimed at least one step off-branch");
    assert.equal(free.cuts, 0);
  });

  test("it costs nothing on geometry that was never near a branch", () => {
    const shapes = (limit?: number) => {
      const mesh = new Mesh();

      mesh.switchSplineType(SPowerClothoid, limit === undefined ? SPowerSolver : solverWith({ branchLimit: limit }));

      const verts = [
        [40, 40],
        [160, 90],
        [260, 40],
        [330, 150],
      ].map((p) => mesh.makeVertex([p[0], p[1], 0]));

      for (let i = 0; i < verts.length - 1; i++) {
        mesh.makeEdge(verts[i], verts[i + 1]);
      }

      mesh.solve();

      return [...mesh.edges].map((e) => e.curve.length);
    };

    const bounded = shapes();
    const unbounded = shapes(Infinity);

    for (let i = 0; i < bounded.length; i++) {
      assert.ok(Math.abs(bounded[i] - unbounded[i]) < 1e-9, `edge ${i}: ${bounded[i]} vs ${unbounded[i]}`);
    }
  });
});

/**
 * A zigzag tight enough that no winding-free fit exists at the level it is authored with.
 *
 * §14's Phase 6 table already has this one degrading rather than converging cleanly. What it
 * adds here is that the degradation is the *branch's* doing: held on-branch it cannot land at
 * all, and the joint it gives up is the one the drifting edge was strained at.
 */
const TIGHT = [
  [0, 0],
  [1, 1.2],
  [2, 0],
  [3, 1.2],
  [4, 0],
  [5, 1.2],
];

function zigzag(options: Partial<ConstructorParameters<typeof SPowerSolver>[1]>) {
  const mesh = new Mesh();

  mesh.CurveCls = SPowerClothoid;
  mesh.SolverCls = SPowerSolver;

  const verts = TIGHT.map((p) => mesh.makeVertex([p[0], p[1], 0]));
  const edges = [];

  for (let i = 0; i + 1 < verts.length; i++) {
    edges.push(mesh.makeEdge(verts[i], verts[i + 1]));
  }

  return { edges, report: new SPowerSolver(mesh, { order: 1, ...options }).solve() };
}

/** Both fits are slow enough to be worth sharing, as with the drag above. */
const repaired = zigzag({});
const unrepaired = zigzag({ breaks: 0 });

describe("branch obstruction", () => {
  test("is raised against a joint, and the joint is lowered for it", () => {
    const found = repaired.report.diagnostics.filter((d) => d.condition === "branch-obstruction");

    // Counted by joint, not by record: a rung records the joint that blocked it, so a joint
    // the ladder retries through is named once per attempt.
    const joints = new Set(found.map((d) => d.at));

    assert.equal(joints.size, 1);

    for (const d of found) {
      assert.equal(d.action, "degraded");

      // The whole reason the guard's bookkeeping is per chain rather than per component: a
      // record that named the chain would have nothing for §8's ladder to act on.
      assert.ok(d.at > 0, `an obstruction must name an interior joint, got ${d.at}`);
    }
  });

  test("and what comes back is a curve, which without the ladder it is not", () => {
    assert.ok(repaired.report.ok);

    for (const e of repaired.edges) {
      assert.ok(Number.isFinite(e.curve.length), `edge length ${e.curve.length}`);
    }

    // Breaking is the only thing between these two. Same points, same guard, no levels to
    // spend: the coefficients run away and the solve says so rather than handing back `NaN`.
    assert.equal(unrepaired.report.ok, false);
    assert.ok(!Number.isFinite(unrepaired.report.maxResidual));
  });

  test("the budget it spends is per joint and not per chain", () => {
    // Four interior joints, and this fixture needs levels from two of them. Spent per chain,
    // the ladder would have stopped after `breaks` passes holding an unconverged rung.
    assert.ok(
      repaired.report.attempts > defaultSPowerSolverOptions.breaks,
      `${repaired.report.attempts} attempts against a per-chain budget of ${defaultSPowerSolverOptions.breaks}`
    );
  });
});

/**
 * `branchLimit` is the one solver option an application is expected to put in front of a
 * person, so it has a route that does not go through a solver subclass: `Mesh.solverTuning`,
 * `Mesh.setSolverTuning`, and `StrokerOptions.branchLimit`.
 *
 * What these check is that route. The bound's own behaviour is the drag above.
 */
describe("the tuning a host can set", () => {
  const bound = (m: Mesh) => (m.solver as SPowerSolver).options.branchLimit;

  function chain() {
    const mesh = new Mesh();

    mesh.switchSplineType(SPowerClothoid, SPowerSolver);

    const verts = [
      [0, 0],
      [100, 0],
      [200, 10],
    ].map((p) => mesh.makeVertex([p[0], p[1], 0]));

    mesh.makeEdge(verts[0], verts[1]);
    mesh.makeEdge(verts[1], verts[2]);
    mesh.solve();

    return mesh;
  }

  test("reaches the solver, and left alone is the solver's own default", () => {
    assert.equal(bound(chain()), defaultSPowerSolverOptions.branchLimit);
  });

  test("a change rebuilds the solver and an unchanged one does not", () => {
    const mesh = chain();
    const first = mesh.solver;

    mesh.setSolverTuning({ branchLimit: Math.PI }).ensureSolve();

    assert.equal(bound(mesh), Math.PI);
    assert.notEqual(mesh.solver, first);

    const second = mesh.solver;

    mesh.setSolverTuning({ branchLimit: Math.PI }).ensureSolve();

    // §8's hysteresis is per instance, so a panel echoing its own state must not cost it.
    assert.equal(mesh.solver, second);
  });

  test("a solver that cannot wind takes it and ignores it", () => {
    const mesh = chain();

    mesh.switchSplineType(Clothoid, ClothoidSolver);
    mesh.setSolverTuning({ branchLimit: Math.PI }).ensureSolve();

    assert.ok(mesh.solver instanceof ClothoidSolver);

    for (const e of mesh.edges) {
      assert.ok(e.curve.length > 0 && Number.isFinite(e.curve.length), `edge length ${e.curve.length}`);
    }
  });

  test("the stroker hands it to the mesh it fits", () => {
    let dabs = 0;

    const stroker = new Stroker(() => dabs++, {
      CurveCls   : SPowerClothoid,
      SolverCls  : SPowerSolver,
      branchLimit: Math.PI,
    });

    for (let i = 0; i < 8; i++) {
      stroker.onInput(i * 20, Math.sin(i * 0.5) * 15, 10, 0.1);
    }

    assert.ok(stroker.lastMesh, "eight inputs should have produced a fit");
    assert.ok(dabs > 0, "and the fit should have been walked");
    assert.equal(bound(stroker.lastMesh), Math.PI);
  });
});

describe("the iteration cap", () => {
  test("a folded chain converges inside it", () => {
    // 100 was the old default and this band needs 105-188; the failure it produced was a
    // spurious `newton-not-converging` on geometry already correct to half a nanoradian.
    assert.equal(guarded.capped, 0);
  });
});
