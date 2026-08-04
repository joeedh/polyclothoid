import assert from "node:assert/strict";
import { describe, test } from "node:test";
import { Mesh, SPowerClothoid, SPowerSolver, type Vertex } from "../src/index.js";

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

describe("the iteration cap", () => {
  test("a folded chain converges inside it", () => {
    // 100 was the old default and this band needs 105-188; the failure it produced was a
    // spurious `newton-not-converging` on geometry already correct to half a nanoradian.
    assert.equal(guarded.capped, 0);
  });
});
