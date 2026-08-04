/**
 * Phase 6 of `docs/plans/spower-solver.md` §12 — authored levels, and §8's stability breaking.
 *
 * Three things are under test, and only the first has an exact answer.
 *
 * **The partitions** (§3) are combinatorics: given a set of authored pairings, which ends share
 * which block entries is a fact, and the interesting case is the transitivity a pairwise lookup
 * would miss. **The layout** is checked through the geometry it produces rather than through
 * indices — a level is a claim about what is continuous, so the test is a tangent gap and a
 * curvature jump, measured on the solved curves.
 *
 * **The breaking** has no exact answer, because a threshold is a judgement. What is pinned is
 * §8's two requirements on any break: that it is *reported* — requested level, delivered level,
 * and the condition that caused it — and that it is *hysteretic and deterministic*, which here
 * means the same input twice delivers the same levels, and a joint that lost one is judged
 * against the stricter restore threshold when it asks for it back.
 */
import assert from "node:assert/strict";
import { describe, it } from "node:test";

import { columnPivotedQR, rankRatio } from "../src/math/index.js";
import {
  Clothoid,
  ClothoidSolver,
  SPowerClothoid,
  SPowerSolver,
  type Diagnostic,
  defaultLevel,
  maxLevel,
  pairingLevel,
  sharing,
  stabilityThresholds,
  vertexPartitions,
} from "../src/curve/index.js";
import { ChainSystem, referenceLengths, defaultSPowerSolverOptions } from "../src/curve/spower_solver.js";
import { chains } from "../src/curve/topology.js";
import { Mesh, type CurveConstructor, type CurveSolverConstructor, type Edge, type Vertex } from "../src/mesh/index.js";
import { DEFAULT_CORNER_THRESHOLD, Stroker } from "../src/stroke.js";

const ZIGZAG = [
  [0, 0],
  [1, 0.35],
  [2, 0.2],
  [3, 0.9],
  [4, 0.7],
];

function polyline(
  points: number[][],
  CurveCls: CurveConstructor = SPowerClothoid,
  SolverCls: CurveSolverConstructor = SPowerSolver
) {
  const mesh = new Mesh();

  mesh.CurveCls = CurveCls;
  mesh.SolverCls = SolverCls;

  const verts = points.map((p) => mesh.makeVertex(p));
  const edges: Edge[] = [];

  for (let i = 0; i + 1 < verts.length; i++) {
    edges.push(mesh.makeEdge(verts[i], verts[i + 1]));
  }

  return { mesh, verts, edges };
}

/** Signed angle in `(−π, π]`, so a gap of `2π` does not read as a discontinuity. */
function wrap(a: number) {
  return a - Math.PI * 2.0 * Math.round(a / (Math.PI * 2.0));
}

/** The tangent-angle gap and the curvature jump across `verts[i]`, on the solved curves. */
function joint(edges: Edge[], i: number) {
  const before = edges[i - 1].curve as SPowerClothoid;
  const after = edges[i].curve as SPowerClothoid;

  const a = before.derivative(before.length);
  const b = after.derivative(0);

  return {
    tangent  : Math.abs(wrap(Math.atan2(b[1], b[0]) - Math.atan2(a[1], a[0]))),
    curvature: Math.abs(before.curvature(before.length) - after.curvature(0)),
  };
}

/** A stroker whose dabs go nowhere — {@link Stroker.markCorners} does not emit any. */
function corners(options = {}) {
  return new Stroker(() => undefined, options);
}

function pick(diagnostics: Diagnostic[], condition: string) {
  return diagnostics.filter((d) => d.condition === condition);
}

/** Everything about a record except the element references, which are identity-compared. */
function shape(diagnostics: Diagnostic[]) {
  return diagnostics.map((d) => [d.condition, d.severity, d.action, d.at, d.requested, d.delivered]);
}

describe("Phase 6: levels and partitions", () => {
  it("defaults valence 2 to the ceiling and a junction to fully split", () => {
    const { mesh, verts } = polyline(ZIGZAG);

    mesh.makeEdge(verts[2], mesh.makeVertex([2, 1.5]));

    assert.equal(defaultLevel(verts[1], 1), maxLevel(1));
    assert.equal(defaultLevel(verts[1], 1), 3);

    // No canonical pairing exists at a junction, and §3 forbids the solver guessing one.
    assert.equal(verts[2].edges.length, 3);
    assert.equal(defaultLevel(verts[2], 1), 0);
  });

  it("clamps an authored level to what the order can represent, and calls it no loss", () => {
    const { verts, edges } = polyline(ZIGZAG);

    verts[1].setPairing(edges[0], edges[1], 9);

    // `p = 1` has no fourth derivative to be continuous in, so nothing was asked for and lost.
    assert.equal(pairingLevel(verts[1], edges[0], edges[1], 1), 3);
    assert.equal(pairingLevel(verts[1], edges[0], edges[1], 3), 5);

    verts[1].setPairing(edges[0], edges[1], -2);
    assert.equal(pairingLevel(verts[1], edges[0], edges[1], 1), 0);
  });

  it("turns a level into shared entries, one per order above the first", () => {
    const { verts, edges } = polyline(ZIGZAG);
    const at = (level: number) => {
      verts[1].setPairing(edges[0], edges[1], level);

      return sharing(verts[1], edges[0], edges[1], 1);
    };

    assert.deepEqual(at(0), { entries: 0, tangent: false });
    assert.deepEqual(at(1), { entries: 0, tangent: true });
    assert.deepEqual(at(2), { entries: 1, tangent: true });
    assert.deepEqual(at(3), { entries: 2, tangent: true });
  });

  /*
    The reason §3 asks for union-finds rather than pairwise lookups. A and C share κ through
    B without anyone having authored A–C, and no pairwise reading of the pairings sees it.
  */
  it("groups two ends nobody paired, through a third they both share", () => {
    const mesh = new Mesh();
    const hub = mesh.makeVertex([0, 0]);

    const a = mesh.makeEdge(hub, mesh.makeVertex([1, 0]));
    const b = mesh.makeEdge(hub, mesh.makeVertex([0, 1]));
    const c = mesh.makeEdge(hub, mesh.makeVertex([-1, 0]));

    hub.setPairing(a, b, 3);
    hub.setPairing(b, c, 2);

    assert.equal(pairingLevel(hub, a, c, 1), 0, "A–C is unauthored, so it defaults to split");
    assert.equal(sharing(hub, a, c, 1).entries, 1, "and shares κ anyway, through B");
    assert.equal(sharing(hub, a, c, 1).tangent, true);

    // The partitions coarsen: order 1 keeps C out, order 0 and the tangent row take it in.
    const parts = vertexPartitions(hub, 1);

    assert.equal(parts[0][0], parts[0][2], "tangent");
    assert.equal(parts[1][0], parts[1][2], "order 0");
    assert.notEqual(parts[2][0], parts[2][2], "order 1 — only A–B reaches level 3");
  });
});

describe("Phase 6: what a level delivers", () => {
  /*
    A level is a claim about what is continuous, so it is measured on the solved curves. The
    corner is the one that also changes the *size* of the system, since a level-0 joint has no
    multiplier and no shared entries.
  */
  it("delivers exactly the continuity that was authored", () => {
    const expected = [
      { level: 0, tangent: 1.0, curvature: 0.0 },
      { level: 1, tangent: 1e-9, curvature: 0.1 },
      { level: 2, tangent: 1e-9, curvature: 1e-9 },
      { level: 3, tangent: 1e-9, curvature: 1e-9 },
    ];

    for (const { level, tangent, curvature } of expected) {
      const { mesh, verts, edges } = polyline(ZIGZAG);

      verts[2].setPairing(edges[1], edges[2], level);

      const report = new SPowerSolver(mesh).solve();
      const gap = joint(edges, 2);

      assert.equal(report.ok, true, `level ${level}`);

      if (level === 0) {
        assert.ok(gap.tangent > tangent, `level 0 kept a corner: ${gap.tangent}`);
      } else {
        assert.ok(gap.tangent < tangent, `level ${level} tangent gap ${gap.tangent}`);
      }

      if (level === 1) {
        assert.ok(gap.curvature > curvature, `level 1 kept a curvature jump: ${gap.curvature}`);
      } else if (level > 1) {
        assert.ok(gap.curvature < curvature, `level ${level} curvature jump ${gap.curvature}`);
      }

      // The joints nobody touched still meet at the ceiling.
      assert.ok(joint(edges, 1).curvature < 1e-9, "an untouched joint lost continuity");
    }
  });

  it("spends fewer unknowns on a corner than on a G3 joint", () => {
    const block = defaultSPowerSolverOptions.dof.blockLength(defaultSPowerSolverOptions.order);
    const size = (level: number) => {
      const { mesh, verts, edges } = polyline(ZIGZAG);

      verts[2].setPairing(edges[1], edges[2], level);

      const system = new ChainSystem([...chains(mesh)][0], referenceLengths(mesh), defaultSPowerSolverOptions);

      return { n: system.n, rows: system.rows.length, shared: [...system.shared] };
    };

    const ceiling = size(3);
    const corner = size(0);

    // A corner splits the block in two and drops the multiplier: `+block`, `−1`.
    assert.equal(corner.n, ceiling.n + block - 1);
    assert.equal(corner.rows, ceiling.rows - 1);
    assert.equal(corner.shared[2], 0);
    assert.equal(ceiling.shared[2], block);

    // G1 shares nothing either, but keeps the row that makes the tangents agree.
    assert.equal(size(1).rows, ceiling.rows);
    assert.equal(size(1).shared[2], 0);
  });
});

describe("Phase 6: the rank test", () => {
  it("counts an exactly dependent column out", () => {
    // Third column is the sum of the first two, so the rank is 2 whatever the pivoting order.
    const a = Float64Array.from([1, 0, 1, 0, 1, 1, 1, 1, 2]);
    const { rank, diagonal, pivots } = columnPivotedQR(a, 3, 3);

    assert.equal(rank, 2);
    assert.equal(pivots.length, 3);
    assert.ok(diagonal[0] >= diagonal[1], "pivoting must leave the diagonal non-increasing");
    assert.ok(diagonal[2] < 1e-14, `dependent pivot ${diagonal[2]}`);
  });

  it("calls an all-zero block rank 0 rather than rank 1", () => {
    assert.equal(columnPivotedQR(new Float64Array(6), 2, 3).rank, 0);
    assert.equal(rankRatio(new Float64Array(6), 2, 3), 0.0);
  });

  it("is invariant under a uniform scaling of the block", () => {
    const rows = [1, 2, 3, 2, 4.5, 6];
    const plain = rankRatio(Float64Array.from(rows), 2, 3);
    const scaled = rankRatio(
      Float64Array.from(rows, (v) => v * 1e6),
      2,
      3
    );

    assert.ok(Math.abs(plain - scaled) < 1e-12, `${plain} vs ${scaled}`);
  });

  /*
    §8's criterion is vacuous at valence 2 and has to survive being so: a chain joint's local
    block is a single row, whose only way to be deficient is to vanish outright. That it can
    vanish at all is the frozen Jacobian's doing, which is what the second half checks.
  */
  it("reads a chain joint as full rank, and a vanished row as not", () => {
    const { mesh } = polyline(ZIGZAG);
    const refs = referenceLengths(mesh);
    const [chain] = [...chains(mesh)];

    const system = new ChainSystem(chain, refs, defaultSPowerSolverOptions);

    system.run();

    for (let i = 1; i < system.rows.length + 1; i++) {
      assert.ok(system.rankAt(i) > stabilityThresholds.rankBreak, `joint ${i} rank ${system.rankAt(i)}`);
    }

    const frozen = new ChainSystem(chain, refs, { ...defaultSPowerSolverOptions, jacobian: "frozen" });

    frozen.run();

    assert.equal(frozen.rankAt(0), 0.0, "a chain end has no row to rank");
  });
});

describe("Phase 6: stability breaking", () => {
  /*
    The frozen Jacobian on a four-edge chain is §14's stall — the one fixture in the suite
    that genuinely fails to converge. Phase 6 is what is allowed to do something about it.
  */
  it("lowers a joint rather than handing back a stalled chain", () => {
    const { mesh } = polyline(ZIGZAG);
    const report = new SPowerSolver(mesh, { jacobian: "frozen" }).solve();

    assert.ok(report.attempts > 1, "nothing was retried");
    assert.equal(report.degraded, 1);
    assert.ok(report.maxResidual < 1e-8, `still stalled at ${report.maxResidual}`);

    const broken = report.diagnostics.filter((d) => d.action === "degraded");

    assert.ok(broken.length > 0);

    for (const d of broken) {
      assert.equal(d.condition, "newton-not-converging");
      assert.equal(d.requested, 3, "the authored level is what was asked for");
      assert.ok(d.delivered! < d.requested!, `${d.requested} -> ${d.delivered}`);
      assert.ok(d.vertex, "a break has to say where");
    }

    // Walked down one level at a time, never further than it had to.
    assert.deepEqual(
      broken.map((d) => d.delivered),
      [2, 1]
    );
  });

  it("refuses instead of degrading when the caller asked for engineering mode", () => {
    const { mesh } = polyline(ZIGZAG);
    const report = new SPowerSolver(mesh, { jacobian: "frozen", mode: "engineering" }).solve();

    assert.equal(report.ok, false, "a refusal is not a successful solve");
    assert.equal(report.degraded, 0, "engineering mode may not break anything");
    assert.equal(report.attempts, 1, "and has nothing to retry");

    const refused = report.diagnostics.filter((d) => d.action === "refused");

    assert.ok(refused.length > 0, "a refusal that says nothing is worse than a degradation");

    for (const d of refused) {
      assert.equal(d.severity, "error");
      assert.ok(d.vertex);
    }
  });

  it("leaves a healthy chain alone", () => {
    const { mesh } = polyline(ZIGZAG);
    const report = new SPowerSolver(mesh).solve();

    assert.equal(report.attempts, 1);
    assert.equal(report.degraded, 0);
    assert.equal(report.diagnostics.filter((d) => d.action !== "none").length, 0);
  });

  /*
    A budget is not a divergence. §8 asks for a geometric fit over at least three iterations,
    and a run too short to produce one must not cost anyone a level — otherwise every
    interactive solve with a tight iteration cap would quietly shed continuity.
  */
  it("does not break a solve that merely ran out of iterations", () => {
    const { mesh } = polyline(ZIGZAG);
    const report = new SPowerSolver(mesh, { iterations: 2 }).solve();

    assert.equal(report.degraded, 0);
    assert.equal(report.attempts, 1);
  });

  it("obeys a breaks budget of zero by leaving the fault in place", () => {
    const { mesh } = polyline(ZIGZAG);
    const report = new SPowerSolver(mesh, { jacobian: "frozen", breaks: 0 }).solve();

    assert.equal(report.degraded, 0);
    assert.ok(report.maxResidual > 1e-3, "the stall is what it was told not to fix");
  });
});

describe("Phase 6: hysteresis", () => {
  /*
    §8's determinism requirement, which is the stronger half: identical input gives identical
    records, whatever order of solves reached it. The first attempt always asks for the
    authored levels, so history moves the *thresholds* and never the starting point.
  */
  it("delivers the same levels for the same input, twice", () => {
    const first = new SPowerSolver(polyline(ZIGZAG).mesh, { jacobian: "frozen" }).solve();
    const second = new SPowerSolver(polyline(ZIGZAG).mesh, { jacobian: "frozen" }).solve();

    assert.ok(first.degraded > 0, "fixture stopped breaking anything");
    assert.deepEqual(shape(first.diagnostics), shape(second.diagnostics));
    assert.equal(first.maxResidual, second.maxResidual);
  });

  it("holds a broken joint to the stricter threshold when it asks for its level back", () => {
    const { mesh } = polyline(ZIGZAG);
    const solver = new SPowerSolver(mesh, { jacobian: "frozen" });

    const first = solver.solve();
    const second = solver.solve();

    assert.equal(second.degraded, first.degraded, "the joint came back");
    assert.equal(second.maxResidual, first.maxResidual);

    // Same outcome, different news: the second pass is retaining a break, not making one.
    assert.equal(pick(first.diagnostics, "level-lowered").length, 0);
    assert.ok(pick(second.diagnostics, "level-lowered").length > 0);

    for (const d of pick(second.diagnostics, "level-lowered")) {
      assert.equal(d.severity, "info");
      assert.equal(d.action, "degraded");
      assert.equal(d.threshold, stabilityThresholds.rateRestore, "the restore threshold, not the break one");
    }
  });

  it("forgets what it withheld when the curve type changes underneath it", () => {
    const { mesh } = polyline(ZIGZAG);

    mesh.solve();
    mesh.switchSplineType(SPowerClothoid, SPowerSolver);

    assert.equal(mesh.solver, undefined);
  });
});

describe("Phase 6: the multiplier as corner evidence", () => {
  /*
    §8 is explicit that `|λ_v|` is not a fault: it says the joint is expensive to hold, which
    is a judgement about whether a corner was wanted there. So it is reported and never acted
    on, and the test is that it stays out of the break path.
  */
  it("reports an outlying multiplier as information, against the chain's own median", () => {
    const { mesh } = polyline(ZIGZAG);
    const refs = referenceLengths(mesh);
    const [chain] = [...chains(mesh)];

    const system = new ChainSystem(chain, refs, defaultSPowerSolverOptions);
    const run = system.run();

    assert.equal(system.rows.length, 3);

    system.z[system.rows[0]] = 1.0;
    system.z[system.rows[1]] = 1.0;
    system.z[system.rows[2]] = 100.0;

    const found: Diagnostic[] = [];

    system.diagnose(0, run, found);

    const large = pick(found, "large-multiplier");

    assert.equal(large.length, 1);
    assert.equal(large[0].severity, "info");
    assert.equal(large[0].action, "none", "the solver may not act on a multiplier");
    assert.equal(large[0].at, 3);
    assert.equal(large[0].measured, 100.0);
    assert.equal(large[0].threshold, stabilityThresholds.multiplierRatio);

    // And it is not a fault, so it can never cause a break.
    assert.equal(system.faults(run).length, 0);
  });
});

describe("Phase 6: cornerThreshold as client-side inference", () => {
  it("authors a corner where the input turned sharply, and nothing where it did not", () => {
    const { mesh, verts, edges } = polyline([
      [0, 0],
      [1, 0],
      [2, 0.05],
      [1.2, 0.15],
    ]);

    corners().markCorners(mesh);

    assert.equal(verts[1].pairings.length, 0, "a gentle turn is not a corner");
    assert.equal(verts[2].pairings.length, 1);
    assert.equal(verts[2].pairing(edges[1], edges[2])?.level, 0);

    // Either order, and the level is what the solver then reads.
    assert.equal(verts[2].pairing(edges[2], edges[1])?.level, 0);
    assert.equal(pairingLevel(verts[2], edges[1], edges[2], 1), 0);
  });

  it("authors none at all when the threshold is turned off", () => {
    const { mesh, verts } = polyline([
      [0, 0],
      [1, 0],
      [1, 0.02],
    ]);

    corners({ cornerThreshold: 0 }).markCorners(mesh);

    assert.equal(verts[1].pairings.length, 0);
    assert.ok(DEFAULT_CORNER_THRESHOLD > Math.PI * 0.3);
  });

  /*
    The port has to leave `ClothoidSolver` behaving as it did, which means the level it now
    reads has to reach the same corner-zeroing pass the angle test used to feed.
  */
  it("still leaves an authored corner as a corner in the clothoid solver", () => {
    const { mesh, verts, edges } = polyline(
      [
        [0, 0],
        [1, 0],
        [0.3, 0.1],
      ],
      Clothoid,
      ClothoidSolver
    );

    corners({ CurveCls: Clothoid, SolverCls: ClothoidSolver }).markCorners(mesh);

    assert.equal(verts[1].pairing(edges[0], edges[1])?.level, 0);

    mesh.solve();

    for (const e of edges) {
      const curve = e.curve as Clothoid;
      const end = (verts[1] as Vertex) === e.v1 ? 0 : curve.order - 1;

      assert.equal(curve.ks[end], 0.0, "the sharp end kept its curvature");
    }
  });
});
