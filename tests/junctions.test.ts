/**
 * Phase 8 of `docs/plans/spower-solver.md` §12 — general valence: the node/chain decomposition,
 * the per-order partitions read as an interface, and the Schur complement that eliminates onto it.
 *
 * Two halves, and they fail differently. The **decomposition** is combinatorics with an exact
 * answer: given the authored pairings, which chains have to be solved together, and which
 * unknowns is it that makes them. The **elimination** has no exact answer to compare against, so
 * it is measured the way §5 says a coupling either works or does not — the geometry at the
 * junction. A shared entry that never reached the other chain shows up as a tangent gap or a
 * curvature jump straight through the node.
 *
 * The three cases §5 singles out each get a test: a plain junction, which must produce *no*
 * interface at all and fall back to Phase 2's independent banded solves; a level-1 pairing,
 * which shares no block and couples through a G1 row alone; and a closed chain, cut at an
 * arbitrary vertex and reconciled through the same machinery as two chains meeting at a node.
 */
import assert from "node:assert/strict";
import { describe, it } from "node:test";

import { Vector2 } from "../src/math/index.js";
import {
  SPowerClothoid,
  SPowerSolver,
  type Component,
  chainEnds,
  chains,
  components,
  cutOpen,
  interfaceSlots,
  maxLevel,
} from "../src/curve/index.js";
import { Mesh, type Edge, type Vertex } from "../src/mesh/index.js";

/** Signed angle in `(−π, π]`. */
function wrap(a: number) {
  return a - Math.PI * 2.0 * Math.round(a / (Math.PI * 2.0));
}

/** Direction of travel at `v`, walking *into* it along `e` if `arriving`, out of it if not. */
function travelAngle(e: Edge, v: Vertex, arriving: boolean) {
  const c = e.curve as SPowerClothoid;
  const atV2 = v === e.v2;

  const t = new Vector2(atV2 ? c.derivative(c.length) : c.derivative(0.0));

  if (atV2 !== arriving) {
    t.negate();
  }

  return Math.atan2(t[1], t[0]);
}

/** Signed curvature at `v` in the same direction of travel {@link travelAngle} uses. */
function travelCurvature(e: Edge, v: Vertex, arriving: boolean) {
  const c = e.curve as SPowerClothoid;
  const atV2 = v === e.v2;
  const k = atV2 ? c.curvature(c.length) : c.curvature(0.0);

  return atV2 === arriving ? k : -k;
}

/** World-space G1 gap at `v` on the path that arrives along `e1` and leaves along `e2`. */
function jointGap(e1: Edge, e2: Edge, v: Vertex) {
  return Math.abs(wrap(travelAngle(e2, v, false) - travelAngle(e1, v, true)));
}

/** World-space G2 gap on the same path. */
function curvatureGap(e1: Edge, e2: Edge, v: Vertex) {
  return Math.abs(travelCurvature(e2, v, false) - travelCurvature(e1, v, true));
}

/** The component holding `e`, which is unique because every edge lands in exactly one chain. */
function componentOf(list: Component[], e: Edge) {
  const found = list.find((c) => c.chains.some((chain) => (chain.edges as unknown[]).includes(e)));

  assert.ok(found, "every edge belongs to a component");

  return found;
}

/** Every chain of a mesh, cut open — what the solver decomposes. */
function paths(mesh: Mesh) {
  return chains(mesh).map(cutOpen);
}

/** A hub of three two-edge arms, laid out so A and B pass roughly straight through it. */
function hub() {
  const mesh = new Mesh();

  mesh.CurveCls = SPowerClothoid;
  mesh.SolverCls = SPowerSolver;

  const v = mesh.makeVertex([0, 0]);
  const arms = [
    [
      [1, 0.2],
      [2, 0.55],
    ],
    [
      [-1, 0.15],
      [-2, 0.5],
    ],
    [
      [0.1, -1],
      [0.4, -2],
    ],
  ].map((points) => {
    const verts = points.map((p) => mesh.makeVertex(p));

    return {
      edges: [mesh.makeEdge(v, verts[0]), mesh.makeEdge(verts[0], verts[1])],

      /** The arm's own interior joint, which stays chain-private however the hub is paired. */
      mid: verts[0],
    };
  });

  return { mesh, v, a: arms[0].edges, b: arms[1].edges, c: arms[2].edges, mids: arms.map((x) => x.mid) };
}

/** A regular `n`-gon on the unit circle, as a closed all-valence-2 loop. */
function ring(n: number) {
  const mesh = new Mesh();

  mesh.CurveCls = SPowerClothoid;
  mesh.SolverCls = SPowerSolver;

  const verts: Vertex[] = [];

  for (let i = 0; i < n; i++) {
    const t = (i / n) * Math.PI * 2.0;

    verts.push(mesh.makeVertex([Math.cos(t), Math.sin(t)]));
  }

  const edges = verts.map((v, i) => mesh.makeEdge(v, verts[(i + 1) % n]));

  return { mesh, verts, edges };
}

describe("Phase 8: the interface", () => {
  it("promotes nothing at a junction nobody paired, which is §3's default", () => {
    const { mesh, v } = hub();

    assert.equal(v.edges.length, 3);
    assert.deepEqual(interfaceSlots(paths(mesh), 1), []);

    // The point of the default: with an empty interface every chain is its own component, and
    // the Schur complement is 0 × 0 — exactly the independent banded solves Phase 2 ran.
    const found = components(paths(mesh), 1);

    assert.equal(found.length, 3);
    assert.deepEqual(
      found.map((c) => c.chains.length),
      [1, 1, 1]
    );
    assert.deepEqual(
      found.map((c) => c.gamma.length),
      [0, 0, 0]
    );
  });

  /*
    §5's first coupling caveat: a level-1 pairing shares no block entry and still couples, so
    the interface has to be the union of shared blocks *and* node G1 multipliers.
  */
  it("couples through a G1 row alone when the pairing shares no block", () => {
    const { mesh, v, a, b, c } = hub();

    v.setPairing(a[0], b[0], 1);

    const slots = interfaceSlots(paths(mesh), 1);

    assert.equal(slots.length, 1);
    assert.equal(slots[0].kind, "row", "level 1 is a row and nothing else");

    const found = components(paths(mesh), 1);

    assert.equal(found.length, 2);
    assert.equal(componentOf(found, a[0]), componentOf(found, b[0]), "A and B solve together");
    assert.notEqual(componentOf(found, a[0]), componentOf(found, c[0]), "C is untouched");
  });

  it("promotes one entry per shared order, plus the row that comes with it", () => {
    const { mesh, v, a, b } = hub();

    const at = (level: number) => {
      v.setPairing(a[0], b[0], level);

      const slots = interfaceSlots(paths(mesh), 1);

      return {
        entries: slots.filter((s) => s.kind === "entry").map((s) => (s.kind === "entry" ? s.order : -1)),
        rows   : slots.filter((s) => s.kind === "row").length,
      };
    };

    assert.deepEqual(at(0), { entries: [], rows: 0 });
    assert.deepEqual(at(1), { entries: [], rows: 1 });
    assert.deepEqual(at(2), { entries: [0], rows: 1 });
    assert.deepEqual(at(3), { entries: [0, 1], rows: 1 });
  });

  /*
    A group of m ends needs a spanning tree of G1 rows, m − 1 of them and not m(m−1)/2. The
    extra rows would be exactly dependent, and a KKT block with dependent constraint rows is
    what §6 spends the quasi-definite shift keeping factorable — it should not have to.
  */
  it("spends m − 1 rows on a tangent group of m, not every pair", () => {
    const { mesh, v, a, b, c } = hub();

    v.setPairing(a[0], b[0], 1);
    v.setPairing(b[0], c[0], 1);

    const slots = interfaceSlots(paths(mesh), 1);

    assert.equal(slots.filter((s) => s.kind === "row").length, 2, "three ends, two rows");
    assert.equal(components(paths(mesh), 1).length, 1, "and all three chains in one component");
  });

  /*
    The grouping is p-invariant even though the slot list is not, which is what lets §9's
    ladder re-derive only the interface per rung. Lowering p drops the *higher* shared orders,
    and tangent grouping survives every clamp down to level 1.
  */
  it("groups the same chains at every order, however the interface shrinks", () => {
    const { mesh, v, a, b, c } = hub();

    v.setPairing(a[0], b[0], 3);
    v.setPairing(b[0], c[0], 2);

    const sizes = new Set<string>();

    for (let p = 0; p <= 3; p++) {
      const found = components(paths(mesh), p);

      sizes.add(found.map((x) => x.chains.length).join(","));

      assert.ok(interfaceSlots(paths(mesh), p).length >= 2, `p = ${p} keeps both rows`);
    }

    assert.deepEqual([...sizes], ["3"], "one component of three chains at every order");
    assert.equal(maxLevel(0), 2, "and at p = 0 the level-3 pairing has clamped down to the level-2 one");
  });

  /*
    §5's closed loops. A cyclic band is not a band and not a rank-one correction of one, so an
    arbitrary vertex is promoted to a cut point and the wrap-around goes through the junction
    machinery. The cut vertex is valence 2, so it keeps the ceiling level and shares everything.
  */
  it("makes a closed chain couple to itself through the cut", () => {
    const { mesh } = ring(6);
    const all = paths(mesh);

    assert.equal(all.length, 1);
    assert.equal(all[0].closed, true, "the cut survives, so a caller knows the vertex is doubled");
    assert.equal(all[0].verts.length, all[0].edges.length + 1);
    assert.equal(all[0].verts[0], all[0].verts[all[0].verts.length - 1]);

    const ends = chainEnds(all).get(all[0].verts[0]);

    assert.equal(ends?.length, 2, "both ends land on the cut vertex");
    assert.notEqual(ends?.[0].edge, ends?.[1].edge);

    const found = components(all, 1);

    assert.equal(found.length, 1);
    assert.equal(found[0].chains.length, 1, "one chain — and an interface all the same");
    assert.deepEqual(
      found[0].gamma.map((s) => s.kind),
      ["entry", "entry", "row"]
    );
  });
});

describe("Phase 8: the coupled solve", () => {
  it("reaches tangent continuity straight through a paired junction", () => {
    const { mesh, v, a, b, c, mids } = hub();

    v.setPairing(a[0], b[0], 1);

    const report = new SPowerSolver(mesh).solve();

    assert.ok(report.ok, `did not converge: ${report.maxResidual}`);
    assert.equal(report.unenforced, 0);

    assert.ok(jointGap(a[0], b[0], v) < 1e-8, `tangent gap through the junction: ${jointGap(a[0], b[0], v)}`);

    // The arms' own interior joints are chain-private and still enforced, and C — which nobody
    // paired — is still free to leave at whatever angle its own energy prefers.
    assert.ok(jointGap(a[0], a[1], mids[0]) < 1e-8);
    assert.ok(jointGap(a[0], c[0], v) > 1e-3, "an unpaired end is not smoothed by accident");
  });

  it("reaches curvature continuity when the pairing asks for it", () => {
    const { mesh, v, a, b } = hub();

    v.setPairing(a[0], b[0], 2);

    const report = new SPowerSolver(mesh).solve();

    assert.ok(report.ok, `did not converge: ${report.maxResidual}`);
    assert.ok(jointGap(a[0], b[0], v) < 1e-8);
    assert.ok(curvatureGap(a[0], b[0], v) < 1e-8, `curvature jump: ${curvatureGap(a[0], b[0], v)}`);
  });

  /*
    Orientation, which §5 does not warn about. A vertex block is written along the direction its
    chain travels *through* the node, so an end running the other way mirrors the odd-order
    entries. "The other way" cannot be read off a pair — three arms all leaving the hub cannot
    pairwise be smooth continuations — so the signs follow the same spanning tree the rows do.
    Reading them against a single reference end instead converges just as happily and leaves a
    curvature jump of order 1 behind, which is why this measures G2 and not just `report.ok`.
  */
  it("orients a tangent group of three, not just of two", () => {
    const { mesh, v, a, b, c } = hub();

    v.setPairing(a[0], b[0], 2);
    v.setPairing(b[0], c[0], 2);

    const report = new SPowerSolver(mesh).solve();

    assert.ok(report.ok, `did not converge: ${report.maxResidual}`);
    assert.equal(components(paths(mesh), 1).length, 1, "one component of three chains");

    for (const [x, y, name] of [
      [a, b, "A–B"],
      [b, c, "B–C"],
    ] as const) {
      assert.ok(jointGap(x[0], y[0], v) < 1e-7, `${name} tangent gap: ${jointGap(x[0], y[0], v)}`);
      assert.ok(curvatureGap(x[0], y[0], v) < 1e-7, `${name} curvature jump: ${curvatureGap(x[0], y[0], v)}`);
    }
  });

  /*
    The coupling has to be *live*: an interface that assembled and then never reached the other
    chain would still converge, and would give each arm the answer it had before the pairing.
  */
  it("changes both arms, and not only the one whose row it is", () => {
    const alone = hub();
    const paired = hub();

    paired.v.setPairing(paired.a[0], paired.b[0], 2);

    new SPowerSolver(alone.mesh).solve();
    new SPowerSolver(paired.mesh).solve();

    for (const arm of ["a", "b"] as const) {
      const before = alone[arm][0].curve as SPowerClothoid;
      const after = paired[arm][0].curve as SPowerClothoid;

      assert.ok(Math.abs(before.ks[0] - after.ks[0]) > 1e-4, `arm ${arm} was left where it was`);
    }
  });

  /*
    The cut vertex is not a real one, so nothing about the solved loop may distinguish it. On a
    regular polygon symmetry makes that measurable exactly: every edge has to come out the same
    length with the same profile, and a mishandled Schur complement shows up as the two edges
    either side of the cut disagreeing with the rest.
  */
  it("leaves no trace of where a closed chain was cut", () => {
    const { mesh, verts, edges } = ring(8);
    const report = new SPowerSolver(mesh).solve();

    assert.ok(report.ok, `did not converge: ${report.maxResidual}`);
    assert.equal(report.cutClosed, 1);
    assert.equal(report.unenforced, 0);

    const curves = edges.map((e) => e.curve as SPowerClothoid);
    const spread = (xs: number[]) => Math.max(...xs) - Math.min(...xs);
    const lengths = curves.map((c) => c.length);

    assert.ok(spread(lengths) / lengths[0] < 1e-6, `arclength varies by ${spread(lengths)} around the ring`);

    // Absolutely, not relatively: the higher coefficients of a circle are zero, and a zero has
    // no scale of its own to be measured against.
    for (let k = 0; k < curves[0].ks.length; k++) {
      const ks = curves.map((c) => c.ks[k]);

      assert.ok(spread(ks) < 1e-6, `coefficient ${k} varies by ${spread(ks)} around the ring`);
    }

    // Including at the cut itself, which is the joint the band could not have carried.
    assert.ok(jointGap(edges[edges.length - 1], edges[0], verts[0]) < 1e-8);
    assert.ok(curvatureGap(edges[edges.length - 1], edges[0], verts[0]) < 1e-7);
  });
});
