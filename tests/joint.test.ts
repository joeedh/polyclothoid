/**
 * Phase 1 of `docs/plans/spower-solver.md` §12 — settling `curvatureConstraint`'s
 * orientation at runtime.
 *
 * A curve's signed curvature depends on which way you walk it, so a joint constraint has
 * to know whether the two incident edges are parameterized *through* the shared vertex or
 * *away from* it. `tangentConstraint` decides that with `isV1e1 === isV1e2`;
 * `curvatureConstraint` decides it with `isV1e1 !== isV1e2`. Exactly one of them can be
 * right, and `clothoids.md` §8 has carried the discrepancy as an open defect.
 *
 * The check here does not trust either. It builds the same three-point polyline four
 * times, once per assignment of edge directions, and asks two questions of the solved
 * result that are stated entirely in world space:
 *
 * 1. Do the two segments agree on the curvature at the joint, in a single consistent
 *    direction of travel?
 * 2. Does the solved *geometry* come out the same regardless of how the edges happen to be
 *    oriented? Edge direction is bookkeeping, not shape.
 *
 * The two tests are exact negations of each other, so whichever is wrong is wrong in all
 * four configurations rather than two. Question 1 answered it: `curvatureConstraint` is
 * the inverted one. Question 2 turned out to be nearly blind to the sign rule and is
 * documented on its own test.
 */
import assert from "node:assert/strict";
import { describe, it } from "node:test";

import { type Clothoid, ClothoidSolver } from "../src/curve/index.js";
import { Mesh, type Edge, type Vertex } from "../src/mesh/index.js";

const A = [0.0, 0.0];
const V = [1.0, 0.4];
const B = [2.0, 0.2];

/**
 * Curvature at `v`, signed for travel *away from* `v` along `e`.
 *
 * The edge's own parameterization runs `v1 -> v2`, so reading the far end means reading it
 * against the direction of travel, hence the negation.
 */
function outwardCurvature(e: Edge, v: Vertex) {
  const curve = e.curve as Clothoid;

  return v === e.v1 ? curve.curvature(0.0) : -curve.curvature(curve.length);
}

/**
 * G2 residual at the joint, in world units.
 *
 * Walking the chain enters along one edge and leaves along the other, so the two outward
 * curvatures are read in opposite directions of travel and must *sum* to zero.
 */
function jointG2(e1: Edge, e2: Edge, v: Vertex) {
  return outwardCurvature(e1, v) + outwardCurvature(e2, v);
}

/** Tangent residual at the joint, as an angle. Same sign reasoning as {@link jointG2}. */
function jointG1(e1: Edge, e2: Edge, v: Vertex) {
  const c1 = e1.curve as Clothoid;
  const c2 = e2.curve as Clothoid;

  const t1 = v === e1.v1 ? c1.derivative(0.0) : c1.derivative(c1.length);
  const t2 = v === e2.v1 ? c2.derivative(0.0) : c2.derivative(c2.length);

  const dot = (v === e1.v1) === (v === e2.v1) ? -t1.dot(t2) : t1.dot(t2);

  return Math.acos(Math.max(-1.0, Math.min(1.0, dot)));
}

const CONFIGS = [
  { name: "e1 = A->V, e2 = V->B", flip1: false, flip2: false },
  { name: "e1 = V->A, e2 = V->B", flip1: true, flip2: false },
  { name: "e1 = A->V, e2 = B->V", flip1: false, flip2: true },
  { name: "e1 = V->A, e2 = B->V", flip1: true, flip2: true },
];

function solveConfig(flip1: boolean, flip2: boolean) {
  const mesh = new Mesh();

  const a = mesh.makeVertex(A);
  const v = mesh.makeVertex(V);
  const b = mesh.makeVertex(B);

  const e1 = flip1 ? mesh.makeEdge(v, a) : mesh.makeEdge(a, v);
  const e2 = flip2 ? mesh.makeEdge(b, v) : mesh.makeEdge(v, b);

  new ClothoidSolver(mesh, { enableG2: true }).solve();

  return { mesh, v, e1, e2 };
}

/**
 * Sample the chain A -> V -> B in world space, always walking from `A` to `B` whatever the
 * edges' own directions are. Two configurations that describe the same curve must produce
 * the same samples.
 */
function samplePath(e1: Edge, e2: Edge, v: Vertex, count: number) {
  const out: number[][] = [];

  for (const e of [e1, e2]) {
    const curve = e.curve as Clothoid;

    // e1 is walked toward `v`, e2 away from it.
    const reversed = e === e1 ? v === e.v1 : v !== e.v1;

    for (let i = 0; i <= count; i++) {
      const t = i / count;
      const p = curve.evaluate((reversed ? 1.0 - t : t) * curve.length);

      out.push([p[0], p[1]]);
    }
  }

  return out;
}

function maxPathDistance(a: number[][], b: number[][]) {
  let worst = 0.0;

  for (let i = 0; i < a.length; i++) {
    worst = Math.max(worst, Math.hypot(a[i][0] - b[i][0], a[i][1] - b[i][1]));
  }

  return worst;
}

describe("Phase 1: joint orientation", () => {
  it("reaches G1 and G2 at the joint in all four edge orientations", (t) => {
    for (const { name, flip1, flip2 } of CONFIGS) {
      const { v, e1, e2 } = solveConfig(flip1, flip2);

      const g1 = jointG1(e1, e2, v);
      const g2 = jointG2(e1, e2, v);

      t.diagnostic(`  ${name.padEnd(22)} G1 ${g1.toExponential(2)} rad   G2 ${g2.toExponential(2)}`);

      assert.ok(g1 < 1e-3, `${name}: tangent residual ${g1} should be near zero`);
      assert.ok(Math.abs(g2) < 1e-3, `${name}: curvature residual ${g2} should be near zero`);
    }
  });

  /**
   * A property check on the solver, not on the sign rule — and worth being explicit that
   * it did *not* detect the sign bug: the worst deviation was 1.68e-4 with the wrong rule
   * and 1.65e-4 with the right one, because the G2 projection barely moves the shape at
   * all in this solver. What it does bound is orientation invariance, and only to about
   * the solver's own convergence level (~1e-4, next to a G1 residual of 1.25e-4). Exact
   * invariance is not a property of the Kaczmarz sweep: `derivative(0)` depends on the
   * curvature samples only through `KTH` while `derivative(length)` depends on all of
   * them, so reversing an edge genuinely changes the sensitivity structure the descent
   * sees. The §5 Gauss-Newton solver should tighten this.
   */
  it("solves the same geometry whichever way the edges point", (t) => {
    const paths = CONFIGS.map(({ flip1, flip2 }) => {
      const { v, e1, e2 } = solveConfig(flip1, flip2);

      return samplePath(e1, e2, v, 16);
    });

    let worst = 0.0;

    for (let i = 1; i < paths.length; i++) {
      const d = maxPathDistance(paths[0], paths[i]);
      worst = Math.max(worst, d);

      t.diagnostic(`  ${CONFIGS[i].name.padEnd(22)} deviates from the reference by ${d.toExponential(2)}`);
    }

    assert.ok(worst < 1e-3, `geometry should not depend on edge direction, worst deviation ${worst}`);
  });
});
