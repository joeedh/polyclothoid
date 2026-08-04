/**
 * Phase 2 of `docs/plans/spower-solver.md` §12 — vertex blocks and one Newton step.
 *
 * Four layers, tested separately because they fail differently:
 *
 * 1. {@link chains} — the node/chain decomposition of §5. Pure graph work, no arithmetic.
 * 2. {@link SPowerClothoid} — one segment against closed-form geometry. A constant s-power
 *    profile is a circular arc, and the chord-radius relation `C = 2R·sin(Φ/2)` pins the
 *    arclength, the curvature and the sagitta without reference to anything in the solver.
 * 3. {@link ChainSystem} — the assembled KKT matrix. This is the layer with no independent
 *    oracle in the solve itself: a transposed ordering and a flipped sign both present as
 *    slow convergence. So the assembled matrix is checked directly against quantities the
 *    curves compute for themselves.
 * 4. {@link SPowerSolver} — the whole loop. Every assertion about a joint is made in world
 *    space, from `derivative()`, so an error inside the solver's own residual cannot cancel
 *    itself out.
 *
 * Both Jacobians are exercised here. The exact one is the default and carries most of the
 * assertions; the frozen one keeps the two tests that pin what §5's wrong sparsity pattern
 * actually costs, so §14's claims stay measurable rather than remembered. The row itself is
 * differenced in `tests/spower_jacobian.test.ts`, which is the only oracle it has.
 */
import assert from "node:assert/strict";
import { describe, it } from "node:test";

import { Vector2 } from "../src/math/index.js";
import { SPowerClothoid, SPowerSolver, chains, edgeHessian } from "../src/curve/index.js";
import {
  ChainSystem,
  type SPowerSolverOptions,
  defaultSPowerSolverOptions,
  referenceLengths,
} from "../src/curve/spower_solver.js";
import { Mesh, type Edge, type Vertex } from "../src/mesh/index.js";

const TAU = Math.PI * 2.0;

function wrap(a: number) {
  return a - TAU * Math.round(a / TAU);
}

/** Deterministic uniform `[-1, 1)`. Tests never call `Math.random()`. */
function lcg(seed: number) {
  let state = seed >>> 0;

  return () => {
    state = (Math.imul(state, 1664525) + 1013904223) >>> 0;

    return (state / 0x100000000) * 2.0 - 1.0;
  };
}

/** A mesh whose edges carry s-power curves, built from a list of points. */
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

/** Direction of travel at `v`, walking *into* it along `e` if `arriving`, out of it if not. */
function travelAngle(e: Edge, v: Vertex, arriving: boolean) {
  const c = e.curve as SPowerClothoid;
  const atV2 = v === e.v2;

  const t = new Vector2(atV2 ? c.derivative(c.length) : c.derivative(0.0));

  // Arriving at v1, or leaving from v2, means walking the edge backwards.
  if (atV2 !== arriving) {
    t.negate();
  }

  return Math.atan2(t[1], t[0]);
}

/** World-space G1 gap at `v` between the edge arriving along `e1` and the one leaving along `e2`. */
function jointGap(e1: Edge, e2: Edge, v: Vertex) {
  return wrap(travelAngle(e2, v, false) - travelAngle(e1, v, true));
}

function worstGap(edges: Edge[], verts: Vertex[]) {
  let worst = 0.0;

  for (let i = 1; i < edges.length; i++) {
    worst = Math.max(worst, Math.abs(jointGap(edges[i - 1], edges[i], verts[i])));
  }

  return worst;
}

const ZIGZAG = [
  [0, 0],
  [1, 0.35],
  [2, 0.2],
  [3, 0.9],
  [4, 0.7],
];

/** One interior joint — as much as the frozen Jacobian can solve. See §14. */
const ELBOW = ZIGZAG.slice(0, 3);

describe("Phase 2: chain decomposition", () => {
  it("walks an open chain end to end", () => {
    const { mesh, verts, edges } = polyline(ZIGZAG);

    const found = chains(mesh);

    assert.equal(found.length, 1);
    assert.equal(found[0].closed, false);
    assert.deepEqual(found[0].verts, verts);
    assert.deepEqual(found[0].edges, edges);
  });

  it("walks the chain in chain order whatever the edges' own directions are", () => {
    const { mesh, verts } = polyline(ZIGZAG, [false, true, true, false]);

    const found = chains(mesh);

    assert.equal(found.length, 1);
    assert.deepEqual(found[0].verts, verts);
  });

  it("splits at a valence-3 node into one chain per incident edge", () => {
    const mesh = new Mesh();
    mesh.CurveCls = SPowerClothoid;

    const hub = mesh.makeVertex([0, 0]);
    const arms = [mesh.makeVertex([1, 0]), mesh.makeVertex([0, 1]), mesh.makeVertex([-1, 0])];

    for (const a of arms) {
      mesh.makeEdge(hub, a);
    }

    const found = chains(mesh);

    assert.equal(found.length, 3);

    for (const c of found) {
      assert.equal(c.edges.length, 1);
      assert.equal(c.verts[0], hub);
      assert.equal(c.closed, false);
    }
  });

  it("reports an all-valence-2 cycle as closed, with one vertex per edge", () => {
    const { mesh, verts, edges } = polyline([
      [0, 0],
      [1, 0],
      [1, 1],
      [0, 1],
    ]);

    mesh.makeEdge(verts[3], verts[0]);

    const found = chains(mesh);

    assert.equal(found.length, 1);
    assert.equal(found[0].closed, true);
    assert.equal(found[0].verts.length, 4);
    assert.equal(found[0].edges.length, 4);
    assert.equal(found[0].edges[0], edges[0]);
  });

  it("assigns every edge to exactly one chain, and runs a whole valence-2 hub through", () => {
    const mesh = new Mesh();
    mesh.CurveCls = SPowerClothoid;

    const hub = mesh.makeVertex([0, 0]);
    const arms = [
      [mesh.makeVertex([1, 0]), mesh.makeVertex([2, 0])],
      [mesh.makeVertex([0, 1]), mesh.makeVertex([0, 2])],
      [mesh.makeVertex([-1, 0]), mesh.makeVertex([-2, 0])],
    ];

    for (const [near, far] of arms) {
      mesh.makeEdge(hub, near);
      mesh.makeEdge(near, far);
    }

    const found = chains(mesh);
    const seen = new Set<Edge>();

    for (const chain of found) {
      assert.equal(chain.edges.length, 2);
      assert.equal(chain.verts[0], hub);

      for (const e of chain.edges) {
        assert.equal(seen.has(e as Edge), false);
        seen.add(e as Edge);
      }
    }

    assert.equal(found.length, 3);
    assert.equal(seen.size, 6);
  });
});

describe("Phase 2: SPowerClothoid", () => {
  it("is a straight line at zero curvature, with the endpoints hit exactly", () => {
    const { edges } = polyline([
      [0, 0],
      [3, 4],
    ]);

    const c = edges[0].curve as SPowerClothoid;

    assert.ok(Math.abs(c.length - 5.0) < 1e-12);

    const start = new Vector2(c.evaluate(0.0));
    const end = new Vector2(c.evaluate(c.length));

    assert.ok(start.vectorDistance(new Vector2([0, 0])) < 1e-15);
    assert.ok(end.vectorDistance(new Vector2([3, 4])) < 1e-12);

    for (let i = 0; i <= 8; i++) {
      assert.ok(Math.abs(c.curvature((i / 8) * c.length)) < 1e-15);
    }
  });

  /*
    A constant profile q = Φ is a circular arc of canonical radius 1/Φ and turning Φ, so
    C = 2R sin(Φ/2) and L = RΦ pin both outputs of the similarity transform. The tolerance
    is the 19-step quadrature's, not the algebra's — §14's Phase 0a table.
  */
  it("reproduces a circular arc from a constant profile", () => {
    const chord = 2.0;

    for (const phi of [0.3, 1.0, 2.5]) {
      const { edges } = polyline([
        [0, 0],
        [chord, 0],
      ]);

      const c = edges[0].curve as SPowerClothoid;

      c.ks.set([phi, phi, 0.0, 0.0]);
      c.recalc = 1;

      const radius = chord / (2.0 * Math.sin(phi * 0.5));

      assert.ok(Math.abs(c.turning - phi) < 1e-15, `turning at Φ=${phi}`);
      assert.ok(Math.abs(c.length - radius * phi) / (radius * phi) < 1e-5, `arclength at Φ=${phi}`);

      for (let i = 0; i <= 8; i++) {
        const k = c.curvature((i / 8) * c.length);

        assert.ok(Math.abs(k * radius - 1.0) < 1e-5, `curvature at Φ=${phi}, i=${i}`);
      }

      const end = new Vector2(c.evaluate(c.length));

      assert.ok(end.vectorDistance(new Vector2([chord, 0])) < 1e-5, `endpoint at Φ=${phi}`);

      /*
        The arc bulges away from its centre of curvature, so a left-turning arc sits to the
        right of a chord laid along +x. Getting this backwards is a very easy sign error.
      */
      const mid = new Vector2(c.evaluate(c.length * 0.5));
      const sagitta = radius * (1.0 - Math.cos(phi * 0.5));

      assert.ok(Math.abs(mid[0] - chord * 0.5) < 1e-5, `sagitta foot at Φ=${phi}`);
      assert.ok(Math.abs(mid[1] + sagitta) < 1e-5, `sagitta at Φ=${phi}`);
    }
  });

  it("mirrors exactly when the profile is negated", () => {
    const { edges } = polyline([
      [0, 0],
      [2, 0],
    ]);

    const c = edges[0].curve as SPowerClothoid;

    c.ks.set([1.0, 0.4, 0.3, 0.0]);
    c.recalc = 1;

    const up = [0, 1, 2, 3, 4].map((i) => new Vector2(c.evaluate((i / 4) * c.length)));
    const length = c.length;

    c.ks.set([-1.0, -0.4, -0.3, 0.0]);
    c.recalc = 1;

    assert.ok(Math.abs(c.length - length) < 1e-15);

    for (let i = 0; i <= 4; i++) {
      const down = new Vector2(c.evaluate((i / 4) * c.length));

      assert.ok(Math.abs(up[i][0] - down[0]) < 1e-14, `x at ${i}`);
      assert.ok(Math.abs(up[i][1] + down[1]) < 1e-14, `y at ${i}`);
    }
  });
});

/** The chain system for one polyline, with the DOF filled deterministically and `λ` left zero. */
function loadedSystem(
  points: number[][],
  flip: boolean[] = [],
  options: Partial<SPowerSolverOptions> = {},
  seed = 7919
) {
  const { mesh, verts, edges } = polyline(points, flip);
  const chain = chains(mesh)[0];

  const system = new ChainSystem(chain, referenceLengths(mesh), { ...defaultSPowerSolverOptions, ...options });
  const rand = lcg(seed);

  for (let i = 0; i < chain.verts.length; i++) {
    for (let k = 0; k < system.block; k++) {
      system.z[system.blockAt[i] + k] = rand();
    }
  }

  system.write();
  system.measure();
  system.assemble();

  return { mesh, verts, edges, system, product: system.kkt.apply(system.z) };
}

describe("Phase 2: KKT assembly", () => {
  /*
    With every multiplier zero, z·Kz is xᵀHx and nothing else, so this checks the scatter
    and the DOF ordering against the element matrices they came from.
  */
  it("scatters the energy Hessian to the same quadratic form the element matrices give", () => {
    const { system, product } = loadedSystem(ZIGZAG, [false, true, true, false]);

    const { p, block, frames, z } = system;
    const d = 2 * block;

    let assembled = 0.0;

    for (let i = 0; i < z.length; i++) {
      assembled += z[i] * product[i];
    }

    let expected = 0.0;

    for (let e = 0; e < frames.length; e++) {
      const h = edgeHessian(p, frames[e], defaultSPowerSolverOptions.alpha);

      for (let r = 0; r < d; r++) {
        for (let c = 0; c < d; c++) {
          expected += h[r * d + c] * z[system.at(e, r)] * z[system.at(e, c)];
        }
      }
    }

    assert.ok(expected > 1e-6, "the test would pass trivially on a zero form");
    assert.ok(Math.abs(assembled - expected) / expected < 1e-12, `${assembled} vs ${expected}`);
  });

  /*
    Each frozen G1 row is a signed sum of the turning integrals of whichever incident edges
    reach the joint at their own far end. `SPowerClothoid.turning` computes that integral in
    closed form from the coefficients, independently of `edgeTurningRow`, so this pins the
    row's signs, its scatter and the transform behind it in one comparison. The exact row
    has no such closed form to check against — `tests/spower_jacobian.test.ts` differences it.
  */
  it("builds each G1 row out of the turning integrals the curves report", () => {
    const { system, product } = loadedSystem(ZIGZAG, [false, true, true, false], { jacobian: "frozen" });

    const { frames, curves, lambdaAt } = system;
    let checked = 0;

    for (let i = 1; i < frames.length; i++) {
      if (lambdaAt[i] < 0) {
        continue;
      }

      const incoming = frames[i - 1].forward ? curves[i - 1].turning : 0.0;
      const outgoing = frames[i].forward ? 0.0 : curves[i].turning;
      const expected = outgoing - incoming;

      assert.ok(Math.abs(expected) > 1e-6, `row ${i} would pass trivially`);
      assert.ok(Math.abs(product[lambdaAt[i]] - expected) < 1e-12, `row ${i}: ${product[lambdaAt[i]]} vs ${expected}`);

      checked++;
    }

    assert.ok(checked >= 2, `only ${checked} rows exercised`);
  });

  it("keeps every assembled entry inside the declared band", () => {
    const { system } = loadedSystem(ZIGZAG);

    assert.ok(system.kkt.bandwidth < system.n);
    assert.ok(system.kkt.bandwidth >= 2 * system.block - 1);
  });
});

describe("Phase 2: SPowerSolver", () => {
  it("leaves a straight polyline straight", () => {
    const { mesh, edges } = polyline([
      [0, 0],
      [1, 0],
      [2, 0],
      [3, 0],
    ]);

    const report = new SPowerSolver(mesh).solve();

    assert.equal(report.ok, true);
    assert.equal(report.unenforced, 0);
    assert.ok(report.maxResidual < 1e-12);

    for (const e of edges) {
      const c = e.curve as SPowerClothoid;

      for (let i = 0; i <= 4; i++) {
        assert.ok(Math.abs(c.curvature((i / 4) * c.length)) < 1e-9);
      }
    }
  });

  it("reaches tangent continuity at a single interior joint", () => {
    const { mesh, verts, edges } = polyline(ELBOW);

    const report = new SPowerSolver(mesh).solve();

    assert.equal(report.ok, true);
    assert.equal(report.chains, 1);
    assert.equal(report.unenforced, 0);

    const gap = worstGap(edges, verts);

    assert.ok(gap < 1e-8, `world-space G1 gap ${gap}`);
    assert.ok(report.maxResidual < 1e-8, `reported residual ${report.maxResidual}`);
  });

  it("converges at a joint both edges reach at their own far ends", () => {
    const { mesh, verts, edges } = polyline(ELBOW, [false, true]);

    const report = new SPowerSolver(mesh).solve();

    assert.equal(report.unenforced, 0);
    assert.ok(worstGap(edges, verts) < 1e-8);
  });

  /*
    Phase 3's whole point, and the one thing the frozen Jacobian could not do at any step
    length: three interior joints coupled through the blocks they share. Eight steps against
    the frozen version's hundred-and-stall — §12 asked for the counts side by side.
  */
  it("couples neighbouring joints, which is what the frozen Jacobian could not", () => {
    const { mesh, verts, edges } = polyline(ZIGZAG);

    const report = new SPowerSolver(mesh).solve();

    assert.equal(report.ok, true);
    assert.equal(report.unenforced, 0);
    assert.ok(report.steps < 20, `${report.steps} steps`);

    const gap = worstGap(edges, verts);

    assert.ok(gap < 1e-8, `world-space G1 gap ${gap}`);
  });

  /*
    The row is built per edge from that edge's own frame, so nothing about it should depend
    on which way the edge was constructed. Identical step counts, not merely similar ones,
    is the assertion — the frozen version varies from 6 to 100 across the same three.
  */
  it("takes the same steps whichever way the edges were built", () => {
    const counts = new Set<number>();

    for (const flip of [[], [true, true, true, true], [false, true, true, false]]) {
      const { mesh, verts, edges } = polyline(ZIGZAG, flip);
      const report = new SPowerSolver(mesh).solve();

      assert.equal(report.unenforced, 0);
      assert.ok(worstGap(edges, verts) < 1e-8, `gap with flip ${JSON.stringify(flip)}`);

      counts.add(report.steps);
    }

    assert.equal(counts.size, 1, `step counts differed: ${[...counts].join(", ")}`);
  });

  /*
    Neither the energy nor the constraint depends on which way an edge was built, so with the
    exact Jacobian the two runs are the same problem and land on the same shape — 8e-7 here,
    which is the solve tolerance rather than an asymmetry. Frozen they differ by 5e-4,
    because each orientation stops at a fixed point of `Hz + Jᵀλ = 0` with a different wrong
    `J`. A sign error would show up as O(1) either way.
  */
  it("agrees with itself when the whole chain's edges are reversed", () => {
    const points = ZIGZAG;

    const forward = polyline(points);
    const reversed = polyline(points, [true, true, true, true]);

    const forwardReport = new SPowerSolver(forward.mesh).solve();
    const report = new SPowerSolver(reversed.mesh).solve();

    assert.equal(report.unenforced, 0);
    assert.ok(forwardReport.maxResidual < 1e-8, `forward residual ${forwardReport.maxResidual}`);
    assert.ok(report.maxResidual < 1e-8, `reversed residual ${report.maxResidual}`);

    let worst = 0.0;

    for (let i = 0; i < forward.edges.length; i++) {
      const a = forward.edges[i].curve as SPowerClothoid;
      const b = reversed.edges[i].curve as SPowerClothoid;

      for (let j = 0; j <= 8; j++) {
        const t = j / 8;
        const pa = new Vector2(a.evaluate(t * a.length));
        const pb = new Vector2(b.evaluate((1.0 - t) * b.length));

        worst = Math.max(worst, pa.vectorDistance(pb));
      }
    }

    assert.ok(worst < 1e-5, `orientation mismatch ${worst}`);
  });

  /*
    §5's wrong sparsity pattern at its sharpest. With KTH frozen an edge reaches a joint
    only through its turning integral, which lives at its own far end — so this one
    arrangement leaves the row identically zero. Phase 3 is what fixes it.
  */
  it("reports the joint it cannot constrain with the Jacobian frozen", () => {
    const { mesh } = polyline(
      [
        [0, 0],
        [1, 0.35],
        [2, 0.2],
      ],
      [true, false]
    );

    const report = new SPowerSolver(mesh, { jacobian: "frozen" }).solve();

    assert.equal(report.unenforced, 1);
    assert.equal(report.steps, 0);

    // The same joint with the exact Jacobian: both edges reach it through ∂KTH/∂x.
    const exact = new SPowerSolver(mesh).solve();

    assert.equal(exact.unenforced, 0);
    assert.ok(exact.maxResidual < 1e-8, `exact residual ${exact.maxResidual}`);
  });

  /*
    The bound on Phase 2, kept as the record of what Phase 3 moved. A second interior joint
    needs the coupling the frozen row does not have, and no step length rescues it — but it
    has to fail bounded, not at 1e+70. See §14.
  */
  it("stalls rather than diverging on a chain the frozen Jacobian cannot solve", () => {
    const { mesh, edges } = polyline(ZIGZAG);

    // `breaks: 0` because Phase 6 would otherwise rescue it by lowering a joint, which is the
    // right answer to a stall and the wrong thing to measure here.
    const report = new SPowerSolver(mesh, { jacobian: "frozen", breaks: 0 }).solve();

    assert.equal(report.ok, true);
    assert.equal(report.unenforced, 0);
    assert.ok(report.maxResidual > 1e-3, `unexpectedly converged to ${report.maxResidual}`);

    for (const e of edges) {
      const c = e.curve as SPowerClothoid;

      assert.ok(Number.isFinite(c.length), "arclength went non-finite");

      for (const k of c.ks) {
        assert.ok(Math.abs(k) < 1e2, `coefficient ran away to ${k}`);
      }
    }
  });

  it("skips a closed chain rather than solving a cyclic band", () => {
    const { mesh, verts } = polyline([
      [0, 0],
      [1, 0],
      [1, 1],
      [0, 1],
    ]);

    mesh.makeEdge(verts[3], verts[0]);

    const report = new SPowerSolver(mesh).solve();

    assert.equal(report.chains, 1);
    assert.equal(report.skippedClosed, 1);
    assert.equal(report.steps, 0);
  });

  it("solves each chain of a junction independently", () => {
    const mesh = new Mesh();

    mesh.CurveCls = SPowerClothoid;
    mesh.SolverCls = SPowerSolver;

    const hub = mesh.makeVertex([0, 0]);
    const a = mesh.makeVertex([1, 0.2]);
    const b = mesh.makeVertex([2, 0.6]);
    const c = mesh.makeVertex([-1, 0.2]);
    const d = mesh.makeVertex([-2, 0.6]);
    const e = mesh.makeVertex([0, -1]);

    const armA = [mesh.makeEdge(hub, a), mesh.makeEdge(a, b)];
    const armB = [mesh.makeEdge(hub, c), mesh.makeEdge(c, d)];

    mesh.makeEdge(hub, e);

    const report = new SPowerSolver(mesh).solve();

    assert.equal(report.chains, 3);
    assert.equal(report.unenforced, 0);
    assert.ok(Math.abs(jointGap(armA[0], armA[1], a)) < 1e-8);
    assert.ok(Math.abs(jointGap(armB[0], armB[1], c)) < 1e-8);
  });

  it("is invariant under a uniform zoom of the input", () => {
    const scale = 37.0;

    const small = polyline(ZIGZAG);
    const large = polyline(ZIGZAG.map(([x, y]) => [x * scale, y * scale]));

    new SPowerSolver(small.mesh).solve();
    new SPowerSolver(large.mesh).solve();

    for (let i = 0; i < small.edges.length; i++) {
      const a = small.edges[i].curve as SPowerClothoid;
      const b = large.edges[i].curve as SPowerClothoid;

      assert.ok(Math.abs(b.length / a.length - scale) / scale < 1e-9, `arclength on edge ${i}`);

      for (let k = 0; k < a.ks.length; k++) {
        assert.ok(Math.abs(a.ks[k] - b.ks[k]) < 1e-9, `coefficient ${k} on edge ${i}`);
      }
    }
  });

  it("rejects a non-positive alpha rather than silently losing definiteness", () => {
    const { mesh } = polyline([
      [0, 0],
      [1, 0],
    ]);

    assert.throws(() => new SPowerSolver(mesh, { alpha: 0.0 }), /alpha/);
  });
});
