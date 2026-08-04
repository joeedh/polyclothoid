/**
 * Phase 9 of `docs/plans/spower-solver.md` §12 — width as a second profile, and §10's three
 * differences from `κ`.
 *
 * The structural half is shared with everything Phases 2–8 already cover, so it is not retested
 * here: the same basis, the same blocks, the same bands, the same Schur complement. What is
 * tested is only what is *different*, and §10 names exactly three things.
 *
 * 1. **The transform.** Width is unsigned and carries no leading `L_e`, so its parity table and
 *    its scale are not `κ`'s. The consequence that matters downstream is that order 0 comes
 *    through a reversed end with sign `+1`, which is what lets a lower bound survive mirroring
 *    across a junction — so that is asserted directly rather than inferred from a solve.
 * 2. **The bound.** A QP is not a linear solve, and the property to check is not "close to the
 *    data" but "never under the floor, and off the floor wherever the data does not need it".
 *    The second half is the one an over-eager active set fails.
 * 3. **The channel.** §10 says a width break is independent of the level `κ` carries at the same
 *    joint. That is a claim about two partitions of the same vertex, and it is cheap to state
 *    exactly.
 *
 * The `h·|κ| < 1` check is a diagnostic and is tested as one: it fires where the offset really
 * does cusp and stays quiet where it does not, and it changes nothing about the solve.
 */
import assert from "node:assert/strict";
import { describe, it } from "node:test";

import {
  type EdgeFrame,
  SPOWER_ORDER,
  SPowerClothoid,
  SPowerSolver,
  WidthSolver,
  type WidthSample,
  chains,
  components,
  curvatureKind,
  cutOpen,
  evalSPower,
  interfaceSlots,
  sPowerLength,
  sPowerValueWeights,
  scalarEntry,
  sharing,
  vertexPartitions,
  widthChannel,
  widthCoefficients,
  widthEnergy,
  widthKind,
} from "../src/curve/index.js";
import { Mesh, type Edge, type Vertex } from "../src/mesh/index.js";

const P = SPOWER_ORDER;

function frame(forward: boolean, arclength = 1.0): EdgeFrame {
  return { chord: arclength, arclength, rEarlier: 1.0, rLater: 1.0, forward };
}

/** The canonical width profile of an edge whose chain-order DOF are `dof`. */
function profileOf(f: EdgeFrame, dof: number[]) {
  const a = widthCoefficients(P, f, Float64Array.from(dof));

  return (u: number) => evalSPower(a, a.length, u);
}

/** An open chain of `n` unit edges along `x`, with the s-power curve type installed. */
function strip(n: number, bend = 0.0) {
  const mesh = new Mesh();

  mesh.CurveCls = SPowerClothoid;
  mesh.SolverCls = SPowerSolver;

  const verts: Vertex[] = [];

  for (let i = 0; i <= n; i++) {
    verts.push(mesh.makeVertex([i, bend * i * (n - i)]));
  }

  const edges = verts.slice(1).map((v, i) => mesh.makeEdge(verts[i], v));

  return { mesh, verts, edges };
}

/** A regular `n`-gon on a circle of radius `r` — a closed chain of near-constant curvature `1/r`. */
function ring(n: number, r: number) {
  const mesh = new Mesh();

  mesh.CurveCls = SPowerClothoid;
  mesh.SolverCls = SPowerSolver;

  const verts: Vertex[] = [];

  for (let i = 0; i < n; i++) {
    const t = (i / n) * Math.PI * 2.0;

    verts.push(mesh.makeVertex([r * Math.cos(t), r * Math.sin(t)]));
  }

  const edges = verts.map((v, i) => mesh.makeEdge(v, verts[(i + 1) % n]));

  return { mesh, verts, edges };
}

/** A hub of three two-edge arms, as `junctions.test.ts` builds it. */
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

    return [mesh.makeEdge(v, verts[0]), mesh.makeEdge(verts[0], verts[1])];
  });

  return { mesh, v, a: arms[0], b: arms[1], c: arms[2] };
}

/** Uniform samples of a target width function over every edge of a mesh. */
function sample(edges: Edge[], per: number, target: (e: Edge, u: number) => number): WidthSample[] {
  const out: WidthSample[] = [];

  for (const edge of edges) {
    for (let i = 0; i < per; i++) {
      const u = i / (per - 1);

      out.push({ edge, u, width: target(edge, u) });
    }
  }

  return out;
}

/** Width at the two ends of an edge, in its own `v1 → v2` direction. */
function ends(e: Edge) {
  const c = e.curve as SPowerClothoid;

  return [c.width(0.0), c.width(c.length)];
}

describe("Phase 9: the width transform", () => {
  it("carries order 0 through unchanged whichever way the edge runs", () => {
    /*
      The load-bearing difference from `κ`. A width does not change sign when you walk the other
      way, so its reversal factor is `(−1)ⁿ` and order 0 is even — `w ≥ w_min` means the same
      thing at both ends of the mirroring, which is what makes the bound well posed on the
      interface. `κ` has the opposite parity here and flips.
    */
    for (const forward of [true, false]) {
      for (const earlier of [true, false]) {
        const f = frame(forward);

        assert.equal(scalarEntry(f, earlier, 0, widthKind).weight, 1.0, `width ${forward} ${earlier}`);
      }
    }

    assert.equal(scalarEntry(frame(false), true, 0, curvatureKind).weight, -1.0);
    assert.equal(scalarEntry(frame(true), true, 0, curvatureKind).weight, 1.0);
  });

  it("scales by no power of the arclength, where curvature scales by one", () => {
    const f = frame(true, 3.0);

    assert.equal(scalarEntry(f, true, 0, widthKind).weight, 1.0);
    assert.equal(scalarEntry(f, true, 0, curvatureKind).weight, 3.0);

    // Order `n` still picks up the `Rᵥ → L_e` rescale; that part is the same for both.
    assert.equal(scalarEntry(f, true, 1, widthKind).weight, 3.0);
    assert.equal(scalarEntry(f, true, 1, curvatureKind).weight, 9.0);
  });

  it("reproduces the DOF as the profile's endpoint values and slopes", () => {
    const dof = [2.0, 0.5, 3.0, -0.25];
    const w = profileOf(frame(true), dof);
    const h = 1e-6;

    assert.ok(Math.abs(w(0.0) - dof[0]) < 1e-12);
    assert.ok(Math.abs(w(1.0) - dof[2]) < 1e-12);

    assert.ok(Math.abs((w(h) - w(0.0)) / h - dof[1]) < 1e-5);
    assert.ok(Math.abs((w(1.0) - w(1.0 - h)) / h - dof[3]) < 1e-5);
  });

  it("gives the same world profile when the edge is built the other way round", () => {
    const dof = [2.0, 0.5, 3.0, -0.25];

    const forward = profileOf(frame(true), dof);
    const backward = profileOf(frame(false), dof);

    for (let i = 0; i <= 8; i++) {
      const u = i / 8;

      // Edge-local `u` runs against chain order on a reversed edge, and nothing else changes.
      assert.ok(Math.abs(forward(u) - backward(1.0 - u)) < 1e-12, `at ${u}`);
    }
  });

  it("reads the profile's value at a sample the same way the evaluator does", () => {
    const n = sPowerLength(P);
    const a = Float64Array.from([1.5, -0.5, 2.0, 0.25]);

    for (let i = 0; i <= 10; i++) {
      const u = i / 10;
      const row = sPowerValueWeights(P, u);

      let dot = 0.0;

      for (let k = 0; k < n; k++) {
        dot += row[k] * a[k];
      }

      assert.ok(Math.abs(dot - evalSPower(a, n, u)) < 1e-14, `at ${u}`);
    }
  });
});

describe("Phase 9: the width energy", () => {
  it("is the arclength integral of the slope, so refining the mesh does not change it", () => {
    /*
      Chord power 1, not `κ`'s 3: `q = L_e·κ` carries a leading `L_e` that `ŵ` does not, and §6
      measures per unit `u` where this measures per unit `s`. The test of that exponent is the
      one property it exists to give — a ramp cut into more pieces has the same energy.
    */
    const total = 4.0;
    const drop = 3.0;

    const energyOf = (pieces: number) => {
      const h = total / pieces;
      const e = widthEnergy(P, h, 0.0);
      const n = sPowerLength(P);

      let sum = 0.0;

      for (let i = 0; i < pieces; i++) {
        const w0 = (drop * i) / pieces;
        const w1 = (drop * (i + 1)) / pieces;
        /* The entry is `Rᵥⁿ·dⁿw/dsⁿ` and `Rᵥ = 1` here, so this is the world slope outright. */
        const slope = (w1 - w0) / h;

        const a = widthCoefficients(P, frame(true, h), Float64Array.from([w0, slope, w1, slope]));

        for (let r = 0; r < n; r++) {
          for (let c = 0; c < n; c++) {
            sum += a[r] * e[r * n + c] * a[c];
          }
        }
      }

      return sum;
    };

    const exact = (drop * drop) / total;

    for (const pieces of [1, 2, 4, 8]) {
      assert.ok(Math.abs(energyOf(pieces) - exact) < 1e-9, `${pieces}: ${energyOf(pieces)} vs ${exact}`);
    }
  });

  it("charges a constant profile only through the foundation term", () => {
    const n = sPowerLength(P);
    const chord = 2.0;
    const alpha = 0.25;
    const e = widthEnergy(P, chord, alpha);

    const a = new Float64Array(n);
    a[0] = 1.0;
    a[1] = 1.0;

    let sum = 0.0;

    for (let r = 0; r < n; r++) {
      for (let c = 0; c < n; c++) {
        sum += a[r] * e[r * n + c] * a[c];
      }
    }

    assert.ok(Math.abs(sum - (alpha * alpha) / chord) < 1e-12);
  });
});

describe("Phase 9: the width channel", () => {
  it("shares one more entry per level than curvature, there being no row to buy", () => {
    const { mesh, verts, edges } = strip(2);
    const v = verts[1];

    for (let level = 0; level <= 3; level++) {
      v.setPairing(edges[0], edges[1], level, "width");

      assert.equal(sharing(v, edges[0], edges[1], P, widthChannel).entries, Math.min(level, P + 1), `L${level}`);
    }

    assert.ok(mesh.edges.length === 2);
  });

  it("never groups the tangents, which have no width meaning", () => {
    const { verts, edges } = strip(2);
    const v = verts[1];

    v.setPairing(edges[0], edges[1], 99, "width");

    const parts = vertexPartitions(v, P, widthChannel);

    assert.notEqual(parts[0][0], parts[0][1]);
    assert.equal(parts[1][0], parts[1][1]);
  });

  it("breaks width at a joint without touching the curvature there, and the reverse", () => {
    const { verts, edges } = strip(2);
    const v = verts[1];

    const width = () => sharing(v, edges[0], edges[1], P, widthChannel).entries;
    const curvature = () => sharing(v, edges[0], edges[1], P);

    assert.equal(width(), P + 1);
    assert.equal(curvature().entries, P + 1);

    v.setPairing(edges[0], edges[1], 0, "width");

    assert.equal(width(), 0);
    assert.equal(curvature().entries, P + 1, "a width break is not a curvature break");
    assert.equal(curvature().tangent, true);

    v.setPairing(edges[0], edges[1], 99, "width");
    v.setPairing(edges[0], edges[1], 0, "curvature");

    assert.equal(width(), P + 1, "a curvature corner is not a width break");
    assert.equal(curvature().tangent, false);
  });

  it("puts no multiplier rows on the interface, width continuity being sharing alone", () => {
    const { mesh, v, a, b } = hub();

    v.setPairing(a[0], b[0], 99, "width");
    v.setPairing(a[0], b[0], 99, "curvature");

    const paths = chains(mesh).map(cutOpen);

    const width = interfaceSlots(paths, P, widthChannel);
    const curvature = interfaceSlots(paths, P);

    assert.ok(width.length > 0);
    assert.deepEqual([...new Set(width.map((s) => s.kind))], ["entry"]);
    assert.ok(
      curvature.some((s) => s.kind === "row"),
      "curvature does need a G1 row"
    );
  });
});

describe("Phase 9: the bounded fit", () => {
  it("recovers a constant width from constant samples", () => {
    const { mesh, edges } = strip(4);

    new SPowerSolver(mesh).solve();

    const report = new WidthSolver(
      mesh,
      sample(edges, 5, () => 2.0),
      { alpha: 1e-8 }
    ).solve();

    assert.equal(report.ok, true);
    assert.equal(report.components, 1);
    assert.ok(report.rms < 1e-6, `rms ${report.rms}`);

    for (const e of edges) {
      const c = e.curve as SPowerClothoid;

      for (let i = 0; i <= 8; i++) {
        assert.ok(Math.abs(c.width((i / 8) * c.length) - 2.0) < 1e-6);
      }
    }
  });

  it("tracks a ramp the samples describe", () => {
    const { mesh, edges } = strip(4);

    new SPowerSolver(mesh).solve();

    const at = (e: Edge, u: number) => 1.0 + (e.v1[0] + u * (e.v2[0] - e.v1[0])) * 0.5;
    const report = new WidthSolver(mesh, sample(edges, 5, at), { alpha: 1e-8, smoothing: 1e-4 }).solve();

    assert.ok(report.rms < 1e-4, `rms ${report.rms}`);
    assert.deepEqual(report.undershoot, []);

    const first = ends(edges[0]);
    const last = ends(edges[edges.length - 1]);

    assert.ok(Math.abs(first[0] - 1.0) < 1e-3);
    assert.ok(Math.abs(last[1] - 3.0) < 1e-3);
  });

  it("keeps width continuous across an interior joint, and lets a level-0 pairing break it", () => {
    const { mesh, verts, edges } = strip(2);

    new SPowerSolver(mesh).solve();

    const at = (e: Edge) => (e === edges[0] ? 1.0 : 4.0);
    const opts = { alpha: 1e-8, smoothing: 1e-3 };

    new WidthSolver(mesh, sample(edges, 6, at), opts).solve();

    const joined = Math.abs(ends(edges[0])[1] - ends(edges[1])[0]);
    assert.ok(joined < 1e-9, `shared entry, so identical: ${joined}`);

    verts[1].setPairing(edges[0], edges[1], 0, "width");

    const report = new WidthSolver(mesh, sample(edges, 6, at), opts).solve();

    assert.equal(report.ok, true);

    const split = ends(edges[0])[1] - ends(edges[1])[0];
    assert.ok(split < -2.0, `a hard width step, ${split}`);

    // §10's independence, measured rather than asserted: the geometry did not move.
    const c = edges[0].curve as SPowerClothoid;
    assert.ok(Math.abs(c.curvature(0.0)) < 1e-6);
  });

  it("clamps at the floor only where the data pushes under it", () => {
    const { mesh, edges } = strip(4);

    new SPowerSolver(mesh).solve();

    const at = (e: Edge, u: number) => 2.0 - (e.v1[0] + u * (e.v2[0] - e.v1[0]));
    const report = new WidthSolver(mesh, sample(edges, 5, at), {
      alpha    : 1e-8,
      smoothing: 1e-4,
      minimum  : 0.25,
    }).solve();

    assert.equal(report.ok, true);

    // Three vertices sit where the data is `≤ 0`, and the set has to settle on exactly those. An
    // active set that oscillates lands on the iteration cap instead, which is a silent failure.
    assert.equal(report.active, 3);
    assert.ok(report.iterations <= 3, `converged in ${report.iterations}`);

    for (const e of edges) {
      for (const w of ends(e)) {
        assert.ok(w >= 0.25 - 1e-9, `${w} under the floor`);
      }
    }

    // The other half: an active set that clamps everything also satisfies the bound. The wide
    // end is nowhere near the floor and must not have been dragged onto it.
    assert.ok(ends(edges[0])[0] > 1.9, `${ends(edges[0])[0]}`);

    // And the bound is a collocation at the vertices, so the interior is *not* covered — see the
    // module doc. It is reported instead, and this fixture is one that reports.
    assert.ok(report.undershoot.length > 0);
    assert.ok(report.undershoot.every((u) => u.width < 0.25));
  });

  it("leaves the bound alone when nothing violates it", () => {
    const { mesh, edges } = strip(3);

    new SPowerSolver(mesh).solve();

    const report = new WidthSolver(
      mesh,
      sample(edges, 5, () => 1.0),
      { alpha: 1e-8, minimum: 0.1 }
    ).solve();

    assert.equal(report.active, 0);
    assert.equal(report.iterations, 1, "one solve, nothing to add and nothing to release");
  });

  it("solves a chain that carries no samples at all rather than going singular", () => {
    const { mesh, edges } = strip(3);

    new SPowerSolver(mesh).solve();

    const report = new WidthSolver(mesh, [], { alpha: 0.5 }).solve();

    assert.equal(report.ok, true);
    assert.equal(report.samples, 0);
    assert.ok(Number.isFinite(report.residual));

    for (const e of edges) {
      for (const w of ends(e)) {
        assert.ok(Math.abs(w) < 1e-9, `the foundation pulls an unsampled profile to zero, got ${w}`);
      }
    }
  });
});

describe("Phase 9: width across a junction", () => {
  it("agrees at a node two chains reach from the same side", () => {
    const { mesh, v, a, b } = hub();

    v.setPairing(a[0], b[0], 99, "width");

    new SPowerSolver(mesh).solve();

    const edges = [...a, ...b];
    const at = (e: Edge) => (a.includes(e) ? 1.0 : 3.0);

    const report = new WidthSolver(mesh, sample(edges, 6, at), { alpha: 1e-8, smoothing: 1e-2 }).solve();

    assert.equal(report.ok, true);

    // Both arms leave `v` along their own `v1`, so the two ends mirror each other — and order 0
    // mirrors with `+1`, which is the whole point of the unsigned parity.
    const gap = Math.abs(ends(a[0])[0] - ends(b[0])[0]);
    assert.ok(gap < 1e-9, `shared order-0 entry, ${gap}`);

    // The third arm was not paired, so it is a component of its own and saw no samples.
    const found = components(chains(mesh).map(cutOpen), P, widthChannel);
    assert.equal(found.length, 2);
    assert.deepEqual(found.map((c) => c.chains.length).sort(), [1, 2]);
  });

  it("holds the floor through the interface, where the mirroring could have flipped it", () => {
    const { mesh, v, a, b } = hub();

    v.setPairing(a[0], b[0], 99, "width");

    new SPowerSolver(mesh).solve();

    const edges = [...a, ...b];
    const report = new WidthSolver(
      mesh,
      sample(edges, 6, () => -1.0),
      {
        alpha    : 1e-8,
        smoothing: 1e-3,
        minimum  : 0.5,
      }
    ).solve();

    assert.equal(report.ok, true);
    assert.ok(report.active > 0);

    for (const e of edges) {
      for (const w of ends(e)) {
        assert.ok(w >= 0.5 - 1e-9, `${w}`);
      }
    }
  });
});

describe("Phase 9: the h·|k| < 1 check", () => {
  it("finds the cusp where the offset turns inside out, and stays quiet where it does not", () => {
    const tight = ring(8, 1.0);
    new SPowerSolver(tight.mesh).solve();

    const wide = new WidthSolver(
      tight.mesh,
      sample(tight.edges, 5, () => 4.0),
      { alpha: 1e-8 }
    ).solve();

    assert.ok(wide.cusps.length > 0, "half-width 2 on a unit circle");
    assert.ok(wide.cusps.every((c) => c.product >= 1.0));
    assert.ok(wide.cusps.every((c) => Math.abs(c.product - Math.abs(c.curvature) * c.width * 0.5) < 1e-12));

    const thin = new WidthSolver(
      tight.mesh,
      sample(tight.edges, 5, () => 0.5),
      { alpha: 1e-8 }
    ).solve();

    assert.deepEqual(thin.cusps, []);
  });

  it("reports and does not act — the same width comes out either way", () => {
    const tight = ring(8, 1.0);
    new SPowerSolver(tight.mesh).solve();

    const samples = sample(tight.edges, 5, () => 4.0);
    const report = new WidthSolver(tight.mesh, samples, { alpha: 1e-8 }).solve();

    assert.ok(report.cusps.length > 0);

    for (const e of tight.edges) {
      for (const w of ends(e)) {
        assert.ok(Math.abs(w - 4.0) < 1e-4, `${w}: the cusp scan is a diagnostic, not a clamp`);
      }
    }
  });
});
