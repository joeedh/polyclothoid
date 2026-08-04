/**
 * Phase 3 of `docs/plans/spower-solver.md` §12 — the exact constraint Jacobian.
 *
 * §12: "`∂C(1)/∂a_j` in the quadrature loop, verified against finite differences on a single
 * edge before anything else." That ordering is the point of this file. Phase 2 shipped a row
 * that was wrong by a factor of three and structurally missing a whole block, and *nothing in
 * the solve said so* — it presented as slow convergence. Central differences are the only
 * oracle for a Jacobian that says otherwise.
 *
 * Three layers, narrowest first: the raw `∂C(1)/∂a_j` against a differenced endpoint, the
 * transform derivatives `∂KTH/∂a` and `∂L/∂a` built on top of it, and finally the assembled
 * constraint row against a differenced residual of the actual chain system.
 */
import assert from "node:assert/strict";
import { describe, it } from "node:test";

import { Vector2 } from "../src/math/index.js";
import {
  SPowerClothoid,
  SPowerSolver,
  chains,
  defaultSPowerQuadrature,
  integrateProfile,
  integrateProfileJacobian,
  sPowerLength,
  sPowerProfile,
  transformJacobian,
} from "../src/curve/index.js";
import { ChainSystem, defaultSPowerSolverOptions, referenceLengths } from "../src/curve/spower_solver.js";
import { Mesh, type Edge } from "../src/mesh/index.js";

const KLEN = sPowerLength(1);

/** Coefficient vectors chosen to span the interesting cases: flat, arc, front-loaded, S. */
const PROFILES: [string, number[]][] = [
  ["straight", [0.0, 0.0, 0.0, 0.0]],
  ["arc", [1.2, 1.2, 0.0, 0.0]],
  ["front-loaded", [2.0, 0.1, 0.6, 0.0]],
  ["s-curve", [1.5, -1.5, 0.4, -0.3]],
  ["strong", [2.6, 1.1, -0.9, 0.7]],
];

function endpoint(ks: Float64Array) {
  return new Vector2(integrateProfile(sPowerProfile, ks, KLEN, 0.0, 1.0, defaultSPowerQuadrature));
}

/** Central difference of `f` in coefficient `j`, restoring `ks` afterwards. */
function central(ks: Float64Array, j: number, h: number, f: (ks: Float64Array) => number) {
  const keep = ks[j];

  ks[j] = keep + h;
  const plus = f(ks);

  ks[j] = keep - h;
  const minus = f(ks);

  ks[j] = keep;

  return (plus - minus) / (2.0 * h);
}

describe("Phase 3: dC(1)/da", () => {
  /*
    `integrateProfileJacobian` differentiates the discrete sum, not the integral it
    approximates, so the agreement here is to roundoff at any step count. Quadraturing the
    continuous integrand instead leaves 3e-4 relative at nineteen steps — a tolerance chosen
    to accommodate that would also have accepted a genuinely wrong row.
  */
  it("matches a differenced endpoint on every profile", () => {
    const h = 1e-6;

    for (const [name, coefficients] of PROFILES) {
      const ks = new Float64Array(coefficients);
      const { end, d } = integrateProfileJacobian(sPowerProfile, ks, KLEN, defaultSPowerQuadrature);

      assert.ok(end.vectorDistance(endpoint(ks)) < 1e-12, `${name}: endpoint disagrees with integrateProfile`);

      for (let j = 0; j < KLEN; j++) {
        const dx = central(ks, j, h, (k) => endpoint(k)[0]);
        const dy = central(ks, j, h, (k) => endpoint(k)[1]);

        assert.ok(Math.abs(d[j * 2] - dx) < 1e-8, `${name}: dx/da${j} ${d[j * 2]} vs ${dx}`);
        assert.ok(Math.abs(d[j * 2 + 1] - dy) < 1e-8, `${name}: dy/da${j} ${d[j * 2 + 1]} vs ${dy}`);
      }
    }
  });

  it("is not trivially zero on a curved profile", () => {
    const ks = new Float64Array([1.2, 1.2, 0.0, 0.0]);
    const { d } = integrateProfileJacobian(sPowerProfile, ks, KLEN, defaultSPowerQuadrature);

    let peak = 0.0;

    for (const v of d) {
      peak = Math.max(peak, Math.abs(v));
    }

    assert.ok(peak > 0.05, `every entry near zero: ${peak}`);
  });
});

describe("Phase 3: transform derivatives", () => {
  /*
    §5 derives KSCALE and KTH from C(1) as `chord/|C|` and `−arg C`. Differencing them
    through the same closed form is a weak check; differencing them through
    `SPowerClothoid._update`, which is what the solver actually measures, is not.
  */
  it("matches KSCALE and KTH differenced through the curve itself", () => {
    const chord = 1.7;
    const h = 1e-6;

    const mesh = new Mesh();

    mesh.CurveCls = SPowerClothoid;

    const v1 = mesh.makeVertex([0, 0]);
    const v2 = mesh.makeVertex([chord, 0]);
    const curve = mesh.makeEdge(v1, v2).curve as SPowerClothoid;

    const scaleOf = (ks: Float64Array) => {
      curve.ks.set(ks);
      curve.recalc = 1;

      return curve.length;
    };

    // Reading `length` is what forces `_update`, so `th` is only current once it has been read.
    const kthOf = (ks: Float64Array) => {
      scaleOf(ks);

      return curve.th;
    };

    for (const [name, coefficients] of PROFILES) {
      const ks = new Float64Array(coefficients);
      const j = transformJacobian(sPowerProfile, ks, KLEN, chord, defaultSPowerQuadrature);

      assert.ok(Math.abs(j.scale - scaleOf(ks)) < 1e-12, `${name}: scale`);

      for (let c = 0; c < KLEN; c++) {
        const dScale = central(ks, c, h, scaleOf);
        const dKth = central(ks, c, h, kthOf);

        assert.ok(Math.abs(j.dScale[c] - dScale) < 1e-8, `${name}: dL/da${c} ${j.dScale[c]} vs ${dScale}`);
        assert.ok(Math.abs(j.dKth[c] - dKth) < 1e-8, `${name}: dKTH/da${c} ${j.dKth[c]} vs ${dKth}`);
      }
    }
  });
});

const ZIGZAG = [
  [0, 0],
  [1, 0.35],
  [2, 0.2],
  [3, 0.9],
  [4, 0.7],
];

function system(points: number[][], flip: boolean[] = [], seed = 4093) {
  const mesh = new Mesh();

  mesh.CurveCls = SPowerClothoid;
  mesh.SolverCls = SPowerSolver;

  const verts = points.map((p) => mesh.makeVertex(p));
  const edges: Edge[] = [];

  for (let i = 0; i + 1 < verts.length; i++) {
    edges.push(flip[i] ? mesh.makeEdge(verts[i + 1], verts[i]) : mesh.makeEdge(verts[i], verts[i + 1]));
  }

  const chain = chains(mesh)[0];
  const out = new ChainSystem(chain, referenceLengths(mesh), { ...defaultSPowerSolverOptions });

  // A nonzero state, so the row is differentiated somewhere the curvature actually varies.
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

describe("Phase 3: the assembled constraint row", () => {
  /*
    The whole chain, not one edge: this is where the pullback through `M_e`, the sign
    convention against `measure`, and the scatter into the band all get their only check.
    `L_e` is held fixed on both sides — it is an outer-iteration quantity, so the row is the
    exact derivative *at frozen arclength*, and differencing through `write` holds the same
    thing fixed.
  */
  for (const [name, flip] of [
    ["all forward", []],
    ["all reversed", [true, true, true, true]],
    ["mixed", [false, true, true, false]],
    ["alternating", [true, false, true, false]],
  ] as [string, boolean[]][]) {
    it(`matches a differenced residual, ${name}`, () => {
      const { system: s } = system(ZIGZAG, flip);
      const h = 1e-6;

      s.assemble();

      let checked = 0;
      let peak = 0.0;

      for (let i = 1; i < s.frames.length; i++) {
        const at = s.lambdaAt[i];

        assert.ok(at >= 0, `joint ${i} has no row with the exact Jacobian`);

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

          assert.ok(Math.abs(row - fd) < 1e-8, `joint ${i}, dof ${c}: ${row} vs ${fd}`);
          checked++;
        }
      }

      s.write();
      s.measure();

      assert.ok(checked >= 30, `only ${checked} entries compared`);
      assert.ok(peak > 0.1, `the whole Jacobian is near zero, peak ${peak}`);
    });
  }

  /*
    Phase 2's row, kept as the record of what it got wrong. The failure is not a scale
    factor — three of the four blocks in a joint's stencil are structurally absent — so no
    amount of damping recovers it. See §14.
  */
  it("is what the frozen row is not", () => {
    const { system: s } = system(ZIGZAG);
    const frozen = system(ZIGZAG);

    frozen.system.options.jacobian = "frozen";

    s.assemble();
    frozen.system.assemble();

    let worst = 0.0;

    for (let i = 1; i < s.frames.length; i++) {
      for (let c = 0; c < s.n; c++) {
        if (s.rows.includes(c)) {
          continue;
        }

        worst = Math.max(
          worst,
          Math.abs(s.kkt.get(s.lambdaAt[i], c) - frozen.system.kkt.get(frozen.system.lambdaAt[i], c))
        );
      }
    }

    assert.ok(worst > 0.1, `the two Jacobians agree to ${worst}, so one of them is not what it claims`);
  });
});
