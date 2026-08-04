import { CacheRing, Constraint, Solver, Vector2, clamp, fract } from "../math/index.js";
import { Curve, type CurveEdge } from "./curve.js";
import { emptySolveReport } from "./diagnostics.js";
import { type CurveSolver, type SolvableMesh, type SolvableVertex } from "./mesh_types.js";
import { activeProfile } from "./profile.js";
import * as nstructjs from "nstructjs";

/*
  A clothoid segment here is not defined by control points. It is defined by a sampled
  curvature profile k(s), and the geometry is recovered by integrating it:

      th(s) = integral of k(s) ds          tangent angle
      x(s)  = integral of cos th(s) ds     position
      y(s)  = integral of sin th(s) ds

  Because s is arclength, the unit tangent is exactly (cos th, sin th) and needs no
  normalization. That is the whole payoff of this parameterization, and it is why
  derivative() below is analytic and cheap while evaluate() is the expensive one.

  Classic Euler spirals are the special case k(s) = a*s. Here k is an arbitrary
  piecewise-linear function over `order` samples, so one segment can be a straight line,
  a circular arc, a true clothoid, or anything smoother.

  Original reduce-algebra derivation of the quadrature terms, kept verbatim:

      operator k, th, dk;
      operator isin, icos;

      x := icos(th(s));
      y := isin(th(s));

      forall s let df(icos(th(s)), s) = cos(th(s));
      forall s let df(isin(th(s)), s) = sin(th(s));

      let df(th(s), s) = k(s);
      let df(k(s), s)  = dk(s);
      let df(dk(s), s) = 0;

      dx2 := df(x, s, 2);   dy2 := df(y, s, 2);
      dx3 := df(x, s, 3);   dy3 := df(y, s, 3);
*/

/** Number of curvature samples per segment — the per-segment degrees of freedom. */
export const KORDER = 12;

/*
  The `ks` array does double duty: the first `order` slots are the curvature samples the
  solver is free to move, and these fixed tail slots cache the similarity transform
  derived from them. Keeping both in one Float64Array is what lets `_ks` be a view over
  just the free parameters, which is what the solver is handed.
*/
export const KSCALE = 16;
export const KTH = 17;
export const KOFFX = 18;
export const KOFFY = 19;
export const KARCSCALE = 20;
export const KTOT = 21;

const quadratureRets = new CacheRing(() => new Vector2(), 128);
const evalRets = new CacheRing(() => new Vector2(), 128);
const dvRets = new CacheRing(() => new Vector2(), 128);
const dv2Rets = new CacheRing(() => new Vector2(), 128);

const QUADRATURE_STEPS = 19;

/**
 * Integrate the position integrals from 0 to `s` in canonical space.
 *
 * Second-order Taylor expansion per step:
 *
 *     dx = cos th - k sin th (ds/2) - (cos th k^2 + dk sin th)(ds^2/6)
 *     dy = sin th + k cos th (ds/2) + (cos th dk - k^2 sin th)(ds^2/6)
 *
 * `s` runs over [-0.5, 0.5], so the shape is built centred on the origin and placed
 * afterwards by the similarity transform in {@link Clothoid._update}.
 */
function quadrature(ks: Float64Array, klen: number, s: number) {
  const ds = s / QUADRATURE_STEPS;
  const ret = quadratureRets.next().zero();

  let s2 = 0.0;
  let x = 0.0;
  let y = 0.0;

  for (let i = 0; i < QUADRATURE_STEPS; i++, s2 += ds) {
    const s3 = clamp(s2 + 0.5, 0.0, 1.0);

    const dk = activeProfile.dCurvature(ks, klen, s3);
    const k = activeProfile.curvature(ks, klen, s3);
    const th = activeProfile.integral(ks, klen, s3);

    const cos = Math.cos(th);
    const sin = Math.sin(th);

    x += cos - k * sin * ds * 0.5 - (cos * k * k + dk * sin) * ds * ds * (1.0 / 6.0);
    y += sin + k * cos * ds * 0.5 + (cos * dk - k * k * sin) * ds * ds * (1.0 / 6.0);
  }

  ret[0] = x * ds;
  ret[1] = y * ds;

  return ret;
}

/**
 * A curve segment defined by its curvature profile.
 *
 * The curvature samples control *shape only*. Endpoint interpolation is satisfied by a
 * similarity transform computed in {@link _update}, not by constraining the integral —
 * which is what removes two constraints per segment from the solver's job, and why every
 * cross-segment curvature comparison has to divide through by `KSCALE` first.
 */
export class Clothoid extends Curve {
  static override STRUCT = nstructjs.inlineRegister(
    this,
    `
Clothoid {
  ks    : array(float);
  order : int;
}
`
  );

  order = KORDER;

  ks = new Float64Array(KTOT);

  /** View over just the free curvature samples. This is what the solver is handed. */
  _ks: Float64Array;

  recalc = 1;

  v1!: SolvableVertex;
  v2!: SolvableVertex;

  constructor(v1?: SolvableVertex, v2?: SolvableVertex) {
    super();

    this.ks.fill(0);
    this._ks = new Float64Array(this.ks.buffer, 0, this.order);

    if (v1 !== undefined && v2 !== undefined) {
      this.v1 = v1;
      this.v2 = v2;
    }
  }

  override get length() {
    if (this.recalc) {
      this._update();
    }

    return this.ks[KSCALE];
  }

  override init(e: CurveEdge) {
    this.v1 = e.v1 as SolvableVertex;
    this.v2 = e.v2 as SolvableVertex;
  }

  override update(e?: CurveEdge) {
    if (e) {
      this.v1 = e.v1 as SolvableVertex;
      this.v2 = e.v2 as SolvableVertex;
    }

    this.recalc = 1;

    return this;
  }

  /**
   * Recompute the similarity transform mapping the canonical shape onto `v1 -> v2`.
   *
   * Integrating both endpoints gives the canonical chord; the offset, rotation and
   * uniform scale that carry that chord onto the real one are cached in the `ks` tail.
   */
  _update() {
    this.recalc = 0;

    const s = new Vector2(quadrature(this.ks, this.order, -0.5));
    const e = new Vector2(quadrature(this.ks, this.order, 0.5));

    const ks = this.ks;

    ks[KOFFX] = -s[0];
    ks[KOFFY] = -s[1];

    ks[KSCALE] = this.v1.vectorDistance(this.v2) / s.vectorDistance(e);
    ks[KARCSCALE] = 1.0 / ks[KSCALE];

    e.sub(s);

    const th1 = Math.atan2(this.v2[1] - this.v1[1], this.v2[0] - this.v1[0]);
    const th2 = Math.atan2(e[1], e[0]);

    ks[KTH] = th1 - th2;
  }

  override evaluate(s: number) {
    if (this.recalc) {
      this._update();
    }

    const ks = this.ks;

    s = clamp(s * ks[KARCSCALE], 0.0, 1.0) - 0.5;

    const p = evalRets.next().load(quadrature(ks, this.order, s));

    p[0] += ks[KOFFX];
    p[1] += ks[KOFFY];

    return p.rot2d(ks[KTH]).mulScalar(ks[KSCALE]).add(this.v1);
  }

  /**
   * Analytic tangent. Since `s` is arclength this is exactly `(cos th, sin th)` — no
   * normalization, no finite difference.
   */
  override derivative(s: number) {
    if (this.recalc) {
      this._update();
    }

    const ks = this.ks;

    s = clamp(s * ks[KARCSCALE], 0.0, 1.0);

    const th = activeProfile.integral(this.ks, this.order, s) + ks[KTH];

    const ret = dvRets.next();
    ret[0] = Math.cos(th);
    ret[1] = Math.sin(th);

    return ret;
  }

  /**
   * Analytic second derivative: `dth/ds` is just `k(s)`, so this is `k` times the tangent
   * rotated a quarter turn. The original finite-differenced {@link derivative} here,
   * which fed noise into {@link curvature}.
   */
  override derivative2(s: number) {
    if (this.recalc) {
      this._update();
    }

    const ks = this.ks;
    const sn = clamp(s * ks[KARCSCALE], 0.0, 1.0);

    const th = activeProfile.integral(this.ks, this.order, sn) + ks[KTH];
    const k = activeProfile.curvature(this.ks, this.order, sn) * ks[KARCSCALE];

    const ret = dv2Rets.next();
    ret[0] = -Math.sin(th) * k;
    ret[1] = Math.cos(th) * k;

    return ret;
  }

  /** Curvature is carried directly by the profile; no cross product needed. */
  override curvature(s: number) {
    if (this.recalc) {
      this._update();
    }

    const ks = this.ks;
    const sn = clamp(s * ks[KARCSCALE], 0.0, 1.0);

    return activeProfile.curvature(ks, this.order, sn) * ks[KARCSCALE];
  }

  afterSTRUCT(v1: SolvableVertex, v2: SolvableVertex) {
    this.v1 = v1;
    this.v2 = v2;
    this.recalc = 1;
  }

  loadSTRUCT(reader: (obj: this) => void) {
    reader(this);

    const ks = new Float64Array(KTOT);
    ks.set(this.ks.subarray(0, Math.min(this.ks.length, KTOT)));

    this.ks = ks;
    this._ks = new Float64Array(this.ks.buffer, 0, this.order);
    this.recalc = 1;
  }
}

/** Tuning for {@link ClothoidSolver}. */
export interface ClothoidSolverOptions {
  /**
   * Register the curvature-continuity constraint (G2) alongside tangent continuity (G1).
   *
   * Off by default, matching the original. See {@link curvatureConstraint} for why this
   * constraint is a projection rather than a residual.
   */
  enableG2: boolean;

  /**
   * Solve coarse-to-fine: converge at 2 curvature samples per edge, resample up, repeat
   * to {@link KORDER}. Off by default, matching the original, where the whole path sat
   * behind `if (0)`. A 12-DOF solve seeded from a converged 2-DOF solve should be more
   * robust than the cold start from 0.001, but this is untested.
   */
  progressiveRefinement: boolean;

  /** Corners sharper than this are excluded from the solve and left as corners. */
  cornerThreshold: number;

  iterations: number;
  relaxation: number;
}

export const defaultClothoidSolverOptions: ClothoidSolverOptions = {
  enableG2             : false,
  progressiveRefinement: false,
  cornerThreshold      : Math.PI * 0.4,
  iterations           : 55,
  relaxation           : 0.7,
};

interface JointParams {
  v: SolvableVertex;
  e1: Clothoid;
  e2: Clothoid;
  disabled: boolean;
}

/** Endpoint curvature of `e` at whichever end is shared, in canonical units. */
function endCurvature(e: Clothoid, isV1: boolean) {
  return isV1 ? e.ks[0] : e.ks[e.order - 1];
}

function setEndCurvature(e: Clothoid, isV1: boolean, k: number) {
  if (isV1) {
    e.ks[0] = k;
  } else {
    e.ks[e.order - 1] = k;
  }
}

/**
 * Tangent continuity at a shared vertex. A true residual: the solver finite-differences
 * it and descends.
 */
function tangentConstraint(params: JointParams) {
  const { v, e1, e2 } = params;

  e1.update();
  e2.update();

  const t1 = new Vector2(e1.derivative(v === e1.v1 ? 0 : e1.length)).normalize();
  const t2 = new Vector2(e2.derivative(v === e2.v1 ? 0 : e2.length)).normalize();

  if (!(v === e1.v1) === !(v === e2.v1)) {
    t2.negate();
  }

  return Math.acos(clamp(t1.dot(t2), -1.0, 1.0));
}

/**
 * Curvature continuity at a shared vertex — deliberately *not* a residual.
 *
 * This returns 0.0 unconditionally, which short-circuits the solver twice: `Constraint.
 * evaluate` skips finite-differencing any gradient, and `Solver.solveStep` skips touching
 * any parameter. The constraint is therefore invoked once per iteration purely for its
 * side effects, costing one call instead of 2 x order evaluations.
 *
 * The effect is a direct projection interleaved Gauss-Seidel style with the gradient
 * descent of {@link tangentConstraint}: each iteration the tangent constraint descends,
 * then this pulls the two endpoint curvatures halfway toward their mean. The 0.5
 * under-relaxation is what stops the projection fighting the descent. Registration order
 * therefore matters.
 *
 * Do not "fix" this into returning a residual — see `docs/research/clothoids.md` §7.
 */
function curvatureConstraint(params: JointParams) {
  const { v, e1, e2 } = params;

  if (params.disabled) {
    for (let i = 0; i < e1.order; i++) {
      e1.ks[i] *= 0.98;
      e2.ks[i] *= 0.98;
    }

    return 0.0;
  }

  const isV1e1 = v === e1.v1;
  const isV1e2 = v === e2.v1;

  const scale1 = e1.ks[KSCALE];
  const scale2 = e2.ks[KSCALE];

  /*
    Both edges parameterized away from `v` (or both toward it) means the through-path walks
    one of them backwards, which negates its signed curvature. Matches the rule
    `tangentConstraint` uses. See `tests/joint.test.ts`.
  */
  const flip = isV1e1 === isV1e2;

  let k1 = endCurvature(e1, isV1e1) / scale1;
  let k2 = endCurvature(e2, isV1e2) / scale2;

  if (flip) {
    k2 = -k2;
  }

  const k = (k1 + k2) * 0.5;
  const fac = 0.5;

  k1 += (k - k1) * fac;
  k2 += (k - k2) * fac;

  if (flip) {
    k2 = -k2;
  }

  setEndCurvature(e1, isV1e1, k1 * scale1);
  setEndCurvature(e2, isV1e2, k2 * scale2);

  e1.update();
  e2.update();

  return 0.0;
}

/** Resample every edge's curvature profile to `order` samples, preserving its shape. */
function changeOrder(mesh: SolvableMesh, order: number) {
  const temp = new Float64Array(KTOT);

  for (const e of mesh.edges) {
    const curve = e.curve as Clothoid;
    const ks = curve.ks;

    temp.set(ks.subarray(0, curve.order));

    for (let i = 0; i < order; i++) {
      let i1 = (i * (curve.order - 1)) / (order - 1);
      const t = fract(i1);

      i1 = ~~i1;

      const i2 = i1 + 1;

      ks[i] = i2 < order ? temp[i1] + (temp[i2] - temp[i1]) * t : temp[i1];
    }

    curve.order = order;
    curve._ks = new Float64Array(curve.ks.buffer, 0, order);
  }
}

/**
 * Fits every clothoid in a mesh so the segments agree at shared vertices.
 *
 * Gauss-Seidel: constraints are visited in order and each steps the parameters
 * immediately. See {@link ClothoidSolverOptions} for the parked experiments.
 *
 * The §8 report comes back empty. That is not an omission so much as the shape of the
 * thing: this iteration has no per-step residual to fit a rate to and no line search to
 * count backtracks in, and its one act of breaking — `cornerThreshold`, decided before the
 * solve runs — leaves for `Stroker` in Phase 6 rather than gaining a channel here.
 */
export class ClothoidSolver implements CurveSolver {
  options: ClothoidSolverOptions;

  constructor(
    public mesh: SolvableMesh,
    options: Partial<ClothoidSolverOptions> = {}
  ) {
    this.options = { ...defaultClothoidSolverOptions, ...options };
  }

  solve() {
    const mesh = this.mesh;
    const opts = this.options;

    const solver = new Solver();
    solver.threshold = 0.001;

    for (const e of mesh.edges) {
      const curve = e.curve as Clothoid;

      curve.ks.fill(0.001);
      curve.update(e);
    }

    const corners = new Set<SolvableVertex>();

    for (const v of mesh.verts) {
      if (v.edges.length !== 2) {
        continue;
      }

      const [edge1, edge2] = v.edges;
      const e1 = edge1.curve as Clothoid;
      const e2 = edge2.curve as Clothoid;

      const t1 = new Vector2(edge1.otherVertex(v)).sub(v).normalize();
      const t2 = new Vector2(edge2.otherVertex(v)).sub(v).normalize();

      if (Math.acos(clamp(t1.dot(t2) * 0.99999, -1.0, 1.0)) < opts.cornerThreshold) {
        corners.add(v);
        continue;
      }

      const params: JointParams = { v, e1, e2, disabled: false };
      const klst = [e1._ks, e2._ks];

      // Order matters: the projection runs after the descent within each iteration.
      solver.add(new Constraint("tan_c", tangentConstraint, klst, params, 1.0));

      if (opts.enableG2) {
        solver.add(new Constraint("curv_c", curvatureConstraint, klst, params, 1.0));
      }
    }

    if (opts.progressiveRefinement) {
      for (let order = 2; order <= KORDER; order++) {
        changeOrder(mesh, order);
        solver.solve(opts.iterations, opts.relaxation);
      }
    } else {
      solver.solve(opts.iterations, opts.relaxation);
    }

    for (const e of mesh.edges) {
      e.curve.update();
    }

    /*
      Corners are left as corners: zero the curvature at the sharp end of both edges. The
      original indexed with an undefined `mesh.order`, so this wrote to ks[NaN] — silently
      a no-op on a typed array — and only the v1 branch ever took effect.
    */
    for (const v of corners) {
      for (const edge of v.edges) {
        setEndCurvature(edge.curve as Clothoid, v === edge.v1, 0.0);
      }
    }

    return emptySolveReport();
  }
}
