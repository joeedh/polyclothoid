import { CacheRing, Vector2, clamp } from "../math/index.js";
import { Curve, type Canvas2DLike, type CurveEdge } from "./curve.js";
import { type CurveSolver, type SolvableEdge, type SolvableMesh, type SolvableVertex } from "./mesh_types.js";
import { walk } from "./topology.js";
import * as nstructjs from "nstructjs";

export class BSplinePoint {
  static STRUCT = nstructjs.inlineRegister(
    this,
    `
BSplinePoint {
  co    : vec2;
  k     : float;
  flag  : int;
  index : int;
}
`
  );

  index = 0;
  flag = 0;
  co = new Vector2();

  /** Knot spacing contributed by this point. */
  k = 1.0;
}

/** Guards the zero-length knot spans that repeated end knots produce. */
function safeInv(n: number) {
  return n === 0 ? 100000.0 : 1.0 / n;
}

const evalRets = new CacheRing(() => new Vector2(), 512);
const dvRets = new CacheRing(() => new Vector2(), 512);

/** Cox-de Boor basis function of degree `n` for knot index `i`. */
function basis(ks: number[], s: number, i: number, n: number): number {
  const len = ks.length;

  const kn = clamp(i + 1, 0, len - 1);
  const knn = clamp(i + n, 0, len - 1);
  const knn1 = clamp(i + n + 1, 0, len - 1);
  const ki = clamp(i, 0, len - 1);

  if (n === 0) {
    return s >= ks[ki] && s < ks[kn] ? 1 : 0;
  }

  const a = (s - ks[ki]) * safeInv(ks[knn] - ks[ki] + 0.0001);
  const b = (ks[knn1] - s) * safeInv(ks[knn1] - ks[kn] + 0.0001);

  return a * basis(ks, s, i, n - 1) + b * basis(ks, s, i + 1, n - 1);
}

/**
 * Derivative of {@link basis}, from the standard identity
 *
 *     dB(i,n)/ds = n/(k[i+n] - k[i]) B(i,n-1) - n/(k[i+n+1] - k[i+1]) B(i+1,n-1)
 *
 * The original derived this in reduce as `df(a*basis(s,i,n-1) + b*basis(s,i+1,n-1), s)`
 * and left three competing implementations stacked behind the live `return`. Only this
 * one was reachable; the others are recorded in git history.
 */
export function dbasis(ks: number[], s: number, i: number, n: number) {
  const len = ks.length;

  const kn = clamp(i + 1, 0, len - 1);
  const knn = clamp(i + n, 0, len - 1);
  const knn1 = clamp(i + n + 1, 0, len - 1);
  const ki = clamp(i, 0, len - 1);

  const a = n * safeInv(ks[knn] - ks[ki]);
  const b = n * safeInv(ks[knn1] - ks[kn]);

  return a * basis(ks, s, i, n - 1) - b * basis(ks, s, i + 1, n - 1);
}

export const BSplineRecalc = {
  KNOTS : 1,
  TABLES: 2,
  FULL  : 1 | 2,
};

/**
 * Uniform B-spline segment.
 *
 * Unlike {@link Clothoid} and {@link CubicBezier}, this one is **not** arclength
 * parameterized — {@link length} returns 1.0 and `s` is the knot parameter. It is here as
 * a comparison curve for the curvature plots, not as a production stroke primitive.
 */
export class BSpline extends Curve {
  static override STRUCT = nstructjs.inlineRegister(
    this,
    `
BSpline {
  points : array(BSplinePoint);
  degree : int;
}
`
  );

  v1?: SolvableVertex;
  v2?: SolvableVertex;

  degree: number;
  points: BSplinePoint[] = [];
  knots: number[] = [];

  /** Multiplicity of the clamped end knots. */
  prefix: number;

  regen = BSplineRecalc.FULL;

  constructor(v1?: SolvableVertex, v2?: SolvableVertex, points = 4, degree = 3) {
    super();

    this.v1 = v1;
    this.v2 = v2;
    this.degree = degree;
    this.prefix = degree;

    for (let i = 0; i < points; i++) {
      this.points.push(new BSplinePoint());
    }

    if (v1 && v2) {
      this.points[0].co.load(v1);
      this.points[this.points.length - 1].co.load(v2);
    }
  }

  /** Not arclength — see the class doc. */
  override get length() {
    return 1.0;
  }

  /** Clamped uniform knot vector, normalized to [0, 1]. */
  regenKnots() {
    /* The original forced degree 2 here, overriding the constructor. Preserved: the
       basis code above was only ever exercised quadratic. */
    this.degree = 2;
    this.prefix = this.degree;

    this.regen &= ~BSplineRecalc.KNOTS;

    const knots = this.knots;
    knots.length = 0;

    for (let i = 0; i < this.prefix; i++) {
      knots.push(0.0);
    }

    let k = 0.0;
    let sumk = 0.0;

    for (const p of this.points) {
      knots.push(k);

      k += p.k;
      sumk += p.k;
    }

    for (let i = 0; i < this.prefix; i++) {
      knots.push(sumk);
    }

    const mulk = safeInv(sumk);
    for (let i = 0; i < knots.length; i++) {
      knots[i] *= mulk;
    }
  }

  override evaluate(s: number) {
    if (this.regen & BSplineRecalc.KNOTS) {
      this.regenKnots();
    }

    const ret = evalRets.next().zero();
    const knots = this.knots;
    const ps = this.points;

    for (let i = 0; i < knots.length; i++) {
      const p = ps[clamp(i, 0, ps.length - 1)];

      ret.addFac(p.co, basis(knots, s, i, this.degree));
    }

    return ret;
  }

  /**
   * Hodograph: the derivative of a degree-n spline is a degree-(n-1) spline over the
   * scaled forward differences of the control points.
   */
  override derivative(s: number) {
    if (this.regen & BSplineRecalc.KNOTS) {
      this.regenKnots();
    }

    const knots = this.knots;
    const ps = this.points;
    const degree = this.degree;

    const ret = dvRets.next().zero();
    const dv = dvRets.next().zero();

    for (let i = 1; i < knots.length - 1; i++) {
      const i1 = clamp(i - 1, 0, ps.length - 1);
      const i2 = clamp(i, 0, ps.length - 1);
      const ip1 = clamp(i + degree + 1, 0, ps.length - 1);

      const w = basis(knots, s, i, degree - 1);

      dv.load(ps[i2].co).sub(ps[i1].co);
      dv.mulScalar((degree - 1) * safeInv(knots[ip1] - knots[i1]));

      ret.addFac(dv, w);
    }

    return ret;
  }

  override derivative2(s: number) {
    const df = 0.0001;

    const a = new Vector2(this.derivative(s));
    const b = new Vector2(this.derivative(s + df));

    return b.sub(a).divScalar(df);
  }

  /**
   * Seed the control points from the edge and its two neighbours on each side.
   *
   * Interior points start on the chord, then the two points adjacent to the endpoints are
   * pulled onto averaged tangents, which is the same rule {@link BezierSolver} uses.
   */
  override init(e: CurveEdge) {
    const edge = e as SolvableEdge;
    const ps = this.points;

    if (this.regen & BSplineRecalc.KNOTS) {
      this.regenKnots();
    }

    this.v1 = edge.v1;
    this.v2 = edge.v2;

    ps[0].co.load(edge.v1);
    ps[ps.length - 1].co.load(edge.v2);

    for (let i = 1; i < ps.length - 1; i++) {
      ps[i].co.load(edge.v1).interp(edge.v2, i / (ps.length - 1));
    }

    const pv1 = walk(edge.v1, edge);
    const nv1 = walk(edge.v2, edge);

    const t1 = new Vector2();
    const t2 = new Vector2();

    t1.load(edge.v2).sub(edge.v1);
    t2.load(edge.v1).sub(pv1);
    t1.interp(t2, 0.5);

    ps[1].co.load(edge.v1).addFac(t1, 1.0 / 3.0);

    t1.load(nv1).sub(edge.v2);
    t2.load(edge.v2).sub(edge.v1);
    t1.interp(t2, 0.5);

    ps[ps.length - 2].co.load(edge.v2).addFac(t1, -1.0 / 3.0);
  }

  override update(e?: CurveEdge) {
    if (e) {
      this.v1 = (e as SolvableEdge).v1;
      this.v2 = (e as SolvableEdge).v2;
    }

    this.regen |= BSplineRecalc.FULL;

    return this;
  }

  /** Debug overlay: control points, sampled curve, and the basis functions themselves. */
  override draw(g: Canvas2DLike) {
    const w = 5;
    const ps = this.points;
    const knots = this.knots;

    g.beginPath();
    g.fillStyle = "rgba(255, 175, 55, 0.5)";
    for (const p of ps) {
      g.rect(p.co[0] - w / 2, p.co[1] - w / 2, w, w);
    }
    g.fill();

    const steps = 32;
    const ds = 1.0 / (steps - 1);

    g.fillStyle = "green";
    g.beginPath();
    for (let i = 0, s = 0.0; i < steps; i++, s += ds) {
      const p = this.evaluate(s);

      g.rect(p[0] - w * 0.5, p[1] - w * 0.5, w, w);
    }
    g.fill();

    const p1 = ps[0].co;
    const p2 = ps[ps.length - 1].co;

    const elen = p1.vectorDistance(p2);
    const th = Math.atan2(p2[1] - p1[1], p2[0] - p1[0]);
    const yScale = 0.2;

    const k0 = knots[0];
    const k1 = knots[knots.length - 1];

    /* The basis plot is drawn in curve space, then rotated onto the chord. */
    const co = new Vector2();
    const place = (x: number, y: number) =>
      co
        .loadXY(x * elen, y * elen * yScale)
        .rot2d(th)
        .add(p1);

    g.strokeStyle = "grey";
    g.beginPath();
    for (let knoti = 0; knoti < knots.length; knoti++) {
      for (let i = 0, s = 0.0; i < steps; i++, s += ds) {
        place(s, basis(knots, s * (k1 - k0) + k0, knoti, this.degree));

        if (i === 0) {
          g.moveTo(co[0], co[1]);
        } else {
          g.lineTo(co[0], co[1]);
        }
      }
    }
    g.stroke();

    /* Partition-of-unity check: this should be a flat line at y = 1. */
    g.beginPath();
    g.strokeStyle = "green";
    for (let i = 0, s = 0.0; i < steps; i++, s += ds) {
      let sum = 0.0;

      for (let j = 0; j < knots.length; j++) {
        sum += basis(knots, s, j, this.degree);
      }

      place(s, sum);

      if (i === 0) {
        g.moveTo(co[0], co[1]);
      } else {
        g.lineTo(co[0], co[1]);
      }
    }
    g.stroke();

    g.beginPath();
    g.strokeStyle = "rgba(0,0,0,0.5)";
    place(0.0, 1.0);
    g.moveTo(co[0], co[1]);
    place(1.0, 1.0);
    g.lineTo(co[0], co[1]);
    g.stroke();
  }

  afterSTRUCT(v1: SolvableVertex, v2: SolvableVertex) {
    this.v1 = v1;
    this.v2 = v2;
    this.regen = BSplineRecalc.FULL;

    return this;
  }

  loadSTRUCT(reader: (obj: this) => void) {
    reader(this);

    this.regen = BSplineRecalc.FULL;
  }
}

/** B-splines are placed entirely by {@link BSpline.init}; there is nothing to iterate. */
export class BSplineSolver implements CurveSolver {
  constructor(public mesh: SolvableMesh) {}

  solve() {
    for (const e of this.mesh.edges) {
      e.curve.init(e);
    }
  }
}
