import {Curve} from './curve.js';
import {Vector2, nstructjs, Matrix4, util, math} from './pathux.js';

export const BSplineFlags = {};

export class BSplinePoint {
  constructor() {
    this.index = 0;
    this.flag = 0;
    this.co = new Vector2();
    this.k = 1.0;
  }
}

BSplinePoint.STRUCT = `
BSplinePoint {
  co    : vec2;
  k     : float;
  flag  : int;
  index : int;
}
`;
nstructjs.register(BSplinePoint);

function safe_inv(n) {
  return n === 0 ? 100000.0 : 1.0/n;
}

let bs_evalrets = util.cachering.fromConstructor(Vector2, 512);
let bs_dvrets = util.cachering.fromConstructor(Vector2, 512);
let bs_dv2rets = util.cachering.fromConstructor(Vector2, 512);

function basis(ks, s, i, n) {
  let len = ks.length;
  let kp = Math.min(Math.max(i - 1, 0), len - 1);
  let kn = Math.min(Math.max(i + 1, 0), len - 1);
  let knn = Math.min(Math.max(i + n, 0), len - 1);
  let knn1 = Math.min(Math.max(i + n + 1, 0), len - 1);
  let ki = Math.min(Math.max(i, 0), len - 1);

  if (n === 0) {
    return s >= ks[ki] && s < ks[kn] ? 1 : 0;
  } else {

    let a = (s - ks[ki])*safe_inv(ks[knn] - ks[ki] + 0.0001);
    let b = (ks[knn1] - s)*safe_inv(ks[knn1] - ks[kn] + 0.0001);

    return a*basis(ks, s, i, n - 1) + b*basis(ks, s, i + 1, n - 1);
  }
}

function dbasis(ks, s, i, n) {
  let len = ks.length;
  let kp = Math.min(Math.max(i - 1, 0), len - 1);
  let kn = Math.min(Math.max(i + 1, 0), len - 1);
  let knn = Math.min(Math.max(i + n, 0), len - 1);
  let knn1 = Math.min(Math.max(i + n + 1, 0), len - 1);
  let ki = Math.min(Math.max(i, 0), len - 1);

  let a, b;
  a = n*safe_inv(ks[knn] - ks[ki]);
  b = n*safe_inv(ks[knn1] - ks[kn]);

  return a*basis(ks, s, i, n - 1) - b*basis(ks, s, i + 1, n - 1);

  if (n === 0) {
    return s >= ks[ki] && s < ks[kn] ? 1 : 0;
  } else {

    /*
    operator basis;

    a := (s - ki) / (knn - ki);
    b := (knn1 - s) / (knn1 - kn);

    ff := a*basis(s, i, n-1) + b*basis(s, i+1, n-1);

    df(ff, s);

     */


    a = n*safe_inv(knn - ki);
    b = n*safe_inv(knn1 - kn);

    return dbasis(ks, s, i, n - 1)*a - dbasis(ks, s, i + 1, n - 1)*b;

    a = (ks[ki] - s)*dbasis(ks, s, i, n - 1) - basis(ks, s, i, n - 1);
    a *= (ks[kn] - ks[knn1]);

    b = (ks[knn1] - s)*dbasis(ks, s, i + 1, n - 1) - basis(ks, s, i + 1, n - 1);
    b *= (ks[ki] - ks[knn]);

    return (a - b)*safe_inv((ks[ki] - ks[knn])*(ks[kn] - ks[knn1]));

    //let a = (s - ks[ki])*safe_inv(ks[knn] - ks[ki] + 0.0001);
    //let b = (ks[knn1] - s)*safe_inv(ks[knn1] - ks[kn] + 0.0001);


    //return a*basis(ks, s, i, n - 1) + b*basis(ks, s, i + 1, n - 1);
  }
}

const BSplineRecalc = {
  KNOTS : 1,
  TABLES: 2,
  FULL  : 1 | 2,
};

export class BSpline extends Curve {
  constructor(v1, v2, points = 4, degree = 3) {
    super();

    this.v1 = v1;
    this.v2 = v2;

    this.degree = degree;
    this.points = []; //extra points

    this.knots = [];
    this._points = [];

    for (let i = 0; i < points; i++) {
      this.points.push(new BSplinePoint());
    }

    this.points[0].co.load(v1);
    this.points[this.points.length - 1].co.load(v2);

    this.prefix = this.degree;

    this.regen = BSplineRecalc.FULL;

    this.arcTable = undefined;
    this.table = undefined;
  }

  get length() {
    return 1.0; //should be arc length!
  }

  derivative(s) {
    let knots = this.knots;
    let ps = this.points;
    let degree = this.degree;

    let ret = bs_dvrets.next().zero();
    let dv = bs_dvrets.next().zero();

    /* multiplicity decreases by one */
    for (let i = 1; i < knots.length - 1; i++) {
      let i0 = i - 2;
      let i1 = i - 1;
      let i2 = i;
      let i3 = i + 1;
      let i4 = i + 2;
      let ip1 = i + degree + 1;
      let ip2 = i + degree + 2;
      let ip0 = i + degree;

      i0 = Math.min(Math.max(i0, 0), ps.length - 1);
      i1 = Math.min(Math.max(i1, 0), ps.length - 1);
      i2 = Math.min(Math.max(i2, 0), ps.length - 1);
      i3 = Math.min(Math.max(i3, 0), ps.length - 1);
      i4 = Math.min(Math.max(i4, 0), ps.length - 1);
      ip0 = Math.min(Math.max(ip0, 0), ps.length - 1);
      ip1 = Math.min(Math.max(ip1, 0), ps.length - 1);
      ip2 = Math.min(Math.max(ip2, 0), ps.length - 1);

      let kp = Math.min(Math.max(i - 1, 0), knots.length - 1);
      let k = i;
      let kn = Math.min(Math.max(i + 1, 0), knots.length - 1);

      let w;
      //w = dbasis(knots, s, i2, degree-1);
      let di = 1.0/this.points.length;

      w = basis(knots, s, i, degree -1);

      dv.load(ps[i2].co).sub(ps[i1].co);
      dv.mulScalar((this.degree-1)*safe_inv((knots[ip1] - knots[i1])));
      //dv = ps[i3].co;
      //dv.load(ps[i2].co);

      ret.addFac(dv, w);
    }

    //ret.normalize();

    return ret;
  }

  regenKnots() {
    this.degree = 2;
    this.prefix = this.degree;

    this.regen &= ~BSplineRecalc.KNOTS;

    this.knots.length = 0;
    let knots = this.knots;

    let p = this.points[0];
    let k = -p.k;//*(this.prefix);

    k = 0.0;
    for (let i = 0; i < this.prefix; i++) {
      knots.push(k);
      //k += p.k;
    }

    k = 0.0;
    let sumk = 0.0;
    for (let p of this.points) {
      knots.push(k);

      k += p.k;
      sumk += p.k;
    }

    p = this.points[this.points.length - 1];

    let mulk = 1.0/(sumk);
    let lastk;

    //k += p.k;

    for (let i = 0; i < this.prefix; i++) {
      knots.push(sumk);

      //lastk = k;
      //k += p.k;
    }

    //mulk = 1.0 / (lastk);

    for (let i = 0; i < knots.length; i++) {
      knots[i] *= mulk;
    }
  }

  update() {
    this.flag |= BSplineRecalc.FULL;
    return this;
  }

  init(e) {
    let ps = this.points;

    if (this.regen & BSplineRecalc.KNOTS) {
      this.regenKnots();
    }

    ps[0].co.load(e.v1);
    ps[ps.length - 1].co.load(e.v2);

    let ewalk;

    function walk(v) {
      let e2 = v.otherEdge(ewalk);

      if (e2) {
        ewalk = e2;
        return ewalk.otherVertex(v);
      } else {
        return v;
      }
    }

    ewalk = e;
    let pv1 = walk(e.v1);
    let pv2 = walk(pv1);

    ewalk = e;
    let nv1 = walk(e.v2);
    let nv2 = walk(nv1);

    let elen = e.v1.vectorDistance(e.v2);

    for (let i = 1; i < ps.length - 1; i++) {
      let p = ps[i];
      let s = i/(ps.length - 1);

      p.co.load(e.v1).interp(e.v2, s);
    }

    let t1 = new Vector2();
    let t2 = new Vector2();

    t1.load(e.v2).sub(e.v1);
    t2.load(e.v1).sub(pv1);
    t1.interp(t2, 0.5);

    ps[1].co.load(e.v1).addFac(t1, 1.0/3.0);

    t1.load(nv1).sub(e.v2);
    t2.load(e.v2).sub(e.v1);
    t1.interp(t2, 0.5);

    ps[ps.length - 2].co.load(e.v2).addFac(t1, -1.0/3.0);

  }

  evaluate(s) {
    if (this.regen & BSplineRecalc.KNOTS) {
      this.regenKnots();
    }

    let ret = bs_evalrets.next().zero();
    let knots = this.knots;
    let ps = this.points;
    let prefix = this.prefix;
    let degree = this.degree;

    for (let i = 0; i < knots.length; i++) {
      let pi = Math.min(Math.max(i, 0), ps.length - 1);
      let p = ps[pi];

      let w = basis(knots, s, i, degree);

      ret.addFac(p.co, w);
    }

    return ret;
  }

  draw(g) {
    let w = 5;
    g.beginPath();
    g.fillStyle = "rgba(255, 175, 55, 0.5)";
    for (let p of this.points) {
      g.rect(p.co[0] - w/2, p.co[1] - w/2, w, w);
    }
    g.fill();

    g.fillStyle = "green";
    g.beginPath();
    let steps = 16;
    let s = 0.0, ds = 1.0/(steps - 1);
    for (let i = 0; i < steps; i++, s += ds) {
      let p = this.evaluate(s);

      g.rect(p[0] - w*0.5, p[1] - w*0.5, w, w);

    }
    g.fill();

    const yScale = 0.2;
    let ps = this.points;

    let elen = this.points[0].co.vectorDistance(this.points[this.points.length - 1].co);

    let p1 = this.points[0].co;
    let p2 = this.points[this.points.length - 1].co;
    let dx = p2[0] - p1[0];
    let dy = p2[1] - p1[1];

    let th = Math.atan2(dy, dx);
    let knots = this.knots;
    let prefix = this.prefix;


    function tok(s) {
      return s*(knots[knots.length - 1] - knots[0]) + knots[0];
    }

    function fromk(k) {
      return (k - knots[0])/(knots[knots.length - 1] - knots[0]);
    }

    let drawBasis = (knoti => {
      let steps = 32;
      let s = 0.0, ds = 1.0/(steps - 1);

      let co = new Vector2();

      for (let i = 0; i < steps; i++, s += ds) {
        let x = s*elen;
        let y = elen*yScale;

        let w = basis(this.knots, tok(s), knoti, this.degree);

        y *= w;

        co[0] = x;
        co[1] = y;
        co.rot2d(th);
        co.add(this.points[0].co);

        if (i === 0) {
          g.moveTo(co[0], co[1]);
        } else {
          g.lineTo(co[0], co[1]);
        }
      }

    });

    g.strokeStyle = 'grey'
    g.beginPath();
    let totk = this.knots.length;
    for (let i = 0; i < totk; i++) {
      drawBasis(i);
    }
    g.stroke();

    steps = 32;
    s = 0;
    ds = 1.0/(steps - 1);
    let co = new Vector2();

    //0/1 markers
    g.beginPath();
    co.loadXY(fromk(0.0)*elen, 0.0).rot2d(th).add(ps[0].co);
    g.moveTo(co[0], co[1]);
    co.loadXY(fromk(0.0)*elen, elen*yScale).rot2d(th).add(ps[0].co);
    g.lineTo(co[0], co[1]);
    co.loadXY(fromk(1.0)*elen, 0.0).rot2d(th).add(ps[0].co);
    g.moveTo(co[0], co[1]);
    co.loadXY(fromk(1.0)*elen, elen*yScale).rot2d(th).add(ps[0].co);
    g.lineTo(co[0], co[1]);
    g.stroke();

    g.beginPath();
    g.strokeStyle = "green";
    for (let i = 0; i < steps - 1; i++, s += ds) {
      let w = 0.0;

      let x = s*elen;
      let y = elen*yScale;

      for (let j = 0; j < knots.length; j++) {
        w += basis(knots, s, j, this.degree);
      }

      y *= w;
      co[0] = x;
      co[1] = y;
      co.rot2d(th);
      co.add(this.points[0].co);

      if (i === 0) {
        g.moveTo(co[0], co[1]);
      } else {
        g.lineTo(co[0], co[1]);
      }
    }

    g.stroke();

    g.beginPath();
    g.strokeStyle = "rgba(0,0,0,0.5)";
    s = 0.0;
    for (let i = 0; i < steps - 1; i++, s += ds) {
      let w = 0.0;

      let x = s*elen;
      let y = elen*yScale;

      co[0] = x;
      co[1] = y;
      co.rot2d(th);
      co.add(this.points[0].co);

      if (i === 0) {
        g.moveTo(co[0], co[1]);
      } else {
        g.lineTo(co[0], co[1]);
      }
    }
    g.stroke();
  }

  loadStruct(reader) {
    this.regen |= BSplineFlags.FULL;
    reader(this);
  }

  afterSTRUCT(v1, v2) {
    this.v1 = v1;
    this.v2 = v2;
    this.regen = BSplineRecalc.FULL;

    return this;
  }
}

BSpline.STRUCT = `
BSpline {
  points : array(BSplinePoint);
  degere : int;
}
`;
nstructjs.register(BSpline);

export class BSplineSolver {
  constructor(mesh) {
    this.mesh = mesh;
  }

  solve() {
    let mesh = this.mesh;

    for (let e of mesh.edges) {
      e.curve.init(e);
    }

  }
}