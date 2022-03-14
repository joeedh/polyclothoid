import {Curve} from './curve.js';
import {util, Vector2, nstructjs, math, Constraint, Solver, binomial} from './pathux.js';

export const KORDER = 12;
export const KSCALE = 16;
export const KTH = 17;
export const KOFFX = 18;
export const KOFFY = 19;
export const KARCSCALE = 20;
export const KTOT = 21;

function step(ks, klen, s) {
  //return bstep(ks, klen, s);

  let eps = 0.00001;
  //s = eps + s * (1.0 - eps*2.0);

  let i1 = s*(klen - 1);
  //let i1 = s*klen*0.9999;
  let t = Math.fract(i1);
  i1 = ~~(i1 + 0.00001);

  let i2 = i1 + 1;

  if (i2 < klen - 1) {
    return ks[i1] + (ks[i2] - ks[i1])*t;
  } else {
    return ks[i1];
  }
}

function bernstein(i, n, s) {
  return binomial(n, i)*Math.pow(s, i)*Math.pow(1.0 - s, n - i);
}

function dstep(ks, klen, s) {
  let df = 0.00001;

  return (step(ks, klen, s + df) - step(ks, klen, s))/df;

  return ks[~~(s*(klen - 1) + 0.00001)];
  //return ks[~~(s*(klen*0.9999))];
}

function bstep(ks, klen, s) {
  let sum = 0.0;
  for (let i = 0; i < klen; i++) {
    sum += bernstein(i, klen - 1, s);
  }

  return sum;
}

function imix(a, b, s) {
  return -((s - 2.0)*a - b*s)*s*0.5;
}

function istep2(ks, klen, s) {
  let eps = 0.00001;

  let klen2 = klen - 1;
  //s = eps + s * (1.0 - eps*2.0);

  let i1 = s*(klen - 1);
  //let i1 = s*klen*0.9999;
  let t = Math.fract(i1);
  i1 = ~~(i1 + 0.00001);
  let i2 = i1 + 1;

  let sum = 0.0;
  for (let i = 0; i < i1; i++) {
    sum += imix(ks[i], ks[i + 1], 1.0)/klen2;
  }

  i2 = Math.min(Math.max(i2, 0), klen - 1);
  if (i2 !== i1) {
    sum += imix(ks[i1], ks[i2], t)/klen2;
  }

  return sum;

  i2 = i1 + 1;

  if (i2 < klen - 1) {
    return ks[i1] + (ks[i2] - ks[i1])*t;
  } else {
    return ks[i1];
  }
}

function istep(ks, klen, s) {
  return istep2(ks, klen, s);
  let steps = 32;
  let ds = s/steps;
  let s2 = 0.0;

  let ret = 0.0;

  for (let i = 0; i < steps; i++, s2 += ds) {
    ret += step(ks, klen, s2)*ds + dstep(ks, klen, s2)*ds*ds*0.5;
  }

  return ret;
}


let piecewise_linear = [dstep, step, istep];

let circle_arc = [
  function(ks, klen, s) { //derivative of curvature
    return 0.0;
  },
  function(ks, klen, s) { //curvature
    let si = s * (klen-1);

    return ks[~~si];
  },
  function(ks, klen, s) { //integral of curvature
    let si = s * (klen-1);
    let t = Math.fract(si);
    si = ~~si;

    let sum = 0.0;

    let ds = 1.0 / klen;

    for (let i=0; i<si; i++) {
      sum += ks[i]*ds;
    }

    sum += t * ks[si]*ds;

    return sum;
  }
]

//let funcs = circle_arc;
let funcs = piecewise_linear;

/*

operator k, th, dk;
operator isin, icos;

x := icos(th(s));
y := isin(th(s));

forall s let df(icos(th(s)), s) = cos(th(s));
forall s let df(isin(th(s)), s) = sin(th(s));

let df(th(s), s) = k(s);
let df(k(s), s) = dk(s);
let df(dk(s), s) = 0;

dx2 := df(x, s, 2);
dy2 := df(y, s, 2);

dx3 := df(x, s, 3);
dy3 := df(y, s, 3);

dx := cos(th(s))*ds*0.5 + dx2*ds*ds*(1.0/6.0) + dx3*ds*ds*ds*(1.0/24.0);
dy := sin(th(s))*ds*0.5 + dy2*ds*ds*(1.0/6.0) + dy3*ds*ds*ds*(1.0/24.0);

dx/ds;
dy/ds;

*/

let rets = util.cachering.fromConstructor(Vector2, 128);
let evalrets = util.cachering.fromConstructor(Vector2, 128);
let dvrets = util.cachering.fromConstructor(Vector2, 128);
let dv2rets = util.cachering.fromConstructor(Vector2, 128);

function quadrature(ks, klen, s) {
  let steps = 19;
  let s2 = 0.0, ds = s/steps;

  let ret = rets.next().zero();
  let x = 0.0;
  let y = 0.0;

  for (let i = 0; i < steps; i++, s2 += ds) {
    let s3 = s2 + 0.5;
    s3 = Math.min(Math.max(s3, 0.0), 1.0);

    let dk = funcs[0](ks, klen, s3);
    let k = funcs[1](ks, klen, s3);
    let th = funcs[2](ks, klen, s3);

    let cos = Math.cos(th);
    let sin = Math.sin(th);

    let dx = -cos*k*k*ds*ds + 12*cos - dk*sin*ds*ds - 4*k*sin*ds;
    dx *= 0.25;

    let dy = cos*dk*ds*ds + 4*cos*k*ds - k*k*sin*ds*ds + 12*sin;
    dy /= 24.0;

    dx = cos - k*sin*ds*0.5 - (cos*k*k + dk*sin)*ds*ds*(1.0/6.0);
    dy = sin + k*cos*ds*0.5 + (cos*dk - k*k*sin)*ds*ds*(1.0/6.0);

    x += dx;
    y += dy;
  }

  ret[0] = x*ds;
  ret[1] = y*ds;

  return ret;
}

export class Clothoid extends Curve {
  constructor(v1, v2) {
    super();

    this.order = KORDER;

    this.ks = new Float64Array(KTOT);
    this.ks.fill(0);
    this._ks = new Float64Array(this.ks.buffer, 0, this.order);

    this.v1 = v1;
    this.v2 = v2;

    this.recalc = 1;
  }

  get length() {
    if (this.recalc) {
      this._update();
    }

    return this.ks[KSCALE];
  }

  init(e) {
    this.v1 = e.v1;
    this.v2 = e.v2;
  }

  update(e) {
    if (e) {
      this.v1 = e.v1;
      this.v2 = e.v2;
    }

    this.recalc = 1;
    return this;
  }

  _update() {
    this.recalc = 0;

    let s = quadrature(this.ks, this.order, -0.5);
    let e = quadrature(this.ks, this.order, 0.5);

    let ks = this.ks;

    ks[KOFFX] = -s[0];
    ks[KOFFY] = -s[1];

    ks[KSCALE] = this.v1.vectorDistance(this.v2)/s.vectorDistance(e);
    ks[KARCSCALE] = 1.0/ks[KSCALE];
    e.sub(s);

    let th1 = Math.atan2(this.v2[1] - this.v1[1], this.v2[0] - this.v1[0]);
    let th2 = Math.atan2(e[1], e[0]);
    ks[KTH] = th1 - th2;
  }

  evaluate(s) {
    if (this.recalc) {
      this._update();
    }

    let ks = this.ks;
    s *= ks[KARCSCALE];

    s = Math.min(Math.max(s, 0.0), 1.0);
    s -= 0.5;

    let p = quadrature(ks, this.order, s);

    p[0] += ks[KOFFX];
    p[1] += ks[KOFFY];

    p.rot2d(ks[KTH]).mulScalar(ks[KSCALE]).add(this.v1);

    return p;
  }

  derivative(s) {
    if (this.recalc) {
      this._update();
    }

    if (0) {
      let df = 0.0001;
      let a = this.evaluate(s-df);
      let b = this.evaluate(s+df);
      b.sub(a).divScalar(2.0*df);

      return dvrets.next().load(b);
    }

    let ks = this.ks;
    s *= ks[KARCSCALE];

    s = Math.min(Math.max(s, 0.0), 1.0);

    let th = funcs[2](this.ks, this.order, s);
    th += ks[KTH];

    let ret = dvrets.next();
    ret[0] = Math.cos(th);
    ret[1] = Math.sin(th);

    return ret;
  }

  derivative2(s) {
    let df = 0.0001;
    let a = this.derivative(s);
    let b = this.derivative(s + df);
    b.sub(a).mulScalar(1.0/df);

    return dv2rets.next().load(b);
  }

  curvature(s) {
    if (0) {
      let ks = this.ks;
      s *= ks[KARCSCALE];

      s = Math.min(Math.max(s, 0.0), 1.0);

      return funcs[1](ks, this.order, s)*ks[KARCSCALE];
    }

    let dv1 = this.derivative(s);
    let dv2 = this.derivative2(s);

    return (dv1[0]*dv2[1] - dv1[1]*dv2[0])/Math.pow(dv1.dot(dv1), 3.0/2.0);
  }

  afterSTRUCT(v1, v2) {
    this.v1 = v1;
    this.v2 = v2;
    this.regen = 1;
  }

  loadSTRUCT(reader) {
    reader(this);

    while (this.ks.length < KTOT) {
      this.ks.push(0.0);
    }

    this.ks = new Float64Array(this.ks);
    this._ks = new Float64Array(this.ks.buffer, 0, this.order);
  }

  normal() {
    return dv2rets.next().zero();
  }

  draw(g) {

  }
}

Clothoid.STRUCT = `
Clothoid {
  ks   : array(float); 
}
`;
nstructjs.register(Clothoid);

export class ClothoidSolver {
  constructor(mesh) {
    this.mesh = mesh;
  }

  solve() {
    let mesh = this.mesh;

    let temp = new Array(32);
    let changeOrder = (order) => {
      for (let e of mesh.edges) {
        let ks = e.curve.ks;

        for (let i = 0; i < e.curve.order; i++) {
          temp[i] = ks[i];
        }

        for (let i = 0; i < order; i++) {
          let i1 = i*(e.curve.order - 1)/(order - 1);
          let t = Math.fract(i1);
          i1 = ~~i1;

          let i2 = i1 + 1;

          let k;
          if (i2 < order) {
            k = temp[i1] + (temp[i2] - temp[i])*t;
          } else {
            k = temp[i1];
          }

          ks[i] = k;
        }

        e._ks = new Float32Array(e.curve.ks.buffer, 0, order);
        e.curve.order = order;
      }
    }

    let solver = new Solver();
    solver.threshold = 0.001;

    for (let e of mesh.edges) {
      if (!(e.curve.ks instanceof Float64Array)) {
        e.curve.ks = new Float64Array(e.curve.ks);
        e.curve._ks = new Float64Array(e.curve.ks.buffer, 0, e.order);
      }

      e.curve.ks.fill(0.001);
      e.curve.update(e);
    }

    function tan_c(params) {
      let [v, e1, e2] = params;

      e1.update();
      e2.update();

      let s1 = v === e1.v1 ? 0 : e1.length;
      let s2 = v === e2.v1 ? 0 : e2.length;

      let t1 = e1.derivative(s1);
      let t2 = e2.derivative(s2);

      if (!(v === e1.v1) === !(v === e2.v1)) {
        t2.negate();
      }

      t1.normalize();
      t2.normalize();

      //return t1.vectorDistanceSqr(t2);
      return Math.acos(t1.dot(t2));
    }

    function curv_c(params) {
      let [v, e1, e2, disabled] = params;

      let k1, k2;

      if (disabled) {
        for (let i=0; i<e1.order; i++) {
          e1.ks[i] *= 0.98;
          e2.ks[i] *= 0.98;
        }
        return 0.0;
      }

      for (let i=0; i<e1.order; i++) {
        //e1.ks[i] *= 0.98;
       // e2.ks[i] *= 0.98;
      }

      if (0) {
        let bad = false;

        for (let e of [e1, e2]) {
          let scale = e.ks[KSCALE];
          e._update();

          if (scale > 4.0 || bad) {
            params[3] = true;

            for (let i = 0; i < e.order; i++) {
              let eps = 0.0005;
              //e.ks[i] /= ((scale*0.05 - 0.05*1.0)**1.0)*eps + (1.0 - eps);
              e.ks[i] *= 0.0;
            }

            e._update();
            bad = true;
          }
        }

        if (bad) {
          return 0.0;
        }
      }
      //return 0.0;

      //return 0.0;
      k1 = v === e1.v1 ? e1.ks[0] : e1.ks[e1.order - 1];
      k2 = v === e2.v1 ? e2.ks[0] : e2.ks[e2.order - 1];

      let scale1 = e1.ks[KSCALE];
      let scale2 = e2.ks[KSCALE];

      const flip = (v === e1.v1) !== (v === e2.v1);

      if (flip) {
        k2 = -k2;
      }

      k1 /= scale1;
      k2 /= scale2;

      let k = (k1 + k2)*0.5;

      const fac = 0.5;
      k1 += (k - k1)*fac;
      k2 += (k - k2)*fac;

      if (flip) {
        k2 = -k2;
      }

      if (v === e1.v1) {
        e1.ks[0] = k1*scale1;
      } else {
        e1.ks[e1.order - 1] = k1*scale1;
      }

      if (v === e2.v1) {
        e2.ks[0] = k2*scale2;
      } else {
        e2.ks[e2.order - 1] = k2*scale2;
      }

      e1.update();
      e2.update();

      return 0.0;
    }

    const order = mesh.order;

    let badvs = new Set();

    for (let v of mesh.verts) {
      if (v.edges.length !== 2) {
        continue;
      }

      let [e1, e2] = v.edges;
      let v1 = e1.otherVertex(v);
      let v2 = e2.otherVertex(v);

      e1 = e1.curve;
      e2 = e2.curve;


      let t1 = new Vector2(v1).sub(v).normalize();
      let t2 = new Vector2(v2).sub(v).normalize();
      let th = Math.acos(t1.dot(t2)*0.99999);
      if (th < Math.PI*0.4) {
        badvs.add(v);
        continue;
      }

      let ks1 = e1._ks, ks2 = e2._ks;

      let con;

      con = new Constraint("curv_c", curv_c, [ks1, ks2], [v, e1, e2, false], 1.0);
      //solver.add(con);

      con = new Constraint("tan_c", tan_c, [ks1, ks2], [v, e1, e2], 1.0);
      solver.add(con);
    }

    let err;

    if (0) {
      let goalOrder = KORDER;

      for (let order = 2; order <= goalOrder; order++) {
        changeOrder(order);
        err = solver.solve(25, 0.7);// Math.random() > 0.9);
      }
    } else {
      //changeOrder(3);
      //err = solver.solve(15, 0.7, Math.random() > 0.9);
      //changeOrder(KORDER);
      err = solver.solve(55, 0.7);// Math.random() > 0.9);
    }

    for (let e of mesh.edges) {
      e.curve.update();
    }

    for (let v of badvs) {
      let [e1, e2] = v.edges;
      if (v === e1.v1) {
        e1.curve.ks[0] = 0.0;
      } else {
        e1.curve.ks[order-1] = 0.0;
      }
      if (v === e2.v1) {
        e2.curve.ks[0] = 0.0;
      } else {
        e2.curve.ks[order-1] = 0.0;
      }
    }
    //console.log("error:", err.toFixed(3));
  }
}