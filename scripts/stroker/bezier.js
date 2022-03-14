import {util, math, Vector2, Vector3, nstructjs} from './pathux.js';

const evalrets = util.cachering.fromConstructor(Vector2, 512);

export function cubic(a, b, c, d, t) {
  let k1 = a + (b - a)*t;
  let k2 = b + (c - b)*t;
  let k3 = c + (d - c)*t;

  let ka = k1 + (k2 - k1)*t;
  let kb = k2 + (k3 - k2)*t;

  return ka + (kb - ka)*t;
}

export class CubicBezier {
  constructor(v1, h1, h2, v2) {
    this.v1 = new Vector2(v1);
    this.v2 = new Vector2(v2);
    this.h1 = new Vector2(h1);
    this.h2 = new Vector2(h2);

    this.regen = 1;
    this.stable = new Array(1024);
    this.arcLength = 0;
  }

  get length() {
    if (this.regen) {
      this.genTable();
    }

    return this.arcLength;
  }

  _evaluate(t) {
    let ret = evalrets.next();

    let {v1, h1, h2, v2} = this;

    for (let i = 0; i < 2; i++) {
      ret[i] = cubic(v1[i], h1[i], h2[i], v2[i], t);
    }

    return ret;
  }

  afterSTRUCT() {
    this.regen = 1;
  }

  genTable() {
    this.regen = 0;

    let stable = this.stable;
    stable.fill(undefined);

    let t = 0, dt = 1.0/(stable.length - 1);
    let lastp;
    let s = 0.0;

    let ss = [];

    stable[0] = 0;

    for (let i = 0; i < stable.length; i++, t += dt) {
      let p = this._evaluate(t);

      if (lastp) {
        s += p.vectorDistance(lastp);
      }

      ss.push(s);
      ss.push(t);

      lastp = p;
    }

    let tots = new Array(this.stable.length);
    tots.fill(0);

    let len = this.arcLength = s;

    let ilen = len !== 0.0 ? 1.0/len : 0.0;
    let tsize = this.stable.length;

    for (let i = 0; i < ss.length; i += 2) {
      let s = ss[i], t = ss[i + 1];

      let si = ~~(s*ilen*(tsize - 1));

      if (!tots[si]) {
        stable[si] = t;
      } else {
        stable[si] += t;
      }

      tots[si]++;
    }

    for (let i = 0; i < stable.length; i++) {
      if (tots[i]) {
        stable[i] /= tots[i];
      }
    }

    /* figure out endpoint */
    let si = stable.length - 1;
    while (si >= 0 && stable[si] === undefined) {
      si--;
    }

    s = stable[si] ?? 1.0;
    while (si < stable.length) {
      stable[si] = s;
      si++;
    }

    /*interpolate empty space*/
    for (let i = 0; i < stable.length - 1; i++) {
      if (stable[i + 1] !== undefined) {
        continue;
      }

      let i1 = i;
      let i2 = i + 1;

      while (stable[i2] === undefined) {
        i2++;
      }

      let dt = 1.0/(i2 - i1);
      let t = dt;

      let a = stable[i1];
      let b = stable[i2];

      for (let j = i1 + 1; j <= i2; j++, t += dt) {
        stable[j] = a + (b - a)*t;
      }
    }
  }

  evaluate(s) {
    if (this.regen) {
      this.genTable();
    }

    s = Math.min(Math.max(s, 0.0), this.arcLength);
    s /= this.arcLength;

    let si = s*(this.stable.length - 1);
    let i1 = ~~si;
    let i2 = i1 + 1;

    si = Math.fract(si);

    let t;
    if (i2 < this.stable.length) {
      t = this.stable[i1] + (this.stable[i2] - this.stable[i1])*si;
    } else {
      t = this.stable[i1];
    }

    return this._evaluate(t);
  }

  update() {
    this.regen = 1;
    return this;
  }

  derivative(s) {
    let df = 0.0001;

    let a = this.evaluate(s);
    let b = this.evaluate(s + df);

    b.sub(a).divScalar(df);
    return b;
  }

  derivative2(s) {
    let df = 0.0001;

    let a = this.derivative(s);
    let b = this.derivative(s + df);

    b.sub(a).divScalar(df);
    return b;
  }

  curvature(s) {
    let dv1 = this.derivative(s);
    let dv2 = this.derivative2(s);

    return (dv1[0]*dv2[1] - dv1[1]*dv2[0]) / Math.pow(dv1.dot(dv1), 3.0/2.0);
  }

  normal(s) {
    let dv1 = this.derivative(s);
    let tmp = dv1[0];

    dv1[0] = -dv1[1];
    dv1[1] = tmp;

    dv1.normalize();
    return dv1;
  }

  draw() {

  }
}

CubicBezier.STRUCT = `
CubicBezier {
  v1 : vec2;
  h1 : vec2;
  h2 : vec2;
  v2 : vec2;
}
`;
nstructjs.register(CubicBezier);

export class BezierSolver {
  constructor(mesh) {
    this.mesh = mesh;
  }

  solve() {
    let mesh = this.mesh;

    function walkv(v, e, depth = 0) {
      for (let e2 of v.edges) {
        if (e2 !== e) {
          if (depth > 0) {
            return walkv(e2.otherVertex(v), e2, depth - 1);
          } else {
            return e2.otherVertex(v);
          }
        }
      }

      return v;
    }

    for (let e of mesh.edges) {
      let vp = walkv(e.v1, e);
      let v1 = e.v1;
      let v2 = e.v2;
      let vn = walkv(e.v2, e);

      let h1 = new Vector2(v1).sub(vp);
      let h2 = new Vector2(v2).sub(v1);

      h1.interp(h2, 0.5).mulScalar(1.0/3.0).add(v1);

      let h3 = new Vector2(vn).sub(v2);
      h2.interp(h3, 0.5).mulScalar(-1.0/3.0).add(v2);

      //h1.load(v1).interp(v2, 1.0/3.0);
      //h2.load(v1).interp(v2, 2.0/3.0);

      e.curve.v1.load(v1);
      e.curve.h1.load(h1);
      e.curve.h2.load(h2);
      e.curve.v2.load(v2);

      e.curve.update();
    }
  }
}