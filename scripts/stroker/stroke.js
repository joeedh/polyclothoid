import {util, Vector2, Vector4} from './pathux.js';
import {Mesh} from './mesh.js';
import {BezierSolver, CubicBezier} from './bezier.js';

export class Stroker {
  /** callback(x, y, dx, dy, interp_t) */
  constructor(callback, doFirst, firstX, firstY, radius, spacing) {
    this.mpos = new Vector2();

    this.lag = 1.0;

    this.last = {
      mpos1: new Vector2(),
      mpos2: new Vector2(),
      mpos3: new Vector2(),
      mpos4: new Vector4(),
      mpos5: new Vector4(),
      mpos6: new Vector4(),
      time1: util.time_ms(),
      time2: util.time_ms(),
      time3: util.time_ms(),
      v1 : 0,
      v2 : 0,
      v3 : 0,
      v4 : 0,
    }

    this.lastMpos = new Vector2();
    this.callback = callback;
    this.x = undefined;
    this.y = undefined;
    this.haveXY = false;

    this.first = true;
    this.first2 = true;

    if (firstX !== undefined) {
      //this.x = firstX;
      //this.y = firstY;
      //this.haveXY = true;
      console.log(firstX, firstY, radius, spacing);

      if (doFirst) {
        this.callback(firstX, firstY, 0, 0, 0);
      }

      this.onInput(firstX, firstY, radius, spacing);

    }
  }

  onInput(x, y, radius, spacing) {
    let mpos = new Vector2().loadXY(x, y);
    this.mpos.load(mpos);

    if (this.first) {
      this.lasts1 = 0;

      this.lastMpos.load(mpos);
      this.last.mpos1.load(mpos);

      this.first = false;
      return;
    }


    let dis = mpos.vectorDistance(this.last.mpos1);
    let dt = dis/(radius*2.0);

    if (dt > spacing*this.lag) {
      /*
      let v = mesh.makeVertex(mpos);
      v[2] = 0.0;

      if (this.lastV) {
        mesh.makeEdge(this.lastV, v);
        mesh.regenSolve();
      }

      mesh.ensureSolve();
      */

      let v = true;
      this.lastV = v;

      let ok = this.last.v4;
      ok = ok || (this.first2 && this.last.v3);

      if (ok) {
        let ds = spacing*(2.0*radius);

        let mesh = new Mesh();
        let v1, v2, v3, v4, v5;
        let e;

        if (!this.first2) {
          v1 = mesh.makeVertex(this.last.mpos4);
          v2 = mesh.makeVertex(this.last.mpos3);
          v3 = mesh.makeVertex(this.last.mpos2);
          v4 = mesh.makeVertex(this.last.mpos1);
          v5 = mesh.makeVertex(mpos);

          v1[2] = v2[2] = v3[2] = v4[2] = v5[2] = 0.0;

          mesh.makeEdge(v1, v2);
          e = mesh.makeEdge(v2, v3);
          mesh.makeEdge(v3, v4);
          mesh.makeEdge(v4, v5);
        } else if (this.first2) {
          this.first2 = false;

          v1 = mesh.makeVertex(this.last.mpos4);
          v2 = mesh.makeVertex(this.last.mpos3);
          v3 = mesh.makeVertex(this.last.mpos2);
          v4 = mesh.makeVertex(mpos);

          this.lasts1 = ds;

          v1[2] = v2[2] = v3[2] = v4[2] = 0.0;
          console.log("DIS", v1.vectorDistance(v2));

          e = mesh.makeEdge(v1, v3);
          mesh.makeEdge(v3, v4);
        }

        mesh.solve();

        let s = this.lasts1;
        let elen = e.length;

        if (isNaN(elen) || isNaN(ds)) {
          console.error("NaN!");
          return;
        }

        if (ds === 0.0) {
          console.error("Spacing was zero!");
          return;
        }

        if (!e) {
          debugger;
          console.error("Missing edge");
          return;
        }

        while (s < elen) {
          let p = e.evaluate(s);
          let dv = e.derivative(s);

          let t = s/elen;

          this.callback(p[0], p[1], dv[0], dv[1], t);
          s += ds;
        }

        this.lasts1 = s - elen;
      }

      this.last.v6 = this.last.v5;
      this.last.v5 = this.last.v4;
      this.last.v4 = this.last.v3;
      this.last.v3 = this.last.v2;
      this.last.v2 = this.last.v1;
      this.last.v1 = true;

      this.last.mpos6.load(this.last.mpos5);
      this.last.mpos5.load(this.last.mpos4);
      this.last.mpos4.load(this.last.mpos3);
      this.last.mpos3.load(this.last.mpos2);
      this.last.mpos2.load(this.last.mpos1);
      this.last.mpos1.load(mpos);

      this.last.time3 = this.last.time2;
      this.last.time2 = this.last.time1;
      this.last.time1 = this.time;
    }

    this.time = util.time_ms();
    this.lastMpos.load(mpos);
    window.redraw_all();
  }
}