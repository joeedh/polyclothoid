import {
  util, nstructjs, simple, FloatProperty,
  Vec2Property, ToolOp, Vector2, Vector3, Quat, Vector4,
  Matrix4, BoolProperty
} from '../path.ux/pathux.js';
import {Mesh} from '../stroker/mesh.js';
import {BezierSolver, CubicBezier} from '../stroker/bezier.js';
import {Stroker} from '../stroker/stroke.js';

/** Note: this is a bitmask, appearances are deceiving */
export const ColorOptions = {
  CLOTHOID      : 1,
  BEZIER        : 2,
  BOTH          : 3,
  CONTROL_POINTS: 4
};

export const ColorValues = {
  CLOTHOID      : 0,
  BEZIER        : 1,
  CONTROL_POINTS: 2
};

export const BrushModes = {
  NONE : 0,
  DEBUG: 1,
  BASIC: 2,
};

export class Brush {
  constructor() {
    this.radius = 25;
    this.spacing = 0.1;
  }

  static defineAPI(api, st) {
    st.float("radius", "radius", "Radius").unit("pixel").range(0.25, 500.0).slideSpeed(3.0);
    st.float("spacing", "spacing", "Spacing").noUnits().range(0.01, 3);
  }
}

Brush.STRUCT = `
Brush {
  radius       : float;
  spacing      : float;
}
`;
simple.DataModel.register(Brush);

export class DotSample {
  constructor(radius, spacing, x, y, dx, dy, v, time) {
    this.radius = radius;
    this.spacing = spacing;
    this.x = x;
    this.y = y;
    this.dx = dx;
    this.dy = dy;
    this.v = v;
    this.time = time;
  }
}

export class BrushTestOp extends ToolOp {
  static tooldef() {
    return {
      uiname  : "Test Stroke",
      toolpath: "brush.test_stroke",
      inputs  : {
        x     : new FloatProperty(),
        y     : new FloatProperty(),
        haveXY: new BoolProperty(),
      },
      outputs : {},
      is_modal: true,
    }
  }

  modalStart(ctx) {
    let {x, y, haveXY} = this.getInputs();

    if (haveXY) {
      [x, y] = this.modal_ctx.workspace.getLocalMouse(x, y);
    } else {
      x = undefined;
      y = undefined;
    }

    let ret = super.modalStart(ctx);

    let brush = ctx.brush;
    this.stroker = new Stroker(this.on_point.bind(this), haveXY, x, y, brush.radius, brush.spacing)
    this.stroker.lag = 1;

    return ret;
  }

  on_point(x, y, dx, dy, t) {
    let ctx = this.modal_ctx;
    let workspace = ctx.workspace;
    let toolmode = workspace.toolmode;
    let brush = ctx.brush;

    toolmode.pushDot(x, y, brush.radius);
    window.redraw_all();
  }

  on_pointermove(e) {
    let ctx = this.modal_ctx;
    let workspace = ctx.workspace;
    let toolmode = workspace.toolmode;
    let brush = ctx.brush;

    let [x, y] = workspace.getLocalMouse(e.x, e.y);
    this.stroker.onInput(x, y, brush.radius, brush.spacing);
  }

  on_pointerup(e) {
    this.modalEnd(false);
  }
}

ToolOp.register(BrushTestOp);

export class BrushOp extends ToolOp {
  constructor() {
    super();

    this.samples = [];
    this.lastMpos = new Vector2();
    this.lastTime = util.time_ms();
    this.mpos = new Vector2();

    this.first = true;
    this.first2 = 2;

    this.time = util.time_ms();
    this.t = 0;
    this.lastT = 0;

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
    }
  }

  static tooldef() {
    return {
      uiname  : "Stroke",
      toolpath: "brush.stroke",
      inputs  : {
        x     : new FloatProperty(),
        y     : new FloatProperty(),
        haveXY: new BoolProperty(),
      },
      outputs : {},
      is_modal: true,
    }
  }

  modalStart(ctx) {
    this.first = true;
    this.first2 = 2

    super.modalStart(ctx);
  }

  undoPre(ctx) {
    let workspace = ctx.workspace;
    this._undo = workspace.toolmode.dots.concat([]);
  }

  undo(ctx) {
    let workspace = ctx.workspace;
    let toolmode = workspace.toolmode;

    let tmp = toolmode.dots;
    toolmode.dots = this._undo;
    this._undo = tmp;

    window.redraw_all();
  }

  redo(ctx) {
    this.undo(ctx);
  }

  calcUndoMem(ctx) {
    return 1;
  }

  on_pointermove(e) {
    let ctx = this.modal_ctx;
    let brush = ctx.brush;
    let workspace = ctx.workspace;
    let mesh = ctx.mesh;

    let radius = brush.radius;
    let spacing = brush.spacing;

    let mpos = workspace.getLocalMouse(e.x, e.y);
    this.mpos.load(mpos);


    if (this.first) {
      let {x, y, haveXY} = this.getInputs();

      this.lasts1 = 0; //clothoid
      this.lasts2 = 0; //bezier

      if (haveXY) {
        this.lastMpos.loadXY(x, y);
        this.last.mpos1.loadXY(x, y);
        this.last.mpos2.loadXY(x, y);
        this.last.mpos3.loadXY(x, y);
      } else {
        this.lastMpos.load(mpos);
        this.last.mpos1.load(mpos);
        this.last.mpos2.load(mpos);
        this.last.mpos3.load(mpos);
        return;
      }

      this.first = false;
    }


    const lag = ctx.properties.lag;

    if (Math.random() < ctx.properties.eventDrop) {
      return;
    }

    let dis = mpos.vectorDistance(this.last.mpos1);
    let dt = dis/(radius*2.0);

    if (dt > spacing*lag) {
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
          console.log("DIS", v1.vectorDistance(v2));

          mesh.makeEdge(v1, v2);
          e = mesh.makeEdge(v2, v3);
          mesh.makeEdge(v3, v4);
          mesh.makeEdge(v4, v5);
        } else if (this.first2) {
          this.first2 = false;

          this.lasts1 = spacing*(2.0 * radius);

          v1 = mesh.makeVertex(this.last.mpos3);
          v2 = mesh.makeVertex(this.last.mpos2);
          v3 = mesh.makeVertex(this.last.mpos1);
          v4 = mesh.makeVertex(mpos);

          v1[2] = v2[2] = v3[2] = v4[2] = 0.0;
          console.log("DIS", v1.vectorDistance(v2));

          e = mesh.makeEdge(v1, v2);
          mesh.makeEdge(v3, v4);
        }

        mesh.solve();

        for (let step = 0; step < 2; step++) {
          let ds = spacing*(2.0*radius);
          let s = !step ? this.lasts1 : this.lasts2;
          let elen = e.length;

          if (isNaN(elen) || isNaN(ds)) {
            console.error("NaN!");
            return;
          }

          if (ds === 0.0) {
            console.error("Spacing was zero!");
            return;
          }

          if (e === undefined) {
            debugger;
            console.error("Missing edge");
            return;
          }

          console.log(elen);
          while (s < elen) {
            let p = e.evaluate(s);

            workspace.toolmode.pushDot(p[0], p[1], radius, step);
            s += ds;
          }

          if (!step) {
            this.lasts1 = s - elen;
          } else {
            this.lasts2 = s - elen;
          }

          mesh.switchSplineType(CubicBezier, BezierSolver);
          mesh.solve();
        }

        workspace.toolmode.pushDot(mpos[0], mpos[1], radius, ColorValues.CONTROL_POINTS);
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

  finish() {
    this.modalEnd(false);
  }

  on_pointerup(e) {
    this.finish();
  }
}

ToolOp.register(BrushOp);
