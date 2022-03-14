import {
  simple, nstructjs, util, math, Vector2, UIBase, Icons, KeyMap, haveModal, ToolOp, Menu
} from '../path.ux/pathux.js';
import {getElemColor} from '../stroker/mesh.js';
import {MeshEditor} from './mesh_editor.js';
import {BrushModes} from './brush.js';

export const EditMenu = [];

export class LoadDefaultsOp extends ToolOp {
  static tooldef() {
    return {
      uiname  : "Load Defaults",
      toolpath: "app.load_defaults",
      inputs  : {},
      outputs : {}
    }
  }

  exec(ctx) {
    ctx.state.createNewFile(true);
    window.redraw_all();
  }
}

ToolOp.register(LoadDefaultsOp);

simple.Editor.registerAppMenu(function (ctx, con, menubarEditor) {
  con.menu("File", [
    "app.new"
  ]);

  con.menu("Edit", EditMenu);
});

export class Workspace extends simple.Editor {
  constructor() {
    super();

    this._last_update_key = undefined;

    this.canvas = document.createElement("canvas");
    this.g = this.canvas.getContext("2d");

    this.toolmode = new MeshEditor();
    this.shadow.appendChild(this.canvas);

    this.keymap = new KeyMap();

    let eventBad = (e) => {
      if (haveModal()) {
        return true;
      }

      let elem = this.ctx.screen.pickElement(e.x, e.y);
      return elem && elem !== this && elem !== this.canvas;
    }

    this.addEventListener("pointerdown", (e) => {
      if (eventBad(e)) {
        return;
      }

      let mpos = this.getLocalMouse(e.x, e.y);
      if (this.ctx.properties.brushMode === BrushModes.DEBUG) {
        let [x, y] = this.getLocalMouse(e.x, e.y);

        this.ctx.api.execTool(this.ctx, `brush.stroke(x=${x} y=${y} haveXY=true)`, {
          x, y, haveXY: true
        });
      } else if (this.ctx.properties.brushMode === BrushModes.BASIC) {
        let x = e.x, y = e.y;

        this.ctx.api.execTool(this.ctx, `brush.test_stroke(x=${x} y=${y} haveXY=true)`, {
          x, y, haveXY: true
        });
      } else {
        this.toolmode.on_mousedown(mpos[0], mpos[1], e);
      }
    });

    this.addEventListener("pointermove", (e) => {
      if (eventBad(e)) {
        return;
      }

      let mpos = this.getLocalMouse(e.x, e.y);
      this.toolmode.on_mousemove(mpos[0], mpos[1], e);
    });

    this.addEventListener("pointerup", (e) => {
      if (eventBad(e)) {
        return;
      }

      let mpos = this.getLocalMouse(e.x, e.y);
      this.toolmode.on_mouseup(mpos[0], mpos[1], e);
    });
  }

  static defineAPI(api, st) {

  }

  static define() {
    return {
      tagname : "workspace-editor-x",
      areaname: "workspace-editor-x",
      uiname  : "Workspace",
    }
  }

  getGlobalMouse(x, y) {
    let mpos = new Vector2();
    let r = this.canvas.getBoundingClientRect();

    let dpi = UIBase.getDPI();

    mpos[0] = x/dpi + r.x;
    mpos[1] = y/dpi + r.y;

    return mpos;
  }

  getLocalMouse(x, y) {
    let mpos = new Vector2();
    let r = this.canvas.getBoundingClientRect();

    let dpi = UIBase.getDPI();

    mpos[0] = (x - r.x)*dpi;
    mpos[1] = (y - r.y)*dpi;

    return mpos;
  }

  getKeyMaps() {
    return [this.keymap, this.toolmode.keymap];
  }

  init() {
    super.init();

    this.toolmode.ctx = this.ctx;

    EditMenu.length = 0;
    EditMenu.push(["Undo", () => this.ctx.toolstack.undo(), "CTRL-Z"]);
    EditMenu.push(["Redo", () => this.ctx.toolstack.undo(), "CTRL-SHIFT-Z"]);
    EditMenu.push(Menu.SEP);

    let makeCB = (toolpath) => {
      return () => this.ctx.api.execTool(this.ctx, toolpath);
    }

    for (let keymap of this.getKeyMaps()) {
      for (let hk of keymap) {
        if (typeof hk.action === "string") {
          let tool = this.ctx.api.parseToolPath(hk.action);
          let tdef = tool.tooldef();
          let uiname = hk.uiname ?? tdef.uiname;


          EditMenu.push([uiname, makeCB(hk.action), hk.buildString()]);
        }
      }
    }

    EditMenu.push("mesh.flip_edge");

    //EditMenu.push(
    let sidebar = this.makeSideBar();

    let header = this.header;
    let row;

    row = header.row();
    row.iconbutton(Icons.UNDO, "Undo", () => {
      this.ctx.toolstack.undo();
    });
    row.iconbutton(Icons.REDO, "Redo", () => {
      this.ctx.toolstack.redo();
    });

    row.button("Save Defaults", () => {
      _appstate.saveStartupFile();
    })

    row.tool("app.load_defaults()");
    row.tool("app.reset()");

    row.prop("brush.radius");
    row.prop("brush.spacing");

    row.prop("properties.brushMode");

    let tab;
    tab = sidebar.tab("Options");

    let props = UIBase.createElement("props-bag-editor-x");
    props.setAttribute("datapath", "properties");

    tab.add(props);
  }

  draw() {
    if (!this.ctx) {
      return;
    }

    let canvas = this.canvas;

    let dpi = UIBase.getDPI();
    let w = ~~(this.size[0]*dpi);
    let h = ~~(this.size[1]*dpi) - 50*dpi;

    if (w !== canvas.width || h !== canvas.height) {
      canvas.width = w;
      canvas.height = h;

      canvas.style["width"] = "" + (w/dpi) + "px";
      canvas.style["height"] = "" + (h/dpi) + "px";
    }

    this.g.clearRect(0, 0, canvas.width, canvas.height);
    console.log("draw!");

    this.toolmode.draw(this.ctx, this.canvas, this.g);
  }

  update() {
    let key = "" + this.size[0] + ":" + this.size[1];

    if (key !== this._last_update_key) {
      this._last_update_key = key;
      window.redraw_all();
    }
  }

  setCSS() {
    this.canvas.style["position"] = "absolute";
  }
}

Workspace.STRUCT = nstructjs.inherit(Workspace, simple.Editor, "Workspace") + `
}`;
simple.Editor.register(Workspace);

