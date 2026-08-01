import { BezierSolver, CubicBezier, Mesh, MeshFlags, Stroker, Vector2, type Vertex } from "../index.js";
import { DabKind, DabList } from "./dabs.js";
import { drawScene } from "./render.js";
import { buildPanel, loadSettings, saveSettings } from "./ui.js";

/** Pick radius in CSS pixels. */
const PICK_DIST = 20;

const settings = loadSettings();
const dabs = new DabList();

let mesh = makeStartMesh();
let redrawPending = false;

function makeStartMesh() {
  const m = new Mesh();
  const s = 120;
  const d = 260;

  const v1 = m.makeVertex([s, s, 0]);
  const v2 = m.makeVertex([s, s + d, 0]);
  const v3 = m.makeVertex([s + d, s + d, 0]);
  const v4 = m.makeVertex([s + d, s, 0]);

  m.makeFace([v1, v2, v3, v4]);

  return m;
}

const canvas = document.createElement("canvas");
const g = canvas.getContext("2d")!;

const panel = document.createElement("aside");
const app = document.createElement("div");

app.id = "app";
panel.id = "panel";
canvas.id = "view";

app.append(canvas, panel);
document.body.append(app);

function redraw() {
  if (redrawPending) {
    return;
  }

  redrawPending = true;

  requestAnimationFrame(() => {
    redrawPending = false;
    drawScene(g, mesh, dabs, settings, canvas.width / dpr(), canvas.height / dpr());
  });
}

function dpr() {
  return window.devicePixelRatio || 1;
}

function resize() {
  const rect = canvas.getBoundingClientRect();
  const scale = dpr();

  canvas.width = Math.max(Math.round(rect.width * scale), 1);
  canvas.height = Math.max(Math.round(rect.height * scale), 1);

  g.setTransform(scale, 0, 0, scale, 0, 0);

  redraw();
}

const panelHooks = buildPanel(panel, settings, {
  onChange() {
    saveSettings(settings);
    redraw();
  },
  clearDabs() {
    dabs.clear();
    redraw();
  },
  resetMesh() {
    mesh = makeStartMesh();
    redraw();
  },
});

/* ---------------------------------------------------------------- picking */

function localPos(e: PointerEvent | MouseEvent) {
  const rect = canvas.getBoundingClientRect();

  return new Vector2([e.clientX - rect.left, e.clientY - rect.top]);
}

function pickVertex(pos: Vector2) {
  const p = new Vector2(pos);

  let best: Vertex | undefined;
  let bestDist = PICK_DIST;

  for (const v of mesh.verts.visible) {
    const dist = Math.sqrt((v[0] - p[0]) ** 2 + (v[1] - p[1]) ** 2);

    if (dist < bestDist) {
      bestDist = dist;
      best = v;
    }
  }

  return best;
}

/* ---------------------------------------------------------------- editing */

const dragOffsets = new Map<Vertex, Vector2>();
const lastPos = new Vector2();

let dragging = false;

function startDrag(pos: Vector2, shift: boolean) {
  const v = pickVertex(pos);

  if (!v) {
    if (!shift) {
      mesh.selectNone();
      mesh.setActive(undefined);
      redraw();
    }

    return;
  }

  if (shift) {
    mesh.verts.setSelect(v, !(v.flag & MeshFlags.SELECT));
  } else if (!(v.flag & MeshFlags.SELECT)) {
    mesh.selectNone();
    mesh.verts.setSelect(v, true);
  }

  mesh.setActive(v);

  dragOffsets.clear();
  for (const sv of mesh.verts.selected) {
    dragOffsets.set(sv, new Vector2([sv[0] - pos[0], sv[1] - pos[1]]));
  }

  dragging = true;
  redraw();
}

function moveDrag(pos: Vector2) {
  for (const [v, offset] of dragOffsets) {
    v[0] = pos[0] + offset[0];
    v[1] = pos[1] + offset[1];
  }

  mesh.regenSolve();
  redraw();
}

/* --------------------------------------------------------------- painting */

let strokers: Stroker[] = [];

function startStroke() {
  const push = (kind: number) => (x: number, y: number) => {
    dabs.push(x, y, settings.radius, kind);
    redraw();
  };

  const clothoid = new Stroker(push(DabKind.CLOTHOID));
  const bezier = new Stroker(push(DabKind.BEZIER), {
    CurveCls : CubicBezier,
    SolverCls: BezierSolver,
  });

  strokers = [clothoid, bezier];

  for (const s of strokers) {
    s.lag = settings.lag;
  }
}

function moveStroke(pos: Vector2) {
  for (const s of strokers) {
    s.onInput(pos[0], pos[1], settings.radius, settings.spacing);
  }

  dabs.push(pos[0], pos[1], 2.5, DabKind.INPUT);
  redraw();
}

/* ----------------------------------------------------------------- events */

canvas.addEventListener("pointerdown", (e) => {
  const pos = localPos(e);

  canvas.setPointerCapture(e.pointerId);
  lastPos.load(pos);

  if (settings.mode === "paint") {
    startStroke();
    moveStroke(pos);
  } else {
    startDrag(pos, e.shiftKey);
  }
});

canvas.addEventListener("pointermove", (e) => {
  const pos = localPos(e);

  lastPos.load(pos);

  if (settings.mode === "paint") {
    if (strokers.length) {
      moveStroke(pos);
    }

    return;
  }

  if (dragging) {
    moveDrag(pos);
    return;
  }

  if (mesh.setHighlight(pickVertex(pos))) {
    redraw();
  }
});

function endPointer(e: PointerEvent) {
  if (canvas.hasPointerCapture(e.pointerId)) {
    canvas.releasePointerCapture(e.pointerId);
  }

  strokers = [];
  dragging = false;
  dragOffsets.clear();
}

canvas.addEventListener("pointerup", endPointer);
canvas.addEventListener("pointercancel", endPointer);

window.addEventListener("keydown", (e) => {
  if (e.target instanceof HTMLInputElement || e.target instanceof HTMLSelectElement) {
    return;
  }

  const selected = [...mesh.verts.selected];

  switch (e.key.toLowerCase()) {
    case "a":
      if (mesh.verts.selected.size > 0) {
        mesh.selectNone();
      } else {
        mesh.verts.selectAll();
      }
      break;
    case "e":
      for (const edge of [...mesh.edges]) {
        if (edge.v1.flag & edge.v2.flag & MeshFlags.SELECT) {
          mesh.splitEdge(edge);
        }
      }
      break;
    case "d":
      for (const v of selected) {
        if (v.edges.length === 2) {
          mesh.dissolveVertex(v);
        }
      }
      break;
    case "x":
    case "delete":
      for (const v of selected) {
        mesh.killVertex(v);
      }
      break;
    case "c":
      dabs.clear();
      break;
    case " ":
      settings.mode = settings.mode === "edit" ? "paint" : "edit";
      saveSettings(settings);
      panelHooks.sync();
      break;
    default:
      return;
  }

  e.preventDefault();
  mesh.regenSolve();
  redraw();
});

window.addEventListener("resize", resize);

resize();
