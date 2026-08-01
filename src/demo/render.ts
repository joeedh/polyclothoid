import { Vector2, Vector4, getElemColor, type Mesh } from "../index.js";
import { DabKind, type DabList } from "./dabs.js";
import { type Settings } from "./ui.js";

const DAB_COLORS = ["rgba(255, 130, 60, 0.16)", "rgba(120, 120, 130, 0.16)", "rgba(30, 80, 240, 0.9)"];

const VERT_SIZE = 6;
const CURVE_STEPS = 48;

function css(c: Vector4) {
  const to255 = (x: number) => Math.round(Math.min(Math.max(x, 0), 1) * 255);

  return `rgba(${to255(c[0])}, ${to255(c[1])}, ${to255(c[2])}, ${c[3]})`;
}

function dabVisible(kind: number, settings: Settings) {
  switch (kind) {
    case DabKind.CLOTHOID:
      return settings.showClothoid;
    case DabKind.BEZIER:
      return settings.showBezier;
    default:
      return settings.showInput;
  }
}

function drawDabs(g: CanvasRenderingContext2D, dabs: DabList, settings: Settings) {
  let paint = -1;

  dabs.forEach((x, y, radius, kind) => {
    if (!dabVisible(kind, settings)) {
      return;
    }

    /* Batching by colour keeps this to one fill per kind instead of one per dab. */
    if (kind !== paint) {
      if (paint >= 0) {
        g.fill();
      }

      paint = kind;
      g.fillStyle = DAB_COLORS[kind];
      g.beginPath();
    }

    g.moveTo(x + radius, y);
    g.arc(x, y, radius, -Math.PI, Math.PI);
  });

  if (paint >= 0) {
    g.fill();
  }
}

function drawCurves(g: CanvasRenderingContext2D, mesh: Mesh) {
  g.strokeStyle = "rgb(40, 90, 200)";
  g.lineWidth = 2;
  g.beginPath();

  for (const e of mesh.edges.visible) {
    const length = e.length;

    if (!isFinite(length) || length === 0.0) {
      continue;
    }

    const ds = length / (CURVE_STEPS - 1);

    for (let i = 0, s = 0.0; i < CURVE_STEPS; i++, s += ds) {
      const p = e.evaluate(s);

      if (i === 0) {
        g.moveTo(p[0], p[1]);
      } else {
        g.lineTo(p[0], p[1]);
      }
    }
  }

  g.stroke();
  g.lineWidth = 1;
}

/** Curvature comb: a normal at each sample, scaled by the signed curvature there. */
function drawCombs(g: CanvasRenderingContext2D, mesh: Mesh, settings: Settings) {
  const p = new Vector2();

  g.strokeStyle = "rgb(40, 160, 60)";
  g.beginPath();

  for (const e of mesh.edges.visible) {
    const length = e.length;

    if (!isFinite(length) || length === 0.0) {
      continue;
    }

    const steps = Math.min(Math.max(Math.floor(length / 5), 2), 512);
    const ds = length / (steps - 1);

    for (let i = 0, s = 0.0; i < steps; i++, s += ds) {
      const n = e
        .normal(s)
        .normalize()
        .mulScalar(e.curvature(s) * settings.normalScale);

      p.load(e.evaluate(s));

      g.moveTo(p[0], p[1]);
      g.lineTo(p[0] + n[0], p[1] + n[1]);
    }
  }

  g.stroke();
}

function drawMesh(g: CanvasRenderingContext2D, mesh: Mesh, settings: Settings) {
  mesh.ensureSolve();

  drawCurves(g, mesh);

  if (settings.drawNormals) {
    drawCombs(g, mesh, settings);
  }

  for (const e of mesh.edges.visible) {
    if (settings.drawExtra) {
      e.curve.draw(g);
    }

    g.strokeStyle = css(getElemColor(mesh.edges, e));
    g.beginPath();
    g.moveTo(e.v1[0], e.v1[1]);
    g.lineTo(e.v2[0], e.v2[1]);
    g.stroke();
  }

  for (const f of mesh.faces.visible) {
    const color = new Vector4(getElemColor(mesh.faces, f));
    color[3] = 0.15;

    g.fillStyle = css(color);
    g.beginPath();

    for (const list of f.lists) {
      let first = true;

      for (const l of list) {
        if (first) {
          first = false;
          g.moveTo(l.v[0], l.v[1]);
        } else {
          g.lineTo(l.v[0], l.v[1]);
        }
      }

      g.closePath();
    }

    g.fill();
  }

  for (const v of mesh.verts.visible) {
    g.fillStyle = css(getElemColor(mesh.verts, v));
    g.beginPath();
    g.rect(v[0] - VERT_SIZE * 0.5, v[1] - VERT_SIZE * 0.5, VERT_SIZE, VERT_SIZE);
    g.fill();
  }
}

export function drawScene(
  g: CanvasRenderingContext2D,
  mesh: Mesh,
  dabs: DabList,
  settings: Settings,
  width: number,
  height: number
) {
  g.clearRect(0, 0, width, height);

  drawDabs(g, dabs, settings);

  if (settings.drawMesh) {
    drawMesh(g, mesh, settings);
  }
}
