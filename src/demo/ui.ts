/** Which curve type the mesh solves with. See {@link Mesh.switchSplineType}. */
export type SplineType = "spower" | "clothoid" | "bezier" | "bspline";

/** Everything the control panel can change. Persisted to localStorage between reloads. */
export interface Settings {
  mode: "edit" | "paint";
  spline: SplineType;

  radius: number;
  spacing: number;
  lag: number;

  showClothoid: boolean;
  showBezier: boolean;
  showInput: boolean;

  drawMesh: boolean;
  drawExtra: boolean;
  drawNormals: boolean;
  normalScale: number;
}

export const defaultSettings: Settings = {
  mode  : "edit",
  spline: "spower",

  radius : 25,
  spacing: 0.1,
  lag    : 1.0,

  showClothoid: true,
  showBezier  : true,
  showInput   : true,

  drawMesh   : true,
  drawExtra  : true,
  drawNormals: false,
  normalScale: 4000,
};

const STORAGE_KEY = "polyclothoid.demo.settings";

export function loadSettings(): Settings {
  const saved = localStorage.getItem(STORAGE_KEY);

  if (!saved) {
    return { ...defaultSettings };
  }

  try {
    return { ...defaultSettings, ...(JSON.parse(saved) as Partial<Settings>) };
  } catch {
    console.warn("demo: could not parse saved settings");
    return { ...defaultSettings };
  }
}

export function saveSettings(settings: Settings) {
  localStorage.setItem(STORAGE_KEY, JSON.stringify(settings));
}

/** Keys of `Settings` whose value has type `T`, so the control helpers stay typed. */
type KeysOf<T> = { [K in keyof Settings]: Settings[K] extends T ? K : never }[keyof Settings];

function row(parent: HTMLElement, label: string) {
  const wrap = document.createElement("label");
  const text = document.createElement("span");

  text.textContent = label;
  wrap.className = "row";
  wrap.append(text);
  parent.append(wrap);

  return wrap;
}

function checkbox(parent: HTMLElement, settings: Settings, key: KeysOf<boolean>, label: string, onChange: () => void) {
  const input = document.createElement("input");

  input.type = "checkbox";
  input.checked = settings[key];
  input.addEventListener("change", () => {
    settings[key] = input.checked;
    onChange();
  });

  row(parent, label).append(input);
}

function slider(
  parent: HTMLElement,
  settings: Settings,
  key: KeysOf<number>,
  label: string,
  min: number,
  max: number,
  step: number,
  onChange: () => void
) {
  const input = document.createElement("input");
  const readout = document.createElement("output");

  input.type = "range";
  input.min = String(min);
  input.max = String(max);
  input.step = String(step);
  input.value = String(settings[key]);
  readout.textContent = String(settings[key]);

  input.addEventListener("input", () => {
    settings[key] = Number(input.value);
    readout.textContent = input.value;
    onChange();
  });

  row(parent, label).append(input, readout);
}

function button(parent: HTMLElement, label: string, onClick: () => void) {
  const el = document.createElement("button");

  el.type = "button";
  el.textContent = label;
  el.addEventListener("click", onClick);
  parent.append(el);
}

function group(parent: HTMLElement, title: string) {
  const box = document.createElement("section");
  const head = document.createElement("h2");

  head.textContent = title;
  box.append(head);
  parent.append(box);

  return box;
}

export interface PanelActions {
  onChange: () => void;
  clearDabs: () => void;
  resetMesh: () => void;
  switchSpline: () => void;
}

/** Returns a `sync` hook so code that changes settings outside the panel can refresh it. */
export function buildPanel(root: HTMLElement, settings: Settings, actions: PanelActions) {
  const { onChange } = actions;

  const modes = group(root, "Mode");
  const modeSelect = document.createElement("select");

  for (const mode of ["edit", "paint"] as const) {
    const option = document.createElement("option");

    option.value = mode;
    option.textContent = mode === "edit" ? "Edit mesh" : "Paint stroke";
    option.selected = settings.mode === mode;

    modeSelect.append(option);
  }

  modeSelect.addEventListener("change", () => {
    settings.mode = modeSelect.value as Settings["mode"];
    onChange();
  });

  row(modes, "Tool").append(modeSelect);

  const brush = group(root, "Brush");
  slider(brush, settings, "radius", "Radius", 1, 120, 0.5, onChange);
  slider(brush, settings, "spacing", "Spacing", 0.01, 1.0, 0.01, onChange);
  slider(brush, settings, "lag", "Input lag", 1, 40, 0.5, onChange);

  const dabs = group(root, "Dabs");
  checkbox(dabs, settings, "showClothoid", "Clothoid", onChange);
  checkbox(dabs, settings, "showBezier", "Bezier", onChange);
  checkbox(dabs, settings, "showInput", "Input points", onChange);
  button(dabs, "Clear dabs", actions.clearDabs);

  const mesh = group(root, "Mesh");
  const splineSelect = document.createElement("select");

  const splines: [SplineType, string][] = [
    ["spower", "Clothoid (s-power)"],
    ["clothoid", "Clothoid (legacy)"],
    ["bezier", "Cubic bezier"],
    ["bspline", "B-spline"],
  ];

  for (const [type, label] of splines) {
    const option = document.createElement("option");

    option.value = type;
    option.textContent = label;
    option.selected = settings.spline === type;

    splineSelect.append(option);
  }

  splineSelect.addEventListener("change", () => {
    settings.spline = splineSelect.value as SplineType;
    actions.switchSpline();
    onChange();
  });

  row(mesh, "Spline").append(splineSelect);

  checkbox(mesh, settings, "drawMesh", "Draw mesh", onChange);
  checkbox(mesh, settings, "drawExtra", "Curve debug overlay", onChange);
  checkbox(mesh, settings, "drawNormals", "Curvature combs", onChange);
  slider(mesh, settings, "normalScale", "Comb scale", 100, 20000, 100, onChange);
  button(mesh, "Reset mesh", actions.resetMesh);

  const help = group(root, "Keys");
  const list = document.createElement("dl");

  const keys: [string, string][] = [
    ["A", "select all / none"],
    ["E", "split selected edges"],
    ["D", "dissolve selected vertices"],
    ["X / Del", "delete selected vertices"],
    ["C", "clear dabs"],
    ["Space", "toggle edit / paint"],
  ];

  for (const [key, what] of keys) {
    const dt = document.createElement("dt");
    const dd = document.createElement("dd");

    dt.textContent = key;
    dd.textContent = what;
    list.append(dt, dd);
  }

  help.append(list);

  return {
    sync() {
      modeSelect.value = settings.mode;
    },
  };
}
