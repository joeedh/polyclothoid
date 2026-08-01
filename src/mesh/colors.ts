import { Vector4 } from "../math/index.js";
import { MeshFlags, type Element, type ElementArray } from "./mesh.js";

const mix = (a: Vector4, b: Vector4, fac: number) => new Vector4(a).interp(b, fac);

const sel = new Vector4([1, 0.8, 0, 1]);
const high = new Vector4([1, 0.8, 0.7, 1]);
const act = new Vector4([0, 0.3, 0.8, 1]);
const actsel = new Vector4([0.5, 0.3, 0.8, 1]);

/** Indexed by the bit mask select|active<<1|highlight<<2. */
export const ElemColors = [
  new Vector4([0, 0, 0, 1]),
  sel,
  act,
  mix(sel, actsel, 0.25),
  high,
  mix(high, sel, 0.5),
  mix(high, actsel, 0.5),
  new Vector4(high)
    .add(sel)
    .add(actsel)
    .mulScalar(1.0 / 3.0),
];

export function getElemColor(list: ElementArray, e: Element) {
  let mask = 0;

  if (e.flag & MeshFlags.SELECT) {
    mask |= 1;
  }

  if (e === list.active) {
    mask |= 2;
  }

  if (e === list.highlight) {
    mask |= 4;
  }

  return ElemColors[mask];
}
