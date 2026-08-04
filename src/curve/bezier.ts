import { CacheRing, Vector2, clamp, fract, type VecLike } from "../math/index.js";
import { Curve, type CurveEdge } from "./curve.js";
import { emptySolveReport } from "./diagnostics.js";
import { type CurveSolver, type SolvableMesh } from "./mesh_types.js";
import { walk } from "./topology.js";
import * as nstructjs from "nstructjs";

const evalRets = new CacheRing(() => new Vector2(), 512);

/** de Casteljau evaluation of a scalar cubic. */
export function cubic(a: number, b: number, c: number, d: number, t: number) {
  const k1 = a + (b - a) * t;
  const k2 = b + (c - b) * t;
  const k3 = c + (d - c) * t;

  const ka = k1 + (k2 - k1) * t;
  const kb = k2 + (k3 - k2) * t;

  return ka + (kb - ka) * t;
}

const ARC_TABLE_SIZE = 1024;

/**
 * Cubic bezier with an arclength reparameterization table.
 *
 * A bezier's natural parameter `t` is not arclength, so `evaluate(s)` looks `s` up in a
 * table built by {@link genTable} and evaluates at the corresponding `t`. That keeps this
 * class interchangeable with {@link Clothoid}, whose parameter is arclength natively.
 */
export class CubicBezier extends Curve {
  static override STRUCT = nstructjs.inlineRegister(
    this,
    `
CubicBezier {
  v1 : vec2;
  h1 : vec2;
  h2 : vec2;
  v2 : vec2;
}
`
  );

  v1: Vector2;
  v2: Vector2;
  h1: Vector2;
  h2: Vector2;

  regen = 1;
  arcLength = 0;

  /** s -> t lookup, uniformly sampled in s. Holes are filled in by {@link genTable}. */
  stable: number[] = new Array<number>(ARC_TABLE_SIZE);

  constructor(v1?: VecLike, h1?: VecLike, h2?: VecLike, v2?: VecLike) {
    super();

    this.v1 = new Vector2(v1);
    this.v2 = new Vector2(v2);
    this.h1 = new Vector2(h1);
    this.h2 = new Vector2(h2);
  }

  override get length() {
    if (this.regen) {
      this.genTable();
    }

    return this.arcLength;
  }

  /** Evaluate at the natural bezier parameter, not arclength. */
  _evaluate(t: number) {
    const ret = evalRets.next();
    const { v1, h1, h2, v2 } = this;

    for (let i = 0; i < 2; i++) {
      ret[i] = cubic(v1[i], h1[i], h2[i], v2[i], t);
    }

    return ret;
  }

  /**
   * Build the arclength table by walking `t` uniformly, binning each sample by its
   * accumulated arclength, averaging collisions, then interpolating across empty bins.
   */
  genTable() {
    this.regen = 0;

    const stable = this.stable;
    const size = stable.length;
    const samples: number[] = [];

    let lastp: Vector2 | undefined;
    let s = 0.0;
    let t = 0.0;

    const dt = 1.0 / (size - 1);

    for (let i = 0; i < size; i++, t += dt) {
      const p = this._evaluate(t);

      if (lastp) {
        s += p.vectorDistance(lastp);
      }

      samples.push(s, t);
      lastp = new Vector2(p);
    }

    this.arcLength = s;

    const filled: number[] = new Array<number>(size).fill(0);
    const ilen = s !== 0.0 ? 1.0 / s : 0.0;

    stable.fill(0);

    for (let i = 0; i < samples.length; i += 2) {
      const si = clamp(~~(samples[i] * ilen * (size - 1)), 0, size - 1);

      stable[si] += samples[i + 1];
      filled[si]++;
    }

    for (let i = 0; i < size; i++) {
      if (filled[i]) {
        stable[i] /= filled[i];
      }
    }

    /* Flat-fill the tail past the last populated bin. */
    let last = size - 1;
    while (last >= 0 && !filled[last]) {
      last--;
    }

    const tail = last >= 0 ? stable[last] : 1.0;
    for (let i = Math.max(last, 0); i < size; i++) {
      stable[i] = tail;
      filled[i] = 1;
    }

    /* Interpolate across interior bins that no sample landed in. */
    for (let i = 0; i < size - 1; i++) {
      if (filled[i + 1]) {
        continue;
      }

      let i2 = i + 1;
      while (i2 < size && !filled[i2]) {
        i2++;
      }

      const a = stable[i];
      const b = stable[i2];
      const step = 1.0 / (i2 - i);

      let f = step;
      for (let j = i + 1; j <= i2; j++, f += step) {
        stable[j] = a + (b - a) * f;
        filled[j] = 1;
      }
    }
  }

  override evaluate(s: number) {
    if (this.regen) {
      this.genTable();
    }

    if (this.arcLength === 0.0) {
      return this._evaluate(0.0);
    }

    const stable = this.stable;

    let si = (clamp(s, 0.0, this.arcLength) / this.arcLength) * (stable.length - 1);

    const i1 = ~~si;
    const i2 = i1 + 1;

    si = fract(si);

    const t = i2 < stable.length ? stable[i1] + (stable[i2] - stable[i1]) * si : stable[i1];

    return this._evaluate(t);
  }

  override derivative(s: number) {
    const df = 0.0001;

    const a = new Vector2(this.evaluate(s));
    const b = new Vector2(this.evaluate(s + df));

    return b.sub(a).divScalar(df);
  }

  override derivative2(s: number) {
    const df = 0.0001;

    const a = new Vector2(this.derivative(s));
    const b = new Vector2(this.derivative(s + df));

    return b.sub(a).divScalar(df);
  }

  override update(e?: CurveEdge) {
    if (e) {
      this.v1.load(e.v1);
      this.v2.load(e.v2);
    }

    this.regen = 1;

    return this;
  }

  afterSTRUCT() {
    this.regen = 1;
  }
}

/**
 * Places bezier handles by a Catmull-Rom style rule.
 *
 * Each handle is a third of the way along the average of the incoming and outgoing chord
 * directions, which gives tangent continuity but says nothing about curvature — the
 * difference this project exists to measure.
 */
export class BezierSolver implements CurveSolver {
  constructor(public mesh: SolvableMesh) {}

  solve() {
    for (const e of this.mesh.edges) {
      const v1 = e.v1;
      const v2 = e.v2;
      const vp = walk(v1, e);
      const vn = walk(v2, e);

      const h1 = new Vector2(v1).sub(vp);
      const h2 = new Vector2(v2).sub(v1);

      h1.interp(h2, 0.5)
        .mulScalar(1.0 / 3.0)
        .add(v1);

      const h3 = new Vector2(vn).sub(v2);
      h2.interp(h3, 0.5)
        .mulScalar(-1.0 / 3.0)
        .add(v2);

      const curve = e.curve as CubicBezier;

      curve.v1.load(v1);
      curve.h1.load(h1);
      curve.h2.load(h2);
      curve.v2.load(v2);

      curve.update();
    }

    return emptySolveReport();
  }
}
