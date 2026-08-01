import { type Vector2, type Vec3Mixin } from "../math/index.js";
import * as nstructjs from "nstructjs";

/**
 * The subset of `CanvasRenderingContext2D` the debug `draw` hooks use.
 *
 * Declared structurally so the library carries no DOM dependency — a real canvas context
 * satisfies this, and so does any stub in a test. {@link CanvasPaint} is wide enough to
 * accept the gradient and pattern values a real context allows.
 */
export type CanvasPaint = string | object;

export interface Canvas2DLike {
  fillStyle: CanvasPaint;
  strokeStyle: CanvasPaint;
  lineWidth: number;

  beginPath(): void;
  moveTo(x: number, y: number): void;
  lineTo(x: number, y: number): void;
  rect(x: number, y: number, w: number, h: number): void;
  arc(x: number, y: number, r: number, a0: number, a1: number): void;
  fill(): void;
  stroke(): void;
}

/**
 * A mesh vertex, from a curve's point of view: an xyz-indexable point that knows its
 * incident edges.
 */
export interface CurveVertex extends Vec3Mixin {
  edges: CurveEdge[];
  otherEdge(e: CurveEdge): CurveEdge | undefined;
}

/**
 * Anything a curve needs from the edge that owns it.
 *
 * Curves are handed their owning edge on `init`/`update` so they can look at neighbours
 * — the bezier and b-spline solvers walk outward from an edge to place handles.
 */
export interface CurveEdge {
  v1: CurveVertex;
  v2: CurveVertex;

  otherVertex(v: CurveVertex): CurveVertex;
}

/**
 * Base class for every curve primitive.
 *
 * **`s` is arclength, not a normalized parameter.** `evaluate(0.5)` means half a unit
 * along the curve, not halfway along it. Implementations that work internally in a
 * normalized parameter convert on entry.
 *
 * Returned vectors are borrowed from a {@link CacheRing} and are only valid for a bounded
 * number of subsequent calls. Copy anything you intend to keep.
 */
export abstract class Curve {
  /**
   * Empty on purpose: subclasses carry all the state. Registering the base lets
   * `Edge.curve` be serialized as `abstract(Curve)`, so a mesh round-trips whichever
   * curve type it was built with. The original hardcoded the clothoid's struct name.
   */
  static STRUCT = nstructjs.inlineRegister(
    this,
    `
Curve {
}
`
  );

  abstract get length(): number;

  abstract evaluate(s: number): Vector2;
  abstract derivative(s: number): Vector2;
  abstract derivative2(s: number): Vector2;

  /** Signed curvature at arclength `s`. */
  curvature(s: number) {
    const dv1 = this.derivative(s);
    const dv2 = this.derivative2(s);

    return (dv1[0] * dv2[1] - dv1[1] * dv2[0]) / Math.pow(dv1.dot(dv1), 3.0 / 2.0);
  }

  /** Unit normal: the tangent rotated a quarter turn. */
  normal(s: number) {
    const dv = this.derivative(s);
    const tmp = dv[0];

    dv[0] = -dv[1];
    dv[1] = tmp;

    return dv.normalize();
  }

  /** Mark the curve dirty; recomputation is deferred to the next query. */
  abstract update(e?: CurveEdge): this;

  /** One-time setup against the owning edge. Defaults to `update`. */
  init(e: CurveEdge): void {
    this.update(e);
  }

  /** Optional debug overlay. Most curves have nothing extra worth drawing. */
  draw(_g: Canvas2DLike): void {
    /* no-op */
  }
}
