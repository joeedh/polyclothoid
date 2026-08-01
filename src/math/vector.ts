/**
 * Vector types for the stroker.
 *
 * These replace path.ux's `vectormath.js`, which generated its component-wise
 * operations with `eval` at module load. The set of methods here is the subset the
 * library actually uses, written out by hand and unrolled per component — the stroker's
 * inner loop runs these per brush dab.
 *
 * Vectors extend `Array<number>`, which is what lets {@link Vertex} and {@link Handle}
 * carry vector behaviour through a mixin while still extending the mesh `Element` base.
 */

/** Anything indexable by number — a vector, a plain array, or a mesh vertex. */
export type VecLike = ArrayLike<number>;

const NORMALIZE_EPS = 1e-8;

abstract class BaseVector extends Array<number> {
  /**
   * Array methods such as `map` and `filter` build a plain Array rather than trying to
   * construct a fixed-size vector from a length.
   */
  static override get [Symbol.species]() {
    return Array;
  }

  abstract dot(b: VecLike): number;

  vectorLength() {
    return Math.sqrt(this.dot(this as unknown as VecLike));
  }

  vectorLengthSqr() {
    return this.dot(this as unknown as VecLike);
  }

  normalize() {
    const l = this.vectorLength();

    if (l > NORMALIZE_EPS) {
      this.mulScalar(1.0 / l);
    }

    return this;
  }

  abstract mulScalar(s: number): this;
}

export class Vector2 extends BaseVector {
  constructor(data?: VecLike) {
    super(2);

    this[0] = this[1] = 0.0;

    if (data !== undefined) {
      this.load(data);
    }
  }

  copy() {
    return new Vector2(this);
  }

  load(b: VecLike) {
    this[0] = b[0];
    this[1] = b[1];

    return this;
  }

  loadXY(x: number, y: number) {
    this[0] = x;
    this[1] = y;

    return this;
  }

  toJSON() {
    return [this[0], this[1]];
  }

  loadJSON(obj: VecLike) {
    return this.load(obj);
  }

  zero() {
    this[0] = this[1] = 0.0;

    return this;
  }

  add(b: VecLike) {
    this[0] += b[0];
    this[1] += b[1];

    return this;
  }

  sub(b: VecLike) {
    this[0] -= b[0];
    this[1] -= b[1];

    return this;
  }

  addFac(b: VecLike, fac: number) {
    this[0] += b[0] * fac;
    this[1] += b[1] * fac;

    return this;
  }

  interp(b: VecLike, t: number) {
    this[0] += (b[0] - this[0]) * t;
    this[1] += (b[1] - this[1]) * t;

    return this;
  }

  mulScalar(s: number) {
    this[0] *= s;
    this[1] *= s;

    return this;
  }

  divScalar(s: number) {
    this[0] /= s;
    this[1] /= s;

    return this;
  }

  negate() {
    this[0] = -this[0];
    this[1] = -this[1];

    return this;
  }

  dot(b: VecLike) {
    return this[0] * b[0] + this[1] * b[1];
  }

  /** 2D cross product magnitude — the z of the 3D cross. Signed. */
  cross(b: VecLike) {
    return this[0] * b[1] - this[1] * b[0];
  }

  vectorDistance(b: VecLike) {
    const d0 = this[0] - b[0];
    const d1 = this[1] - b[1];

    return Math.sqrt(d0 * d0 + d1 * d1);
  }

  vectorDistanceSqr(b: VecLike) {
    const d0 = this[0] - b[0];
    const d1 = this[1] - b[1];

    return d0 * d0 + d1 * d1;
  }

  /** Rotate about the origin by `a` radians. `axis === 1` rotates the other way. */
  rot2d(a: number, axis?: number) {
    const x = this[0];
    const y = this[1];
    const cos = Math.cos(a);
    const sin = Math.sin(a);

    if (axis === 1) {
      this[0] = x * cos + y * sin;
      this[1] = y * cos - x * sin;
    } else {
      this[0] = x * cos - y * sin;
      this[1] = y * cos + x * sin;
    }

    return this;
  }
}

export class Vector3 extends BaseVector {
  constructor(data?: VecLike) {
    super(3);

    this[0] = this[1] = this[2] = 0.0;

    if (data !== undefined) {
      this.load(data);
    }
  }

  copy() {
    return new Vector3(this);
  }

  load(b: VecLike) {
    this[0] = b[0];
    this[1] = b[1];
    this[2] = b[2] ?? 0.0;

    return this;
  }

  loadXYZ(x: number, y: number, z: number) {
    this[0] = x;
    this[1] = y;
    this[2] = z;

    return this;
  }

  toJSON() {
    return [this[0], this[1], this[2]];
  }

  loadJSON(obj: VecLike) {
    return this.load(obj);
  }

  zero() {
    this[0] = this[1] = this[2] = 0.0;

    return this;
  }

  add(b: VecLike) {
    this[0] += b[0];
    this[1] += b[1];
    this[2] += b[2];

    return this;
  }

  sub(b: VecLike) {
    this[0] -= b[0];
    this[1] -= b[1];
    this[2] -= b[2];

    return this;
  }

  addFac(b: VecLike, fac: number) {
    this[0] += b[0] * fac;
    this[1] += b[1] * fac;
    this[2] += b[2] * fac;

    return this;
  }

  interp(b: VecLike, t: number) {
    this[0] += (b[0] - this[0]) * t;
    this[1] += (b[1] - this[1]) * t;
    this[2] += (b[2] - this[2]) * t;

    return this;
  }

  mulScalar(s: number) {
    this[0] *= s;
    this[1] *= s;
    this[2] *= s;

    return this;
  }

  divScalar(s: number) {
    this[0] /= s;
    this[1] /= s;
    this[2] /= s;

    return this;
  }

  negate() {
    this[0] = -this[0];
    this[1] = -this[1];
    this[2] = -this[2];

    return this;
  }

  dot(b: VecLike) {
    return this[0] * b[0] + this[1] * b[1] + this[2] * b[2];
  }

  cross(b: VecLike) {
    const x = this[1] * b[2] - this[2] * b[1];
    const y = this[2] * b[0] - this[0] * b[2];
    const z = this[0] * b[1] - this[1] * b[0];

    this[0] = x;
    this[1] = y;
    this[2] = z;

    return this;
  }

  vectorDistance(b: VecLike) {
    const d0 = this[0] - b[0];
    const d1 = this[1] - b[1];
    const d2 = this[2] - (b[2] ?? 0.0);

    return Math.sqrt(d0 * d0 + d1 * d1 + d2 * d2);
  }

  vectorDistanceSqr(b: VecLike) {
    const d0 = this[0] - b[0];
    const d1 = this[1] - b[1];
    const d2 = this[2] - (b[2] ?? 0.0);

    return d0 * d0 + d1 * d1 + d2 * d2;
  }

  rot2d(a: number, axis?: number) {
    const x = this[0];
    const y = this[1];
    const cos = Math.cos(a);
    const sin = Math.sin(a);

    if (axis === 1) {
      this[0] = x * cos + y * sin;
      this[1] = y * cos - x * sin;
    } else {
      this[0] = x * cos - y * sin;
      this[1] = y * cos + x * sin;
    }

    return this;
  }
}

export class Vector4 extends BaseVector {
  constructor(data?: VecLike) {
    super(4);

    this[0] = this[1] = this[2] = this[3] = 0.0;

    if (data !== undefined) {
      this.load(data);
    }
  }

  copy() {
    return new Vector4(this);
  }

  load(b: VecLike) {
    this[0] = b[0];
    this[1] = b[1];
    this[2] = b[2] ?? 0.0;
    this[3] = b[3] ?? 0.0;

    return this;
  }

  toJSON() {
    return [this[0], this[1], this[2], this[3]];
  }

  loadJSON(obj: VecLike) {
    return this.load(obj);
  }

  zero() {
    this[0] = this[1] = this[2] = this[3] = 0.0;

    return this;
  }

  add(b: VecLike) {
    this[0] += b[0];
    this[1] += b[1];
    this[2] += b[2];
    this[3] += b[3];

    return this;
  }

  sub(b: VecLike) {
    this[0] -= b[0];
    this[1] -= b[1];
    this[2] -= b[2];
    this[3] -= b[3];

    return this;
  }

  addFac(b: VecLike, fac: number) {
    this[0] += b[0] * fac;
    this[1] += b[1] * fac;
    this[2] += b[2] * fac;
    this[3] += b[3] * fac;

    return this;
  }

  interp(b: VecLike, t: number) {
    this[0] += (b[0] - this[0]) * t;
    this[1] += (b[1] - this[1]) * t;
    this[2] += (b[2] - this[2]) * t;
    this[3] += (b[3] - this[3]) * t;

    return this;
  }

  mulScalar(s: number) {
    this[0] *= s;
    this[1] *= s;
    this[2] *= s;
    this[3] *= s;

    return this;
  }

  divScalar(s: number) {
    this[0] /= s;
    this[1] /= s;
    this[2] /= s;
    this[3] /= s;

    return this;
  }

  dot(b: VecLike) {
    return this[0] * b[0] + this[1] * b[1] + this[2] * b[2] + this[3] * b[3];
  }

  vectorDistance(b: VecLike) {
    const d0 = this[0] - b[0];
    const d1 = this[1] - b[1];
    const d2 = this[2] - (b[2] ?? 0.0);
    const d3 = this[3] - (b[3] ?? 0.0);

    return Math.sqrt(d0 * d0 + d1 * d1 + d2 * d2 + d3 * d3);
  }
}
