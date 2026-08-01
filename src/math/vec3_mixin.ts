import { type VecLike } from "./vector.js";

/**
 * Vector behaviour for classes that cannot extend {@link Vector3} because they already
 * extend something else — {@link Vertex} and {@link Handle}, which extend the mesh
 * `Element` base.
 *
 * The original walked `Vector3`'s whole prototype chain with `Reflect.ownKeys` and copied
 * every key it found, including `Array.prototype` methods, onto a non-array object. This
 * version states the methods explicitly instead: same behaviour, checkable by the
 * compiler, and no array methods landing on things that are not arrays.
 *
 * Apply with {@link applyVec3Mixin} at class definition time, and declare the types by
 * merging {@link Vec3Mixin} into the class interface:
 *
 * ```ts
 * export class Vertex extends Element { ... }
 * export interface Vertex extends Vec3Mixin {}
 * applyVec3Mixin(Vertex);
 * ```
 */
export interface Vec3Mixin {
  /* The index signature plus `length` is what makes a mixed-in class an `ArrayLike`, so
     it can be passed anywhere a `Vector2`/`Vector3` is expected. */
  [index: number]: number;
  0: number;
  1: number;
  2: number;
  length: number;

  initVector3(): this;
  load(b: VecLike): this;
  loadXYZ(x: number, y: number, z: number): this;
  zero(): this;
  copyTo(b: { 0: number; 1: number; 2: number }): void;

  add(b: VecLike): this;
  sub(b: VecLike): this;
  addFac(b: VecLike, fac: number): this;
  interp(b: VecLike, t: number): this;
  mulScalar(s: number): this;
  divScalar(s: number): this;
  negate(): this;

  dot(b: VecLike): number;
  vectorLength(): number;
  vectorLengthSqr(): number;
  vectorDistance(b: VecLike): number;
  vectorDistanceSqr(b: VecLike): number;
  normalize(): this;
}

const NORMALIZE_EPS = 1e-8;

const impl: Vec3Mixin = {
  0     : 0.0,
  1     : 0.0,
  2     : 0.0,
  length: 3,

  initVector3(this: Vec3Mixin) {
    this.length = 3;
    this[0] = this[1] = this[2] = 0.0;

    return this;
  },

  load(this: Vec3Mixin, b: VecLike) {
    this[0] = b[0];
    this[1] = b[1];
    this[2] = b[2] ?? 0.0;

    return this;
  },

  loadXYZ(this: Vec3Mixin, x: number, y: number, z: number) {
    this[0] = x;
    this[1] = y;
    this[2] = z;

    return this;
  },

  zero(this: Vec3Mixin) {
    this[0] = this[1] = this[2] = 0.0;

    return this;
  },

  copyTo(this: Vec3Mixin, b: { 0: number; 1: number; 2: number }) {
    b[0] = this[0];
    b[1] = this[1];
    b[2] = this[2];
  },

  add(this: Vec3Mixin, b: VecLike) {
    this[0] += b[0];
    this[1] += b[1];
    this[2] += b[2];

    return this;
  },

  sub(this: Vec3Mixin, b: VecLike) {
    this[0] -= b[0];
    this[1] -= b[1];
    this[2] -= b[2];

    return this;
  },

  addFac(this: Vec3Mixin, b: VecLike, fac: number) {
    this[0] += b[0] * fac;
    this[1] += b[1] * fac;
    this[2] += b[2] * fac;

    return this;
  },

  interp(this: Vec3Mixin, b: VecLike, t: number) {
    this[0] += (b[0] - this[0]) * t;
    this[1] += (b[1] - this[1]) * t;
    this[2] += (b[2] - this[2]) * t;

    return this;
  },

  mulScalar(this: Vec3Mixin, s: number) {
    this[0] *= s;
    this[1] *= s;
    this[2] *= s;

    return this;
  },

  divScalar(this: Vec3Mixin, s: number) {
    this[0] /= s;
    this[1] /= s;
    this[2] /= s;

    return this;
  },

  negate(this: Vec3Mixin) {
    this[0] = -this[0];
    this[1] = -this[1];
    this[2] = -this[2];

    return this;
  },

  dot(this: Vec3Mixin, b: VecLike) {
    return this[0] * b[0] + this[1] * b[1] + this[2] * b[2];
  },

  vectorLength(this: Vec3Mixin) {
    return Math.sqrt(this[0] * this[0] + this[1] * this[1] + this[2] * this[2]);
  },

  vectorLengthSqr(this: Vec3Mixin) {
    return this[0] * this[0] + this[1] * this[1] + this[2] * this[2];
  },

  vectorDistance(this: Vec3Mixin, b: VecLike) {
    const d0 = this[0] - b[0];
    const d1 = this[1] - b[1];
    const d2 = this[2] - (b[2] ?? 0.0);

    return Math.sqrt(d0 * d0 + d1 * d1 + d2 * d2);
  },

  vectorDistanceSqr(this: Vec3Mixin, b: VecLike) {
    const d0 = this[0] - b[0];
    const d1 = this[1] - b[1];
    const d2 = this[2] - (b[2] ?? 0.0);

    return d0 * d0 + d1 * d1 + d2 * d2;
  },

  normalize(this: Vec3Mixin) {
    const l = this.vectorLength();

    if (l > NORMALIZE_EPS) {
      this.mulScalar(1.0 / l);
    }

    return this;
  },
};

/** Methods only — the numeric seeds and `length` are set per instance by `initVector3`. */
const METHOD_KEYS = Object.keys(impl).filter((k) => typeof impl[k as keyof Vec3Mixin] === "function");

export function applyVec3Mixin(cls: { prototype: object }) {
  for (const key of METHOD_KEYS) {
    Object.defineProperty(cls.prototype, key, {
      value       : impl[key as keyof Vec3Mixin],
      writable    : true,
      configurable: true,
      enumerable  : false,
    });
  }
}
