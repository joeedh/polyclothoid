/**
 * Small utilities vendored from path.ux, plus the global patches it used to install.
 *
 * path.ux defined `Math.fract` and `Array.prototype.remove` on the built-ins. A library
 * has no business doing that, so they are plain exports here.
 */

/**
 * A fixed-size ring of reusable objects.
 *
 * Curve `evaluate`/`derivative` calls return vectors borrowed from one of these rather
 * than allocating, because the stroker calls them once per brush dab. **A borrowed value
 * is only valid until the ring wraps** — roughly `size` further calls. Callers that need
 * to keep a value must copy it.
 *
 * The original extended `Array`; this holds its storage in a field instead, which is
 * what lets the ring be generic in its element type.
 */
export class CacheRing<T> {
  private items: T[] = [];
  private cur = 0;

  constructor(factory: () => T, size: number) {
    for (let i = 0; i < size; i++) {
      this.items.push(factory());
    }
  }

  static fromConstructor<T>(cls: new () => T, size: number) {
    return new CacheRing<T>(() => new cls(), size);
  }

  get length() {
    return this.items.length;
  }

  next(): T {
    const ret = this.items[this.cur];
    this.cur = (this.cur + 1) % this.items.length;

    return ret;
  }
}

/** Monotonic element id source. Mesh elements are keyed by these. */
export class IDGen {
  cur = 1;

  next() {
    return this.cur++;
  }

  toJSON() {
    return { cur: this.cur };
  }

  loadJSON(obj: { cur: number }) {
    this.cur = obj.cur;

    return this;
  }

  static fromJSON(obj: { cur: number }) {
    return new IDGen().loadJSON(obj);
  }
}

/** Fractional part, always in `[0, 1)` — `fract(-0.25)` is `0.75`, not `-0.25`. */
export function fract(f: number) {
  return f - Math.floor(f);
}

export function clamp(f: number, min: number, max: number) {
  return Math.min(Math.max(f, min), max);
}

export function time_ms() {
  return globalThis.performance?.now() ?? Date.now();
}

/** Remove the first occurrence of `item`, in place. Returns whether it was found. */
export function listRemove<T>(list: T[], item: T) {
  const i = list.indexOf(item);

  if (i < 0) {
    return false;
  }

  list.splice(i, 1);

  return true;
}

const binomialCache = new Map<number, number>();

/**
 * Binomial coefficient `n choose i`.
 *
 * path.ux's version recursed into an undefined `bin()` and threw for every input that
 * was not a trivial base case. This one is correct.
 */
export function binomial(n: number, i: number): number {
  if (i > n || i < 0) {
    throw new Error(`bad call to binomial(${n}, ${i})`);
  }

  if (i === 0 || i === n) {
    return 1;
  }

  const key = n * 4096 + i;
  const hit = binomialCache.get(key);

  if (hit !== undefined) {
    return hit;
  }

  const ret = binomial(n - 1, i - 1) + binomial(n - 1, i);
  binomialCache.set(key, ret);

  return ret;
}
