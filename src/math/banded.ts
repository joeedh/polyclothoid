/**
 * Symmetric banded `LDLᵀ` with quasi-definite regularization, Ruiz equilibration and
 * iterative refinement — the linear algebra `docs/plans/spower-solver.md` §5 calls for.
 *
 * The KKT systems this exists to solve are *indefinite*: the (2,2) block is structurally
 * zero, which puts zeros on the diagonal exactly where unpivoted `LDLᵀ` breaks down.
 * Bunch-Kaufman would fix that but chooses pivots dynamically, which destroys both the
 * static ordering and the band — the two properties that make the solve `O(n·b²)` instead
 * of `O(n³)`. The standard alternative, and the one used here, is Vanderbei's:
 * perturb the (2,2) block to `−δI` so the matrix becomes *quasi-definite* (symmetric, with
 * a positive definite (1,1) block and a negative definite (2,2) block), for which
 * unpivoted `LDLᵀ` exists and is stable under **any** symmetric permutation. The accuracy
 * lost to `δ` comes back from iterative refinement against the unperturbed matrix.
 */

/**
 * A symmetric matrix stored as its lower band.
 *
 * Column-major within each column: entry `(i, j)` with `j <= i <= j + bandwidth` lives at
 * `data[j * (bandwidth + 1) + (i - j)]`. Everything outside the band is structurally zero
 * and is neither stored nor touched.
 */
export class BandedSymmetric {
  data: Float64Array;

  constructor(
    public n: number,
    public bandwidth: number
  ) {
    this.data = new Float64Array(n * (bandwidth + 1));
  }

  /** Index of `(i, j)` in {@link data}, or `-1` if the entry lies outside the band. */
  index(i: number, j: number) {
    if (i < j) {
      const t = i;
      i = j;
      j = t;
    }

    return i - j > this.bandwidth ? -1 : j * (this.bandwidth + 1) + (i - j);
  }

  get(i: number, j: number) {
    const at = this.index(i, j);

    return at === -1 ? 0.0 : this.data[at];
  }

  set(i: number, j: number, v: number) {
    const at = this.index(i, j);

    if (at === -1) {
      throw new Error(`(${i}, ${j}) lies outside a bandwidth of ${this.bandwidth}`);
    }

    this.data[at] = v;
  }

  /** Accumulate into `(i, j)`. Assembly is a sum of element contributions, so this is the hot path. */
  add(i: number, j: number, v: number) {
    const at = this.index(i, j);

    if (at === -1) {
      throw new Error(`(${i}, ${j}) lies outside a bandwidth of ${this.bandwidth}`);
    }

    this.data[at] += v;
  }

  zero() {
    this.data.fill(0.0);

    return this;
  }

  copyTo(out = new BandedSymmetric(this.n, this.bandwidth)) {
    out.data.set(this.data);

    return out;
  }

  /*
    `out` is annotated rather than inferred from its default: inference would pin it to
    `Float64Array<ArrayBuffer>` and reject any caller holding a plainly declared one.
  */
  /** `out = A·x`, reading the stored lower band as symmetric. */
  apply(x: Float64Array, out: Float64Array = new Float64Array(this.n)) {
    const { n, bandwidth, data } = this;

    out.fill(0.0);

    for (let j = 0; j < n; j++) {
      const base = j * (bandwidth + 1);
      const hi = Math.min(n - 1, j + bandwidth);

      out[j] += data[base] * x[j];

      for (let i = j + 1; i <= hi; i++) {
        const v = data[base + (i - j)];

        out[i] += v * x[j];
        out[j] += v * x[i];
      }
    }

    return out;
  }
}

/**
 * An `LDLᵀ` factorization in place of the band it was built from.
 *
 * `L` is unit lower triangular and shares the band; `D` is diagonal and may carry negative
 * entries, which is the whole point — a quasi-definite matrix has as many negative `D`
 * entries as it has constraint rows.
 */
export class BandedLDL {
  d: Float64Array;

  constructor(public a: BandedSymmetric) {
    this.d = new Float64Array(a.n);
  }

  /**
   * Factor in place. Overwrites `a`'s band with the strict lower triangle of `L`.
   *
   * Returns `false` if a pivot underflowed, which for a genuinely quasi-definite matrix
   * means `δ` was too small rather than that the matrix was singular.
   */
  factor() {
    const { n, bandwidth, data } = this.a;
    const d = this.d;

    for (let j = 0; j < n; j++) {
      const base = j * (bandwidth + 1);
      const lo = Math.max(0, j - bandwidth);

      let dj = data[base];

      for (let k = lo; k < j; k++) {
        const ljk = data[k * (bandwidth + 1) + (j - k)];

        dj -= ljk * ljk * d[k];
      }

      if (Math.abs(dj) < 1e-300) {
        return false;
      }

      d[j] = dj;
      data[base] = 1.0;

      const hi = Math.min(n - 1, j + bandwidth);

      for (let i = j + 1; i <= hi; i++) {
        let sum = data[base + (i - j)];

        // k >= i - bandwidth is the binding lower limit, and it implies k >= j - bandwidth.
        for (let k = Math.max(0, i - bandwidth); k < j; k++) {
          const kb = k * (bandwidth + 1);

          sum -= data[kb + (i - k)] * data[kb + (j - k)] * d[k];
        }

        data[base + (i - j)] = sum / dj;
      }
    }

    return true;
  }

  /** Solve `A y = b` in place on `b`, using the factorization. */
  solveInPlace(b: Float64Array) {
    const { n, bandwidth, data } = this.a;
    const d = this.d;

    for (let i = 0; i < n; i++) {
      let sum = b[i];

      for (let k = Math.max(0, i - bandwidth); k < i; k++) {
        sum -= data[k * (bandwidth + 1) + (i - k)] * b[k];
      }

      b[i] = sum;
    }

    for (let i = 0; i < n; i++) {
      b[i] /= d[i];
    }

    for (let i = n - 1; i >= 0; i--) {
      let sum = b[i];
      const hi = Math.min(n - 1, i + bandwidth);
      const base = i * (bandwidth + 1);

      for (let k = i + 1; k <= hi; k++) {
        sum -= data[base + (k - i)] * b[k];
      }

      b[i] = sum;
    }

    return b;
  }
}

/**
 * Ruiz symmetric diagonal equilibration: find `s` such that every row of `S A S` has
 * infinity norm near 1, with `S = diag(s)`.
 *
 * §3 needs this because no single reference length `Rᵥ` can equalize column scales when
 * incident chord lengths are disparate — the `Rᵥⁿ/n!` form gets the `n`-dependence right
 * and equilibration cleans up the rest. It is part of the solver, not a polish step.
 *
 * The scaling is applied to `a` in place and `s` is returned; a system `A x = b` becomes
 * `(S A S)(S⁻¹x) = S b`, so callers scale the right-hand side by `s` and the solution by
 * `s` again.
 */
export function ruizEquilibrate(a: BandedSymmetric, iterations = 8, out = new Float64Array(a.n)) {
  const { n, bandwidth, data } = a;

  out.fill(1.0);

  const row = new Float64Array(n);

  for (let iter = 0; iter < iterations; iter++) {
    row.fill(0.0);

    for (let j = 0; j < n; j++) {
      const base = j * (bandwidth + 1);
      const hi = Math.min(n - 1, j + bandwidth);

      for (let i = j; i <= hi; i++) {
        const v = Math.abs(data[base + (i - j)]);

        if (v > row[i]) {
          row[i] = v;
        }

        if (v > row[j]) {
          row[j] = v;
        }
      }
    }

    let done = true;

    for (let i = 0; i < n; i++) {
      // An empty row cannot be scaled into range; leave it alone rather than dividing by 0.
      const r = row[i] > 0.0 ? 1.0 / Math.sqrt(row[i]) : 1.0;

      if (Math.abs(r - 1.0) > 1e-12) {
        done = false;
      }

      row[i] = r;
      out[i] *= r;
    }

    for (let j = 0; j < n; j++) {
      const base = j * (bandwidth + 1);
      const hi = Math.min(n - 1, j + bandwidth);

      for (let i = j; i <= hi; i++) {
        data[base + (i - j)] *= row[i] * row[j];
      }
    }

    if (done) {
      break;
    }
  }

  return out;
}

export interface KKTSolveOptions {
  /**
   * Quasi-definite shift applied to the constraint block. Small enough not to perturb the
   * answer much, large enough to keep pivots away from zero; refinement removes what is
   * left.
   */
  delta: number;

  /** Iterative refinement passes against the *unshifted* matrix. */
  refinement: number;

  /** Ruiz sweeps. Zero disables equilibration. */
  equilibration: number;
}

export const defaultKKTSolveOptions: KKTSolveOptions = {
  delta        : 1e-8,
  refinement   : 2,
  equilibration: 8,
};

export interface KKTSolveResult {
  x: Float64Array;

  /** Infinity norm of the final residual `‖A x − b‖∞`, in equilibrated units. */
  residual: number;

  /** False if the factorization underflowed — `delta` is too small for this system. */
  ok: boolean;
}

/**
 * One equilibrated, shifted `LDLᵀ` of a KKT matrix, reusable across right-hand sides.
 *
 * Separate from {@link solveKKT} because §5's substructuring needs `A⁻¹B` as well as
 * `A⁻¹b` — one solve per interface column, all against the same matrix. Factoring per
 * column would put the `O(n·b²)` factorization inside a loop it does not belong in.
 *
 * `a` is not modified: the equilibration and the shift are applied to working copies, so
 * refinement has an unperturbed matrix to measure against.
 */
export class KKTFactorization {
  options: KKTSolveOptions;

  /** Equilibrated but unshifted — what refinement measures the residual against. */
  scaled: BandedSymmetric;

  /** Ruiz scaling, applied to the right-hand side going in and the solution coming out. */
  scaling: Float64Array;

  ldl: BandedLDL;

  /** False if a pivot underflowed, which means `delta` is too small for this system. */
  ok: boolean;

  private work: Float64Array;
  private rhs: Float64Array;

  constructor(a: BandedSymmetric, constraintRows: Iterable<number>, options: Partial<KKTSolveOptions> = {}) {
    const opts = { ...defaultKKTSolveOptions, ...options };
    const n = a.n;

    this.options = opts;
    this.scaled = a.copyTo();
    this.scaling =
      opts.equilibration > 0 ? ruizEquilibrate(this.scaled, opts.equilibration) : new Float64Array(n).fill(1.0);

    const shifted = this.scaled.copyTo();

    for (const i of constraintRows) {
      shifted.add(i, i, -opts.delta);
    }

    this.ldl = new BandedLDL(shifted);
    this.ok = this.ldl.factor();

    this.work = new Float64Array(n);
    this.rhs = new Float64Array(n);
  }

  /** Solve `A x = b`, refining against the unshifted matrix. `b` is not modified. */
  solve(b: Float64Array, out?: Float64Array): KKTSolveResult {
    const { scaled, scaling, ldl, options: opts } = this;
    const n = scaled.n;
    const x = out ?? new Float64Array(n);

    if (!this.ok) {
      x.fill(0.0);

      return { x, residual: Infinity, ok: false };
    }

    const rhs = this.rhs;

    for (let i = 0; i < n; i++) {
      rhs[i] = b[i] * scaling[i];
      x[i] = rhs[i];
    }

    ldl.solveInPlace(x);

    const work = this.work;
    let residual = Infinity;

    for (let pass = 0; pass <= opts.refinement; pass++) {
      scaled.apply(x, work);

      residual = 0.0;

      for (let i = 0; i < n; i++) {
        work[i] = rhs[i] - work[i];
        residual = Math.max(residual, Math.abs(work[i]));
      }

      if (pass === opts.refinement) {
        break;
      }

      ldl.solveInPlace(work);

      for (let i = 0; i < n; i++) {
        x[i] += work[i];
      }
    }

    for (let i = 0; i < n; i++) {
      x[i] *= scaling[i];
    }

    return { x, residual, ok: true };
  }
}

/**
 * Solve the symmetric indefinite system `A x = b`, where `A` is quasi-definite once the
 * rows listed in `constraintRows` are shifted by `−δ`.
 *
 * `a` is not modified. Factor once and solve many with {@link KKTFactorization} instead when
 * there is more than one right-hand side.
 */
export function solveKKT(
  a: BandedSymmetric,
  b: Float64Array,
  constraintRows: Iterable<number>,
  options: Partial<KKTSolveOptions> = {}
): KKTSolveResult {
  return new KKTFactorization(a, constraintRows, options).solve(b);
}
