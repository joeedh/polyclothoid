/**
 * A small dense `LU` with partial pivoting, for the Schur complement of
 * `docs/plans/spower-solver.md` §5.
 *
 * Everything else in this directory is banded because it is solving over a chain. The
 * interface between chains is not a chain — it is a handful of node DOF and multipliers with
 * no useful ordering — so it gets the general method instead, on the grounds that it is
 * `O(#shared groups × (p+1))` on a side and cubing a small number is free.
 *
 * Partial pivoting rather than the quasi-definite trick `banded.ts` uses: a Schur complement
 * is not quasi-definite even when both of its parents are, and the multiplier rows arrive
 * with structurally zero diagonals, which is exactly the case a pivot search handles and an
 * unpivoted factorization does not.
 */

/** `PA = LU` with `L` unit lower triangular, both stored in one `n × n` row-major array. */
export class DenseLU {
  lu: Float64Array;
  pivot: Int32Array;

  /** False if a column was numerically empty, which for a Schur complement means singular. */
  ok = true;

  constructor(
    a: ArrayLike<number>,
    public n: number
  ) {
    const lu = new Float64Array(n * n);

    for (let i = 0; i < n * n; i++) {
      lu[i] = a[i];
    }

    this.lu = lu;
    this.pivot = new Int32Array(n);

    for (let k = 0; k < n; k++) {
      let best = k;
      let mag = Math.abs(lu[k * n + k]);

      for (let i = k + 1; i < n; i++) {
        const v = Math.abs(lu[i * n + k]);

        if (v > mag) {
          mag = v;
          best = i;
        }
      }

      this.pivot[k] = best;

      if (mag < 1e-300) {
        this.ok = false;

        return;
      }

      if (best !== k) {
        for (let c = 0; c < n; c++) {
          const t = lu[k * n + c];

          lu[k * n + c] = lu[best * n + c];
          lu[best * n + c] = t;
        }
      }

      const d = lu[k * n + k];

      for (let i = k + 1; i < n; i++) {
        const f = (lu[i * n + k] /= d);

        if (f === 0.0) {
          continue;
        }

        for (let c = k + 1; c < n; c++) {
          lu[i * n + c] -= f * lu[k * n + c];
        }
      }
    }
  }

  /** Solve `A x = b` in place on `b`. Leaves `b` untouched when the factorization failed. */
  solveInPlace(b: Float64Array) {
    const { lu, n, pivot } = this;

    if (!this.ok) {
      return b;
    }

    for (let k = 0; k < n; k++) {
      const p = pivot[k];

      if (p !== k) {
        const t = b[k];

        b[k] = b[p];
        b[p] = t;
      }

      for (let i = k + 1; i < n; i++) {
        b[i] -= lu[i * n + k] * b[k];
      }
    }

    for (let i = n - 1; i >= 0; i--) {
      let sum = b[i];

      for (let c = i + 1; c < n; c++) {
        sum -= lu[i * n + c] * b[c];
      }

      b[i] = sum / lu[i * n + i];
    }

    return b;
  }
}
