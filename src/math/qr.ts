/**
 * Column-pivoted Householder QR, for rank on small dense blocks.
 *
 * The one caller is the G1 rank test of `docs/plans/spower-solver.md` §8, which is a few rows
 * by `2(p+1)` columns per node group. §8 is specific about why it is *this* and not something
 * cheaper: pivot magnitude during a factorization is not a rank test — Kahan's matrix is the
 * standard counterexample — so the rank question is answered on its own block, away from the
 * banded solve, where pivoting is free to choose its own column order.
 *
 * Rank-revealing in the usual practical sense rather than the guaranteed one. Column-pivoted
 * QR can be fooled (Kahan again), but only by matrices whose leading columns hide the
 * deficiency, and a G1 block is a handful of nearly-parallel tangent rows — the case pivoting
 * is good at. Anything stronger wants an SVD, which is not worth carrying for `2 × 4` blocks.
 */

/** What {@link columnPivotedQR} measured. `diagonal` is `|R_kk|`, non-increasing. */
export interface QRRank {
  rank: number;

  /** `|R_kk|` for `k = 0 … min(rows, cols) − 1`. `diagonal[0]` estimates the largest singular value. */
  diagonal: number[];

  /** Column order pivoting chose, longest-remaining-norm first. */
  pivots: number[];
}

/**
 * Factor `a` (row-major `rows × cols`, overwritten) and count pivots above `tol · |R₀₀|`.
 *
 * The threshold is relative by construction, so the answer is invariant under a uniform
 * scaling of the block — which matters here, because the rows are energy-weighted and their
 * absolute size says nothing about whether they are independent. An all-zero block has rank
 * `0`, not rank `1` with a zero pivot: `|R₀₀| = 0` makes every comparison false.
 */
export function columnPivotedQR(a: Float64Array, rows: number, cols: number, tol = 1e-10): QRRank {
  const pivots = Array.from({ length: cols }, (_, i) => i);
  const norms = new Float64Array(cols);
  const diagonal: number[] = [];
  const v = new Float64Array(rows);

  for (let c = 0; c < cols; c++) {
    let sum = 0.0;

    for (let r = 0; r < rows; r++) {
      sum += a[r * cols + c] * a[r * cols + c];
    }

    norms[c] = sum;
  }

  const steps = Math.min(rows, cols);

  for (let k = 0; k < steps; k++) {
    let best = k;

    for (let c = k + 1; c < cols; c++) {
      if (norms[c] > norms[best]) {
        best = c;
      }
    }

    if (best !== k) {
      for (let r = 0; r < rows; r++) {
        const t = a[r * cols + k];

        a[r * cols + k] = a[r * cols + best];
        a[r * cols + best] = t;
      }

      [pivots[k], pivots[best]] = [pivots[best], pivots[k]];
      [norms[k], norms[best]] = [norms[best], norms[k]];
    }

    let alpha = 0.0;

    for (let r = k; r < rows; r++) {
      v[r] = a[r * cols + k];
      alpha += v[r] * v[r];
    }

    alpha = Math.sqrt(alpha);
    diagonal.push(alpha);

    if (!(alpha > 0.0)) {
      continue;
    }

    // Sign chosen away from v[k] so the reflector never cancels catastrophically.
    const beta = v[k] >= 0.0 ? -alpha : alpha;

    v[k] -= beta;

    let vv = 0.0;

    for (let r = k; r < rows; r++) {
      vv += v[r] * v[r];
    }

    if (vv > 0.0) {
      for (let c = k; c < cols; c++) {
        let dot = 0.0;

        for (let r = k; r < rows; r++) {
          dot += v[r] * a[r * cols + c];
        }

        dot = (2.0 * dot) / vv;

        for (let r = k; r < rows; r++) {
          a[r * cols + c] -= dot * v[r];
        }
      }
    }

    for (let c = k + 1; c < cols; c++) {
      norms[c] -= a[k * cols + c] * a[k * cols + c];
      norms[c] = Math.max(norms[c], 0.0);
    }
  }

  const limit = diagonal.length === 0 ? 0.0 : diagonal[0] * tol;
  let rank = 0;

  for (const d of diagonal) {
    if (d > limit && d > 0.0) {
      rank++;
    }
  }

  return { rank, diagonal, pivots };
}

/**
 * `σ_min / σ_max` as the pivoted-QR diagonal sees it — the number §8's rank test compares.
 *
 * Not the true singular-value ratio; `|R_kk|` bounds it from above. It is what the threshold
 * is calibrated against, and it is `0` when the block is short of full column rank.
 */
export function rankRatio(a: Float64Array, rows: number, cols: number) {
  const { diagonal } = columnPivotedQR(a, rows, cols);

  if (diagonal.length === 0 || !(diagonal[0] > 0.0)) {
    return 0.0;
  }

  return diagonal[diagonal.length - 1] / diagonal[0];
}
