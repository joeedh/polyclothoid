/**
 * The banded `LDLᵀ` / KKT layer of `docs/plans/spower-solver.md` §5.
 *
 * Everything here is checked against dense linear algebra written independently in the
 * test — Gaussian elimination with partial pivoting — rather than against the band solver's
 * own inverse. The interesting claim is not "it solves systems" but the specific one §5
 * rests on: that unpivoted `LDLᵀ` on a *quasi-definite* matrix is stable with no pivoting
 * at all, so the static band ordering survives.
 */
import assert from "node:assert/strict";
import { describe, it } from "node:test";

import { BandedLDL, BandedSymmetric, ruizEquilibrate, solveKKT } from "../src/math/banded.js";

/** A deterministic generator — `Math.random()` would make a failure impossible to rerun. */
function lcg(seed: number) {
  let state = seed >>> 0;

  return () => {
    state = (Math.imul(state, 1664525) + 1013904223) >>> 0;

    return state / 0x100000000;
  };
}

function toDense(a: BandedSymmetric) {
  const m: number[][] = [];

  for (let i = 0; i < a.n; i++) {
    m.push([]);

    for (let j = 0; j < a.n; j++) {
      m[i].push(a.get(i, j));
    }
  }

  return m;
}

/** Dense solve by Gaussian elimination with partial pivoting — the independent reference. */
function denseSolve(m: number[][], b: number[]) {
  const n = b.length;
  const a = m.map((row, i) => [...row, b[i]]);

  for (let col = 0; col < n; col++) {
    let best = col;

    for (let i = col + 1; i < n; i++) {
      if (Math.abs(a[i][col]) > Math.abs(a[best][col])) {
        best = i;
      }
    }

    const t = a[col];
    a[col] = a[best];
    a[best] = t;

    for (let i = col + 1; i < n; i++) {
      const f = a[i][col] / a[col][col];

      for (let j = col; j <= n; j++) {
        a[i][j] -= f * a[col][j];
      }
    }
  }

  const x = new Array<number>(n).fill(0);

  for (let i = n - 1; i >= 0; i--) {
    let sum = a[i][n];

    for (let j = i + 1; j < n; j++) {
      sum -= a[i][j] * x[j];
    }

    x[i] = sum / a[i][i];
  }

  return x;
}

function maxAbsDiff(a: ArrayLike<number>, b: ArrayLike<number>) {
  let worst = 0.0;

  for (let i = 0; i < a.length; i++) {
    worst = Math.max(worst, Math.abs(a[i] - b[i]));
  }

  return worst;
}

/** Diagonally dominant, hence positive definite, with a random band. */
function randomSPD(n: number, bandwidth: number, rand: () => number) {
  const a = new BandedSymmetric(n, bandwidth);

  for (let j = 0; j < n; j++) {
    let sum = 0.0;

    for (let i = j + 1; i <= Math.min(n - 1, j + bandwidth); i++) {
      const v = rand() * 2.0 - 1.0;

      a.set(i, j, v);
      sum += Math.abs(v);
    }

    a.set(j, j, sum + 1.0);
  }

  // The dominance check above only counted each off-diagonal against one of its two rows.
  for (let i = 0; i < n; i++) {
    let sum = 0.0;

    for (let j = 0; j < n; j++) {
      if (i !== j) {
        sum += Math.abs(a.get(i, j));
      }
    }

    a.set(i, i, sum + 1.0);
  }

  return a;
}

/**
 * A quasi-definite KKT matrix in the interleaved ordering §5 describes: blocks of `dof`
 * primal variables separated by single multiplier rows, with the multiplier at vertex `i`
 * coupling the blocks on either side of it.
 */
function randomKKT(blocks: number, dof: number, rand: () => number) {
  const stride = dof + 1;
  const n = blocks * stride - 1;
  const bandwidth = 2 * stride - 1;

  const a = new BandedSymmetric(n, bandwidth);
  const constraintRows: number[] = [];

  for (let bi = 0; bi < blocks; bi++) {
    const base = bi * stride;

    for (let i = 0; i < dof; i++) {
      for (let j = 0; j <= i; j++) {
        a.add(base + i, base + j, (rand() * 2.0 - 1.0) * 0.25);
      }

      a.add(base + i, base + i, 2.0 + rand());

      if (bi + 1 < blocks) {
        for (let j = 0; j < dof; j++) {
          a.add(base + stride + i, base + j, (rand() * 2.0 - 1.0) * 0.2);
        }
      }
    }

    const lambda = base + dof;

    if (lambda < n) {
      constraintRows.push(lambda);

      for (let i = 0; i < dof; i++) {
        a.add(lambda, base + i, rand() * 2.0 - 1.0);
        a.add(lambda, base + stride + i, rand() * 2.0 - 1.0);
      }
    }
  }

  return { a, constraintRows };
}

describe("banded symmetric algebra", () => {
  it("multiplies as a symmetric matrix from a lower-band store", () => {
    const rand = lcg(7);
    const a = randomSPD(9, 3, rand);
    const dense = toDense(a);

    const x = Float64Array.from({ length: 9 }, () => rand() * 2.0 - 1.0);
    const got = a.apply(x);

    for (let i = 0; i < 9; i++) {
      let sum = 0.0;

      for (let j = 0; j < 9; j++) {
        sum += dense[i][j] * x[j];
      }

      assert.ok(Math.abs(sum - got[i]) < 1e-12, `row ${i}: ${sum} vs ${got[i]}`);
    }
  });

  it("refuses to store outside the band", () => {
    const a = new BandedSymmetric(6, 2);

    assert.throws(() => a.set(5, 0, 1.0), /outside a bandwidth/);
    assert.equal(a.get(5, 0), 0.0);
  });

  it("factors and solves a positive definite band against a dense reference", () => {
    const rand = lcg(11);

    for (const [n, bandwidth] of [
      [8, 1],
      [12, 3],
      [20, 5],
    ]) {
      const a = randomSPD(n, bandwidth, rand);
      const dense = toDense(a);
      const b = Array.from({ length: n }, () => rand() * 2.0 - 1.0);

      const ldl = new BandedLDL(a.copyTo());
      assert.ok(ldl.factor(), `n=${n}, b=${bandwidth}: factorization failed`);

      const got = ldl.solveInPlace(Float64Array.from(b));
      const want = denseSolve(dense, b);

      assert.ok(maxAbsDiff(got, want) < 1e-9, `n=${n}, b=${bandwidth}: ${maxAbsDiff(got, want)}`);
    }
  });

  it("has all-positive D on a positive definite matrix", () => {
    const a = randomSPD(14, 4, lcg(3));
    const ldl = new BandedLDL(a);

    assert.ok(ldl.factor());

    for (let i = 0; i < ldl.d.length; i++) {
      assert.ok(ldl.d[i] > 0.0, `d[${i}] = ${ldl.d[i]} should be positive`);
    }
  });
});

describe("Ruiz equilibration", () => {
  it("drives every row norm to 1 and preserves the solution under rescaling", (t) => {
    const rand = lcg(29);
    const a = randomSPD(16, 3, rand);

    // Wreck the scaling first, by 12 orders of magnitude, so there is something to fix.
    for (let i = 0; i < a.n; i++) {
      for (let j = Math.max(0, i - 3); j <= i; j++) {
        a.set(i, j, a.get(i, j) * Math.pow(10.0, (i % 5) + (j % 5) - 4));
      }
    }

    const before = toDense(a);
    const s = ruizEquilibrate(a.copyTo(a), 20);

    let worst = 0.0;

    for (let i = 0; i < a.n; i++) {
      let rowMax = 0.0;

      for (let j = 0; j < a.n; j++) {
        rowMax = Math.max(rowMax, Math.abs(a.get(i, j)));
      }

      worst = Math.max(worst, Math.abs(rowMax - 1.0));
    }

    t.diagnostic(`worst row-norm deviation from 1: ${worst.toExponential(2)}`);
    assert.ok(worst < 1e-6, `row norms should equilibrate to 1, worst deviation ${worst}`);

    // S A S must be exactly the scaled original, entry by entry.
    for (let i = 0; i < a.n; i++) {
      for (let j = 0; j < a.n; j++) {
        const want = before[i][j] * s[i] * s[j];
        const got = a.get(i, j);

        assert.ok(Math.abs(got - want) < 1e-12 * Math.max(1.0, Math.abs(want)), `(${i},${j}): ${got} vs ${want}`);
      }
    }
  });
});

describe("quasi-definite KKT solve", () => {
  it("matches a dense pivoted solve on an indefinite system", (t) => {
    const rand = lcg(101);

    for (const [blocks, dof] of [
      [4, 1],
      [6, 3],
      [9, 4],
    ]) {
      const { a, constraintRows } = randomKKT(blocks, dof, rand);
      const dense = toDense(a);
      const b = Array.from({ length: a.n }, () => rand() * 2.0 - 1.0);

      const got = solveKKT(a, Float64Array.from(b), constraintRows);
      const want = denseSolve(dense, b);

      const err = maxAbsDiff(got.x, want);

      t.diagnostic(
        `blocks=${blocks} dof=${dof} n=${a.n}   error ${err.toExponential(2)}   residual ${got.residual.toExponential(2)}`
      );

      assert.ok(got.ok, `blocks=${blocks} dof=${dof}: factorization failed`);
      assert.ok(err < 1e-8, `blocks=${blocks} dof=${dof}: ${err}`);
    }
  });

  it("has exactly one negative D entry per constraint row", () => {
    const { a, constraintRows } = randomKKT(6, 3, lcg(5));

    const shifted = a.copyTo();

    for (const i of constraintRows) {
      shifted.add(i, i, -1e-8);
    }

    const ldl = new BandedLDL(shifted);
    assert.ok(ldl.factor());

    let negatives = 0;

    for (const d of ldl.d) {
      if (d < 0.0) {
        negatives++;
      }
    }

    // Sylvester's law of inertia: the signature is a property of the matrix, not the
    // factorization, so this is what makes "quasi-definite" checkable rather than asserted.
    assert.equal(negatives, constraintRows.length);
  });

  /**
   * Worth stating precisely, because §5 overstates it slightly: a structurally zero (2,2)
   * diagonal does *not* make unpivoted `LDLᵀ` break down generically. The pivot at a
   * multiplier row is `0 − Σ L²d`, and the accumulated update is nonzero for a generic
   * constraint row — the random systems above factor unshifted with a smallest `|d|` of
   * 1.3e-1. What the shift buys is a *guarantee* that holds under any symmetric ordering,
   * and it becomes load-bearing rather than insurance exactly when a constraint row goes
   * rank deficient. That is not a hypothetical here: §5 notes that with the transform
   * frozen, `e2`'s degrees of freedom drop out of its G1 row entirely.
   */
  it("is what keeps a rank-deficient constraint row factorable", () => {
    const { a, constraintRows } = randomKKT(6, 3, lcg(17));
    const dead = constraintRows[2];

    for (let j = Math.max(0, dead - a.bandwidth); j <= Math.min(a.n - 1, dead + a.bandwidth); j++) {
      a.set(dead, j, 0.0);
    }

    assert.equal(new BandedLDL(a.copyTo()).factor(), false, "a zero constraint row should give a zero pivot");

    const shifted = solveKKT(a, new Float64Array(a.n).fill(1.0), constraintRows);

    assert.ok(shifted.ok, "the quasi-definite shift should make it factorable again");
  });

  it("recovers with refinement the accuracy the shift costs", (t) => {
    const { a, constraintRows } = randomKKT(6, 3, lcg(17));
    const b = Array.from({ length: a.n }, (_, i) => Math.sin(i * 1.7));
    const want = denseSolve(toDense(a), b);

    const coarse = solveKKT(a, Float64Array.from(b), constraintRows, { delta: 1e-4, refinement: 0 });
    const refined = solveKKT(a, Float64Array.from(b), constraintRows, { delta: 1e-4, refinement: 3 });

    t.diagnostic(`delta=1e-4  no refinement ${maxAbsDiff(coarse.x, want).toExponential(2)}`);
    t.diagnostic(`delta=1e-4  3 refinements ${maxAbsDiff(refined.x, want).toExponential(2)}`);

    assert.ok(
      maxAbsDiff(refined.x, want) < maxAbsDiff(coarse.x, want) * 1e-3,
      "refinement should recover at least three digits lost to a deliberately large delta"
    );
    assert.ok(maxAbsDiff(refined.x, want) < 1e-9);
  });
});
