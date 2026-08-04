/**
 * The polynomial s-power basis of Sánchez-Reyes.
 *
 * Background and derivations are in `docs/research/spower.md`; the design that uses this
 * is `docs/plans/spower-solver.md`. This module is pure scalar math — no curves, no mesh,
 * no solver — and is the Phase 0b kernel of that plan.
 *
 * **Naming.** `docs/research/spower.md` and both papers call the symmetric parameter `s`.
 * This codebase's `s` is arclength (`CLAUDE.md` invariant 2), so the symmetric parameter
 * is `sym` here and the normalized segment parameter is `u`:
 *
 * ```
 * sym = (1 - u) * u        sym in [0, 1/4], max at u = 1/2
 * ```
 *
 * **Layout.** An order-`p` polynomial has degree `2p + 1` and `2p + 2` coefficients,
 * arranged as `p + 1` *pairs*. Pairs are stored interleaved,
 * `[a_0^0, a_0^1, a_1^0, a_1^1, ...]`, so pair `k` lives at `2k` and `2k + 1`:
 *
 * ```
 * a(u) = sum_k [ (1 - u) * a_k^0 + u * a_k^1 ] * sym^k
 * ```
 *
 * Every function takes the coefficient count `len` explicitly rather than reading
 * `a.length`, so a longer backing array can be used as scratch — the same trick
 * `Clothoid.ks` plays.
 */
import { binomial } from "../math/index.js";

/**
 * Highest order this module will handle, and the size of its scratch buffers.
 *
 * Far above anything the solver ships — `spower-solver.md` §9 puts the useful ladder at
 * `p = 0..2`. The cap exists so the scratch buffers can be module-level, which is what
 * keeps {@link sPowerIntegral} allocation-free inside a quadrature loop.
 */
export const MAX_SPOWER_ORDER = 12;

/** Coefficient count of an order-`p` polynomial: `2p + 2`. */
export function sPowerLength(p: number) {
  return 2 * p + 2;
}

/** Inverse of {@link sPowerLength}. */
export function sPowerOrder(len: number) {
  return len / 2 - 1;
}

/**
 * Evaluate at `u` by Horner in `sym`, with a linear blend per term.
 *
 * `sym <= 1/4` over the whole interval, so successive terms shrink geometrically — the
 * conditioning argument for the basis (`spower.md` §4). Note that this is about
 * *evaluation*; the Gram matrices below condition the other way (`spower-solver.md` §13).
 */
export function evalSPower(a: ArrayLike<number>, len: number, u: number) {
  const sym = (1.0 - u) * u;
  const u1 = 1.0 - u;

  let acc = 0.0;

  for (let k = len - 2; k >= 0; k -= 2) {
    acc = acc * sym + (u1 * a[k] + u * a[k + 1]);
  }

  return acc;
}

/**
 * Derivative, in place of a new array: `spower.md` §6's rule, exact.
 *
 * ```
 * mean_c_k  = (2k+1) * diff_a_k - (k+1) * diff_a_{k+1} / 2
 * diff_c_k  = -2(k+1) * mean_a_{k+1}
 * ```
 *
 * The derivative of a degree-`2p+1` polynomial has degree `2p`, which sits inside the
 * order-`p` basis with `diff_c_p = 0` — so `out` has the same length as `a`, and the
 * top pair comes out symmetric. `out` may alias `a`.
 */
export function differentiateSPower(a: ArrayLike<number>, len: number, out: Float64Array) {
  const p = sPowerOrder(len);

  for (let k = 0; k <= p; k++) {
    const diffA = a[2 * k + 1] - a[2 * k];

    const meanA1 = k < p ? (a[2 * k + 3] + a[2 * k + 2]) * 0.5 : 0.0;
    const diffA1 = k < p ? a[2 * k + 3] - a[2 * k + 2] : 0.0;

    const meanC = (2 * k + 1) * diffA - (k + 1) * diffA1 * 0.5;
    const diffC = -2 * (k + 1) * meanA1;

    out[2 * k] = meanC - diffC * 0.5;
    out[2 * k + 1] = meanC + diffC * 0.5;
  }

  return out;
}

const meanScratch = new Float64Array(MAX_SPOWER_ORDER + 2);
const diffScratch = new Float64Array(MAX_SPOWER_ORDER + 2);
const antiScratch = new Float64Array(sPowerLength(MAX_SPOWER_ORDER) + 2);

/**
 * Antiderivative with `out(0) = value0`: `spower.md` §6's rule, exact.
 *
 * ```
 * mean_a_k = -diff_c_{k-1} / (2k)                            k = 1 .. p+1
 * diff_a_k = (mean_c_k + (k+1) * diff_a_{k+1} / 2) / (2k+1)  k = p .. 0, seeded 0
 * ```
 *
 * Integrating raises the order by one, so `out` must hold `len + 2` coefficients. The
 * result has `diff_a_{p+1} = 0` — degree `2p + 2`, sitting inside the order-`p+1` basis
 * the same way the derivative sits inside order `p`. `out` must not alias `a`.
 */
export function integrateSPower(a: ArrayLike<number>, len: number, out: Float64Array, value0 = 0.0) {
  const p = sPowerOrder(len);

  if (p > MAX_SPOWER_ORDER) {
    throw new Error(`s-power order ${p} exceeds MAX_SPOWER_ORDER`);
  }

  const means = meanScratch;
  const diffs = diffScratch;

  means.fill(0.0);
  diffs.fill(0.0);

  for (let k = 1; k <= p + 1; k++) {
    means[k] = -(a[2 * k - 1] - a[2 * k - 2]) / (2 * k);
  }

  for (let k = p; k >= 0; k--) {
    const meanC = (a[2 * k + 1] + a[2 * k]) * 0.5;

    diffs[k] = (meanC + (k + 1) * diffs[k + 1] * 0.5) / (2 * k + 1);
  }

  // The constant of integration is the one free parameter; pin it by out(0).
  means[0] = value0 + diffs[0] * 0.5;

  for (let k = 0; k <= p + 1; k++) {
    out[2 * k] = means[k] - diffs[k] * 0.5;
    out[2 * k + 1] = means[k] + diffs[k] * 0.5;
  }

  return out;
}

/**
 * `da/du` at `u`, without materializing the derivative coefficients.
 *
 * With `A_k(u) = (1-u) a_k^0 + u a_k^1` and `sym' = 1 - 2u`, the product rule gives
 * `d/du [A_k sym^k] = diff_a_k sym^k + k A_k sym^{k-1} sym'`. Same value as
 * {@link differentiateSPower} followed by {@link evalSPower}, one pass and no scratch —
 * which is what a quadrature inner loop wants.
 */
export function sPowerDerivative(a: ArrayLike<number>, len: number, u: number) {
  const sym = (1.0 - u) * u;
  const dsym = 1.0 - 2.0 * u;
  const u1 = 1.0 - u;

  let acc = 0.0;
  let symK = 1.0;
  let symK1 = 0.0;

  for (let k = 0, i = 0; i < len; k++, i += 2) {
    const value = u1 * a[i] + u * a[i + 1];

    acc += (a[i + 1] - a[i]) * symK;

    if (k > 0) {
      acc += value * k * symK1 * dsym;
    }

    symK1 = symK;
    symK *= sym;
  }

  return acc;
}

/**
 * `d2a/du2` at `u`, likewise in one pass.
 *
 * ```
 * d2/du2 [A_k sym^k] = 2k diff_a_k sym^{k-1} sym'
 *                    + A_k [ k(k-1) sym^{k-2} sym'^2 - 2k sym^{k-1} ]
 * ```
 *
 * (`sym'' = -2`, which is where the last term comes from.) This is what the fourth-order
 * quadrature terms of `spower-solver.md` §7 need, and having it exact rather than
 * finite-differenced is what makes them usable as an error estimator.
 */
export function sPowerCurvature2(a: ArrayLike<number>, len: number, u: number) {
  const sym = (1.0 - u) * u;
  const dsym = 1.0 - 2.0 * u;
  const u1 = 1.0 - u;

  let acc = 0.0;
  let symK = 1.0;
  let symK1 = 0.0;
  let symK2 = 0.0;

  for (let k = 0, i = 0; i < len; k++, i += 2) {
    const value = u1 * a[i] + u * a[i + 1];

    if (k > 0) {
      acc += (2 * k * (a[i + 1] - a[i]) * dsym - 2 * k * value) * symK1;
    }

    if (k > 1) {
      acc += value * k * (k - 1) * symK2 * dsym * dsym;
    }

    symK2 = symK1;
    symK1 = symK;
    symK *= sym;
  }

  return acc;
}

/**
 * `integral of a from 0 to u`, evaluated through a module-scratch antiderivative.
 *
 * Exact — it is {@link integrateSPower} followed by {@link evalSPower}, not a quadrature.
 * The scratch buffer means this is not reentrant, which is fine for the single-threaded
 * uses here but is worth knowing before it is called from a worker.
 */
export function sPowerIntegral(a: ArrayLike<number>, len: number, u: number) {
  integrateSPower(a, len, antiScratch, 0.0);

  return evalSPower(antiScratch, len + 2, u);
}

/**
 * Reverse the parameter, `out(u) = a(1 - u)`.
 *
 * `sym` is invariant under `u -> 1 - u`, so reversal is exactly a swap within each pair.
 * **Signed curvature reverses with a negation on top** — see {@link reverseCurvature}.
 * `out` may alias `a`.
 */
export function reverseSPower(a: ArrayLike<number>, len: number, out: Float64Array) {
  for (let k = 0; k < len; k += 2) {
    const t = a[k];

    out[k] = a[k + 1];
    out[k + 1] = t;
  }

  return out;
}

/**
 * Reverse a *signed curvature* profile: swap each pair and negate it.
 *
 * Under arclength reversal the tangent reverses and signed curvature flips with it, so on
 * top of {@link reverseSPower}'s swap the values change sign. Getting this backwards is
 * easy and quiet — `spower-solver.md` §3 and §13.
 */
export function reverseCurvature(a: ArrayLike<number>, len: number, out: Float64Array) {
  for (let k = 0; k < len; k += 2) {
    const t = a[k];

    out[k] = -a[k + 1];
    out[k + 1] = -t;
  }

  return out;
}

/**
 * The `h(i, k)` pair of `spower.md` §5 — the triangular map from endpoint Taylor data to
 * s-power pairs.
 *
 * ```
 * h(i, k) = ( C(2k-i-1, k-i),  -C(2k-i-1, k-i-1) )        h(k, k) = (1, 0)
 * ```
 *
 * The `h(k, k)` case is a guard for `k = 0`, where `C(-1, 0)` is undefined; for `k >= 1`
 * it is what the binomials already give.
 */
export function hik(i: number, k: number): [number, number] {
  if (i === k) {
    return [1.0, 0.0];
  }

  return [binomial(2 * k - i - 1, k - i), -binomial(2 * k - i - 1, k - i - 1)];
}

/**
 * Build order-`p` pairs from endpoint Taylor data (`spower.md` §5).
 *
 * `f[i]` is `d^i f/du^i (0) / i!` and `g[i]` is `(-1)^i * d^i f/du^i (1) / i!` — note the
 * alternating sign, which is what makes `g` the Taylor data of `u -> f(1 - u)` at 0.
 * Both must carry at least `p + 1` entries.
 *
 * ```
 * a_k^0 = sum_{i=0..k} [ f_i * h(i,k)_0 + g_i * h(i,k)_1 ]
 * a_k^1 = sum_{i=0..k} [ g_i * h(i,k)_0 + f_i * h(i,k)_1 ]
 * ```
 *
 * The two halves of a pair are the same expression with the endpoints exchanged. Truncating
 * at pair `k` leaves endpoint derivatives up to order `k` untouched, which is the whole
 * point of the basis: truncation *is* two-point Hermite interpolation.
 */
export function taylorToPairs(f: ArrayLike<number>, g: ArrayLike<number>, p: number, out: Float64Array) {
  for (let k = 0; k <= p; k++) {
    let a0 = 0.0;
    let a1 = 0.0;

    for (let i = 0; i <= k; i++) {
      const [h0, h1] = hik(i, k);

      a0 += f[i] * h0 + g[i] * h1;
      a1 += g[i] * h0 + f[i] * h1;
    }

    out[2 * k] = a0;
    out[2 * k + 1] = a1;
  }

  return out;
}

/**
 * Endpoint Taylor data from order-`p` pairs — the inverse of {@link taylorToPairs}.
 *
 * `h` is unit lower-triangular in the sense that pair `k` contributes to Taylor
 * coefficient `k` with weight 1 (`h(k,k) = (1, 0)`), so the inversion is a forward
 * substitution and needs no factorization. Writes `p + 1` entries into each of `f`, `g`.
 */
export function pairsToTaylor(a: ArrayLike<number>, p: number, f: Float64Array, g: Float64Array) {
  for (let k = 0; k <= p; k++) {
    let a0 = a[2 * k];
    let a1 = a[2 * k + 1];

    for (let i = 0; i < k; i++) {
      const [h0, h1] = hik(i, k);

      a0 -= f[i] * h0 + g[i] * h1;
      a1 -= g[i] * h0 + f[i] * h1;
    }

    f[k] = a0;
    g[k] = a1;
  }

  return { f, g };
}

const taylorScratch = new Float64Array(sPowerLength(MAX_SPOWER_ORDER));

/**
 * The order-`order` endpoint Taylor coefficients of an order-`p` polynomial, both ends.
 *
 * {@link pairsToTaylor} stops at `p`, which is all the *basis* carries as data. But the
 * polynomial has degree `2p + 1`, so its higher endpoint derivatives exist and are determined
 * — they are simply not free. Reading one of them is what `spower-solver.md` §9's degree
 * continuation needs: the order-`p+1` block entry that would reproduce the order-`p` curve.
 *
 * Returned as `[f, g]` in {@link taylorToPairs}'s convention, so `g` carries the `(-1)^order`.
 * Differentiating in the basis is exact (see {@link differentiateSPower}), and `order > 2p+1`
 * gives zero rather than an error, which is the true answer.
 */
export function endpointTaylor(a: ArrayLike<number>, len: number, order: number): [number, number] {
  const d = taylorScratch;

  for (let i = 0; i < len; i++) {
    d[i] = a[i];
  }

  let factorial = 1.0;

  for (let k = 1; k <= order; k++) {
    differentiateSPower(d, len, d);
    factorial *= k;
  }

  const f = evalSPower(d, len, 0.0) / factorial;
  const g = evalSPower(d, len, 1.0) / factorial;

  return [f, (order & 1) === 0 ? g : -g];
}

/**
 * `(2p+2) x (2p+2)` differentiation matrix, row-major: `D * a` is `a'` in the same basis.
 *
 * The same rule as {@link differentiateSPower}, materialized. This is what makes the
 * stiffness matrix exact — see {@link stiffnessMatrix}.
 */
export function differentiationMatrix(p: number) {
  const n = sPowerLength(p);
  const m = new Float64Array(n * n);

  const col = new Float64Array(n);
  const dcol = new Float64Array(n);

  for (let j = 0; j < n; j++) {
    col.fill(0.0);
    col[j] = 1.0;

    differentiateSPower(col, n, dcol);

    for (let i = 0; i < n; i++) {
      m[i * n + j] = dcol[i];
    }
  }

  return m;
}

/*
  Exact rationals, only as deep as the Gram matrices need: numerators and denominators as
  BigInt, reduced on construction, converted to double exactly once at the end.
*/
function gcd(a: bigint, b: bigint): bigint {
  while (b) {
    [a, b] = [b, a % b];
  }

  return a < 0n ? -a : a;
}

function factorial(n: number) {
  let r = 1n;

  for (let i = 2; i <= n; i++) {
    r *= BigInt(i);
  }

  return r;
}

/** `integral of u^a (1-u)^b du` over `[0, 1]`, exactly: `a! b! / (a + b + 1)!`. */
function betaExact(a: number, b: number) {
  const num = factorial(a) * factorial(b);
  const den = factorial(a + b + 1);
  const d = gcd(num, den);

  return Number(num / d) / Number(den / d);
}

const massCache = new Map<number, Float64Array>();

/**
 * Mass matrix `M_jl = integral of phi_j(u) phi_l(u) du`, row-major `(2p+2)^2`.
 *
 * Basis function `2k` is `u^k (1-u)^{k+1}` and `2k + 1` is `u^{k+1} (1-u)^k`, so every
 * entry is a Beta integral with integer arguments and is computed in exact rational
 * arithmetic before rounding once to double.
 *
 * `spower-solver.md` §6 asks for these as checked-in literals. Building them from exact
 * rationals at load time is the same numbers with a shorter path to being wrong, and the
 * cost is a handful of BigInt multiplies per order, memoized.
 */
export function massMatrix(p: number) {
  const hit = massCache.get(p);

  if (hit !== undefined) {
    return hit;
  }

  const n = sPowerLength(p);
  const m = new Float64Array(n * n);

  for (let j = 0; j < n; j++) {
    const kj = j >> 1;
    const ej = j & 1;

    for (let l = 0; l < n; l++) {
      const kl = l >> 1;
      const el = l & 1;

      m[j * n + l] = betaExact(kj + kl + ej + el, kj + kl + 2 - ej - el);
    }
  }

  massCache.set(p, m);

  return m;
}

const stiffnessCache = new Map<number, Float64Array>();

/**
 * Stiffness matrix `K_jl = integral of phi_j'(u) phi_l'(u) du`, row-major `(2p+2)^2`.
 *
 * Formed as `D^T M D`. That is exact rather than a numerical shortcut: differentiation
 * maps the order-`p` basis into itself (`differentiateSPower`), so `D` carries no
 * truncation and the product is the true stiffness matrix of the basis.
 */
export function stiffnessMatrix(p: number) {
  const hit = stiffnessCache.get(p);

  if (hit !== undefined) {
    return hit;
  }

  const n = sPowerLength(p);
  const d = differentiationMatrix(p);
  const mass = massMatrix(p);
  const k = new Float64Array(n * n);

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let sum = 0.0;

      for (let r = 0; r < n; r++) {
        for (let c = 0; c < n; c++) {
          sum += d[r * n + i] * mass[r * n + c] * d[c * n + j];
        }
      }

      k[i * n + j] = sum;
    }
  }

  stiffnessCache.set(p, k);

  return k;
}

/**
 * Definite integral over `[0, 1]` as a row vector against the coefficients.
 *
 * `integral of phi_j` is the `j`th entry, so `dot(integralWeights(p), a)` is the total
 * turning of a curvature profile. Used by the G1 residual, which needs `integral of kappa`
 * as a *linear functional* of the coefficients rather than as a number
 * (`spower-solver.md` §4).
 */
export function integralWeights(p: number) {
  const n = sPowerLength(p);
  const w = new Float64Array(n);

  for (let j = 0; j < n; j++) {
    const k = j >> 1;
    const e = j & 1;

    w[j] = betaExact(k + e, k + 1 - e);
  }

  return w;
}
