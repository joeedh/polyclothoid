# s-power series

Notes on the s-power basis and s-power series — a two-point ("bidirectional") analogue of
the Taylor expansion, due to Javier Sánchez-Reyes. This is background reading, not a
description of anything currently in `src/`. §9 covers why it is filed here at all.

Three companion files sit in this directory:

| File | What it is |
|---|---|
| `spower-orig-paper.pdf` | the 1997 original |
| `spower-practical.pdf` | the 2000 follow-up, where the algorithms live |
| `s-power.reduce` | a reduce-algebra implementation of the expansion — see §5 |

Formulas below were re-derived and checked rather than transcribed: both PDFs' text layers
silently drop every minus sign.

## 1. Provenance

| Paper | Where | Role |
|---|---|---|
| J. Sánchez-Reyes, *The symmetric analogue of the polynomial power basis* | ACM TOG 16(3), 1997, 319–357. [doi:10.1145/256157.256162](https://doi.org/10.1145/256157.256162) | **The original.** Introduces the basis, its conditioning, evaluation, and the geometric meaning of its coefficients. Local copy: `spower-orig-paper.pdf`. |
| J. Sánchez-Reyes, *Applications of the polynomial s-power basis in geometry processing* | ACM TOG 19(1), 2000, 27–55. [doi:10.1145/343002.343018](https://doi.org/10.1145/343002.343018) | The practical companion: arithmetic, division, square roots, series, convergence. Local copy: `spower-practical.pdf`. |
| J. Sánchez-Reyes, *Inversion approximations for functions via s-power series* | CAGD 18(6), 2001 | Functional inversion. |
| J. Sánchez-Reyes & J. M. Chacón, *Polynomial approximation to clothoids via s-power series* | CAD 35(3), 2003, 263–271. [doi:10.1016/S0010-4485(03)00045-9](https://doi.org/10.1016/S0010-4485\(03\)00045-9) | Directly on this project's topic. See §8. |
| J. Sánchez-Reyes & J. M. Chacón, *s-power series: an alternative to Poisson expansions* | CAGD 22(2), 2005, 103–119 | Representing general analytic functions. |

Independent lineage, five years later and apparently unaware of the CAGD work: López &
Temme, *Two-point Taylor expansions of analytic functions*, Stud. Appl. Math. 109(4),
2002. Same construction from the asymptotics side.

## 2. The basis

Work over `u ∈ [0, 1]` and define the **symmetric parameter**

```
s = (1 − u)·u          s ∈ [0, 1/4],  max at u = 1/2
```

The order-`p` s-power basis (degree `n = 2p + 1`) is the `2p + 2` functions

```
s^k·(1 − u),  s^k·u          k = 0 … p
```

`s^k` is the scaled central Bernstein polynomial of degree `2k`. It is symmetric about
`u = 1/2` and is a `(k+1)`-fold zero at *both* endpoints — every derivative below order `k`
vanishes at `u = 0` and `u = 1`. That is the whole trick: where `u^k` is the degree-`k`
polynomial vanishing to order `k` at one point, the pair `(s^k(1−u), s^k u)` spans the
degree-`(2k+1)` polynomials vanishing to order `k` at two points.

Two notational traps when reading the 1997 paper against the 2000 one. The original writes
the parameter as `t`, and writes the basis as `P̃⁰_k = t^k(1−t)^{k+1}`, `P̃¹_k = t^{k+1}(1−t)^k`
— the same pair, unfactored. It also carries an even-degree case: for `n = 2q` a single
extra function `P̃_q = t^q(1−t)^q = s^q` sits on the symmetry axis, and that term is *lost*
again at `n = 2q+1`. The 2000 paper drops the even case and uses the odd-degree scaled basis
throughout, because it is the one suited to algebraic manipulation. This document follows
the 2000 convention.

## 3. Expansion, and what a "coefficient" is

A polynomial in this basis is written as a power series in `s` whose coefficients are
**linear functions of `u`** rather than constants:

```
a(u) = Σ_{k=0..p} a_k(u)·s^k          a_k(u) = (1 − u)·a_k⁰ + u·a_k¹
```

so each coefficient is a *pair* `a_k = (a_k⁰, a_k¹)` — the Bézier ordinates of that linear
function. Equivalently the whole thing is a convex blend of two ordinary power series in
`s`, one anchored at each end:

```
a(u) = (1 − u)·A⁰(s) + u·A¹(s)        A⁰(s) = Σ a_k⁰ s^k,  A¹(s) = Σ a_k¹ s^k
```

The second normal form, used for every derivation below, splits the pair into a mean and a
difference (basis `{1, u − 1/2}` instead of Bernstein):

```
ā_k = (a_k⁰ + a_k¹)/2        symmetric part
Δa_k = a_k¹ − a_k⁰           skew-symmetric part
a_k(u) = ā_k + (u − 1/2)·Δa_k
```

Sign conventions vary between write-ups of this material; everything below is stated in
this one and was checked against worked examples.

## 4. Why it matters: truncation is Hermite interpolation

Because `s^k` kills the first `k` derivatives at both ends, **truncating the series at term
`k` changes nothing about the endpoint derivatives up to order `k`**. So

```
H_k(a; u) = Σ_{i=0..k} a_i(u)·s^i
```

is exactly the order-`k` two-point Hermite interpolant of `a` — the degree-`(2k+1)`
polynomial matching `a, a', … , a^(k)` at `u = 0` and `u = 1`. Degree reduction is a
truncation, degree elevation is appending zeros, and the true minimum degree is read off
the highest nonzero pair. In the power basis none of these are free; in the Bernstein basis
degree reduction is a least-squares problem.

The order-0 case is the linear function through the two endpoint values, `H_0(a) = (a⁰, a¹)`,
and it is the only "how do I get the coefficients from derivatives" case the 2000 paper
needs explicitly. The general formula is §5.

Evaluation is Horner in `s`, with a linear blend per term:

```
a(u) = a_0(u) + s·(a_1(u) + s·(a_2(u) + s·(…)))
```

The 1997 paper counts this at `n+1` additions and `n+2` multiplications for degree `n` —
about as fast as monomial Horner, and far better conditioned, since `s ≤ 1/4` on the whole
interval so terms shrink geometrically instead of growing.

## 5. Getting the coefficients from endpoint derivatives

Neither paper spells out how to compute `a_k` from the derivatives of an arbitrary
function; the 2000 paper cites Davis, *Interpolation and Approximation*, p. 37 and moves
on, because its own algorithms build coefficients recursively instead. `s-power.reduce` is
that missing piece, implemented in reduce algebra. The construction:

```
f_i = f⁽ⁱ⁾(0)/i!          Taylor coefficients of f at 0
g_i = g⁽ⁱ⁾(0)/i!          same for g(u) = f(1 − u), i.e. g_i = (−1)ⁱ·f⁽ⁱ⁾(1)/i!

h(i,k) = ( C(2k−i−1, k−i),  −C(2k−i−1, k−i−1) )        h(k,k) = (1, 0)

a_k⁰ = Σ_{i=0..k} [ f_i·h(i,k)₀ + g_i·h(i,k)₁ ]
a_k¹ = Σ_{i=0..k} [ g_i·h(i,k)₀ + f_i·h(i,k)₁ ]        (f and g swapped)
```

The two coefficients of a pair are the same expression with the endpoints exchanged — the
symmetry the whole basis is named for. The `h(k,k) = (1,0)` case is only a guard for `k = 0`
(where `C(−1,0)` is undefined); for `k ≥ 1` it is what the binomials already give. At `k = 0`
this collapses to `a_0 = (f(0), f(1))`, agreeing with `H_0` above.

Checked, not assumed. Reimplementing the script's `hik`/`ak`/`sps` in exact rational
arithmetic confirms both properties that matter:

- **Reproduction.** A degree-7 polynomial expanded to `terms = 3` comes back exactly.
- **Hermite contact.** For a degree-10 input, truncation at `k` matches derivatives `0..k`
  at both `u = 0` and `u = 1`, and — as it should — fails to match at order `k+1`.

Notes for anyone running the script: `bernstein` and `chain` are defined but unused, the
`dfsign` parameter of `ak` is accepted and never read, `sfac` is assigned twice, and the
second `sps` term substitutes `u = 1 − t` into `sfac`, which is a no-op because `s` is
symmetric. None of that affects the result. The `return` inside the `for … sum` block is
load-bearing REDUCE block syntax, not a stray statement.

## 6. Arithmetic

Represent a polynomial as two coefficient lists, `a = (A⁰, A¹)`.

**Addition** is componentwise, and unlike Bernstein there is no degree-matching step.

**Multiplication.** For two linear coefficients the product is *not* linear; it picks up one
extra `s`:

```
a(u)·b(u) = ((a⁰b⁰, a¹b¹)) − Δa·Δb·s
```

(the endpoint values interpolate exactly; the deficit at `u = 1/2` is `−ΔaΔb/4 = −ΔaΔb·s`).
Lifting to full lists, with `*` denoting convolution and `shift₁` a right shift by one:

```
C⁰ = A⁰*B⁰ − shift₁(ΔA*ΔB)
C¹ = A¹*B¹ − shift₁(ΔA*ΔB)          ΔA = A¹ − A⁰
```

Three convolutions instead of one, but on lists of half the length — so the cost lands in
the same place as the power basis.

**Conversion from Bernstein** falls out of multiplication: run de Casteljau symbolically,
reading `u` as the pair `(0, 1)` and `1 − u` as `(1, 0)`. The conversion matrices are *not*
ill-conditioned, which is the practical reason this basis is usable at high degree.

**Composition / subdivision.** Substituting `u(v)` needs only the multiply and add above.
Splitting `[0,1]` at `λ` is composition with `u(v) = (0, λ)` and `u(v) = (λ, 1)`.

**Differentiation.** With `ds/du = 1 − 2u = −2(u − 1/2)` and `t² = 1/4 − s`:

```
c(u) = a'(u) = Σ c_k(u)·s^k

c̄_k  = (2k+1)·Δa_k − (k+1)·Δa_{k+1}/2
Δc_k = −2(k+1)·ā_{k+1}
```

Note the cross-coupling: the skew part of the derivative comes from the symmetric part of
`a`, and vice versa. `Δa_{p+1} = 0` forces `Δc_p = 0`, i.e. the derivative really does drop
to degree `2p`.

**Integration** inverts that — a forward assignment plus a backward recursion, with `ā_0`
the constant of integration:

```
ā_k  = −Δc_{k-1} / (2k)                              k = 1 … p+1
Δa_k = (c̄_k + (k+1)·Δa_{k+1}/2) / (2k+1)             k = p … 0,  seeded Δa_{p+1} = 0
```

**Division and square root** are *online* algorithms in the sense of Knuth — they emit
`c_0, c_1, …` one at a time and can run until the caller stops asking, without knowing the
truncation order up front. Division, given `b⁰_0 ≠ 0 ≠ b¹_0`:

```
c := 0;  r := a
for i := 0 … k:
    c_i := H_0(r / b_0)          # componentwise: (r⁰_{i,0}/b⁰_0, r¹_{i,0}/b¹_0)
    c   := c + c_i·s^i
    r   := r − c_i·b·s^i
```

Square root is the same shape, seeded with `c_0 = H_0(√a)` and updating
`r := r − (2c + c_i s^i)·c_i s^i`. Nth roots generalize directly.

## 7. Series and convergence

Nothing above requires `a` to be a polynomial. Letting `p → ∞` gives the **s-power series**,
the two-point analogue of a Taylor series, whose `k`th partial sum is the order-`k` Hermite
interpolant. Coefficients for transcendental functions come from undetermined coefficients
against a differential equation — the 2000 paper works `e^{z(u)}` for linear `z(u)` this way
and gets a two-term recursion on the pairs; logarithms and inverse trig follow by
integrating a quotient or a square root.

Convergence is the Taylor story with **disks replaced by lemniscates**. A Taylor series
about `z₀` converges inside the largest singularity-free disk centered on `z₀`; a two-point
Hermite series over foci `z₀, z₁` converges inside the largest singularity-free *lemniscate*
`|z − z₀|·|z − z₁| = r²`. For foci `0, 1` that condition is just `|s(z)| < r²`, so:

```
the s-power series converges on [0,1]  ⟺  |s(z)| > 1/4 for every singularity z of a
```

i.e. every singularity's image under `z ↦ (1−z)z` lies outside the disk of radius `1/4`. The
limiting case `r = 1/2` is the lemniscate of Bernoulli. **Subdivision always fixes
divergence** unless a singularity sits on `[0,1]` itself: the limiting lemniscates over a
subinterval are strictly smaller and nested inside the original.

The error term is Cauchy's remainder, with `n = 2k + 1`:

```
R_k(a; u) = a^(n+1)(ξ)/(n+1)! · s^{k+1}          ξ ∈ (0,1)
```

Since `s ≤ 1/4` on the interval, each extra order buys a factor of ~4 on top of the
factorial — the practical reason low-order truncations are usually good enough.

Truncation is *not* the minimax approximation. Sánchez-Reyes is explicit about this: the
Hermite interpolant is not optimal in error, it is optimal in endpoint contact, which is
what a piecewise representation actually wants.

## 8. The clothoid application

The 2003 CAD paper is the one to read for this project. Clothoids are transcendental —
Fresnel integrals — so they cannot be represented exactly by anything in a CAD kernel.
Approximating a clothoid segment with the order-`k` s-power truncation gives a
degree-`(2k+1)` polynomial curve reproducing derivatives up to order `k` at both endpoints;
chaining those over subdivided arcs gives a **Hermite spline with C^k joints and a
near-arclength parameterization**. The mechanism is §7's exponential recursion applied to
`e^{iθ(u)}` with `θ` quadratic, plus the integration rule of §6.

Near-arclength parameterization is not incidental here. It is the property `Clothoid`
already has by construction (`CLAUDE.md` invariant 2: `s` is arclength), and it is the
property that any polynomial replacement would have to keep.

## 9. Relevance to this codebase

Nothing here is implemented. Recording the connections that look real, none of them
verified against our solver:

- **Evaluation.** `Clothoid.evaluate` runs `QUADRATURE_STEPS = 19` fixed steps of a
  2nd-order Taylor integrator per query (`clothoids.md` §4). An s-power truncation would
  replace that with a fixed-degree polynomial evaluated by Horner — same object, no per-query
  quadrature. Whether the accuracy trade is favourable at our degrees is untested.
- **Continuity for free.** Endpoint derivative matching is structural in a truncation, not
  something a solver has to converge to. Our G1/G2 continuity is currently the *output* of
  `ClothoidSolver`'s constraints; an s-power formulation moves part of that into the
  representation.
- **Export.** Conversion to Bernstein is well conditioned in both directions, so a Bézier
  approximation of a segment (for SVG, for a CAD kernel, for anything that will not accept a
  curvature profile) is a basis change rather than a fitting problem.
- **Curvature profiles.** Our `k(s)` is piecewise linear over `KORDER = 12` samples, so `θ`
  is piecewise quadratic and subdivision at sample boundaries is required anyway. That is
  the same subdivision §7 wants for convergence, which is convenient.
- **The `s` collision.** This document's `s` is Sánchez-Reyes' symmetric parameter
  `(1−u)u` on a normalized interval. Everywhere else in this repo `s` is arclength. Any
  implementation needs to pick different names before it starts, or it will be unreadable.

## 10. Reading order

1. The 2000 paper (`spower-practical.pdf`) §§2–3 for the basis and the arithmetic — it
   re-derives enough of the 1997 paper to stand alone.
2. `s-power.reduce` alongside it, for the one thing that paper leaves out: how to get the
   coefficients out of a function in the first place.
3. The 2003 CAD paper for the clothoid specialization.
4. The 1997 original (`spower-orig-paper.pdf`) for the parts the later papers cite rather
   than repeat — §§3.5–3.6 on the conditioning of the Bernstein conversion, §4 on
   evaluation and root-finding condition numbers, and §5 on the geometric reading of the
   coefficients and their use as shape handles.
