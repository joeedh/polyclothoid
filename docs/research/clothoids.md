# Clothoid curves

Notes on the generalized-clothoid curve primitive and the mesh solver that fits chains
of them. Written against `scripts/stroker/clothoid.js` as of commit `25a68ae`; that file
has since been ported to `src/curve/clothoid.ts`, and the defect table in §8 records
which of the problems below were fixed on the way across. The math is unchanged.

## 1. The formulation

A clothoid segment here is not defined by control points. It is defined by a sampled
**curvature profile** `k(s)`, and geometry is recovered by integration:

```
θ(s) = ∫ k(s) ds                    tangent angle
x(s) = ∫ cos θ(s) ds                position
y(s) = ∫ sin θ(s) ds
```

Because `s` is arclength, the unit tangent is exactly `(cos θ, sin θ)` — no
normalization needed. That is the main payoff of this parameterization, and it is why
`derivative()` is analytic and cheap while `evaluate()` is the expensive one.

Classic Euler spirals are the special case `k(s) = a·s`. Here `k` is an arbitrary
piecewise-linear function with `order` (= `KORDER` = 12) samples, so a segment can be a
line, an arc, a true clothoid, or anything smoother.

## 2. Data layout

`Clothoid.ks` is a single `Float64Array` of `KTOT` = 21 slots doing double duty:

| Index range | Contents |
|---|---|
| `0 .. order-1` | curvature samples (the actual degrees of freedom) |
| `KSCALE` (16) | arclength scale factor |
| `KTH` (17) | rotation angle |
| `KOFFX`, `KOFFY` (18, 19) | translation |
| `KARCSCALE` (20) | `1 / KSCALE` |

`this._ks` is a `Float64Array` **view** over just the first `order` slots. The solver
hands `_ks` to `Constraint` as the parameter array, which is what keeps the solver from
perturbing the cached transform fields. This aliasing is load-bearing — the view and the
backing array must stay the same element type.

## 3. Curvature-profile samplers

Three parallel functions supply `dk/ds`, `k`, and `∫k` over normalized `s ∈ [0, 1]`.
Two complete sets exist; `funcs` selects between them:

- **`piecewise_linear` = `[dstep, step, istep]`** — active. `step` linearly interpolates
  between samples; `istep`/`istep2` compute the *exact* integral of that piecewise-linear
  function by trapezoid accumulation; `dstep` is a finite difference of `step`.
- **`circle_arc`** — piecewise-constant curvature, i.e. a chain of circular arcs. Present
  and correct, currently unused.
- **`bstep` / `bernstein`** — an abandoned Bernstein-basis profile. Note `bstep` sums the
  basis functions without weighting by `ks[i]`, so it is incomplete, not merely unused.

The Readme's idea of a B-spline curvature function would slot in here as a fourth set.

## 4. Integration — `quadrature()`

19 fixed steps, 2nd-order Taylor expansion per step:

```
dx = cos θ − k·sin θ·(ds/2) − (cos θ·k² + dk·sin θ)·(ds²/6)
dy = sin θ + k·cos θ·(ds/2) + (cos θ·dk − k²·sin θ)·(ds²/6)
```

Integration runs over `s ∈ [−0.5, 0.5]`, so the canonical shape is built centered on the
origin and placed afterwards.

## 5. Endpoint interpolation by similarity transform

`_update()` is the key trick. It integrates the canonical shape once for each endpoint,
then computes the offset / rotation / uniform scale that maps that shape's endpoints onto
`v1 → v2`:

```js
ks[KOFFX]   = -s[0];  ks[KOFFY] = -s[1];
ks[KSCALE]  = v1.vectorDistance(v2) / s.vectorDistance(e);
ks[KTH]     = atan2(v2 - v1) - atan2(e - s);
```

So the curvature samples control **shape only**; endpoint interpolation is satisfied by a
similarity transform rather than by constraining the integral. This removes two
constraints per segment from the solver's job, and it is why `KSCALE` shows up
everywhere in the solver — a curvature value in canonical space must be divided by
`KSCALE` to be comparable across segments of different length.

`length` is `KSCALE`, and `evaluate(s)` maps world arclength back to canonical `s` via
`KARCSCALE` before quadrature.

## 6. Derivatives

- `derivative(s)` — analytic, `(cos(θ + KTH), sin(θ + KTH))`.
- `derivative2(s)` — forward finite difference of `derivative`, `df = 1e-4`.
- `curvature(s)` — the standard cross-product formula on those two. An analytic path
  (`funcs[1](...) * KARCSCALE`) exists but is disabled behind `if (0)`.

Using a finite difference for the second derivative when the first is analytic is odd —
`dθ/ds` is just `k(s)`, so `derivative2` could be `k·(−sin θ, cos θ)` exactly. Worth
testing whether the FD noise here is what makes the analytic `curvature` path unusable.

## 7. `ClothoidSolver`

Fits an entire mesh of clothoid edges simultaneously via pathux `Solver` / `Constraint`.

1. Coerce every edge's `ks` to `Float64Array`, reset all curvatures to `0.001`, update.
2. For each vertex with exactly two edges, add constraints (see below).
3. Skip and flag "bad" vertices where the corner is sharper than ~72° (`th < PI*0.4`).
4. `solver.solve(55, 0.7)` — 55 iterations, 0.7 relaxation.
5. Re-update all edges, then force endpoint curvature to zero at the bad vertices.

### Constraints

- **`tan_c`** (registered) — returns `acos(t1 · t2)` between the tangents of the two
  edges meeting at a vertex, handling the sign flip when the shared vertex is `v1` on one
  edge and `v2` on the other. This is the G1 continuity driver and is the only thing
  actually solving right now.
- **`curv_c`** (written, **`solver.add` commented out at line 573**) — G2 continuity, and
  a deliberately different kind of constraint. It does not descend a gradient. It averages
  the two endpoint curvatures in `1/KSCALE`-normalized space (`fac = 0.5`), writes them
  straight back into `ks`, and returns `0.0`.

  This works *because* of how `Solver.solveStep` is written, not in spite of it. Returning
  `0.0` short-circuits the whole descent machinery twice over: `Constraint.evaluate`
  returns at `solver.js:46` before finite-differencing any gradients, and `solveStep`
  hits `if (r1 === 0.0) continue` at `solver.js:105` before touching the parameters. So
  the constraint is invoked once per iteration purely for its side effects, and costs one
  function call instead of `2 × order` evaluations.

  The effect is a **direct projection interleaved with the gradient-descent constraints**
  — Gauss-Seidel style. Each iteration `tan_c` descends toward tangent continuity, then
  `curv_c` projects the curvature endpoints halfway toward agreement. The `fac = 0.5`
  under-relaxation is what keeps the projection from fighting the descent. Registration
  order matters here.

  Its `disabled` flag (`params[3]`) is the other half of the design: when tripped, it
  decays both edges' curvature profiles by `0.98` per iteration instead of projecting —
  a bail-out that relaxes a segment toward straight rather than letting it blow up.

  Since it is structurally sound, `solver.add` being commented out reads as a G2 on/off
  experiment toggle, not as a disabled-because-broken.

### Corner handling

Corners sharper than ~72° are deliberately left as corners: excluded from the solve, then
zeroed. This is the right call for brush strokes, but the threshold is a hard-coded magic
number with no way to tune it per-stroke.

### Progressive refinement (dormant)

`changeOrder(order)` resamples every edge's curvature profile to a new sample count,
intended for a coarse-to-fine solve (`order = 2 .. KORDER`, solving at each level). The
whole path sits behind `if (0)` at line 581. It is the most promising dormant idea in the
file — a 12-DOF-per-edge solve seeded from a converged 2-DOF solve should be far more
robust than the current cold start from `0.001`.

## 8. Known defects

Ordered by whether they can bite today.

**Live:**

- **`mesh.order` does not exist.** Line 543 `const order = mesh.order` evaluates to
  `undefined` — the string `order` appears nowhere in `mesh.js`. The bad-vertex cleanup
  at lines 599–611 therefore does `ks[NaN] = 0.0`, which typed arrays silently ignore. So
  corner-curvature zeroing only works on the `v === e.v1` branch; the `v2` branch is a
  no-op. Should be `e1.curve.order`.

**In dormant code (would bite if re-enabled):**

- `changeOrder` line 405: `temp[i1] + (temp[i2] - temp[i])*t` — the third index should be
  `i1`, not `i`.
- `changeOrder` line 413: builds a `Float32Array` view over a `Float64Array` buffer, and
  assigns it to `e._ks` rather than `e.curve._ks`. Both wrong; see §2 on why the view
  type matters.
- `bstep` ignores `ks` entirely (§3).

**Cosmetic:**

- `quadrature` lines 193–197 compute a `dx`/`dy` pair that lines 199–200 immediately
  overwrite.
- `dstep`, `istep`, `istep2` all have unreachable code after an early `return`.

## 9. State of the file

This is a research file, not settled code. Multiple formulations are kept side by side
and switched with `if (0)` blocks, module-level `funcs` swaps, and commented-out
`solver.add` calls. That is reasonable for exploration, but it means the version that
actually runs is: piecewise-linear curvature, order 12, G1-only, cold-start, 55
iterations. Everything else in the file is either an alternative under evaluation or
dead.

Next things worth trying, roughly in order of expected payoff:

1. Fix `mesh.order` (§8) — real, silent, active.
2. Re-enable progressive refinement (§7) as a proper coarse-to-fine solve.
3. Register `curv_c` and evaluate G2 (§7). It is a projection by design, not a residual —
   the open question is whether the projection and `tan_c`'s descent converge together,
   and whether `fac = 0.5` is the right under-relaxation. Registration order matters.
4. Make `derivative2` analytic and re-test the analytic `curvature` path.
5. The Readme's B-spline curvature profile — drops in at §3 and would subsume
   `circle_arc` (degree 0) and `piecewise_linear` (degree 1) as special cases.
