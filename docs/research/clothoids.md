# Clothoid curves

Notes on the generalized-clothoid curve primitive and the mesh solver that fits chains of
them. Written against `src/curve/clothoid.ts` and `src/math/solver.ts`. The original was
`scripts/stroker/clothoid.js`, deleted once the TypeScript port landed; git history has it
if a behaviour needs checking against it. The math is unchanged by the port — §8 records
what changed around it.

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

The reduce-algebra derivation of the quadrature terms in §4 is kept verbatim in the header
comment of `clothoid.ts`.

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
backing array must stay the same element type. `loadSTRUCT` rebuilds both after
deserialization for the same reason.

All the `K*` indices are exported, since a caller writing its own constraint needs
`KSCALE` to compare curvatures across segments.

## 3. Curvature-profile samplers

A profile is a `CurvatureProfile`: three mutually consistent members supplying `dk/ds`,
`k`, and `∫k` over normalized `s ∈ [0, 1]`. `activeProfile` selects the one the integrator
uses, and `setCurvatureProfile()` swaps it — a module-level switch, not per-curve state.

- **`piecewiseLinear`** — the default and the only one the solver is exercised against.
  `curvature` interpolates linearly between samples; `integral` computes the *exact*
  integral of that piecewise-linear function by trapezoid accumulation, so it carries no
  step error; `dCurvature` is a finite difference of `curvature` (`df = 1e-5`).
- **`circleArc`** — piecewise-constant curvature, i.e. a chain of circular arcs. Complete
  and switchable, currently unused. Kept because comparing formulations on identical input
  is half the point of the project.
- **`bernsteinCurvature`** — a bare function, not a `CurvatureProfile`, because its
  integral has never been derived and without one it cannot drive `evaluate`. The original
  additionally summed the basis functions without weighting by `ks[i]`, making it a
  constant-1 profile regardless of input; that is fixed here, but it is still parked.

Getting `integral` wrong bends the curve without producing any visible error in the
curvature plot, which is the argument for exact integrals over numeric ones.

The Readme's idea of a B-spline curvature profile would slot in here as a fourth member,
and would subsume `circleArc` (degree 0) and `piecewiseLinear` (degree 1) as special
cases.

## 4. Integration — `quadrature()`

`QUADRATURE_STEPS` = 19 fixed steps, 2nd-order Taylor expansion per step:

```
dx = cos θ − k·sin θ·(ds/2) − (cos θ·k² + dk·sin θ)·(ds²/6)
dy = sin θ + k·cos θ·(ds/2) + (cos θ·dk − k²·sin θ)·(ds²/6)
```

Integration runs over `s ∈ [−0.5, 0.5]`, so the canonical shape is built centered on the
origin and placed afterwards.

Expanding θ to second order (θ, `k`, `dk`) makes the per-step integral third-order
accurate, so the scheme is `O(ds³)` globally. Measured against a dense reference, error in
canonical units on a unit chord:

| N | error | | N | error |
|---|---|---|---|---|
| 3 | 2.1e-3 | | 19 | 2.1e-5 |
| 6 | 1.0e-3 | | 38 | 1.0e-5 |
| 11 | 2.2e-4 | | 76 | 1.2e-6 |

The observed order oscillates between ~1 and ~4 rather than sitting at 3, and at low `N`
the error is not even monotonic (N=3 beats N=4). That is knot aliasing: `order` = 12
samples give 11 intervals and so 10 interior knots, where κ is only C⁰ and θ therefore only
C¹. Below roughly `N = 11` the error is
dominated by where step boundaries fall relative to them, not by `ds³`. Aligning steps to
knots cleans up the mid-range but does not rescue the low end. Dropping to the 3–4 steps
a smooth curvature polynomial would allow requires removing the knots, i.e. §3's
open question about the profile basis — not a change to the integrator.

Until the two defects at the bottom of §8 were fixed, this was `O(ds)` instead: `N = 19`
carried 1.3e-3 of error, sixty times the current figure, and adding steps barely helped.
The step count was compensating for a broken third-order term.

## 5. Endpoint interpolation by similarity transform

`_update()` is the key trick. It integrates the canonical shape once for each endpoint,
then computes the offset / rotation / uniform scale that maps that shape's endpoints onto
`v1 → v2`:

```ts
ks[KOFFX]     = -s[0];  ks[KOFFY] = -s[1];
ks[KSCALE]    = this.v1.vectorDistance(this.v2) / s.vectorDistance(e);
ks[KARCSCALE] = 1.0 / ks[KSCALE];
ks[KTH]       = atan2(v2 - v1) - atan2(e - s);
```

So the curvature samples control **shape only**; endpoint interpolation is satisfied by a
similarity transform rather than by constraining the integral. This removes two
constraints per segment from the solver's job, and it is why `KSCALE` shows up
everywhere in the solver — a curvature value in canonical space must be divided by
`KSCALE` to be comparable across segments of different length.

`length` is `KSCALE`, and `evaluate(s)` maps world arclength back to canonical `s` via
`KARCSCALE` before quadrature. Recomputation is lazy: `update()` sets `recalc`, and every
query calls `_update()` first if it is set.

## 6. Derivatives

All three are analytic. There is no finite differencing left on this path.

- `derivative(s)` — `(cos(θ + KTH), sin(θ + KTH))`.
- `derivative2(s)` — `dθ/ds` is just `k(s)`, so this is `k` times the tangent rotated a
  quarter turn: `k·(−sin θ, cos θ)`, with `k` scaled by `KARCSCALE` into world units.
- `curvature(s)` — read straight off the profile, `activeProfile.curvature(...) *
  KARCSCALE`. No cross product.

`Curve.curvature` in the base class still carries the general cross-product formula
`(x'y'' − y'x'') / |r'|³`; that is what bezier and b-spline use. `Clothoid` overrides it
because the profile *is* the curvature, and going through the derivatives to recover a
value it already holds is both slower and less accurate.

The original finite-differenced `derivative` to get `derivative2` (`df = 1e-4`) and then
ran the cross-product formula on the result, which is why its analytic curvature path sat
disabled behind `if (0)` — the FD noise was the reason the two disagreed.

## 7. `ClothoidSolver`

Fits an entire mesh of clothoid edges simultaneously via the vendored `Solver` /
`Constraint` in `src/math/solver.ts`. Behaviour is configured by `ClothoidSolverOptions`,
whose defaults reproduce the original exactly:

| Option | Default | Effect |
|---|---|---|
| `enableG2` | `false` | register `curv_c` alongside `tan_c` |
| `progressiveRefinement` | `false` | coarse-to-fine solve, `order = 2 .. KORDER` |
| `iterations` | `55` | passed to `solver.solve` |
| `relaxation` | `0.7` | ditto |

`solve()` then:

1. Coerces every edge's `ks` to `0.001` and marks it dirty.
2. For each vertex with exactly two edges, adds constraints (see below).
3. Collects "bad" vertices, those whose authored pairing level is `0`, into a `corners` set
   and skips them.
4. Runs `solver.solve(iterations, relaxation)` — once, or once per order level under
   progressive refinement.
5. Re-updates all edges, then forces endpoint curvature to zero on **both** edges at every
   corner vertex.

### Constraints

- **`tan_c`** — returns `acos(t1 · t2)` between the tangents of the two edges meeting at a
  vertex, handling the sign flip when the shared vertex is `v1` on one edge and `v2` on the
  other. This is the G1 continuity driver and, at default options, the only thing actually
  solving.
- **`curv_c`** (registered only when `enableG2`) — G2 continuity, and a deliberately
  different kind of constraint. It does not descend a gradient. It averages the two
  endpoint curvatures in `1/KSCALE`-normalized space (`fac = 0.5`), writes them straight
  back into `ks`, and returns `0.0`.

  This works *because* of how `Solver.solveStep` is written, not in spite of it. Returning
  `0.0` short-circuits the whole descent machinery twice over: `Constraint.evaluate`
  returns early before finite-differencing any gradients, and `solveStep` hits
  `if (r1 === 0.0) continue` before touching the parameters. So the constraint is invoked
  once per iteration purely for its side effects, and costs one function call instead of
  `2 × order` evaluations. Both call sites carry a comment saying so; see also
  `CLAUDE.md` invariant 1.

  The effect is a **direct projection interleaved with the gradient-descent constraints**
  — Gauss-Seidel style. Each iteration `tan_c` descends toward tangent continuity, then
  `curv_c` projects the curvature endpoints halfway toward agreement. The `fac = 0.5`
  under-relaxation is what keeps the projection from fighting the descent. Registration
  order matters, and `solve()` registers `tan_c` first on purpose.

  Its `disabled` flag on `JointParams` is the other half of the design: when tripped, it
  decays both edges' curvature profiles by `0.98` per iteration instead of projecting — a
  bail-out that relaxes a segment toward straight rather than letting it blow up. Nothing
  currently sets it (see §8).

### Corner handling

Vertices paired at level `0` are deliberately left as corners: excluded from the solve, then
zeroed at both edges' shared end. Which vertices those are is no longer this solver's
decision. The angle test used to live here as `cornerThreshold`, a hard-coded magic number
that became an option; it is now `Stroker.markCorners`, which writes the level onto the
vertex before `Mesh.solve()` runs. `enableG2` still gates G2 globally, so the level a bare
mesh gets is `pairingLevel(v, e1, e2, 0)` clamped by that flag — see
`docs/plans/spower-solver.md` §3 for what the levels mean.

### Progressive refinement (off by default)

`changeOrder(mesh, order)` resamples every edge's curvature profile to a new sample count
and rebuilds `_ks` over the first `order` slots. With `progressiveRefinement` on, `solve()`
walks `order = 2 .. KORDER`, calling `solver.solve` at each level, so a 12-DOF solve is
seeded from a converged 2-DOF one. That should be far more robust than the cold start from
`0.001`, but it is untested — in the original the whole path sat behind `if (0)`.

One wrinkle if it is turned on: the constraints capture `e._ks` at registration time, when
`order` is still `KORDER`, so they keep stepping all 12 slots even while the curve reads
only the first `order` of them. The views alias the same buffer so nothing is corrupted,
and the dead slots finite-difference to a zero gradient, but each iteration pays for
evaluations that cannot move anything.

## 8. Defect status

The original's defects, and what the port did with each. All but the last two are fixed;
the fixed entries are kept because several describe traps that are easy to reintroduce.

| Defect | Status |
|---|---|
| `mesh.order` was undefined, so corner-curvature zeroing did `ks[NaN] = 0` — silently a no-op on a typed array — and only the `v === e.v1` branch ever took effect | Fixed. `setEndCurvature` indexes off `e.order`, and both branches run. |
| `changeOrder` computed `temp[i1] + (temp[i2] - temp[i])*t`, meaning `temp[i1]` for the third index | Fixed. |
| `changeOrder` built a `Float32Array` view over a `Float64Array` buffer and assigned it to `e._ks` rather than `e.curve._ks` | Fixed. Both matter — see §2 on why the view type is load-bearing. |
| `bstep` summed the Bernstein basis without weighting by `ks[i]`, making it a constant-1 profile | Weighting fixed, but the profile is still incomplete: no integral (§3). |
| `quadrature` computed a `dx`/`dy` pair that the next two lines overwrote | Deleted. |
| `dstep`, `istep`, `istep2` had unreachable code after an early `return` | Deleted, keeping the live branch. |
| `derivative2` finite-differenced `derivative`, feeding noise into `curvature` | Fixed — analytic (§6). |
| `linearCurvature` guarded its interpolation with `i2 < klen - 1`, so the last interval was flat at `ks[klen-2]` and **`ks[klen-1]` was never read**, while `integral` ramped to it. The two profile members described different functions | Fixed — `i2 <= klen - 1`. This is the defect that cost the integrator two orders of convergence (§4). |
| Both members computed `t = fract(i1)` from the *unrounded* index and then indexed with `~~(i1 + 1e-5)`, pairing a `t` near 1 with an `i1` already past the knot — a spurious jump of one full sample on a 1e-5-wide window below every knot | Fixed. `t` is now derived from the same `i1` used to index, which makes the boundary continuous and removes the need for the snapping epsilon at all. |
| `circleArc.curvature` indexes with `~~(s·(klen−1))`, i.e. intervals of width `1/(klen−1)`, but `circleArc.integral` accumulates with `ds = 1/klen`. The two disagree by exactly `(klen−1)/klen` — a systematic 8.3% error in θ at `klen = 12` | **Won't fix.** Measured, not fixed: it is an unused code path, and the profile is slated for removal along with the b-spline curve type. §3 describes it as complete and switchable; it is neither, and it should not be used as a reference shape for quadrature measurements in the meantime. |
| `curvatureConstraint` sets `flip = isV1e1 !== isV1e2` (`clothoid.ts:530`) while `tangentConstraint` negates when `isV1e1 === isV1e2` (`clothoid.ts:490`) — opposite tests for the same question, whether the two edges traverse the vertex in the same direction. When both edges have `v` as their `v1` the through-path traverses `e1` backwards and `κ` must flip, which the tangent version does and the curvature version does not | **Fixed** — `const flip = isV1e1 === isV1e2;`. Confirmed at runtime by `tests/joint.test.ts`, which solves the same three-point chain in all four edge orientations and reads the joint curvature in world space: the G2 residual went from `−2.5e-2` to `−3.1e-4` and became uniform across configurations. The two tests are exact negations, so the rule was inverted in every configuration, not two of them. Only reachable when `enableG2` is `true`, which is off by default, so it was unexercised; the fix inverts behaviour for anyone who had enabled G2. |

The two fixed `piecewiseLinear` defects survived the TypeScript port unchanged; they are
original defects, not port regressions. `piecewiseLinear.integral` is now the exact integral of
`piecewiseLinear.curvature` to 1e-12, which is the contract `CurvatureProfile`'s doc
comment states and the thing §3 warns is invisible in a curvature plot when broken.

The first of those also made `enableG2` partly fictional: `endCurvature(e, false)` reads
`ks[order - 1]`, so at every `v2` end the G2 projection and the corner-zeroing pass were
both writing a value that had no effect on the geometry. Together with the `flip` entry,
that is two independent reasons the G2 path had never done what it claims. Both are now
fixed, and `enableG2` reaches G2 at a joint for the first time — though only to about
`3e-4` in world curvature, which is the Kaczmarz sweep's own convergence level rather than
anything about the constraint. `docs/plans/spower-solver.md` retires `enableG2` as a flag
entirely: G2 is the floor of its continuity ladder rather than an option.

One thing `tests/joint.test.ts` measures in passing is worth recording, because it bears on
that plan's §5. The solved geometry is only orientation-invariant to about `1.7e-4` —
reversing an edge changes the answer by roughly the same amount as the solver's residual.
That is inherent to the sweep rather than a defect: `derivative(0)` depends on the
curvature samples only through `KTH` while `derivative(length)` depends on all of them, so
reversing an edge genuinely changes the sensitivity structure the descent sees. It is also
why that test does *not* detect the `flip` bug — the deviation was `1.68e-4` with the wrong
sign rule and `1.65e-4` with the right one, the G2 projection barely moving the shape.

Live gaps, as opposed to defects:

- **Nothing ever sets `JointParams.disabled`.** The blow-up bail-out in `curv_c` is
  written and correct but unreachable; there is no divergence detector to trip it. If G2 is
  turned on for real, that detector is the missing piece.
- **No tests.** Every claim above about convergence is from eyeballing the canvas demo.

## 9. State of the file

The port turned the original's `if (0)` blocks, module-level `funcs` swaps, and
commented-out `solver.add` calls into named options and exported profiles, so the
alternatives are now selectable from a caller rather than by editing the file. What
changed is the switching mechanism, not which branch runs: the default configuration is
still **piecewise-linear curvature, order 12, G1-only, cold-start from 0.001, 55
iterations at 0.7 relaxation**.

Next things worth trying, roughly in order of expected payoff:

1. Turn on `progressiveRefinement` and see whether coarse-to-fine actually beats the cold
   start. Fix the stale-`_ks` waste noted in §7 if it does.
2. Turn on `enableG2` and evaluate. `curv_c` is a projection by design, not a residual —
   the open question is whether the projection and `tan_c`'s descent converge together,
   and whether `fac = 0.5` is the right under-relaxation. Registration order matters.
3. Add the divergence detector that trips `disabled`, without which G2 has no bail-out.
4. Try `circleArc` against `piecewiseLinear` on the same stroke input; the switch is one
   `setCurvatureProfile` call now.
5. The Readme's B-spline curvature profile — drops in at §3 and would subsume `circleArc`
   (degree 0) and `piecewiseLinear` (degree 1) as special cases.
