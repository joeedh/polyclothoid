# An s-power clothoid and a two-level solver

A design for a second curve/solver pair alongside the existing `Clothoid` /
`ClothoidSolver`: the curvature profile represented in the s-power basis instead of
piecewise-linear samples, and a solver built around what that representation makes
possible. Math background is `docs/research/spower.md` (the basis) and
`docs/research/clothoids.md` (the current formulation and its measurements). Nothing here
is implemented.

The existing clothoid stays. Comparing formulations on identical input is half the point
of the project, and the piecewise-linear profile is the baseline this has to beat.

## 1. Premise

Three problems with the current system have the same root, and one change addresses all
three.

| Problem | Current cause |
|---|---|
| Curvature continuity is a self-referential projection that has never been trusted enough to enable by default | `curv_c` normalizes by `KSCALE`, writes `ks`, then calls `update()` — which recomputes `KSCALE`. It projects onto a set that moves because of the projection. |
| Convergence is by Gauss-Seidel sweep with a fixed iteration count, and slows as `O(N²)` in chain length | One scalar constraint per joint, gradient-descended, information propagating one joint per sweep. |
| The Taylor quadrature can't reach its asymptotic order at low step counts | The profile has `order − 1` = 11 knots where θ drops to C². Below ~11 steps the error is dominated by knot aliasing (`clothoids.md` §4). |

The root is the representation. Twelve piecewise-linear samples per edge give a
non-smooth profile with far more degrees of freedom than there are constraints, no
structural relationship between adjacent edges, and no basis-level notion of endpoint
derivative data.

An s-power profile is a smooth polynomial whose coefficients *are* endpoint derivative
data. That makes `Gⁿ` continuity a property of the data layout rather than something a
solver converges to, removes the knots, and drops the DOF count to something commensurate
with the constraint count.

## 2. Naming

`spower.md` §9 flags this as the thing to settle before writing any code, and it is:

| Symbol | Meaning | In code |
|---|---|---|
| `s` | arclength — unchanged, `CLAUDE.md` invariant 2 | `s` |
| `u` | normalized segment parameter, `u ∈ [0, 1]` | `u` |
| `σ` | Sánchez-Reyes' symmetric parameter, `σ = (1 − u)·u ∈ [0, 1/4]` | `sym` |
| `κ` | curvature | `k` (matches existing code) |

Every citation of the papers uses their `s` for our `σ`. Translate on the way in.

## 3. Representation

### The basis

Order `p`, degree `2p + 1`, `2p + 2` coefficients arranged as pairs
`a_k = (a_k⁰, a_k¹)`, `k = 0..p`:

```
κ(u) = Σ_{k=0..p} [ (1 − u)·a_k⁰ + u·a_k¹ ] · σ^k
```

Because `σ^k` is a `(k+1)`-fold zero at both ends, truncating at term `k` leaves the
endpoint derivatives up to order `k` untouched (`spower.md` §4). The consequences that
matter here:

| Order `p` | Degree of κ | Endpoint data matched | Continuity if shared |
|---|---|---|---|
| 0 | 1 | κ | G2 |
| 1 | 3 | κ, κ′ | G3 |
| 2 | 5 | κ, κ′, κ″ | G4 |

Degrees are odd. The 1997 paper's even-degree case is dropped in the 2000 convention this
follows (`spower.md` §2), so the ladder is 1, 3, 5, 7 with no degree 4.

**G2 is the floor of this ladder, not a target on it.** The lowest rung, `p = 0`, is one
scalar per vertex, and sharing it *is* G2. There is no rung below. So `enableG2` does not
survive the port as a flag — the option is `p`, and with it goes the entire `curv_c`
apparatus: `endCurvature`, its orientation flip, and the normalize-by-`KSCALE`-then-
`update()` self-reference that makes row 1 of §1's table a problem. That row is deleted
rather than fixed.

### Degrees of freedom live on vertices, not edges

This is the structural change, and it is what eliminates the continuity constraints rather
than merely linearizing them.

Each vertex `v` owns a block of `p + 1` scalars: the world-space curvature derivatives at
that point, scaled to be dimensionally commensurate,

```
block(v)_n = (Lᵥⁿ / n!) · dⁿκ/dsⁿ |_v          n = 0 .. p
```

with `Lᵥ` a reference length for the vertex (mean of its incident edge lengths). World
units, because curvature continuity is a geometric condition — the current code already
divides by `KSCALE` for exactly this reason.

An edge's own coefficient vector is then a constant linear function of its two endpoint
blocks:

```
a_e = M_e · [ block(v1) ; block(v2) ]
```

`M_e` composes three things, all constant once the edge's arclength `L_e` is fixed: the
rescale from `Lᵥ` to `L_e`, an orientation sign, and the triangular `h(i,k)` map from
Taylor data to s-power pairs (`spower.md` §5). Note that `h` mixes both endpoints, so
`a_k⁰` depends on both blocks — but never on a third, so each edge couples exactly two
adjacent vertex blocks.

The orientation sign is `(−1)ⁿ⁺¹`, i.e. **even** orders flip. Under arclength reversal
`s̃ = L − s` the tangent reverses and signed curvature flips with it, but its derivative
does not:

```
κ̃ = −κ        dκ̃/ds̃ = +κ′        d²κ̃/ds̃² = −κ″        dⁿκ̃/ds̃ⁿ = (−1)ⁿ⁺¹ dⁿκ/dsⁿ
```

Getting this backwards is easy and quiet. The existing `curvatureConstraint` and
`tangentConstraint` currently disagree about it — see §12.

### Valence, pairings, and continuity levels

The mesh already carries arbitrary valence (`Vertex.edges` is unbounded); it is the
solvers that assume chains, skipping anything with `edges.length !== 2`. Generalizing the
block model covers the general graph without a second mechanism.

The authored primitive is a **pairing of two edge-ends at a vertex, carrying a continuity
level**:

| Level | Meaning | Structure |
|---|---|---|
| 0 | position only — a corner | no G1 row, no shared block |
| 1 | G1 — tangent continuous, curvature free | G1 row, no shared block |
| 2 | G2 | G1 row + share block entry 0 (`κ`) |
| `k` | G`k` | G1 row + share block entries `0 .. k−2` |

The ceiling is `p + 2`. Unpaired edge-ends share nothing but position.

A single partition is not enough, because "broken tangent" and "broken curvature" are
different amounts of breaking and both need expressing. With mixed levels around one
vertex — say A–B at level 3 and B–C at level 2 — there is no single partition, but there
is a clean reading: **one partition per order, each a coarsening of the next.** The group
at order `n` is a connected component of the pairing graph restricted to pairings with
level ≥ `n + 2`; find them with a union-find per order. Graded breaking (§8) is then not a
separate mechanism, just an intermediate level.

A group of size `m` needs `m − 1` G1 rows, a spanning tree of the group — not `m(m−1)/2`.
Groups of size > 2 mean three or more edge-ends mutually tangent at a point:
representable, geometrically degenerate, rarely wanted.

Valence 1 is a cap. It has no pairing partner, so its block is governed by energy alone and
KKT stationarity supplies natural boundary conditions for free. Valence 0 is skipped.

**Within a group there is no way to express a `Gⁿ` discontinuity for `n ≤ p + 1`.** Not a
constraint that converges to zero; not a row in a Jacobian. It does not exist in the state
space. Discontinuity is expressed by lowering the pairing level, which is a change of
structure, not of value.

### Default policy, and who is allowed to decide

**The core solver never raises continuity, only lowers it.** The authored specification is
a ceiling and the solve result is ≤ ceiling. No inference lives in `curve/` — not pairing
inference, not corner detection. This is what makes a load/solve/save round-trip safe:
re-solving cannot drift a document toward smoothness nobody asked for.

The default ceiling is valence-dependent, and honestly so:

| Valence | Default | Why |
|---|---|---|
| ≤ 2 | maximum continuity (`p + 2`) | the pairing is unambiguous, so there is nothing to guess |
| ≥ 3 | fully split, no pairings | no canonical pairing exists; any default is a guess |

"Highest possible continuity" is well defined at valence ≤ 2 and *undefined* at valence
≥ 3 — which is exactly why pairing must be authored there. Fully split is the honest
default at a junction, not a degraded one.

Client mesh code is where inference belongs: deriving pairings from incoming geometry,
promoting sharp input angles to level-0 pairings, and acting on the solver's diagnostics
(§8). The core takes the resulting specification as given.

DOF accounting for an open chain of `V` vertices, against the current system:

| | Current | s-power, `p = 1` |
|---|---|---|
| Parameters | `12·(V−1)` | `2V` |
| Constraint rows | `V−2` (G1), plus a G2 projection | `V−2` (G1) |
| Continuity constraints | G2 projected, G3+ inexpressible | G2, G3 structural |

## 4. What stays a constraint

Tangent continuity. It is not a coefficient identity — it is a condition on the *integral*
of κ across the joint plus the segment's placement angle:

```
r_v = wrap( θ_e1(1) − θ_e2(0) )       θ_e(u) = ∫₀^u κ · L_e du + KTH_e
```

`∫κ` is a linear functional of the coefficients, and `spower.md` §6 gives it exactly in
closed form on the pairs. So the residual is *linear in the DOF plus one frozen scalar per
edge*. Endpoint position stays satisfied by construction via the similarity transform, as
today.

Use the wrapped signed angle difference, not `acos(t1·t2)`. The current residual
(`clothoid.ts:493`) has `acos(1−ε) ≈ √(2ε)`, so its gradient diverges exactly at the
solution; that singularity is currently masked by `Constraint.threshold` short-circuiting
to `0.0`, which means the accuracy floor is set by a workaround. The angle is analytic and
already available.

## 5. The remaining nonlinearity, and the solver that follows

With κ in this basis, everything is linear in the DOF except two scalars per edge, both
computed in `_update()`:

```
KTH_e   = atan2(v2 − v1) − atan2(canonical chord)
KSCALE_e = |v2 − v1| / |canonical chord|
```

Both are transcendental functions of the coefficients through `∫cos θ`, `∫sin θ`.
Clothoids are Fresnel integrals; no change of curvature basis removes this. But it is
confined to `2E` scalars, which suggests the solver directly:

**Inner solve — transforms frozen.** `M_e` and `KTH_e` are constants, so the problem is a
quadratic energy subject to linear G1 rows. The KKT system

```
[ H   Jᵀ ] [ x ]   [ −g ]
[ J   0  ] [ λ ] = [ −c ]
```

has `H` block tridiagonal in vertex index (each edge couples adjacent blocks) and `J` one
row per joint touching two adjacent blocks. The whole thing is banded: direct `LDLᵀ` in
`O(V·(p+1)³)`. No iteration, no relaxation parameter, no sweep-order dependence, no
`O(N²)` chain slowdown. A closed loop makes it cyclic banded — a rank-one correction.

**Outer loop — recompute the transforms.** A fixed-point iteration on `2E` scalars. Its
contraction factor is one measurable number per iteration rather than the emergent
behaviour of coupled projections.

That contraction factor is also the divergence detector. It degrades precisely when
`∂KTH/∂a` and `∂KSCALE/∂a` blow up, which happens when the canonical chord shrinks toward
zero — total turning approaching closure. That is a hard geometric bound with a clear
meaning, and it is what should trip the bail-out that `JointParams.disabled` was written
for and that nothing currently sets (`clothoids.md` §8).

### Chain decomposition and substructuring

The banded structure above is stated for a chain, and the vertex groups of §3 are what
make it hold on a general graph. Decompose the mesh into **nodes** (valence ≠ 2) and
**chains** (maximal paths of valence-2 vertices between them). Then:

- With every node fully split into singleton groups — §3's default at valence ≥ 3, so this
  is the common case rather than the fallback — no block is shared across a chain boundary
  and the chains are *completely decoupled*. Each is an independent banded solve of
  `O(len)`, embarrassingly parallel, with no global assembly at all.
- With some node groups shared, the shared blocks are the only cross-chain coupling.
  Eliminate each chain's interior onto its endpoints and solve the resulting Schur
  complement on node DOF only — standard substructuring / nested dissection. The
  complement is small: `O(#shared groups × (p+1))`.

An isolated all-valence-2 cycle has no nodes at all. It needs either an arbitrary cut
point promoted to a node, or the cyclic-banded rank-one correction already noted above.

**Two decompositions, and they are not the same one.** Keeping them separate matters:

| | Basis | Lifetime | Used for |
|---|---|---|---|
| Topological chains | maximal valence-2 runs | static, from the graph | stroking: width profiles, cap/join placement, offsetting |
| Solve blocks | maximal runs of *coupled* DOF | dynamic — re-derived when a level changes | factorization structure |

Solve blocks are topological chains further split at authored corners, and further still
whenever §8's stability breaking lowers a level mid-solve. Breaking a corner splits a solve
block but *not* a topological chain: position is still shared there, only the DOF coupling
parts. Conflating the two will fail exactly when a level changes mid-solve — rare under
§8's policy, which is precisely why it will not be caught by ordinary use.

### Width is the same machinery

The reason to make chains first-class rather than an implementation detail: a width
profile `w(s)` along a chain is another scalar profile, with the same structure as `κ`.
Same s-power basis, same continuity-by-shared-block at nodes, same `∫(w′)²` energy for
smooth taper, same banded solve, and the same authored pairing levels — a hard width
discontinuity is a level-0 pairing on the width profile, independent of the level `κ`
carries at that joint.

It is strictly easier than curvature in one respect — `w` does not feed the similarity
transform, so there is no nonlinear coupling and no outer loop. Width is inner-solve only.

So the solver is better understood as operating on *scalar profiles over chains*, with `κ`
as the first one and the geometric coupling in §5 as `κ`'s peculiarity rather than the
general case. That is the abstraction to build against; variable-width stroking then
reuses it rather than reimplementing it.

## 6. Energy

Nothing in the current system minimizes energy; the Kaczmarz step's minimum-norm property
keeps it near the all-`0.001` seed, which is an accidental `‖Δk‖₂` regularizer. Make it
explicit:

```
E = Σ_e [ (1/L_e)·∫₀¹ κ′(u)² du  +  ε·L_e·∫₀¹ κ(u)² du ]
```

The first term is minimum-variation, which for brush strokes is generally preferable to
the elastica `∫κ²` — it does not bias toward circular arcs. The `ε` term is there because
`∫(κ′)²` alone has a one-dimensional null space (globally constant curvature), which makes
`H` singular.

That weighting is an *artistic* default, not a universal one (§9). Engineering fairing
often wants the opposite bias — arcs and lines where they will do — which is elastica, and
classical CAD fairing frequently targets curvature-plot monotonicity, a different objective
again. Keep the two coefficients as mode parameters rather than constants; the Gram
matrices below are shared either way, so this costs nothing structurally.

Both Gram matrices are constant `(2p+2)²` matrices per order, computable once in exact
rational arithmetic — the entries are `∫₀¹ σ^{j+k}(1−u)^α u^β du`, i.e. Beta functions.
For `p ≤ 3` these are at most 8×8. Compute them offline, check them in as literals, and
record `κ(H)` for each order.

Note the tension with orthogonality: s-power is not an orthogonal basis, so unlike shifted
Legendre it does not diagonalize `∫κ²`. It buys structural endpoint continuity instead.
The two bases are related by a constant invertible matrix per order, so if the local blocks
ever condition badly, that map is available as a preconditioner. At `p ≤ 2` this is very
unlikely to matter; measure before adding it.

The `Lⁿ/n!` scaling in §3 is not cosmetic. Raw `dⁿκ/dsⁿ` values differ by orders of
magnitude across `n`, which wrecks the column scaling of `J` and `H` alike.

## 7. Evaluation and quadrature

The geometry integrals stay transcendental — `κ` polynomial does not make `∫cos(∫κ)`
polynomial — so `quadrature()`'s Taylor stepping is retained unchanged in structure. What
changes is that it can finally reach its asymptotic order.

`clothoids.md` §4 records the measurement: with the piecewise-linear profile the observed
order oscillates between ~1 and ~4 and the error is not monotonic below `N ≈ 11` (three
steps beats four), because 11 knots dominate. A polynomial κ has no knots, so `O(ds³)`
should hold from the first step. Target 3–4 steps for graphic-design quality, configurable
upward.

The fourth-order terms become worth adding at that point — they need `κ″`, which this
basis supplies exactly via §6's differentiation rule:

```
dx += ( (κ³ − κ″)·sin θ − 3κ·κ′·cos θ ) · ds³/24
dy += ( (κ″ − κ³)·cos θ − 3κ·κ′·sin θ ) · ds³/24
```

(These share `(κ³ − κ″)` and `3κκ′`, so two extra scalars per step.)

They also double as an **error estimator**, which is what the engineering mode needs. A
fixed 3–4 steps is an artistic default; CAD wants a tolerance, and the magnitude of the
`ds³` terms bounds the local truncation error of the scheme without them. So adaptive step
control costs only the terms already worth computing: evaluate both, refine where the
difference exceeds tolerance. Exact `κ″` from §6's differentiation rule is what makes the
estimate trustworthy rather than a finite-difference guess.

All three `CurvatureProfile` members are exact and cheap: `κ` by Horner in `σ`, `κ′` by
§6's differentiation rule, `∫κ` by §6's integration rule. That matters — `clothoids.md` §3
notes that an inconsistent `integral` bends the curve with no visible error in the
curvature plot, which is precisely the defect that was found and fixed in the
piecewise-linear profile.

## 8. Constraint breaking

The current mechanism is `cornerThreshold` on the input polyline angle, decided before the
solve runs, blind to segment lengths and to what the rest of the chain is doing. Under the
policy in §3 it does not get a principled replacement inside the solver — it gets moved
out. Breaking splits into two jobs with different owners.

### What the solver may break: stability only

Lowering a pairing level is mechanical — decrement it, re-derive the per-order groups,
refactor. What needs a criterion is *when*, and under §3's invariant the only admissible
reasons are numerical:

1. **Rank deficiency in `J`.** Conflicting or redundant G1 rows around a junction, showing
   up as a negligible pivot in the KKT `(2,2)` block during factorization. Genuine
   infeasibility, detected by the solve itself, with no threshold to tune.
2. **Chord degeneracy.** The §5 bound — the outer fixed point failing to contract. Also a
   hard geometric limit rather than a tuning knob.

Note what is *not* on this list: `|λ_v|`. A large multiplier means the joint is expensive
to hold, which is a judgement about whether a corner was wanted there — modelling, not
stability. Under this policy the solver may not act on it.

**In engineering mode the solver may not break at all.** Detection is the same; the
response is not. The value of asking for G4 is the guarantee, so quietly delivering G2 at
one joint produces a model that looks right and is not — worse than no answer. Degrade and
report in artistic mode, fail and report in engineering mode (§9). This is a mode flag on
the solve, not a per-joint decision.

Two requirements on any break the solver does make:

- **It must be reported.** Otherwise the client cannot distinguish an authored break from a
  solver bail-out, and a load/solve/save round-trip silently promotes solver failures into
  document state. The result needs a per-vertex channel saying "requested level `k`,
  delivered level `j`, because *reason*."
- **It must be hysteretic.** A joint that flickers between broken and smooth while the user
  drags is worse than either state. Break and restore thresholds must differ, and the
  decision must be deterministic given identical input.

### What the client breaks: authoring, with better information

Everything aesthetic — corner detection, pairing inference at junctions, taste — belongs in
client mesh code, run at input time and recorded as authored levels.

The solver's contribution is a better signal than the client can compute on its own.
`λ_v` falls out of the KKT solve for free: the force required to hold tangent continuity at
that joint, in energy units, accounting for segment lengths and for what the rest of the
chain is doing. That is strictly better corner evidence than an input-polyline angle. It is
reported as a diagnostic, and the client may act on it by authoring a level-0 pairing on
the next edit. Suggestion, not decision.

Same for graded breaking, if it is wanted: `Constraint.wlst` already exists, is allocated,
is filled with `1.0` and never touched — the natural hook for IRLS/Huber weighting on the
G1 rows, with the scale set from `δ = 1.345·MAD(residuals)` rather than by eye. But that is
softening a constraint the client asked for, so the client sets the policy.

### Diagnostics

Engineering mode makes this load-bearing: if the answer to a degenerate configuration is
"refuse," the user has to learn *what* and *where* well enough to fix it. But the same
records are the debugging surface for development, which is the stronger reason to build
them properly — every convergence claim in this document is currently checked by eyeballing
the demo (§12).

A diagnostic is a located, typed record, not a boolean: **where** (vertex, edge, or chain),
**what condition**, the **measured value against its threshold**, and **what the solver did**
— nothing, degraded from level `k` to `j`, or refused.

| Condition | Measured from | Meaning |
|---|---|---|
| chord degeneracy | `\|canonical chord\|` in `_update()` | total turning approaching closure; the §5 hard bound |
| outer loop not contracting | ratio of successive outer iterate deltas | the §5 fixed point failing; ≥ 1 is divergence |
| G1 rank deficiency | negligible pivot in the KKT `(2,2)` block | conflicting or redundant tangent rows at a junction |
| ill-conditioned `H` | smallest pivot of the `(1,1)` block | energy Hessian near-singular; usually a degenerate `L_e` |
| level lowered | the active set | artistic mode only; engineering mode refuses instead |
| large multiplier | `λ_v` from the solve | not a fault — corner evidence for the client (above) |

Three rules make these useful rather than decorative:

1. **Measure unconditionally, not only on failure.** A contraction factor of 0.9 is not a
   fault but is exactly what you want to see when a stroke feels sluggish. Severity gates
   what a client *surfaces*; it must not gate what the solver *records*.
2. **It is nearly free, which is what makes (1) affordable.** Every row above is a byproduct
   of a computation already happening — pivots from the factorization, the contraction ratio
   from iterates already held, the chord from `_update()`. Nothing here needs its own pass.
3. **It must be deterministic.** Identical input, identical diagnostics. Same requirement as
   the hysteresis rule above and for the same reason: a report that varies run to run is
   worse than none when bisecting.

Per-iteration history — contraction factor and largest multipliers per outer step — is what
actually diagnoses non-convergence, but it allocates, so keep it opt-in behind a trace flag
rather than always-on.

This also gives `JointParams.disabled` a real replacement. `clothoids.md` §8 records it as
written, correct, and unreachable because nothing detects divergence; row 2 above is that
detector.

## 9. Degree continuation

`p = 0 → 1 → 2`. Because truncation is Hermite interpolation, the order-`p` solution is
exactly the order-`p+1` solution with the new pair zeroed — continuation is *appending
zeros*, with no resampling and no perturbation of converged coefficients. Compare
`changeOrder`, which interpolates and moves every value.

Each rung raises smoothness as well as degree, so the ladder is G2 → G3 → G4. Reported
experience with this class of solver puts the sweet spot at cubic-to-quartic curvature;
`p = 1` (cubic, G3) is the expected default, `p = 0` the continuation seed and the fallback
when the outer loop struggles, `p = 2` (quintic, G4) opt-in for quality work.

### Order is a domain setting, not a quality dial

The two rungs serve different applications, and the applications do not overlap:

| | Artistic — graphic design, brush strokes | Engineering — CAD |
|---|---|---|
| Order | `p = 0..1`, G2–G3 | `p = 2`, G4 |
| Solver may lower a level for stability (§8) | yes, degrade and report | **no — fail and report** |
| Energy (§6) | minimum-variation, no arc bias | elastica or fairing-specific weighting |
| Quadrature (§7) | fixed 3–4 steps | tolerance-driven adaptive |
| Edit locality | matters | does not matter |

Raising `p` makes the chain stiffer as well as smoother — at G4, `κ`, `κ′` and `κ″` are
coupled across every joint, so a disturbance travels much further before the energy term
damps it. That would be a problem for interactive handle-dragging, but it lands only in the
engineering case, where locality is not a requirement. It is a non-issue rather than a
tradeoff, and no localization workaround is needed.

What does follow is that `p` is a per-document (or per-chain) setting with a cluster of
other defaults attached, not a single number to turn up for quality. Treat the rows above
as one choice, not five.

## 10. Fitting into the codebase

New `SPowerClothoid implements Curve` and `SPowerSolver implements CurveSolver`, selected
at runtime the same way the others are — `Mesh.CurveCls`/`SolverCls`,
`Mesh.switchSplineType()`, `StrokerOptions`. `curve/` still must not import `mesh/`.

One genuine architectural friction: **the DOF no longer live on the curve.** `Curve` and
`CurveSolver` both assume the curve owns its parameters, and `Clothoid.ks` is the model
for that. Vertex-owned blocks need somewhere to live and some way for a solver in `curve/`
to reach them through the structural interfaces in `curve/mesh_types.ts`. Options:

- extend `SolvableVertex` with an opaque solver-data slot;
- have `SPowerSolver` own a side table keyed by vertex;
- keep blocks on edges and re-impose sharing as an explicit constraint (defeats the point).

This needs deciding before Phase 2, and it is the one place the design pushes back on the
existing module boundaries.

A second, smaller friction: `curve/mesh_types.ts` exposes `otherEdge(e)`, which is a
valence-2 idiom — `Mesh.Vertex.otherEdge` throws outright on anything else, with the
comment that it is "only meaningful on 2-valence vertices, which is what the curve solvers
walk." The general-graph solver needs to walk `v.edges` and the authored pairings instead.
`otherEdge` can stay for the existing solvers, but `SolvableVertex` needs pairing-and-level
access added alongside it, and the chain decomposition of §5 has to be computable through
the structural interface without reaching into `mesh/`.

Third: **`CurveSolver.solve()` returns `void`.** The diagnostics of §8 need a channel, and
that is the interface change — `solve()` yields a report. Small, but it touches the
structural interface in `curve/mesh_types.ts` that the existing solvers also implement, so
either they return an empty report or the signature widens with a default. Prefer the
former; a solver that reports nothing should say so explicitly.

Fourth: **`cornerThreshold` leaves the solver.** It is currently in
`ClothoidSolverOptions` (`clothoid.ts:441`) and tested at `clothoid.ts:628` against the
*input polyline* angle `otherVertex(v) − v` — solver-side inference on pre-solve geometry,
which §3 disallows on both counts. It belongs in `Stroker` and other client mesh code, run
at input time, emitting authored level-0 pairings.

That is a behaviour change, not just a move. Pointer-device polylines have genuine corners
that the current 72° test catches, so with maximum continuity as the valence-2 default the
stroker will round them off until the inference is reimplemented client-side. Port it
rather than deleting it.

## 11. Phasing

Ordered so the cheapest experiment carrying the most information runs first.

**Phase 0 — math kernel.** Basis evaluation, differentiation, integration, the `h(i,k)`
map, Gram matrices. Standalone, no mesh, no solver. Verify against `s-power.reduce` in
exact rational arithmetic — the two properties `spower.md` §5 already checks (reproduction
and Hermite contact) are the acceptance test. *Small.*

**Phase 1 — profile only, existing solver.** Implement `CurvatureProfile` over s-power
coefficients stored in the existing `ks` array and run the *current* `ClothoidSolver`
against it unchanged. This is a genuine drop-in: same array, same constraints, same
everything. It immediately answers the quadrature question (does `QUADRATURE_STEPS` drop
to 3–4?) at almost no cost, and with the pair layout `[a_0⁰, a_0¹, a_1⁰, …]` the endpoint
curvatures are literally `ks[0]` and `ks[1]`, which makes `curv_c` exact and trivial even
before any refactor. *Small, high information.* **Gate: if the measured quadrature order
does not clean up here, the premise in §1 is wrong and the rest should not be built.**

Run this phase with `enableG2: true`, which makes the `curvatureConstraint` orientation bug
(§12) blocking rather than deferred — and makes it cheap to settle. With the pair layout
there is no `KSCALE` normalization and no `order − 1` indexing in the way: build a two-edge
mesh in each of the four `isV1` configurations, read `ks[0]` / `ks[1]`, and check they agree
in the through-orientation. That is a test rather than an argument, and it is the first
runtime evidence that G2 does anything at all.

**Phase 2 — vertex blocks and the KKT inner solve.** The architectural decision from §10,
then the banded direct solve with transforms frozen. Compare against Gauss-Seidel on
identical input. *Large.*

**Phase 3 — outer fixed point and diagnostics.** Contraction measurement, divergence
detector, retire the unreachable `disabled` flag. Land the §8 diagnostics record and the
`solve()` signature change here rather than with Phase 4: the contraction factor is the
first thing that needs reporting, and having the channel in place makes Phase 3 itself
easier to debug. *Medium.*

**Phase 4 — levels and stability breaking.** Authored pairing levels with the per-order
union-find; the two stability criteria, the artistic/engineering response split, and
hysteresis; `λ_v` added to the Phase 3 report. Port `cornerThreshold` out to `Stroker` in
the same change so default behaviour does not regress in between. *Medium.*

**Phase 5 — continuation.** `p = 0 → 2`. Should be nearly free given Phase 2. *Small.*

**Phase 6 — general valence.** Node/chain decomposition, group partitions, the Schur
complement. Phase 2 should be written with the partition already in the data model even
though every group starts at size 2, so this phase is decomposition and elimination only,
not a representation change. *Medium.*

**Phase 7 — width as a second profile.** Reuse the Phase 2 inner solve on `w`, no outer
loop. Mostly a matter of whether Phase 2's machinery was written generically over profiles
or hardcoded to `κ`; decide that in Phase 2. *Small if so, large if not.*

## 12. What this does not fix, and what could go wrong

- **The geometry integrals stay transcendental.** Quadrature is not eliminated, only made
  well-behaved. Replacing it entirely is the 2003 paper's separate contribution
  (`spower.md` §8), and it trades exact arclength for *near*-arclength parameterization —
  which `CLAUDE.md` invariant 2 depends on. Out of scope here; a later decision.
- **The basis is global within a segment.** Changing one coefficient moves curvature over
  the whole segment, which removes the poke-one-of-twelve-samples editing the current
  profile allows. Distinct from the chain-level stiffness in §9 — that one lands only in
  engineering mode, where locality is not wanted; this one lands in artistic mode, where it
  is. At `p ≤ 1` a segment has 2–4 coefficients so the granularity loss is small, but if
  direct curvature-profile editing is a workflow, this is a real regression.
- **The outer fixed point is unproven.** Its contraction is asserted from structure, not
  measured. Phase 3 could find it does not contract on realistic strokes, in which case the
  fallback is a full Newton solve on the coupled system — more work, same structure.
- **Endpoint interpolation by similarity transform is assumed.** The alternative is dropping
  it and constraining position explicitly: two nonlinear rows per segment instead of two
  nonlinear scalars. Roughly a wash in count, worse in conditioning (world-unit position
  residuals mixed with angle residuals), so the transform is kept — but it is a fork worth
  a prototype if Phase 3 disappoints.
- **Conditioning at `p ≥ 3` is unmeasured.** Hermite-type bases degrade as contact order
  rises. s-power is much better behaved than a monomial Hermite basis — `σ ≤ 1/4` so terms
  shrink geometrically, and `spower.md` §6 notes the Bernstein conversions are well
  conditioned — but that is an argument, not a measurement.
- **Maximum continuity by default is untested at scale.** Nothing in the current system
  ever ran with every valence-2 joint smooth — `cornerThreshold` has always been culling
  them. A long stroke held to G3 end to end may be better-conditioned than the current
  broken-up one or considerably worse, and §5's `O(V·(p+1)³)` band solve is what has to
  absorb it. Measure on real stroke input before committing to the default.
- **Solver-side breaking is a failure path, and failure paths rot.** With the aesthetic
  criteria gone it should fire almost never, which means it will be exercised almost never
  and its hysteresis and reporting will be wrong when it finally matters. Needs synthetic
  degenerate cases in the test suite from the day it lands, not real strokes.
- **The solve-block decomposition is dynamic.** Refactorization on every active-set change
  is correct but wasteful; splitting a block is a low-rank change to a factored system, so
  updating rather than refactoring is available if it matters. Measure first — at these
  sizes a refactor may simply be cheaper than the bookkeeping.
- **Sign conventions at junctions are load-bearing and currently inconsistent in `src/`.**
  `tangentConstraint` negates when `isV1e1 === isV1e2` (`clothoid.ts:490`) while
  `curvatureConstraint` sets `flip = isV1e1 !== isV1e2` (`clothoid.ts:530`). Those are
  opposite tests for the same question. Working it through, `curvatureConstraint` looks
  like the inverted one: when both edges have `v` as their `v1`, the through-path traverses
  `e1` backwards and `κ` must flip, which is what the tangent version does and the
  curvature version does not. This is on the `enableG2: false` path, so it is unexercised —
  a third reason that path has never worked, alongside the last-sample defect in
  `clothoids.md` §8. Not yet verified at runtime, and fixing it inverts behaviour for
  anyone who has enabled G2. **Blocking for Phase 1**, since that phase runs with G2 on;
  Phase 1 also describes the cheap runtime check that settles it.
- **There are still no tests.** `clothoids.md` §8 lists this as a live gap and every
  convergence claim in this document is either measured ad hoc or inherited. The
  quadrature-order harness used for the §4 table should land as a real test with Phase 1,
  not after.

## 13. References

- `docs/research/spower.md` — the basis, arithmetic, and the coefficient-extraction map.
- `docs/research/clothoids.md` — the current formulation, the solver, and the measured
  quadrature behaviour that motivates §7.
- `docs/research/s-power.reduce` — reference implementation for Phase 0 verification.
- Sánchez-Reyes 2000, *Applications of the polynomial s-power basis in geometry
  processing* (`spower-practical.pdf`) §§2–3 — the working reference.
- Sánchez-Reyes & Chacón 2003, *Polynomial approximation to clothoids via s-power series*
  — the geometry-side approximation deliberately left out of scope in §12.
