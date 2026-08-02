# An s-power clothoid and a Gauss-Newton solver

A design for a second curve/solver pair alongside the existing `Clothoid` /
`ClothoidSolver`: the curvature profile represented in the s-power basis instead of
piecewise-linear samples, and a solver built around what that representation makes
possible. Math background is `docs/research/spower.md` (the basis) and
`docs/research/clothoids.md` (the current formulation and its measurements). Nothing here
is implemented.

The existing clothoid stays. Comparing formulations on identical input is half the point
of the project, and the piecewise-linear profile is the baseline this has to beat.

## 1. Premise

Three problems with the current system motivate this, but they do not have the same root,
and conflating them is how the first draft of this plan built an experiment that could not
test its own hypothesis.

| Problem | Current cause | Kind |
|---|---|---|
| Curvature continuity is a self-referential projection that has never been trusted enough to enable by default | `curv_c` writes `ks`, then marks the curve dirty; `_update()` later recomputes `KSCALE`, which `curv_c` had normalized by. It projects onto a set that moves because of the projection. | **solver** |
| Convergence is by Kaczmarz sweep, and slows as `O(N²)` in chain length | One scalar constraint per joint, projected one at a time, information propagating one joint per sweep. | **solver** |
| The Taylor quadrature can't reach its asymptotic order at low step counts | The profile has 10 interior knots where θ drops to C¹. Below `N ≈ 11` the error is dominated by knot aliasing (`clothoids.md` §4). | **representation** |

Only the third row is about the basis. Rows 1 and 2 are fixable on the *existing*
piecewise-linear profile: share the endpoint sample between adjacent edges instead of
projecting, write the tangent residual analytically, assemble a banded KKT system, and
solve it directly. That would kill the self-reference and the `O(N²)` sweep without
touching the representation.

So the basis change has to be justified on row 3 plus the structural argument below, and
the two effects have to be measured separately. §12 phases them so they are: the
quadrature question is settled with no solver at all, and the solver is written generically
enough that the piecewise-linear profile can be run through it as a control.

The structural argument, independent of the measurements: twelve piecewise-linear samples
per edge give a non-smooth profile with far more degrees of freedom than there are
constraints, no structural relationship between adjacent edges, and no basis-level notion
of endpoint derivative data. An s-power profile is a smooth polynomial whose coefficients
*are* endpoint derivative data. That makes `Gⁿ` continuity a property of the data layout
rather than something a solver converges to, removes the knots, and drops the DOF count to
something commensurate with the constraint count.

## 2. Naming

`spower.md` §9 flags this as the thing to settle before writing any code, and it is:

| Symbol | Meaning | In code |
|---|---|---|
| `s` | arclength — unchanged, `CLAUDE.md` invariant 2 | `s` |
| `u` | normalized segment parameter, `u ∈ [0, 1]` | `u` |
| `σ` | Sánchez-Reyes' symmetric parameter, `σ = (1 − u)·u ∈ [0, 1/4]` | `sym` |
| `κ` | curvature | `k` (matches existing code) |

Every citation of the papers uses their `s` for our `σ`. Translate on the way in.

Two lengths per edge, and they are deliberately different (§5, §6):

| Symbol | Meaning | Status during a solve |
|---|---|---|
| `C_e` | chord length, `\|v2 − v1\|` | **constant** — vertex positions are solver input |
| `L_e` | true arclength, i.e. `KSCALE_e` | **an unknown**, an implicit function of the DOF |

## 3. Representation

### The basis

Order `p`, degree `2p + 1`, `2p + 2` coefficients arranged as pairs
`a_k = (a_k⁰, a_k¹)`, `k = 0..p`:

```
κ(u) = Σ_{k=0..p} [ (1 − u)·a_k⁰ + u·a_k¹ ] · σ^k
```

`σ^k` is a `k`-fold zero at both ends, so `(1−u)σ^k` has a `k`-fold zero at `u = 0` and
`u·σ^k` a `(k+1)`-fold one. Consequently everything from term `k+1` upward has contact of
order `k` at both ends, and truncating at term `k` leaves the endpoint derivatives up to
order `k` untouched (`spower.md` §4). The consequences that matter here:

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
apparatus: `endCurvature`, its orientation flip, and the normalize-then-recompute
self-reference that makes row 1 of §1's table a problem. That row is deleted rather than
fixed.

### Degrees of freedom live on vertices, not edges

This is the structural change, and it is what eliminates the continuity constraints rather
than merely linearizing them.

Each vertex `v` owns a block of `p + 1` scalars: the world-space curvature derivatives at
that point, scaled to be dimensionally commensurate,

```
block(v)_n = (Rᵥⁿ / n!) · dⁿκ/dsⁿ |_v          n = 0 .. p
```

with `Rᵥ` a reference length for the vertex. World units, because curvature continuity is a
geometric condition — the current code already divides by `KSCALE` for exactly this reason.

**`Rᵥ` is built from chord lengths, not arclengths.** Arclength is `KSCALE`, a solve output
(§5); making the DOF scaling depend on it would mean the DOF change meaning between
iterations, so any convergence measured in that metric would be measured in a moving frame.
Chords are input. Take the geometric mean of the incident chord lengths — symmetric, and
better behaved than the arithmetic mean, which is dominated by the longest edge.

No single `Rᵥ` is adequate when incident lengths are disparate: chords of 1 and 100 give
ratios `L_e/Rᵥ` of `0.1` and `10`, which at `p = 2` is a `10⁴` spread in column scale on
the very columns the scaling exists to equalize. The `Rᵥⁿ/n!` form gets the `n`-dependence
right, which is the part a diagonal rescale cannot recover on its own; the residual spread
is then handled by **symmetric diagonal equilibration (Ruiz) of the assembled KKT matrix**,
which is cheap and standard. Treat equilibration as part of the solver, not an optional
polish step.

An edge's own coefficient vector is a linear function of its two endpoint blocks:

```
a_e = M_e · [ block(v1) ; block(v2) ]
```

`M_e` composes three things, all constant because all three depend only on input geometry:
the rescale from `Rᵥ` to `C_e`, an orientation sign, and the triangular `h(i,k)` map from
Taylor data to s-power pairs (`spower.md` §5). Note that `h` mixes both endpoints, so
`a_k⁰` depends on both blocks — but never on a third, so each edge couples exactly two
adjacent vertex blocks.

The orientation sign is `(−1)ⁿ⁺¹`, i.e. **even** orders flip. Under arclength reversal
`s̃ = L − s` the tangent reverses and signed curvature flips with it, but its derivative
does not:

```
κ̃ = −κ        dκ̃/ds̃ = +κ′        d²κ̃/ds̃² = −κ″        dⁿκ̃/ds̃ⁿ = (−1)ⁿ⁺¹ dⁿκ/dsⁿ
```

On the s-power pairs themselves, reversal swaps the pair and negates it, since `σ` is
invariant under `u → 1 − u`:

```
(a_k⁰, a_k¹)  →  (−a_k¹, −a_k⁰)
```

Getting either backwards is easy and quiet. The existing `curvatureConstraint` and
`tangentConstraint` disagree about it today — see §13.

### Valence, pairings, and continuity levels

The mesh already carries arbitrary valence (`Vertex.edges` is unbounded); it is the
solvers that assume chains. Generalizing the block model covers the general graph without
a second mechanism.

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
is a clean reading: **one partition per order, each a coarsening of the next.**

There are `p + 2` union-finds, not `p + 1`. Index them by `n = −1 .. p`:

| `n` | Restrict to pairings with | Components determine |
|---|---|---|
| `−1` | level ≥ 1 | which edge-ends are tangent-coupled — the G1 rows |
| `0 .. p` | level ≥ `n + 2` | which edge-ends share block entry `n` |

The `n = −1` partition is easy to omit and its omission is silent: without it there is no
grouping for level-1 pairings, and the G1 rows at a G1-only junction are simply never
emitted. A group of size `m` in partition `n` needs `m − 1` G1 rows only when `n = −1`
— a spanning tree of that group, not `m(m−1)/2`. Groups of size > 2 mean three or more
edge-ends mutually tangent at a point: representable, geometrically degenerate, rarely
wanted.

Valence 1 is a cap. It has no pairing partner, so its block is governed by energy alone and
KKT stationarity supplies natural boundary conditions (§6 — note these are not what the
current code does). Valence 0 is skipped.

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
r_v = wrap( θ_e1(1) − θ_e2(0) )       θ_e(u) = L_e · ∫₀^u κ du + KTH_e
```

`∫κ` is a linear functional of the coefficients, and `spower.md` §6 gives it exactly in
closed form on the pairs. `L_e` and `KTH_e` are not linear — see §5.

Use the wrapped signed angle difference, not `acos(t1·t2)`. The current residual
(`clothoid.ts:493`) has three distinct problems and only the third is about magnitude:

1. **It is non-differentiable at the solution.** `acos(t1·t2) = |Δθ|`, whose derivative is
   `sign(Δθ)` — finite, but with a kink exactly where the solver is trying to land, so a
   finite difference straddling zero chatters instead of converging. (The gradient does
   *not* diverge: `acos′` blows up like `1/√`, but `d(cos Δθ)/dΔθ = −sin Δθ` vanishes at
   the same rate and the two cancel. An earlier draft of this document claimed divergence;
   the prescription was right and the reason was not.)
2. **The sign is discarded.** `|Δθ|` cannot tell the solver which way to move, which is
   what makes the projection direction guesswork.
3. **It loses about half the mantissa.** Near alignment `1 − t1·t2 ≈ Δθ²/2`, so a `1e-8`
   angle is a `5e-17` perturbation of a quantity near 1 — below double precision. The
   accuracy floor is set by cancellation, and is currently masked by
   `Constraint.threshold` short-circuiting to `0.0`.

The signed angle is analytic, cheap, and already available.

Endpoint position stays satisfied by construction via the similarity transform, as today.

## 5. The remaining nonlinearity, and the solver that follows

With κ in this basis, everything is linear in the DOF except two scalars per edge, both
computed in `_update()`:

```
KTH_e    = atan2(v2 − v1) − atan2(canonical chord)
KSCALE_e = |v2 − v1| / |canonical chord|          — and this is L_e
```

`KSCALE` *is* the arclength: the canonical curve is built with unit arclength, so scaling
it by `KSCALE` gives an edge of arclength `KSCALE`. That means `L_e` in §4's residual is an
unknown, not a constant, and the naive reading — "freeze `KTH` and `KSCALE`, solve a
quadratic program, iterate" — is worse than it looks.

### Why freezing the transforms does not work

`θ_e2(0) = KTH_e2` exactly. So with the transform frozen, edge `e2`'s degrees of freedom do
not appear in the G1 row *at all*. That is not a small-perturbation lag; it is the wrong
sparsity pattern. For a circular arc of turning `Φ`, the canonical chord is
`(e^{iΦ} − 1)/(iΦ)`, whose argument is `Φ/2`, so the true sensitivity of the end tangent to
`Φ` is `1/2` against the frozen row's `1` — contraction `0.5`, about 17 outer
factorizations for five digits. On a front-loaded profile — a hook, a stroke terminus,
anything where the turning happens early — `arg(chord) → Φ` and the true sensitivity → 0,
so the contraction factor → 1 and the iteration stalls. Those shapes are not exotic; a
brush stroker produces them constantly.

The claimed benefit was dimensional reduction, and there isn't any: for a chain `E ≈ V`, so
"confined to `2E` scalars" is `2V` scalars against a state space of `(p+1)V = 2V` at
`p = 1`.

### Gauss-Newton, which costs the same

`KTH_e` and `KSCALE_e` are smooth functions of `a_e`, which depends on exactly the two
adjacent vertex blocks. So `∂KTH_e/∂x` and `∂KSCALE_e/∂x` touch **the same two blocks the
frozen row already touched**. The exact Jacobian has the same sparsity, the same bandwidth,
and the same `O(V·(p+1)³)` factorization cost. It is not a fallback; it is the design.

It is also cheap to form. Writing the canonical curve as `C(u) = ∫₀^u e^{iθ_c(v)} dv` with
`θ_c(v) = ∫₀^v κ`, the derivative with respect to coefficient `a_j` is

```
∂C(1)/∂a_j = ∫₀¹ i · (∂θ_c/∂a_j)(v) · e^{iθ_c(v)} dv
```

which accumulates in the quadrature loop that is already running — `2p + 2` extra complex
accumulators, and `∂θ_c/∂a_j` is the same closed-form `∫κ` from §4 applied to a basis
function. `KTH` and `KSCALE` are then `−arg` and `1/|·|` of `C(1)` scaled by the chord, so
their derivatives follow from `∂C(1)` by two lines of complex arithmetic.

### The linear system

Each Newton step solves

```
[ H   Jᵀ ] [ Δx ]   [ −(g + Jᵀλ) ]
[ J  −δI ] [ Δλ ] = [ −c         ]
```

with `H` block tridiagonal in vertex index (each edge couples adjacent blocks) and `J` one
row per G1 pairing touching two adjacent blocks. What is constant and what is not:

| | Constant across the solve? | Why |
|---|---|---|
| `M_e`, `Rᵥ` | **yes** | built from chord lengths, which are solver input |
| `H` | **yes** | §6's weights use `C_e`, not `L_e` — assemble and Gram-scale once |
| `J` | no | carries `∂KTH/∂x`, `∂L_e/∂x` — reassembled per Newton step |
| the KKT factorization | no | refactored per Newton step, `O(V·(p+1)³)` each |

Using chord length rather than arclength in the energy weights is what buys the `H` row,
and it is not a fudge: the energy is a regularizer, and any equivalent quadratic form will
do. `C_e ≤ L_e` with agreement to `O(Φ²)`, and where they diverge — an edge turning toward
closure — the chord-based weight `1/C_e` grows, so the regularizer *strengthens* on the
degenerate edge. The arclength-based weight does the opposite: `L_e → ∞` as the chord
collapses, so `1/L_e → 0` and the regularizer switches itself off on exactly the edge that
is running away. The sign of that feedback loop is the whole argument.

**`LDLᵀ` needs the `−δI`.** The (2,2) block of a KKT system is structurally zero, and the
interleaved ordering that makes the system banded puts those zeros on the diagonal, where
unpivoted `LDLᵀ` breaks down. Bunch–Kaufman would fix it but destroys the static ordering
and the band. Instead use the standard quasi-definite regularization: perturb to `−δI` with
a small `δ`, which makes unpivoted `LDLᵀ` provably stable with *any* symmetric ordering,
then recover accuracy with one or two steps of iterative refinement against the unperturbed
residual. Band preserved, static symbolic factorization preserved, `O(V·(p+1)³)` preserved.

A closed loop makes the system cyclic banded. This is **not** a rank-one correction — the
wrap-around couples a full block plus its G1 row, so the correction has rank `≈ p + 2` per
profile, and Sherman–Morrison–Woodbury on an indefinite system is numerically fragile
regardless. Promote an arbitrary vertex to a cut point and use bordering, or just treat the
cut vertex as a node (below) and let the Schur complement handle it.

Globalize with a line search on `‖r‖² + energy` — Gauss-Newton on a residual with a genuine
kink risk (§4 removes the worst of it, but wrap-around remains) is not unconditionally
convergent, and the step is cheap enough that a backtracking search costs little.

### Chain decomposition and substructuring

Decompose the mesh into **nodes** (valence ≠ 2) and **chains** (maximal paths of valence-2
vertices between them). Then:

- With every node fully split into singleton groups — §3's default at valence ≥ 3, so this
  is the common case — no block is shared across a chain boundary, and the chains are
  decoupled up to the two caveats below. Each is an independent banded solve of `O(len)`,
  embarrassingly parallel.
- With some node groups shared, the shared blocks are the cross-chain coupling. Eliminate
  each chain's interior onto its endpoints and solve the resulting Schur complement on node
  DOF only — standard substructuring / nested dissection. The complement is
  `O(#shared groups × (p+1))`.

Two things break naive decoupling and both are easy to miss:

1. **Level-1 pairings share no block but do emit a G1 row.** A junction authored G1-only
   couples two chains through a single row that touches both. The interface is the union of
   shared blocks *and* level-1-paired edge-ends, not just the former.
2. **`Rᵥ` is shared.** Every edge incident to a node scales its blocks by that node's `Rᵥ`,
   so the chains are coupled through the scaling even when nothing else is shared. This is
   benign only because §3 pins `Rᵥ` to chord lengths — if it were arclength-derived it would
   be a genuine, and very confusing, cross-chain dependency.

An isolated all-valence-2 cycle has no nodes at all and needs an arbitrary cut point
promoted to one.

**Two decompositions, and they are not the same one:**

| | Basis | Lifetime | Used for |
|---|---|---|---|
| Topological chains | maximal valence-2 runs | static, from the graph | stroking: width profiles, cap/join placement, offsetting |
| Solve blocks | maximal runs of *coupled* DOF | dynamic — re-derived when a level changes | factorization structure |

Solve blocks are topological chains further split at authored corners, and further still
whenever §8's stability breaking lowers a level mid-solve. Breaking a corner splits a solve
block but *not* a topological chain: position is still shared there, only the DOF coupling
parts. Conflating the two will fail exactly when a level changes mid-solve — rare under
§8's policy, which is precisely why it will not be caught by ordinary use.

## 6. Energy

Nothing in the current system minimizes energy; the Kaczmarz step's minimum-norm property
keeps it near the all-`0.001` seed, which is an accidental `‖Δk‖₂` regularizer. Make it
explicit:

```
E = Σ_e [ (1/C_e)·∫₀¹ κ′(u)² du  +  ε·C_e·∫₀¹ κ(u)² du ]
```

with `C_e` the chord length, for the reasons in §5. The first term is minimum-variation,
which for brush strokes is generally preferable to the elastica `∫κ²` — it does not bias
toward circular arcs. That weighting is an *artistic* default, not a universal one (§9).
Engineering fairing often wants the opposite bias — arcs and lines where they will do —
and classical CAD fairing frequently targets curvature-plot monotonicity, a different
objective again. Keep the two coefficients as mode parameters; the Gram matrices are
shared either way, so this costs nothing structurally.

### What `ε` is and is not for

Not well-posedness. The null space of `∫(κ′)²` is constant curvature, and the constant mode
is *not* in `ker J`: perturbing κ by a constant `δ` changes the G1 residual by `L_e1·δ ≠ 0`.
So `H` is already positive definite on `ker J`, which is all the KKT system needs. `ε` earns
its place for three other reasons:

- **Isolated edges and short blocks**, where there is no G1 row to supply that definiteness.
- **Conditioning**, which is a real concern here for a reason §13 previously got backwards.
- **Edit locality along a chain**, which is the one that actually matters in artistic mode.
  `∫(κ′)²` alone is a fourth-order operator whose Green's function decays only
  algebraically; adding `ε∫κ²` makes it a beam on an elastic foundation, whose response
  decays exponentially with rate set by `ε^{1/4}`. **`ε` is the locality knob**, and it is
  the answer to the locality worry in §13 — not a workaround bolted on beside it.

The null space is also not one-dimensional overall: it is one dimension **per coupled solve
block**, and the fully-split default at junctions guarantees many of them.

### `ε` is dimensional

`∫(κ′)² ds` scales as `1/L³` and `∫κ² ds` as `1/L`, so `ε` carries units of `1/L²`. A fixed
numeric `ε` makes the solver scale-dependent: the same stroke drawn ten times larger
converges to a different shape. Parameterize it as `ε = (α/C_e)²` with `α` dimensionless,
or equivalently express the whole energy in the edge's own normalized parameter and let the
`C_e` powers in the weights do the work. Either way the shipped constant must be `α`.

### Natural boundary conditions are not what the code does today

At a free end, stationarity of `∫(κ′)²` gives `κ′ = 0`, not `κ = 0`. The current code
zeroes the endpoint curvature at a cap. That is a visible change in cap shape — a straight
run-out becomes a constant-curvature run-out — and it should be introduced deliberately with
a before/after, not slipped in as a consequence of the energy. If the current behaviour is
wanted, it is a `κ = 0` boundary row, cheap to add, and then it is an authored choice.

### Gram matrices

Both are constant `(2p+2)²` matrices per order, computable once in exact rational
arithmetic — the entries are `∫₀¹ σ^{j+k}(1−u)^α u^β du`, i.e. Beta functions. For `p ≤ 3`
these are at most 8×8. Compute them offline, check them in as literals, and record `κ(H)`
for each order.

s-power is not an orthogonal basis, so unlike shifted Legendre it does not diagonalize
`∫κ²`. It buys structural endpoint continuity instead. The two bases are related by a
constant invertible matrix per order, so if the local blocks condition badly that map is
available as a preconditioner — and §13 now expects them to, so budget for measuring it.

The `Rⁿ/n!` scaling in §3 is not cosmetic. Raw `dⁿκ/dsⁿ` values differ by orders of
magnitude across `n`, which wrecks the column scaling of `J` and `H` alike.

## 7. Evaluation and quadrature

The geometry integrals stay transcendental — `κ` polynomial does not make `∫cos(∫κ)`
polynomial — so `quadrature()`'s Taylor stepping is retained unchanged in structure. What
changes is that it can finally reach its asymptotic order.

`clothoids.md` §4 records the measurement: with the piecewise-linear profile the observed
order oscillates between ~1 and ~4 and the error is not monotonic below `N ≈ 11` (three
steps beats four), because 10 interior knots dominate. A polynomial κ has no knots, so
`O(ds³)` should hold from the first step.

### Clean order is not the same as few steps

`O(ds³)` says the error is `C·ds³`; it says nothing about `C`, and the fourth-order terms
below show `C` growing with the cube of the turning `Φ` over the segment. At `Φ = π` with
`N = 3` that is on the order of `5e-2` — an order of magnitude *worse* than the `2.1e-3`
the current piecewise-linear profile manages, despite the cleaner order. "3–4 steps" is
only a defensible default with a bound on turning per segment.

The right criterion is therefore **turning per step, not steps per segment**, and it is
free: `∫κ` is available in closed form, so the step schedule can be chosen for equal
turning directly rather than by trial refinement.

One catch: a coefficient-dependent step schedule makes the residual discontinuous in `x`,
which breaks the Jacobian in §5. Freeze the schedule at the top of each Newton step —
recompute it from the current iterate, hold it fixed while differentiating. This is the
same lagging trick as a lagged mesh in FEM and costs nothing.

### The fourth-order terms

They need `κ″`, which this basis supplies exactly via §6's differentiation rule:

```
dx += ( (κ³ − κ″)·sin θ − 3κ·κ′·cos θ ) · ds³/24
dy += ( (κ″ − κ³)·cos θ − 3κ·κ′·sin θ ) · ds³/24
```

(These share `(κ³ − κ″)` and `3κκ′`, so two extra scalars per step.)

They also double as an **error estimator**, which is what the engineering mode needs and
what makes the equal-turning schedule verifiable rather than assumed: the magnitude of the
`ds³` terms bounds the local truncation error of the scheme without them. Exact `κ″` from
§6's differentiation rule is what makes the estimate trustworthy rather than a
finite-difference guess.

All three `CurvatureProfile` functions are exact and cheap: `κ` by Horner in `σ`, `κ′` by
§6's differentiation rule, `∫κ` by §6's integration rule. That matters — `clothoids.md` §3
notes that an inconsistent `integral` bends the curve with no visible error in the
curvature plot, which is precisely the defect that was found and fixed in the
piecewise-linear profile, and precisely the defect still open in `circleArc` (§13).

## 8. Constraint breaking and diagnostics

The current mechanism is `cornerThreshold` on the input polyline angle, decided before the
solve runs, blind to segment lengths and to what the rest of the chain is doing. Under the
policy in §3 it does not get a principled replacement inside the solver — it gets moved
out. Breaking splits into two jobs with different owners.

### What the solver may break: stability only

Lowering a pairing level is mechanical — decrement it, re-derive the per-order partitions,
refactor. What needs a criterion is *when*, and under §3's invariant the only admissible
reasons are numerical:

**1. Rank deficiency in `J`.** Note that §3's spanning-tree construction makes the G1 rows
*structurally* full rank by design, so "conflicting or redundant rows" cannot arise — the
earlier draft's stated cause was designed out by its own §3. What remains is **numerical**
near-dependence: rows that become nearly parallel when incident edges collapse in length.

Test it where it is cheap and honest: a column-pivoted QR of the *local* `J` block for each
node group, which is at most a few rows by `2(p+1)` columns, with the threshold relative to
that block's largest singular value. `O(1)` per node, scale-invariant.

Do **not** test it by watching pivots during the global factorization. Pivot magnitude is
not a rank test — Kahan's matrix is the standard counterexample — and the ordering that
would make a pivot interpretable (multipliers last) is incompatible with the interleaved
banded ordering §5 needs. Separating the rank test from the factorization resolves that
conflict; the earlier draft required both and could have neither.

**2. Chord degeneracy.** `|canonical chord| → 0`, total turning approaching closure. This is
a genuine geometric limit: `KSCALE → ∞` and the similarity transform stops being
well-conditioned. Measured directly in `_update()`.

**3. Newton divergence.** Not a one-step iterate ratio — the Gauss-Newton iteration is
highly non-normal, so a single step can grow while the sequence converges, and a single step
can shrink while it diverges. Use monotone decrease of the merit function over a window,
plus a geometric fit of the residual over at least three iterations, plus a line-search
failure count. Three signals, none of which is a single-step test.

Note what is *not* on this list: `|λ_v|`. A large multiplier means the joint is expensive to
hold, which is a judgement about whether a corner was wanted there — modelling, not
stability. Under this policy the solver may not act on it.

**In engineering mode the solver may not break at all.** Detection is the same; the response
is not. The value of asking for G4 is the guarantee, so quietly delivering G2 at one joint
produces a model that looks right and is not — worse than no answer. Degrade and report in
artistic mode, fail and report in engineering mode (§9). This is a mode flag on the solve,
not a per-joint decision.

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

The solver's contribution is a better signal than the client can compute on its own. `λ_v`
falls out of the KKT solve for free: the force required to hold tangent continuity at that
joint, in energy units, accounting for segment lengths and for what the rest of the chain is
doing. That is strictly better corner evidence than an input-polyline angle.

But note the timing, because it constrains how it can be used. `cornerThreshold` runs at
input time; `λ_v` exists only *after* a solve. So acting on it means solving, reporting,
and having the client author a level on the *next* edit — an explicit two-pass loop with a
user-visible step in between. It must not become an automatic re-solve, because that is the
solver-side inference §3 forbids, merely laundered through the client. Suggestion surfaced
in the UI, not decision.

Same for graded breaking, if it is wanted: `Constraint.wlst` already exists and is filled
with `1.0` — the natural hook for IRLS/Huber weighting on the G1 rows, with the scale set
from `δ = 1.345·MAD(residuals)` rather than by eye. But that is softening a constraint the
client asked for, so the client sets the policy.

### Diagnostics

Engineering mode makes this load-bearing: if the answer to a degenerate configuration is
"refuse," the user has to learn *what* and *where* well enough to fix it. But the same
records are the debugging surface for development, which is the stronger reason to build
them properly — every convergence claim in this document is currently checked by eyeballing
the demo (§13).

A diagnostic is a located, typed record, not a boolean: **where** (vertex, edge, or chain),
**what condition**, the **measured value against its threshold**, and **what the solver did**
— nothing, degraded from level `k` to `j`, or refused.

| Condition | Measured from | Meaning |
|---|---|---|
| chord degeneracy | `\|canonical chord\|` in `_update()` | total turning approaching closure |
| Newton not converging | merit-function history + geometric fit | §8 criterion 3 |
| line search failing | backtrack count per step | step direction is not a descent direction |
| G1 rank deficiency | pivoted QR of the local `J` block | near-parallel tangent rows at a junction |
| ill-conditioned `H` | condition estimate of the `(1,1)` block | usually a degenerate chord length |
| refinement not converging | residual after iterative refinement | the `δ` regularization is too large for this system |
| level lowered | the solve record | artistic mode only; engineering mode refuses instead |
| large multiplier | `λ_v` | not a fault — corner evidence for the client (above) |

Three rules make these useful rather than decorative:

1. **Measure unconditionally, not only on failure.** A slow convergence rate is not a fault
   but is exactly what you want to see when a stroke feels sluggish. Severity gates what a
   client *surfaces*; it must not gate what the solver *records*.
2. **It is nearly free, which is what makes (1) affordable.** Nearly every row above is a
   byproduct of a computation already happening — the merit history from the line search,
   the chord from `_update()`, the refinement residual from the refinement. The two that are
   not free (the local QR, the `H` condition estimate) are `O(1)` per node.
3. **It must be deterministic.** Identical input, identical diagnostics. Same requirement as
   the hysteresis rule above and for the same reason: a report that varies run to run is
   worse than none when bisecting.

Per-iteration history — merit, step norm, largest multipliers — is what actually diagnoses
non-convergence, but it allocates, so keep it opt-in behind a trace flag rather than
always-on.

This also gives `JointParams.disabled` a real replacement. `clothoids.md` §8 records it as
written, correct, and unreachable because nothing detects divergence.

## 9. Degree continuation and domain modes

`p = 0 → 1 → 2`. Continuation is a **warm start, not an embedding**, and the distinction is
not pedantic.

In *pair* coordinates, a degree-`2p+1` polynomial has zero for its order-`p+1` pair, so
appending a zero pair leaves κ unchanged. But the DOF are *blocks* — Taylor data — and a
degree-`2p+1` polynomial generally has a nonzero `κ^(p+1)` at its endpoints. So:

- Zeroing the new **block entry** does *not* preserve the curve; it forces
  `κ^(p+1) = 0` at both ends, which is a different degree-`2p+3` interpolant.
- Preserving the curve requires a specific nonzero `f_{p+1}` at each edge-end — and that
  value is *edge-specific*, while at the new default ceiling `p + 3` the two edges at a
  joint must **share** it. In general they want different values, so the order-`p` solution
  is not representable in the order-`p+1` state space at all.

And even where it is representable, a feasible point is not a minimizer: the order-`p+1`
energy will use the new freedom.

So continuation is worth doing — the warm start is much better than the `0.001` seed, and
compare `changeOrder`, which resamples and moves every value — but the claim is "fewer
Newton steps," not "already solved." Report the iteration count on the second rung; if it
is not markedly lower than a cold start, the continuation ladder is not earning its
complexity.

### Order is a domain setting, not a quality dial

The rungs serve different applications, and the applications do not overlap:

| | Artistic — graphic design, brush strokes | Engineering — CAD |
|---|---|---|
| Order | `p = 0..1`, G2–G3 | `p = 2`, G4 |
| Solver may lower a level for stability (§8) | yes, degrade and report | **no — fail and report** |
| Energy (§6) | minimum-variation, no arc bias | elastica or fairing-specific weighting |
| Quadrature (§7) | equal-turning, coarse | tolerance-driven adaptive |
| Edit locality | matters — tune `ε` (§6) | does not matter |

`p = 1` (cubic, G3) is the expected shipping default, `p = 0` the continuation seed and the
fallback when Newton struggles, `p = 2` (quintic, G4) opt-in for engineering work.

Raising `p` makes the chain stiffer as well as smoother — at G4, `κ`, `κ′` and `κ″` are
coupled across every joint, so a disturbance travels much further before the energy term
damps it. That lands only in the engineering case, where locality is not a requirement.

What does follow is that `p` is a per-document (or per-chain) setting with a cluster of
other defaults attached, not a single number to turn up for quality. Treat the rows above as
one choice, not five.

## 10. Scalar profiles, and width

The reason to make chains first-class rather than an implementation detail: a width profile
`w(s)` along a chain is another scalar profile with the same structure as `κ`. Same s-power
basis, same continuity-by-shared-block at nodes, same `∫(w′)²` energy for smooth taper, same
banded solve, same authored pairing levels — a hard width discontinuity is a level-0 pairing
on the width profile, independent of the level `κ` carries at that joint.

Write the solver generically over *scalar profiles on chains*, with a per-profile
description of how DOF map to coefficients. Three things fall out of that: width reuses the
machinery, `κ`'s geometric coupling in §5 becomes `κ`'s peculiarity rather than the general
case, and — the one that matters for §1 — **the existing piecewise-linear profile becomes
just another profile**, which is what makes the control experiment in §12 cheap instead of a
second implementation.

Width is *not* strictly easier than curvature, though. Three real differences:

1. **It has inequality constraints.** Unconstrained `∫(w′)²` minimization overshoots
   between data points, and negative width flips the offset inside out. This is a bounded
   variable QP (`w ≥ w_min`), not the linear KKT solve `κ` uses. Bound constraints are the
   easy case — active-set on bounds only, warm-startable, each iteration the same banded
   solve with some variables fixed — but it is a different solver, not the same one.
2. **It couples to geometry after all.** The offset curve at half-width `h` is regular only
   where `h·|κ| < 1`; at `h = 1/κ` the offset cusps. So `w` and `κ` are not independent, and
   a stroker that solves them separately can produce a self-intersecting outline from two
   individually valid profiles. Detecting it is cheap (`κ` is available in closed form);
   deciding what to do about it is a stroker policy question, not a solver one.
3. **It needs a data term.** Curvature is determined by the constraints and the energy;
   width is *given*, per-sample, by pressure input. So the width objective is
   `Σ‖w(sᵢ) − pᵢ‖² + smoothing`, a least-squares fit, not pure energy minimization. That is a
   different assembly (a diagonal data block added to `H`), even though the smoothing term is
   shared.

Width also has its own `L_e` question: `w`'s domain is arclength, which is `KSCALE`, which is
a `κ` output. Solve `κ` first, then `w` on the resulting arclengths — one-directional, no
outer loop, but an ordering constraint worth writing down.

## 11. Fitting into the codebase

New `SPowerClothoid extends Curve` (`Curve` is an abstract class, not an interface, and
`Edge.curve` serializes as `abstract(Curve)` — the inheritance chain is load-bearing) and
`SPowerSolver implements CurveSolver` (that one genuinely is an interface). Selected at
runtime the same way the others are — `Mesh.CurveCls`/`SolverCls`,
`Mesh.switchSplineType()`, `StrokerOptions`. `curve/` still must not import `mesh/`.

One genuine architectural friction: **the DOF no longer live on the curve.** `Curve` and
`CurveSolver` both assume the curve owns its parameters, and `Clothoid.ks` is the model for
that. Vertex-owned blocks need somewhere to live and some way for a solver in `curve/` to
reach them through the structural interfaces in `curve/mesh_types.ts`. Options:

- extend `SolvableVertex` with an opaque solver-data slot;
- have `SPowerSolver` own a side table keyed by vertex;
- keep blocks on edges and re-impose sharing as an explicit constraint (defeats the point).

This needs deciding before the first solver phase, and it is the one place the design pushes
back on the existing module boundaries.

Second: `SolvableVertex` needs pairing-and-level access, and the chain decomposition of §5
has to be computable through the structural interface without reaching into `mesh/`. The
general-graph solver walks `v.edges` and the authored pairings rather than `otherEdge`.

`curve/topology.ts` is where that lives. It currently holds only `walk`, hoisted out of
`bezier.ts` and `bspline.ts`, which had divergent copies under the same doc comment — the
b-spline one went through `v.otherEdge(e)`, which throws on any valence other than 2, so
`BSplineSolver.solve()` threw on *any* open chain rather than only at junctions. Never
noticed because `BSplineSolver` appears in `src/` only in the barrel export. The three
solvers now behave consistently off valence 2:

| | Behaviour off valence 2 |
|---|---|
| `ClothoidSolver` | skips `edges.length !== 2` (`clothoid.ts:617`) |
| `BezierSolver` | `walk` returns `v` at valence 1; picks an arbitrary incident edge at valence ≥ 3 |
| `BSplineSolver` | same, since both now call the shared `walk` |

`topology.ts` is the module §5's node/chain decomposition belongs in, which is why `walk`
was given a home rather than duplicated a third time. Note that `walk` itself does not
generalize — "the other edge" is undefined at valence ≥ 3 and it resolves the ambiguity by
insertion order. The general-graph solver consults authored pairings instead; `walk` stays
for the bezier solver, which is a chain-only construction by design.

That leaves `otherEdge` with **no callers at all** — only the declarations in
`curve.ts:33` and `mesh_types.ts:13` and the implementation at `mesh.ts:76`. It is also
mis-declared: both interfaces say `| undefined` while the implementation throws, so the
next caller to guard the contract correctly will still get an exception. Either fix the
implementation to return `undefined` or drop the member from the structural interfaces;
this plan needs neither, so it is a cleanup, not a dependency.

Third: **`CurveSolver.solve()` returns `void`.** The diagnostics of §8 need a channel, and
that is the interface change — `solve()` yields a report. Small, but it touches the
structural interface the existing solvers also implement, so either they return an empty
report or the signature widens with a default. Prefer the former; a solver that reports
nothing should say so explicitly.

Fourth: **`cornerThreshold` leaves the solver.** It is currently in `ClothoidSolverOptions`
(`clothoid.ts:443`) and tested at `clothoid.ts:628` against the *input polyline* angle
`otherVertex(v) − v` — solver-side inference on pre-solve geometry, which §3 disallows on
both counts. It belongs in `Stroker` and other client mesh code, run at input time, emitting
authored level-0 pairings.

That is a behaviour change, not just a move. Pointer-device polylines have genuine corners
that the current 72° test catches, so with maximum continuity as the valence-2 default the
stroker will round them off until the inference is reimplemented client-side. Port it rather
than deleting it.

Fifth, minor but a trap for anyone reading the current code as a template:
`QUADRATURE_STEPS` (`clothoid.ts:195`, value 19) is a module-private `const`, not an option.
Varying it is a source edit, which is why §12's quadrature experiment is standalone rather
than a knob turned on the existing solver.

## 12. Phasing

Ordered so the cheapest experiment carrying the most information runs first. The first two
phases are gates and neither of them involves the solver.

**Phase 0a — the quadrature gate.** A standalone harness, no mesh, no solver, no
`Clothoid`: fix a coefficient vector, integrate the Taylor scheme at `N` steps, compare
against `N = 2000`, and read the observed order — for an s-power profile and for a
piecewise-linear one, at matched total turning. The ad-hoc harness that produced §4 of
`clothoids.md` already does exactly this for the piecewise-linear case; land it in the repo
and extend it. *Very small.* **Gate: if the order does not clean up and the constant at low `N` is not
competitive at realistic turning, the representational half of §1 is wrong and the rest
should not be built.**

Doing this without a solver is the point. Routing it through `ClothoidSolver` would confound
it five ways — the Kaczmarz min-norm step scales coefficient columns by `~4^{-k}` so
essentially all correction lands in `a_0` and the high-order coefficients never leave the
seed; the `1e-3`-scale finite difference then sits below the quadrature's own error, so the
experiment would be measuring the thing it is trying to measure through itself; and
`order = 12` slots means `p = 5`, an order §13 explicitly cannot vouch for.

**Phase 0b — math kernel.** Basis evaluation, differentiation, integration, the `h(i,k)`
map, the reversal rule, Gram matrices. Verify against `s-power.reduce` in exact rational
arithmetic — the two properties `spower.md` §5 already checks (reproduction and Hermite
contact) are the acceptance test, plus a round-trip of the pair reversal. *Small.*

**Phase 1 — fix the sign bug and land tests.** Settle `curvatureConstraint`'s orientation
(§13) at runtime: build a two-edge mesh in each of the four `isV1` configurations, enable
G2, read the endpoint curvatures, check they agree in the through-orientation. Land the
Phase 0a harness as a real test at the same time. This is the first executable test in the
repository and it unblocks everything that touches sign conventions. *Small.*

Note that the "profile-only drop-in" of an earlier draft is not available: `endCurvature`
reads `ks[order − 1]`, `curvatureConstraint` divides by `ks[KSCALE]`, `changeOrder` would
resample an s-power vector as if it were samples, and `QUADRATURE_STEPS` is private. Any
version of that experiment is a modification, not a substitution, so Phase 0a takes its
place.

**Phase 2 — vertex blocks and one Newton step.** The architectural decision from §11, the
generic scalar-profile abstraction from §10, `M_e`, `H` assembly, equilibration, the
quasi-definite `LDLᵀ` with refinement. Solve with the Jacobian *frozen* first — not as a
shipping scheme, but because it is the same code minus the derivative accumulators and it
isolates assembly bugs from Jacobian bugs. *Large.*

**Phase 3 — the exact Jacobian.** `∂C(1)/∂a_j` in the quadrature loop, verified against
finite differences on a single edge before anything else. Line search, merit function,
convergence. Compare iteration counts against Phase 2's frozen version — the prediction from
§5 is roughly 17 frozen iterations against 3–4 Newton ones on an arc, and no convergence at
all frozen on a front-loaded profile. *Medium.*

**Phase 4 — the control experiment.** Run the *existing* piecewise-linear profile through
the Phase 2/3 solver, using §10's generic profile interface. This is what separates "the
banded KKT and the energy were the win" from "the basis was the win," and §1 is not settled
without it. Cheap if §10's abstraction held; if it did not, that is worth knowing too.
*Small if Phase 2 was written generically.*

**Phase 5 — diagnostics.** The §8 record and the `solve()` signature change. Later than the
earlier draft placed it, because Phase 3's convergence work is what tells you which
measurements are actually diagnostic. *Medium.*

**Phase 6 — levels and stability breaking.** Authored pairing levels with the `p + 2`
union-finds; the three stability criteria, the artistic/engineering response split, and
hysteresis; `λ_v` in the report. Port `cornerThreshold` out to `Stroker` in the same change
so default behaviour does not regress in between. *Medium.*

**Phase 7 — continuation.** `p = 0 → 2`. Report the iteration-count improvement (§9); the
phase is only worth keeping if there is one. *Small.*

**Phase 8 — general valence.** Node/chain decomposition, the per-order partitions, the Schur
complement including the two coupling caveats in §5. Phase 2 should be written with the
partitions already in the data model even though every group starts at size 2, so this phase
is decomposition and elimination only, not a representation change. *Medium.*

**Phase 9 — width as a second profile.** The bounded-variable QP, the data term, the
`h·|κ| < 1` check (§10). *Medium.*

## 13. What this does not fix, and what could go wrong

- **The geometry integrals stay transcendental.** Quadrature is not eliminated, only made
  well-behaved. Replacing it entirely is the 2003 paper's separate contribution
  (`spower.md` §8), and it trades exact arclength for *near*-arclength parameterization —
  which `CLAUDE.md` invariant 2 depends on. Out of scope here; a later decision.
- **The basis is global within a segment.** Changing one coefficient moves curvature over
  the whole segment, which removes the poke-one-of-twelve-samples editing the current
  profile allows. At `p ≤ 1` a segment has 2–4 coefficients so the granularity loss is
  small, but if direct curvature-profile editing is a workflow, this is a real regression.
  Distinct from — and much less important than — the chain-level locality governed by `ε`
  (§6), which is the one that will actually be felt while dragging.
- **Conditioning at higher `p` is likely to get *worse*, not better.** An earlier draft
  argued the opposite from `σ ≤ 1/4`, which is an argument about *evaluation* stability and
  points the wrong way for the Gram matrix: basis functions decaying like `4^{-k}` give a
  Gram diagonal decaying like `16^{-k}`, a ratio of ~462 at `p = 2` before any off-diagonal
  degeneracy. §6's Bernstein or Legendre conversion is the preconditioner, and at `p ≥ 2` it
  should be assumed necessary until measured otherwise.
- **Gauss-Newton is not unconditionally convergent.** §5's argument establishes that it is
  affordable and has the right sparsity, not that it converges from any start. The line
  search and the `p = 0` fallback are the mitigations; if both prove insufficient, the next
  step is a trust region, which preserves the structure.
- **Closed loops have a turning-number branch problem** that nothing in this design
  addresses. A closed chain admits solutions with different total turning numbers, the
  energy does not select among them, and Newton will land in whichever basin the seed is in.
  Dragging a vertex can jump basins, which reads as a sudden flip. Needs an explicit turning
  number, either authored or carried as solve state — decide before Phase 8.
- **Endpoint interpolation by similarity transform is assumed.** The alternative is dropping
  it and constraining position explicitly: two nonlinear rows per segment instead of two
  nonlinear scalars. Roughly a wash in count, worse in conditioning (world-unit position
  residuals mixed with angle residuals), so the transform is kept — but it is a fork worth a
  prototype if Phase 3 disappoints.
- **Maximum continuity by default is untested at scale.** Nothing in the current system ever
  ran with every valence-2 joint smooth — `cornerThreshold` has always been culling them. A
  long stroke held to G3 end to end may be better-conditioned than the current broken-up one
  or considerably worse, and §5's `O(V·(p+1)³)` band solve is what has to absorb it. Measure
  on real stroke input before committing to the default.
- **Solver-side breaking is a failure path, and failure paths rot.** With the aesthetic
  criteria gone it should fire almost never, which means it will be exercised almost never
  and its hysteresis and reporting will be wrong when it finally matters. Needs synthetic
  degenerate cases in the test suite from the day it lands, not real strokes.
- **The solve-block decomposition is dynamic.** Refactorization on every level change is
  correct but wasteful; splitting a block is a low-rank change to a factored system, so
  updating rather than refactoring is available if it matters. Measure first — at these
  sizes a refactor may simply be cheaper than the bookkeeping.
- **Sign conventions at junctions are load-bearing and currently inconsistent in `src/`.**
  `tangentConstraint` negates when `isV1e1 === isV1e2` (`clothoid.ts:490`) while
  `curvatureConstraint` sets `flip = isV1e1 !== isV1e2` (`clothoid.ts:530`). Those are
  opposite tests for the same question. Worked through all four configurations,
  `curvatureConstraint` is the inverted one: when both edges have `v` as their `v1`, the
  through-path traverses `e1` backwards and `κ` must flip, which is what the tangent version
  does and the curvature version does not. The correct line is
  `const flip = isV1e1 === isV1e2;`. This has been derived independently twice but not yet
  confirmed at runtime, which is Phase 1's job. It is only reachable when `enableG2` is
  `true` — off by default — which is why it has gone unnoticed, and fixing it inverts
  behaviour for anyone who has turned G2 on.
- **`circleArc` is inconsistent and is being removed, so do not measure against it.**
  `curvature`/`dCurvature` use `ds = 1/klen` while `integral` uses `1/(klen−1)`, an 8.3%
  error in θ (`clothoids.md` §8). Since the profile is slated for deletion along with the
  b-spline curve type, it is a won't-fix — but it is exactly the shape someone would reach
  for as a quadrature reference, and it would silently bias Phase 0a. Use an analytically
  integrable profile constructed in the harness instead.
- **There are still no tests.** `clothoids.md` §8 lists this as a live gap and every
  convergence claim in this document is either measured ad hoc or inherited. Phase 1 is
  where that changes.

## 14. References

- `docs/research/spower.md` — the basis, arithmetic, and the coefficient-extraction map.
- `docs/research/clothoids.md` — the current formulation, the solver, and the measured
  quadrature behaviour that motivates §7.
- `docs/research/s-power.reduce` — reference implementation for Phase 0b verification.
- Sánchez-Reyes 2000, *Applications of the polynomial s-power basis in geometry
  processing* (`spower-practical.pdf`) §§2–3 — the working reference.
- Sánchez-Reyes & Chacón 2003, *Polynomial approximation to clothoids via s-power series*
  — the geometry-side approximation deliberately left out of scope in §13.
