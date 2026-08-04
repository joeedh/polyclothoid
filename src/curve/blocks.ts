/**
 * Vertex-owned curvature DOF and the per-edge transform that turns them into s-power
 * coefficients — §3 and §6 of `docs/plans/spower-solver.md`.
 *
 * The unknowns live on *vertices*, not edges. A vertex `v` owns one block of `p + 1`
 * scalars,
 *
 * ```
 *   block(v)_n = (Rᵥⁿ / n!) · dⁿκ/dsⁿ |_v ,      n = 0 … p
 * ```
 *
 * where `κ` is world signed curvature, `s` is world arclength, and `Rᵥ` is a reference
 * length that makes the block dimensionally uniform. Every edge then reads the two blocks
 * at its ends and reconstructs its own coefficient vector, `a_e = M_e · [block(v₁); block(v₂)]`.
 * Because that reconstruction is two-point Hermite interpolation, `Gᵖ` continuity at every
 * vertex is *structural*: two edges sharing a vertex read the same derivatives by
 * construction, so there is nothing left for the solver to enforce beyond G1.
 *
 * ## Which length goes where
 *
 * Two lengths per edge, and they do different jobs. `C_e` is the chord `|v2 − v1|`, a
 * constant of the input polyline. `L_e` is the arclength, which is `KSCALE` and therefore a
 * solve *output*.
 *
 * The coefficients this module produces describe the **canonical** profile
 *
 * ```
 *   q(u) = L_e · κ_world(L_e · u) ,     u ∈ [0, 1]
 * ```
 *
 * — the curvature of the unit-arclength curve whose similarity image is the edge, which is
 * exactly what `Clothoid` integrates. It is dimensionless and invariant under a uniform
 * zoom of the drawing, and `∫₀^u q` is the turning directly, with no conversion factor.
 *
 * That forces `L_e`, not `C_e`, into the rescale:
 *
 * ```
 *   (1/n!)·dⁿq/duⁿ|₀ = (L_eⁿ⁺¹/n!)·dⁿκ/dsⁿ|₀ = L_e·(L_e/Rᵥ)ⁿ·block(v)_n
 * ```
 *
 * **This is a deviation from §3 of the plan**, which specifies "the rescale from `Rᵥ` to
 * `C_e`" and lists `M_e` and `H` as constant across the solve. Substituting `C_e` leaves
 * every reconstructed endpoint derivative wrong by `(C_e/L_e)ⁿ⁺¹`, so two edges with
 * different `L/C` ratios do not meet — which is the entire structural-continuity claim. For
 * a circular arc the error is 2–4% at half a radian of turning and 60–84% at half a turn.
 * `M_e` therefore moves with the outer iteration, and `H` is reassembled with it; assembly
 * is `O(V·(p+1)³)`, the same order as the factorization that already runs every step, so
 * the cost model of §5 survives with a worse constant. See §14.
 *
 * `Rᵥ` stays chord-derived, and that part of §3 is right and matters: `Rᵥ` fixes what the
 * *unknowns mean*, so tying it to a solve output would measure convergence in a moving
 * frame. `M_e` is only the reconstruction map, and it drifting is no worse than `KTH_e`
 * drifting.
 *
 * ## Orientation
 *
 * Signed curvature depends on which way you walk, so a block only means something relative
 * to a chosen direction of travel. Blocks are stored in *chain* orientation; an edge whose
 * own `v1 → v2` parameterization runs backwards along the chain sees
 *
 * ```
 *   dⁿκ_e/dtⁿ = (−1)ⁿ⁺¹ dⁿκ/dsⁿ
 * ```
 *
 * — one sign from reversing the field, `n` more from reversing the parameter. Odd orders
 * survive, even orders flip. {@link edgeTransform} folds that into its column signs.
 */
import { type CurvatureProfile, sPowerProfile } from "./profile.js";
import { type QuadratureOptions, defaultSPowerQuadrature } from "./quadrature.js";
import { endpointTaylor, integralWeights, massMatrix, sPowerLength, stiffnessMatrix, taylorToPairs } from "./spower.js";

/** Number of scalars a single vertex owns at order `p`. */
export function blockLength(p: number) {
  return p + 1;
}

/** Number of DOF an edge reads: the two blocks at its ends. Equals `sPowerLength(p)`. */
export function edgeDOFLength(p: number) {
  return 2 * (p + 1);
}

/**
 * The reference length `Rᵥ`, as the geometric mean of the chord lengths incident on `v`.
 *
 * The geometric mean rather than the arithmetic one because `Rᵥ` enters the block as `Rᵥⁿ`:
 * what needs to be well-behaved is the *ratio* `C_e / Rᵥ` raised to a power, and the
 * geometric mean is the average that bounds that ratio symmetrically above and below.
 * Disparate incident lengths still leave the columns unequally scaled — that residue is
 * what Ruiz equilibration in `math/banded.ts` is there to absorb.
 */
export function referenceLength(chords: Iterable<number>) {
  let sum = 0.0;
  let count = 0;

  for (const chord of chords) {
    if (chord > 0.0) {
      sum += Math.log(chord);
      count++;
    }
  }

  return count === 0 ? 1.0 : Math.exp(sum / count);
}

/** How an edge sits against the chain that owns its two vertex blocks. */
export interface EdgeFrame {
  /** Chord length `C_e`: `|v2.co − v1.co|`. Constant across the solve; weights the energy. */
  chord: number;

  /**
   * Arclength `L_e` — the current outer-iteration estimate of `KSCALE`.
   *
   * A solve output, so this moves between Newton steps and anything built from it has to
   * be rebuilt. Seed it with {@link chord}, which is the exact answer for a straight edge
   * and a good one for a gently turning one.
   */
  arclength: number;

  /** `Rᵥ` of the vertex that comes first in chain order. */
  rEarlier: number;

  /** `Rᵥ` of the vertex that comes second in chain order. */
  rLater: number;

  /** True when the edge's own `v1 → v2` direction agrees with chain order. */
  forward: boolean;
}

/**
 * What distinguishes one scalar profile's transform from another's — §10.
 *
 * Two questions, and they are the only two. **How the canonical profile is normalized**: `κ`'s
 * is `q(u) = L_e·κ(L_e u)`, carrying one power of `L_e` so that `∫q du` is the turning
 * outright; width's is `ŵ(u) = w(L_e u)`, carrying none, because a width is a length and the
 * stroker wants it in world units. **Whether the field flips with the direction of travel**:
 * signed curvature does, a width does not.
 *
 * Everything else — the `(L_e/Rᵥ)ⁿ` rescale, the `(−1)ⁿ` parameter reversal, the `g`-side
 * convention — is common, which is why this is two flags rather than two transforms.
 */
export interface ScalarKind {
  /** Powers of `L_e` the canonical profile carries on top of the `(L_e/Rᵥ)ⁿ` rescale. */
  readonly leading: number;

  /** Whether the field itself changes sign when the direction of travel reverses. */
  readonly signed: boolean;
}

export const curvatureKind: ScalarKind = { leading: 1, signed: true };
export const widthKind: ScalarKind = { leading: 0, signed: false };

/**
 * Where entry `order` of one edge-end's block lands in the edge's Taylor data, and scaled by
 * what — one column of `M_e` before {@link taylorToPairs} spreads it over the pairs.
 *
 * `far` picks the `g`-side, the expansion in `1 − u`; see the orientation note above. The
 * weight is the `Rᵥ → L_e` rescale times the parity sign, which factors into two independent
 * halves: `(−1)ⁿ` on the `g`-side because that series is expanded in `1 − u`, and — on a
 * backwards edge — the reversal, `(−1)ⁿ⁺¹` for a signed field and `(−1)ⁿ` for an unsigned one.
 * Both are here rather than inline in {@link scalarTransform} because §9's continuation has to
 * *invert* them, and a transform whose inverse is written out separately is a transform with
 * two conventions.
 */
export function scalarEntry(frame: EdgeFrame, earlier: boolean, order: number, kind: ScalarKind) {
  const { arclength, rEarlier, rLater, forward } = frame;

  const scale = Math.pow(arclength, kind.leading) * Math.pow(arclength / (earlier ? rEarlier : rLater), order);

  const parity = (order & 1) === 0 ? 1.0 : -1.0;
  const far = earlier !== forward;

  let sign = far ? parity : 1.0;

  if (!forward) {
    sign *= kind.signed ? -parity : parity;
  }

  return { far, weight: sign * scale };
}

/** {@link scalarEntry} for curvature — the `κ` transform every phase before §10 uses. */
export function transformEntry(frame: EdgeFrame, earlier: boolean, order: number) {
  return scalarEntry(frame, earlier, order, curvatureKind);
}

/**
 * The `(2p+2) × (2p+2)` matrix `M_e` mapping vertex DOF to s-power coefficients,
 * row-major, with DOF ordered `[block(earlier)₀…ₚ, block(later)₀…ₚ]` in **chain** order.
 *
 * Each column is one DOF pushed through {@link taylorToPairs}, which is linear, so the
 * matrix is exactly "which Taylor slot does this DOF land in, times what". The slot and the
 * sign come from {@link EdgeFrame.forward}; the magnitude is the `Rᵥ → L_e` rescale
 *
 * ```
 *   f_n = L_e · (L_e / Rᵥ)ⁿ · block(v)_n
 * ```
 *
 * which is just `(Rᵥⁿ/n!)·dⁿκ/dsⁿ` re-expressed as `(1/n!)·dⁿq/duⁿ` under `q = L_e·κ` and
 * `s = L_e·u`. Note `g_n` carries an extra `(−1)ⁿ` because the far-end Taylor series is
 * expanded in `1 − u`, matching {@link taylorToPairs}'s convention.
 */
export function scalarTransform(p: number, frame: EdgeFrame, kind: ScalarKind, out?: Float64Array) {
  const n = sPowerLength(p);
  const m = out ?? new Float64Array(n * n);

  m.fill(0.0);

  const f = new Float64Array(p + 1);
  const g = new Float64Array(p + 1);
  const col = new Float64Array(n);

  for (let c = 0; c < n; c++) {
    const earlier = c <= p;
    const order = earlier ? c : c - (p + 1);
    const { far, weight } = scalarEntry(frame, earlier, order, kind);

    f.fill(0.0);
    g.fill(0.0);
    (far ? g : f)[order] = weight;

    taylorToPairs(f, g, p, col);

    for (let i = 0; i < n; i++) {
      m[i * n + c] = col[i];
    }
  }

  return m;
}

export function edgeTransform(p: number, frame: EdgeFrame, out?: Float64Array) {
  return scalarTransform(p, frame, curvatureKind, out);
}

/**
 * The per-edge energy `(K + α²·M) / C_e³` in coefficient space, row-major `(2p+2)²`.
 *
 * §6's energy is
 *
 * ```
 *   E_e = (1/C_e)∫₀¹ κ′(u)² du  +  ε·C_e·∫₀¹ κ(u)² du ,       ε = (α/C_e)²
 * ```
 *
 * and taking `κ` to be the canonical profile `q` collapses both terms onto the same
 * `1/C_e³`, leaving `α` as the only free knob and dimensionless with it. That is §6's
 * edit-locality parameter: a bending energy on an elastic foundation, whose disturbances
 * decay like `exp(−√(α/2)·u)` in units of the edge.
 *
 * Reading `κ` as `q` rather than as world curvature makes this differ from the physical
 * bending energy by `(C_e/L_e)²`, which is deliberate. The energy is a regularizer, so any
 * equivalent quadratic form will do, and a chord-only weight is the one §5 argues for: as
 * an edge turns toward closure `1/C_e` grows and the regularizer *strengthens* on the
 * degenerate edge, where an arclength weight would switch itself off. It also keeps this
 * matrix constant, which the transform no longer is.
 */
export function edgeEnergy(p: number, chord: number, alpha: number, out?: Float64Array) {
  const n = sPowerLength(p);
  const e = out ?? new Float64Array(n * n);

  const k = stiffnessMatrix(p);
  const mass = massMatrix(p);

  const w = 1.0 / (chord * chord * chord);
  const a2 = alpha * alpha;

  for (let i = 0; i < n * n; i++) {
    e[i] = w * (k[i] + a2 * mass[i]);
  }

  return e;
}

/*
  Row-major throughout, with the shape passed alongside rather than carried. Rectangular
  because §10's profiles need not have as many coefficients as DOF — s-power does, sampling
  does not.
*/

/** `out = A B` for `A` of `rows × inner` and `B` of `inner × cols`. */
export function multiply(
  a: ArrayLike<number>,
  b: ArrayLike<number>,
  rows: number,
  inner: number,
  cols: number,
  out?: Float64Array
) {
  const r = out ?? new Float64Array(rows * cols);

  for (let i = 0; i < rows; i++) {
    for (let j = 0; j < cols; j++) {
      let sum = 0.0;

      for (let k = 0; k < inner; k++) {
        sum += a[i * inner + k] * b[k * cols + j];
      }

      r[i * cols + j] = sum;
    }
  }

  return r;
}

/** `out = wᵀ M` for `w` of length `rows` and `M` of `rows × cols`. */
export function rowTimesMatrix(
  w: ArrayLike<number>,
  m: ArrayLike<number>,
  rows: number,
  cols: number,
  out?: Float64Array
) {
  const row = out ?? new Float64Array(cols);

  for (let c = 0; c < cols; c++) {
    let sum = 0.0;

    for (let i = 0; i < rows; i++) {
      sum += w[i] * m[i * cols + c];
    }

    row[c] = sum;
  }

  return row;
}

/** `out = M v` for `M` of `rows × cols` and `v` of length `cols`. */
export function matrixTimesVector(
  m: ArrayLike<number>,
  v: ArrayLike<number>,
  rows: number,
  cols: number,
  out?: Float64Array
) {
  const r = out ?? new Float64Array(rows);

  for (let i = 0; i < rows; i++) {
    let sum = 0.0;

    for (let c = 0; c < cols; c++) {
      sum += m[i * cols + c] * v[c];
    }

    r[i] = sum;
  }

  return r;
}

/**
 * `out = Bᵀ A B`, `cols × cols`, for `A` of `rows × rows` and `B` of `rows × cols`.
 *
 * Symmetry of `A` is not assumed.
 */
export function congruence(a: ArrayLike<number>, b: ArrayLike<number>, rows: number, cols: number, out?: Float64Array) {
  const ab = multiply(a, b, rows, rows, cols);
  const r = out ?? new Float64Array(cols * cols);

  for (let i = 0; i < cols; i++) {
    for (let j = 0; j < cols; j++) {
      let sum = 0.0;

      for (let k = 0; k < rows; k++) {
        sum += b[k * cols + i] * ab[k * cols + j];
      }

      r[i * cols + j] = sum;
    }
  }

  return r;
}

/**
 * The edge's contribution to the energy Hessian, `M_eᵀ · E_e · M_e`, in DOF space.
 *
 * Symmetric positive semi-definite, `(2p+2)²` row-major, with the same DOF ordering as
 * {@link edgeTransform}. `E_e` is constant but `M_e` is not, so this is rebuilt whenever
 * the arclength estimate moves — see the note on `L_e` at the top of this file.
 */
export function edgeHessian(p: number, frame: EdgeFrame, alpha: number, out?: Float64Array) {
  const n = sPowerLength(p);

  return congruence(edgeEnergy(p, frame.chord, alpha), edgeTransform(p, frame), n, n, out);
}

/**
 * The width analogue of {@link edgeEnergy}: `(K + α²·M) / C_e`, row-major `(2p+2)²`.
 *
 * §10 asks for "the same `∫(w′)²` energy for smooth taper", and the natural reading of that
 * for a width is the world-arclength Dirichlet energy `∫(dw/ds)² ds`, which on the canonical
 * profile `ŵ(u) = w(L_e u)` is `(1/L_e)∫₀¹ ŵ′(u)² du`. Chord for arclength for §6's reason —
 * it keeps the matrix constant across the solve, and it strengthens rather than weakens the
 * regularizer on a degenerate edge.
 *
 * The chord power is `1` where {@link edgeEnergy}'s is `3`, and the gap is not a discrepancy:
 * `q = L_e·κ` carries a factor of `L_e` that `ŵ` does not, and §6 measures `κ` per unit `u`
 * where this measures `w` per unit `s`. `α` plays the same edit-locality role in both.
 *
 * Unlike curvature's, this quantity is *dimensional* — a length, against a data term that is
 * a squared length — which is why the weight in front of it in `width.ts` is scaled by a
 * reference length rather than used raw.
 */
export function widthEnergy(p: number, chord: number, alpha: number, out?: Float64Array) {
  const n = sPowerLength(p);
  const e = out ?? new Float64Array(n * n);

  const k = stiffnessMatrix(p);
  const mass = massMatrix(p);

  const w = 1.0 / chord;
  const a2 = alpha * alpha;

  for (let i = 0; i < n * n; i++) {
    e[i] = w * (k[i] + a2 * mass[i]);
  }

  return e;
}

/** {@link widthEnergy} pulled back onto vertex DOF, `M_eᵀ·E_e·M_e` with the width transform. */
export function widthHessian(p: number, frame: EdgeFrame, alpha: number, out?: Float64Array) {
  const n = sPowerLength(p);

  return congruence(widthEnergy(p, frame.chord, alpha), scalarTransform(p, frame, widthKind), n, n, out);
}

/**
 * Coefficient-space row `w` such that `w · a` is the s-power profile's value at `u`.
 *
 * Read straight off {@link evalSPower}'s Horner loop: pair `k` contributes
 * `((1−u)·a_{2k} + u·a_{2k+1})·((1−u)u)^k`. This is what a data term is written against — §10's
 * `Σ‖w(sᵢ) − pᵢ‖²` needs the value at a sample, where `κ`'s solver only ever needed integrals.
 */
export function sPowerValueWeights(p: number, u: number, out?: Float64Array) {
  const n = sPowerLength(p);
  const row = out ?? new Float64Array(n);

  const sym = (1.0 - u) * u;
  let symK = 1.0;

  for (let k = 0; k <= p; k++) {
    row[2 * k] = (1.0 - u) * symK;
    row[2 * k + 1] = u * symK;

    symK *= sym;
  }

  return row;
}

/**
 * Row `wᵀ M_e` such that `row · dof` is the edge's total turning, length `2p+2`.
 *
 * `∫₀¹ q du = ∫₀^{L_e} κ ds` is the turning outright — no conversion factor, which is what
 * the canonical normalization buys over §4's `L_e · ∫₀^u κ du`. The unknown arclength has
 * not gone away, it has moved inside `M_e`: freezing it there is what makes the Phase 2
 * solver linear, and differentiating through it is Phase 3.
 */
export function edgeTurningRow(p: number, frame: EdgeFrame, out?: Float64Array) {
  return edgePullback(p, frame, integralWeights(p), out);
}

/**
 * `wᵀ M_e`: a row written against the edge's *coefficients* pulled back onto its DOF.
 *
 * Anything the solver differentiates lands in coefficient space first — `∂KTH_e/∂a_j` from
 * §5's quadrature Jacobian, or the turning weights above — and has to cross the transform to
 * become a row of `J`. That crossing is the same for every such row, so it is one function.
 */
export function edgePullback(p: number, frame: EdgeFrame, w: ArrayLike<number>, out?: Float64Array) {
  const n = sPowerLength(p);

  return rowTimesMatrix(w, edgeTransform(p, frame), n, n, out);
}

/** {@link edgePullback} across the width transform. */
export function widthPullback(p: number, frame: EdgeFrame, w: ArrayLike<number>, out?: Float64Array) {
  const n = sPowerLength(p);

  return rowTimesMatrix(w, scalarTransform(p, frame, widthKind), n, n, out);
}

/** The edge's coefficient vector `a_e = M_e · dof`, with `dof` in chain order. */
export function edgeCoefficients(p: number, frame: EdgeFrame, dof: ArrayLike<number>, out?: Float64Array) {
  const n = sPowerLength(p);

  return matrixTimesVector(edgeTransform(p, frame), dof, n, n, out);
}

/** {@link edgeCoefficients} across the width transform. */
export function widthCoefficients(p: number, frame: EdgeFrame, dof: ArrayLike<number>, out?: Float64Array) {
  const n = sPowerLength(p);

  return matrixTimesVector(scalarTransform(p, frame, widthKind), dof, n, n, out);
}

/**
 * The order-`p+1` block entry at one end of an edge that reproduces its order-`p` curvature
 * there — §9's warm start, one number.
 *
 * The order-`p` blocks carry derivatives up to `p`, but the profile they reconstruct is a
 * polynomial of degree `2p + 1`, so `dᵖ⁺¹κ/dsᵖ⁺¹` at the edge's ends is determined by them.
 * That is the value §9 says a *preserving* continuation would need — and, in the same
 * breath, says the two edges at a joint generally disagree about, since at the new ceiling
 * they have to share it. So this is one edge's opinion; whoever seeds a shared entry has to
 * reconcile them, and the reconciliation is where the curve moves.
 *
 * Read off the s-power Hermite polynomial `M_e · dof` whatever profile is in play. For the
 * sampled control that is not the profile being *rendered* — a polyline has no `p+1`th
 * derivative — but it is the polynomial the blocks mean, which is what a DOF-space warm start
 * is entitled to use.
 */
export function continuationEntry(p: number, frame: EdgeFrame, dof: ArrayLike<number>, earlier: boolean) {
  const order = p + 1;
  const a = edgeCoefficients(p, frame, dof);
  const taylor = endpointTaylor(a, sPowerLength(p), order);

  const { far, weight } = transformEntry(frame, earlier, order);

  return weight === 0.0 ? 0.0 : taylor[far ? 1 : 0] / weight;
}

/**
 * Everything a solver needs to know about one *scalar profile on a chain* — §10.
 *
 * The solver's job is the same whatever the profile is: assemble an energy Hessian over
 * vertex-owned DOF, constrain the turning across joints, factor, step. What differs is only
 * how DOF become coefficients and how many of each there are. Naming that boundary is what
 * lets width reuse the machinery (§10), and what makes §12's control experiment — the old
 * piecewise-linear profile through the new solver — a second implementation of *this*
 * interface rather than a second solver.
 *
 * Matrix shapes stay behind the interface deliberately. {@link edgeTransform} happens to be
 * square because an s-power edge has as many coefficients as it has DOF; nothing else here
 * depends on that.
 */
export interface ScalarProfileDOF {
  readonly name: string;

  /** The profile that integrates the coefficients {@link coefficients} produces. */
  readonly profile: CurvatureProfile;

  /**
   * How that profile should be integrated.
   *
   * Not a free choice: a profile without `d2Curvature` cannot drive the fourth-order scheme
   * at all, so the achievable order is a property of the basis rather than of the caller.
   * That is half of what §12's control experiment measures.
   */
  readonly quadrature: QuadratureOptions;

  /** Scalars owned by one vertex at order `p`. */
  blockLength(p: number): number;

  /** Coefficients one edge carries, i.e. the `klen` its profile expects. */
  coefficientLength(p: number): number;

  /** Energy Hessian in DOF space, `(2·blockLength)²` row-major, symmetric PSD. */
  hessian(p: number, frame: EdgeFrame, alpha: number, out?: Float64Array): Float64Array;

  /**
   * A coefficient-space row pulled back to DOF space, length `2·blockLength`.
   *
   * `w` has {@link coefficientLength} entries. This is how a derivative taken against the
   * profile's coefficients — §5's `∂KTH_e/∂a_j` — becomes a row of the constraint Jacobian.
   */
  pullback(p: number, frame: EdgeFrame, w: ArrayLike<number>, out?: Float64Array): Float64Array;

  /** Coefficients for one edge, from its two blocks concatenated in chain order. */
  coefficients(p: number, frame: EdgeFrame, dof: ArrayLike<number>, out?: Float64Array): Float64Array;

  /** Coefficient-space row whose dot with the coefficients is the profile's value at `u`. */
  valueWeights(p: number, u: number, out?: Float64Array): ArrayLike<number>;
}

/**
 * A scalar profile that is also a *curvature* — §5's geometric coupling, which §10 wants
 * isolated as `κ`'s peculiarity rather than carried by every profile.
 *
 * Turning is the whole of that peculiarity at this boundary. A width has an integral too, and
 * nothing in the solve wants it: no constraint says the widths along a chain have to add up to
 * anything.
 */
export interface ProfileDOF extends ScalarProfileDOF {
  /** Row whose dot with the edge's DOF is its total turning, length `2·blockLength`. */
  turningRow(p: number, frame: EdgeFrame, out?: Float64Array): Float64Array;

  /**
   * The same total turning as a *coefficient*-space row, length {@link coefficientLength}.
   *
   * `turningRow` is this pulled back. The solver needs the unpulled form so it can add it to
   * `∂KTH/∂a` before crossing the transform, rather than crossing twice.
   */
  turningWeights(p: number): ArrayLike<number>;
}

/** The s-power profile as a {@link ProfileDOF}: Hermite blocks of `p + 1` derivatives. */
export const sPowerDOF: ProfileDOF = {
  name      : "s-power",
  profile   : sPowerProfile,
  quadrature: defaultSPowerQuadrature,

  blockLength,
  coefficientLength: sPowerLength,

  hessian       : edgeHessian,
  turningRow    : edgeTurningRow,
  turningWeights: integralWeights,
  pullback      : edgePullback,
  coefficients  : edgeCoefficients,
  valueWeights  : sPowerValueWeights,
};

/**
 * Width on the same blocks — {@link sPowerDOF} with the unsigned, `L_e`-free transform of
 * {@link widthKind} and the arclength Dirichlet energy of {@link widthEnergy}.
 *
 * Deliberately a {@link ScalarProfileDOF} and not a {@link ProfileDOF}: it is the same basis
 * and the same blocks, and the thing it does not have is a turning.
 */
export const widthDOF: ScalarProfileDOF = {
  name      : "width",
  profile   : sPowerProfile,
  quadrature: defaultSPowerQuadrature,

  blockLength,
  coefficientLength: sPowerLength,

  hessian     : widthHessian,
  pullback    : widthPullback,
  coefficients: widthCoefficients,
  valueWeights: sPowerValueWeights,
};
