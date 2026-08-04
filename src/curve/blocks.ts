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
import { integralWeights, massMatrix, sPowerLength, stiffnessMatrix, taylorToPairs } from "./spower.js";

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
export function edgeTransform(p: number, frame: EdgeFrame, out?: Float64Array) {
  const n = sPowerLength(p);
  const m = out ?? new Float64Array(n * n);
  const { arclength, rEarlier, rLater, forward } = frame;

  m.fill(0.0);

  const f = new Float64Array(p + 1);
  const g = new Float64Array(p + 1);
  const col = new Float64Array(n);

  for (let c = 0; c < n; c++) {
    const earlier = c <= p;
    const order = earlier ? c : c - (p + 1);
    const scale = arclength * Math.pow(arclength / (earlier ? rEarlier : rLater), order);

    const even = (order & 1) === 0;

    // `far` selects the g-side (the 1−u expansion); see the orientation note above.
    const far = earlier !== forward;
    let sign: number;

    if (forward) {
      sign = earlier ? 1.0 : even ? 1.0 : -1.0;
    } else {
      sign = earlier ? -1.0 : even ? -1.0 : 1.0;
    }

    f.fill(0.0);
    g.fill(0.0);
    (far ? g : f)[order] = sign * scale;

    taylorToPairs(f, g, p, col);

    for (let i = 0; i < n; i++) {
      m[i * n + c] = col[i];
    }
  }

  return m;
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

/** `out = Bᵀ A B` for square row-major `n × n` matrices. Symmetry of `A` is not assumed. */
export function congruence(a: ArrayLike<number>, b: ArrayLike<number>, n: number, out?: Float64Array) {
  const ab = new Float64Array(n * n);

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let sum = 0.0;

      for (let k = 0; k < n; k++) {
        sum += a[i * n + k] * b[k * n + j];
      }

      ab[i * n + j] = sum;
    }
  }

  const r = out ?? new Float64Array(n * n);

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let sum = 0.0;

      for (let k = 0; k < n; k++) {
        sum += b[k * n + i] * ab[k * n + j];
      }

      r[i * n + j] = sum;
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

  return congruence(edgeEnergy(p, frame.chord, alpha), edgeTransform(p, frame), n, out);
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
  const m = edgeTransform(p, frame);
  const row = out ?? new Float64Array(n);

  for (let c = 0; c < n; c++) {
    let sum = 0.0;

    for (let i = 0; i < n; i++) {
      sum += w[i] * m[i * n + c];
    }

    row[c] = sum;
  }

  return row;
}

/** The edge's coefficient vector `a_e = M_e · dof`, with `dof` in chain order. */
export function edgeCoefficients(p: number, frame: EdgeFrame, dof: ArrayLike<number>, out?: Float64Array) {
  const n = sPowerLength(p);
  const m = edgeTransform(p, frame);
  const a = out ?? new Float64Array(n);

  for (let i = 0; i < n; i++) {
    let sum = 0.0;

    for (let c = 0; c < n; c++) {
      sum += m[i * n + c] * dof[c];
    }

    a[i] = sum;
  }

  return a;
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
export interface ProfileDOF {
  readonly name: string;

  /** The profile that integrates the coefficients {@link coefficients} produces. */
  readonly profile: CurvatureProfile;

  /** Scalars owned by one vertex at order `p`. */
  blockLength(p: number): number;

  /** Coefficients one edge carries, i.e. the `klen` its profile expects. */
  coefficientLength(p: number): number;

  /** Energy Hessian in DOF space, `(2·blockLength)²` row-major, symmetric PSD. */
  hessian(p: number, frame: EdgeFrame, alpha: number, out?: Float64Array): Float64Array;

  /** Row whose dot with the edge's DOF is its total turning, length `2·blockLength`. */
  turningRow(p: number, frame: EdgeFrame, out?: Float64Array): Float64Array;

  /**
   * The same total turning as a *coefficient*-space row, length {@link coefficientLength}.
   *
   * `turningRow` is this pulled back. The solver needs the unpulled form so it can add it to
   * `∂KTH/∂a` before crossing the transform, rather than crossing twice.
   */
  turningWeights(p: number): ArrayLike<number>;

  /**
   * A coefficient-space row pulled back to DOF space, length `2·blockLength`.
   *
   * `w` has {@link coefficientLength} entries. This is how a derivative taken against the
   * profile's coefficients — §5's `∂KTH_e/∂a_j` — becomes a row of the constraint Jacobian.
   */
  pullback(p: number, frame: EdgeFrame, w: ArrayLike<number>, out?: Float64Array): Float64Array;

  /** Coefficients for one edge, from its two blocks concatenated in chain order. */
  coefficients(p: number, frame: EdgeFrame, dof: ArrayLike<number>, out?: Float64Array): Float64Array;
}

/** The s-power profile as a {@link ProfileDOF}: Hermite blocks of `p + 1` derivatives. */
export const sPowerDOF: ProfileDOF = {
  name   : "s-power",
  profile: sPowerProfile,

  blockLength,
  coefficientLength: sPowerLength,

  hessian       : edgeHessian,
  turningRow    : edgeTurningRow,
  turningWeights: integralWeights,
  pullback      : edgePullback,
  coefficients  : edgeCoefficients,
};
