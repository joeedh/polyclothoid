/**
 * A clothoid segment whose curvature profile is a polynomial in the s-power basis rather
 * than a piecewise-linear sampling — `docs/plans/spower-solver.md` §1 and §3.
 *
 * Geometrically it is the same construction as {@link Clothoid}: integrate `κ`, integrate
 * the tangent, place the result on `v1 → v2` with a similarity transform. What changes is
 * that `κ` is now smooth and closed-form, so the quadrature has no knots capping its order
 * and `d2Curvature` exists for the fourth-order terms (§7).
 *
 * Two differences from {@link Clothoid} that matter when reading them side by side:
 *
 * 1. **The canonical interval is `[0, 1]`, not `[−0.5, 0.5]`.** Integrating from the start
 *    rather than the middle puts the canonical curve's origin at `v1`, so `evaluate(0)`
 *    returns `v1` exactly and there is no `KOFFX`/`KOFFY` to cache. The centred interval
 *    buys nothing here — it existed so the shape was built around the origin before being
 *    placed, which the offset then undid.
 * 2. **The transform lives in named fields, not in tail slots of `ks`.** `Clothoid`'s
 *    `ks` does double duty so that `_ks` can be a view over just the free parameters, which
 *    is what its solver is handed. This curve's solver does not touch `ks` directly at all —
 *    the degrees of freedom live on vertices (§3) and `ks` is derived from them — so there
 *    is nothing to take a view over.
 *
 * The basis is nonetheless a *field* rather than a hard-coded import, because §12's control
 * experiment runs the piecewise-linear profile through this same curve and solver. Everything
 * below reads `ks` through {@link SPowerClothoid.profile}; the class name records what it is
 * for, not what it is limited to.
 */
import { CacheRing, Vector2, clamp } from "../math/index.js";
import { Curve, type CurveEdge } from "./curve.js";
import { type SolvableVertex } from "./mesh_types.js";
import { type CurvatureProfile, profileByName, sPowerProfile } from "./profile.js";
import { MIN_CANONICAL_CHORD, defaultQuadratureFor, defaultSPowerQuadrature, integrateProfile } from "./quadrature.js";
import { sPowerLength } from "./spower.js";
import * as nstructjs from "nstructjs";

/** Default s-power order `p`. Degree `2p + 1`, so cubic curvature and `p + 1 = 2` DOF per vertex. */
export const SPOWER_ORDER = 1;

const evalRets = new CacheRing(() => new Vector2(), 128);
const dvRets = new CacheRing(() => new Vector2(), 128);
const dv2Rets = new CacheRing(() => new Vector2(), 128);

const scratch = new Vector2();

export class SPowerClothoid extends Curve {
  static override STRUCT = nstructjs.inlineRegister(
    this,
    `
SPowerClothoid {
  ks          : array(float);
  ws          : array(float);
  order       : int;
  profileName : string;
}
`
  );

  /** The s-power order `p`. */
  order = SPOWER_ORDER;

  /** The `2p + 2` s-power coefficients of the canonical profile `q(u) = L·κ(L·u)`. */
  ks = new Float64Array(sPowerLength(SPOWER_ORDER));

  /**
   * Width coefficients in the same basis, for the canonical profile `ŵ(u) = w(L·u)` — §10.
   *
   * Empty until a `WidthSolver` has run, which is a state the stroker has to handle anyway:
   * §10's ordering constraint is that `κ` is solved first, so between the two solves every
   * curve in the mesh is geometrically complete and has no width.
   *
   * No leading `L` here where {@link ks} has one. A width is a length and the stroker wants
   * it in world units, so `ŵ` is `w` reparameterized and not rescaled — see `blocks.ts`'s
   * `widthKind`.
   */
  ws = new Float64Array(0);

  /**
   * The basis `ks` is written in.
   *
   * Almost always {@link sPowerProfile}. It is a field rather than a hard-coded import for
   * §12's control experiment, which runs the piecewise-linear profile through this same
   * curve and solver so that "the basis was the win" can be told apart from "the banded KKT
   * and the energy were the win". Set it through {@link setProfile}, never on its own —
   * {@link quadrature} is not independent of it.
   */
  profile: CurvatureProfile = sPowerProfile;

  /** Arclength `L`, the uniform scale of the similarity transform. `Clothoid`'s `KSCALE`. */
  scale = 1.0;

  /** Placement angle. `Clothoid`'s `KTH`. */
  th = 0.0;

  /**
   * `|C(1)|` of the canonical curve, before the {@link MIN_CANONICAL_CHORD} clamp.
   *
   * Its reciprocal is `L_e / C_e`, so it falls towards zero exactly as the segment closes on
   * itself and the similarity transform stops being well-conditioned. §8's chord-degeneracy
   * measurement, kept unclamped because the clamped value cannot say how far past it went.
   *
   * Bounded above by `1`, the canonical curve having unit arclength — but only for the exact
   * integral. The quadrature's truncation terms carry `k²` and `k³`, so a coefficient vector
   * that has run away produces a value far above `1` rather than a saturated one, which is
   * what makes the upper bound worth checking as well as the lower.
   */
  canonical = 1.0;

  /** How {@link profile} is integrated. Constrained by it — see {@link setProfile}. */
  quadrature = defaultSPowerQuadrature;

  recalc = 1;

  v1!: SolvableVertex;
  v2!: SolvableVertex;

  constructor(v1?: SolvableVertex, v2?: SolvableVertex) {
    super();

    if (v1 !== undefined && v2 !== undefined) {
      this.v1 = v1;
      this.v2 = v2;
    }
  }

  /** The profile's name, which is what {@link STRUCT} stores in place of the profile itself. */
  get profileName() {
    return this.profile.name;
  }

  set profileName(name: string) {
    this.setProfile(profileByName(name));
  }

  /**
   * Adopt a profile and the quadrature that goes with it, together.
   *
   * They travel as a pair because the achievable order is a property of the basis: a profile
   * with no `d2Curvature` cannot drive the fourth-order terms, and `integrateProfile` throws
   * rather than quietly dropping them. Setting {@link profile} alone would leave the two
   * disagreeing until whoever set it remembered the second half.
   */
  setProfile(profile: CurvatureProfile, quadrature = defaultQuadratureFor(profile)) {
    this.profile = profile;
    this.quadrature = quadrature;
    this.recalc = 1;

    return this;
  }

  /** Coefficient count, i.e. the `klen` {@link profile} expects. */
  get klen() {
    return this.ks.length;
  }

  /**
   * Resize `ks` for a new order, keeping whatever coefficients still fit.
   *
   * `count` is the coefficient count, which is `2p + 2` only for the s-power basis: a
   * sampled profile carries as many samples as it likes over the same `p + 1` derivatives
   * per vertex. It is the owning {@link ProfileDOF} that knows which.
   */
  setOrder(p: number, count = sPowerLength(p)) {
    const ks = new Float64Array(count);

    ks.set(this.ks.subarray(0, Math.min(this.ks.length, ks.length)));

    this.ks = ks;
    this.order = p;
    this.recalc = 1;

    return this;
  }

  override get length() {
    if (this.recalc) {
      this._update();
    }

    return this.scale;
  }

  override init(e: CurveEdge) {
    this.v1 = e.v1 as SolvableVertex;
    this.v2 = e.v2 as SolvableVertex;
  }

  override update(e?: CurveEdge) {
    if (e) {
      this.v1 = e.v1 as SolvableVertex;
      this.v2 = e.v2 as SolvableVertex;
    }

    this.recalc = 1;

    return this;
  }

  /**
   * Recompute the similarity transform carrying the canonical shape onto `v1 → v2`.
   *
   * The canonical curve starts at the origin, so its chord is just its endpoint and one
   * integration suffices. `scale` is the ratio of the real chord to the canonical one, which
   * is the arclength because the canonical curve has unit arclength by construction.
   */
  _update() {
    this.recalc = 0;

    const end = integrateProfile(this.profile, this.ks, this.klen, 0.0, 1.0, this.quadrature, scratch);

    const chord = this.v1.vectorDistance(this.v2);

    this.canonical = end.vectorLength();

    this.scale = chord / Math.max(this.canonical, MIN_CANONICAL_CHORD);

    const th1 = Math.atan2(this.v2[1] - this.v1[1], this.v2[0] - this.v1[0]);
    const th2 = Math.atan2(end[1], end[0]);

    this.th = th1 - th2;
  }

  /** Normalized parameter from arclength, clamped to the segment. */
  _param(s: number) {
    return clamp(s / this.scale, 0.0, 1.0);
  }

  override evaluate(s: number) {
    if (this.recalc) {
      this._update();
    }

    const u = this._param(s);
    const p = evalRets.next().zero();

    integrateProfile(this.profile, this.ks, this.klen, 0.0, u, this.quadrature, p);

    return p.rot2d(this.th).mulScalar(this.scale).add(this.v1);
  }

  /** Analytic tangent: `s` is arclength, so this is exactly `(cos θ, sin θ)`. */
  override derivative(s: number) {
    if (this.recalc) {
      this._update();
    }

    const th = this.profile.integral(this.ks, this.klen, this._param(s)) + this.th;

    const ret = dvRets.next();
    ret[0] = Math.cos(th);
    ret[1] = Math.sin(th);

    return ret;
  }

  override derivative2(s: number) {
    if (this.recalc) {
      this._update();
    }

    const u = this._param(s);

    const th = this.profile.integral(this.ks, this.klen, u) + this.th;
    const k = this.profile.curvature(this.ks, this.klen, u) / this.scale;

    const ret = dv2Rets.next();
    ret[0] = -Math.sin(th) * k;
    ret[1] = Math.cos(th) * k;

    return ret;
  }

  /** World curvature. The coefficients carry the canonical profile, hence the `1 / L`. */
  override curvature(s: number) {
    if (this.recalc) {
      this._update();
    }

    return this.profile.curvature(this.ks, this.klen, this._param(s)) / this.scale;
  }

  /** Total turning across the segment, `∫₀¹ q du`. Exact, and linear in `ks`. */
  get turning() {
    return this.profile.integral(this.ks, this.klen, 1.0);
  }

  /**
   * Width at arclength `s`, or `NaN` where no width has been solved.
   *
   * `NaN` rather than a default so that a stroker reading a width nobody set produces visible
   * nonsense instead of a plausible constant — §10 makes the width a solve output, and a curve
   * that has not had one is not a curve with width 1.
   */
  width(s: number) {
    if (this.ws.length === 0) {
      return NaN;
    }

    if (this.recalc) {
      this._update();
    }

    return this.profile.curvature(this.ws, this.ws.length, this._param(s));
  }

  afterSTRUCT(v1: SolvableVertex, v2: SolvableVertex) {
    this.v1 = v1;
    this.v2 = v2;
    this.recalc = 1;
  }

  /*
    `array(float)` reads back as a plain Array, and the coefficient count is whatever was
    written — `2p + 2` for the s-power basis but a free sample count otherwise.
  */
  loadSTRUCT(reader: (obj: this) => void) {
    reader(this);

    this.ks = Float64Array.from(this.ks);
    this.ws = Float64Array.from(this.ws ?? []);
    this.recalc = 1;
  }
}
