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
 */
import { CacheRing, Vector2, clamp } from "../math/index.js";
import { Curve, type CurveEdge } from "./curve.js";
import { type SolvableVertex } from "./mesh_types.js";
import { sPowerProfile } from "./profile.js";
import { MIN_CANONICAL_CHORD, type QuadratureOptions, integrateProfile } from "./quadrature.js";
import { evalSPower, sPowerIntegral, sPowerLength } from "./spower.js";
import * as nstructjs from "nstructjs";

/** Default s-power order `p`. Degree `2p + 1`, so cubic curvature and `p + 1 = 2` DOF per vertex. */
export const SPOWER_ORDER = 1;

export const defaultSPowerQuadrature: QuadratureOptions = {
  steps      : 19,
  fourthOrder: true,
};

const evalRets = new CacheRing(() => new Vector2(), 128);
const dvRets = new CacheRing(() => new Vector2(), 128);
const dv2Rets = new CacheRing(() => new Vector2(), 128);

const scratch = new Vector2();

export class SPowerClothoid extends Curve {
  static override STRUCT = nstructjs.inlineRegister(
    this,
    `
SPowerClothoid {
  ks    : array(float);
  order : int;
}
`
  );

  /** The s-power order `p`. */
  order = SPOWER_ORDER;

  /** The `2p + 2` s-power coefficients of the canonical profile `q(u) = L·κ(L·u)`. */
  ks = new Float64Array(sPowerLength(SPOWER_ORDER));

  /** Arclength `L`, the uniform scale of the similarity transform. `Clothoid`'s `KSCALE`. */
  scale = 1.0;

  /** Placement angle. `Clothoid`'s `KTH`. */
  th = 0.0;

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

  /** Coefficient count, i.e. the `klen` {@link sPowerProfile} expects. */
  get klen() {
    return sPowerLength(this.order);
  }

  /** Resize `ks` for a new order, keeping whatever coefficients still fit. */
  setOrder(p: number) {
    const ks = new Float64Array(sPowerLength(p));

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

    const end = integrateProfile(sPowerProfile, this.ks, this.klen, 0.0, 1.0, this.quadrature, scratch);

    const chord = this.v1.vectorDistance(this.v2);
    const canonical = Math.max(end.vectorLength(), MIN_CANONICAL_CHORD);

    this.scale = chord / canonical;

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

    integrateProfile(sPowerProfile, this.ks, this.klen, 0.0, u, this.quadrature, p);

    return p.rot2d(this.th).mulScalar(this.scale).add(this.v1);
  }

  /** Analytic tangent: `s` is arclength, so this is exactly `(cos θ, sin θ)`. */
  override derivative(s: number) {
    if (this.recalc) {
      this._update();
    }

    const th = sPowerIntegral(this.ks, this.klen, this._param(s)) + this.th;

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

    const th = sPowerIntegral(this.ks, this.klen, u) + this.th;
    const k = evalSPower(this.ks, this.klen, u) / this.scale;

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

    return evalSPower(this.ks, this.klen, this._param(s)) / this.scale;
  }

  /** Total turning across the segment, `∫₀¹ q du`. Exact, and linear in `ks`. */
  get turning() {
    return sPowerIntegral(this.ks, this.klen, 1.0);
  }

  afterSTRUCT(v1: SolvableVertex, v2: SolvableVertex) {
    this.v1 = v1;
    this.v2 = v2;
    this.recalc = 1;
  }

  loadSTRUCT(reader: (obj: this) => void) {
    reader(this);

    const ks = new Float64Array(sPowerLength(this.order));
    ks.set(this.ks.subarray(0, Math.min(this.ks.length, ks.length)));

    this.ks = ks;
    this.recalc = 1;
  }
}
