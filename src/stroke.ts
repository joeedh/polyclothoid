import { Vector2, clamp } from "./math/index.js";
import { Mesh, type CurveConstructor, type CurveSolverConstructor } from "./mesh/index.js";

/** Called once per dab. `dx`/`dy` is the unit tangent; `t` is 0..1 along the segment. */
export type StrokeCallback = (x: number, y: number, dx: number, dy: number, t: number) => void;

/** How many past input positions the fitter keeps. */
const HISTORY = 6;

/**
 * Turns a stream of pointer positions into evenly spaced dabs.
 *
 * Raw pointer input is spaced by however fast the user moved, which is the wrong spacing
 * for a brush. So each accepted input triggers a fit: the last few positions are dropped
 * into a throwaway {@link Mesh}, solved into a smooth curve, and the *second* segment of
 * that curve — the one with real neighbours on both sides — is walked at a fixed
 * arclength step. Trailing distance carries over in {@link lastS} so dabs stay evenly
 * spaced across segment boundaries.
 *
 * The lag of one segment is deliberate: fitting the newest segment would mean fitting
 * with no lookahead, and the curve would wobble as each new sample arrived.
 */
/** Overrides the curve type the fitter solves with; defaults to {@link Clothoid}. */
export interface StrokerOptions {
  CurveCls?: CurveConstructor;
  SolverCls?: CurveSolverConstructor;

  /**
   * Input-polyline turns sharper than this are authored as corners — level-0 pairings.
   *
   * This used to be `ClothoidSolverOptions.cornerThreshold`, and it is here because it is
   * *inference*: a guess, from pre-solve geometry, about what the person drawing meant. The
   * solvers take authored levels and never invent one, so the guess has to be made by whoever
   * owns the input — see `docs/plans/spower-solver.md` §3 and §11. Its being a stroker
   * concern is also why it has always been about pointer devices: fast input genuinely turns
   * a corner between two samples, and rounding those off looks like the fitter lagging.
   *
   * Measured as the interior angle between the two chords, so smaller is sharper and `π` is
   * straight through. Set it to `0` to author no corners at all.
   */
  cornerThreshold?: number;
}

/** ~72°, the value `ClothoidSolver` used before the test moved out here. */
export const DEFAULT_CORNER_THRESHOLD = Math.PI * 0.4;

export class Stroker {
  callback: StrokeCallback;
  options: StrokerOptions;

  /** Multiplier on the spacing gate, in units of `radius * 2`. */
  lag = 1.0;

  mpos = new Vector2();
  history: Vector2[] = [];

  /** Accepted inputs so far; the fit needs at least three. */
  count = 0;

  /** Leftover arclength from the previous segment. */
  lastS = 0.0;

  first = true;
  firstSegment = true;

  /** Last solved mesh, kept so a caller can draw the fitted curve. */
  lastMesh?: Mesh;

  constructor(callback: StrokeCallback, options: StrokerOptions = {}) {
    this.callback = callback;
    this.options = options;

    for (let i = 0; i < HISTORY; i++) {
      this.history.push(new Vector2());
    }
  }

  /**
   * Feed one pointer position.
   *
   * `radius` and `spacing` are the brush radius and its spacing as a fraction of
   * diameter, matching the usual painting-app convention.
   */
  onInput(x: number, y: number, radius: number, spacing: number) {
    const mpos = this.mpos.loadXY(x, y);
    const history = this.history;

    if (this.first) {
      this.first = false;
      this.lastS = 0;

      history[0].load(mpos);

      return;
    }

    /* Gate on distance travelled, not time, so slow strokes do not oversample. */
    if (mpos.vectorDistance(history[0]) / (radius * 2.0) <= spacing * this.lag) {
      return;
    }

    const ready = this.count >= 4 || (this.firstSegment && this.count >= 3);

    if (ready) {
      this.emit(mpos, radius, spacing);
    }

    for (let i = HISTORY - 1; i > 0; i--) {
      history[i].load(history[i - 1]);
    }

    history[0].load(mpos);
    this.count++;
  }

  emit(mpos: Vector2, radius: number, spacing: number) {
    const ds = spacing * (2.0 * radius);

    if (ds === 0.0) {
      console.error("Stroker: spacing was zero");
      return;
    }

    const history = this.history;
    const mesh = new Mesh();

    const { CurveCls, SolverCls } = this.options;
    if (CurveCls && SolverCls) {
      mesh.CurveCls = CurveCls;
      mesh.SolverCls = SolverCls;
    }

    let e;

    if (this.firstSegment) {
      this.firstSegment = false;
      this.lastS = ds;

      /* Only three positions are trustworthy yet, and the original skips history[2]. */
      const v1 = mesh.makeVertex(history[3]);
      const v2 = mesh.makeVertex(history[1]);
      const v3 = mesh.makeVertex(mpos);

      e = mesh.makeEdge(v1, v2);
      mesh.makeEdge(v2, v3);
    } else {
      const v1 = mesh.makeVertex(history[3]);
      const v2 = mesh.makeVertex(history[2]);
      const v3 = mesh.makeVertex(history[1]);
      const v4 = mesh.makeVertex(history[0]);
      const v5 = mesh.makeVertex(mpos);

      mesh.makeEdge(v1, v2);
      e = mesh.makeEdge(v2, v3);
      mesh.makeEdge(v3, v4);
      mesh.makeEdge(v4, v5);
    }

    this.markCorners(mesh);

    mesh.solve();
    this.lastMesh = mesh;

    const elen = e.length;

    if (isNaN(elen)) {
      console.error("Stroker: segment length was NaN");
      return;
    }

    let s = this.lastS;

    while (s < elen) {
      const p = e.evaluate(s);
      const dv = e.derivative(s);

      this.callback(p[0], p[1], dv[0], dv[1], s / elen);

      s += ds;
    }

    this.lastS = s - elen;
  }

  /**
   * Author a level-0 pairing wherever the input turned sharply — {@link
   * StrokerOptions.cornerThreshold}.
   *
   * Only at valence 2, which is every vertex a fitted segment has; a junction has no canonical
   * pairing to author and the stroker does not build one. The `0.99999` scaling on the dot
   * product is inherited from the original and keeps `acos` off its infinite-slope endpoints,
   * where a doubled-back chord would otherwise cost several digits.
   */
  markCorners(mesh: Mesh) {
    const threshold = this.options.cornerThreshold ?? DEFAULT_CORNER_THRESHOLD;

    if (!(threshold > 0.0)) {
      return;
    }

    for (const v of mesh.verts) {
      if (v.edges.length !== 2) {
        continue;
      }

      const [e1, e2] = v.edges;

      const t1 = new Vector2(e1.otherVertex(v)).sub(v).normalize();
      const t2 = new Vector2(e2.otherVertex(v)).sub(v).normalize();

      if (Math.acos(clamp(t1.dot(t2) * 0.99999, -1.0, 1.0)) < threshold) {
        v.setPairing(e1, e2, 0);
      }
    }
  }
}
