import { type Vector2 } from "../math/index.js";
import { type Curve, type CurveEdge, type CurveVertex } from "./curve.js";
import { type SolveReport } from "./diagnostics.js";

/**
 * What the per-curve-type solvers need from a mesh.
 *
 * These are structural rather than imports of the concrete `Mesh` classes, which is what
 * keeps `curve/` from importing `mesh/` — the dependency runs one way, mesh -> curve. See
 * `docs/architecture.md`.
 */
/**
 * An authored pairing of two edge-ends at a vertex, carrying a continuity level — §3.
 *
 * The level is a *ceiling*: the solver may deliver less (§8) and never more. `0` is a corner,
 * `1` is G1, and `k` shares block entries `0 … k−2` on top of the G1 row, up to `p + 2`.
 * Levels are authored by client mesh code; nothing in `curve/` infers one.
 */
export interface Pairing {
  a: SolvableEdge;
  b: SolvableEdge;
  level: number;
}

export interface SolvableVertex extends CurveVertex {
  edges: SolvableEdge[];
  otherEdge(e: SolvableEdge): SolvableEdge | undefined;

  /**
   * Authored pairings at this vertex, if the mesh carries any.
   *
   * An absent or empty list is not "everything is a corner" — it means *unspecified*, and
   * `pairing.ts` fills it in from valence. Only valence ≤ 2 has a default worth having.
   */
  pairings?: Pairing[];
}

export interface SolvableEdge<C extends Curve = Curve> extends CurveEdge {
  v1: SolvableVertex;
  v2: SolvableVertex;
  curve: C;

  otherVertex(v: SolvableVertex): SolvableVertex;

  get length(): number;
  evaluate(s: number): Vector2;
  derivative(s: number): Vector2;
  update(): void;
}

export interface SolvableMesh<C extends Curve = Curve> {
  verts: Iterable<SolvableVertex>;
  edges: Iterable<SolvableEdge<C>>;
}

/** Fits every curve in a mesh so they agree at shared vertices. */
export interface CurveSolver {
  /**
   * Solve, and hand back what happened — `docs/plans/spower-solver.md` §8.
   *
   * The report is the only channel through which a caller can tell an authored break from a
   * solver bail-out, so it is a return value rather than a field: a solver that quietly
   * degrades and is never asked about it produces a model that looks right and is not.
   */
  solve(): SolveReport;
}
