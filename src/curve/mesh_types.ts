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
export interface SolvableVertex extends CurveVertex {
  edges: SolvableEdge[];
  otherEdge(e: SolvableEdge): SolvableEdge | undefined;
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
