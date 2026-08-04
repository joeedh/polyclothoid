/**
 * Mesh traversal for solvers, over the structural interfaces in `mesh_types.ts` rather
 * than the concrete `Mesh` classes — `curve/` must not import `mesh/`. See
 * `docs/architecture.md`.
 */
import { type Curve } from "./curve.js";
import { type SolvableEdge, type SolvableMesh, type SolvableVertex } from "./mesh_types.js";

/**
 * Step from `v` across an edge that is not `e`.
 *
 * Returns `v` itself when there is no such edge, i.e. at a valence-1 chain endpoint, so
 * callers get a degenerate zero-length step rather than an exception.
 *
 * This is a chain idiom and it does not generalize: at valence 3 or above "the other
 * edge" is not well defined, and this picks whichever comes first in `v.edges`, which is
 * insertion order. Solvers that care about junctions must consult authored edge pairings
 * instead — see `docs/plans/spower-solver.md` §3.
 */
export function walk(v: SolvableVertex, e: SolvableEdge): SolvableVertex {
  for (const e2 of v.edges) {
    if (e2 !== e) {
      return e2.otherVertex(v);
    }
  }

  return v;
}

/**
 * A maximal run of valence-2 vertices, with the nodes at either end.
 *
 * `edges[i]` joins `verts[i]` to `verts[i + 1]`, so an open chain has one more vertex than
 * it has edges. A closed one wraps: `verts.length === edges.length` and the last edge joins
 * `verts[last]` back to `verts[0]`.
 *
 * The orientation is the walk's, not the edges' — an edge may run against it, which is what
 * `EdgeFrame.forward` in `blocks.ts` records.
 */
export interface Chain<C extends Curve = Curve> {
  verts: SolvableVertex[];
  edges: SolvableEdge<C>[];

  /** True when the walk returned to where it started. */
  closed: boolean;
}

/** Walk one chain out of `start` along `first`, consuming edges from `seen`. */
function follow<C extends Curve>(start: SolvableVertex, first: SolvableEdge<C>, seen: Set<SolvableEdge<C>>) {
  const verts = [start];
  const edges: SolvableEdge<C>[] = [];

  let v = start;
  let e: SolvableEdge<C> | undefined = first;

  while (e !== undefined && !seen.has(e)) {
    seen.add(e);
    edges.push(e);

    v = e.otherVertex(v);
    verts.push(v);

    if (v === start || v.edges.length !== 2) {
      break;
    }

    e = (v.edges[0] === e ? v.edges[1] : v.edges[0]) as SolvableEdge<C>;
  }

  const closed = edges.length > 0 && verts[verts.length - 1] === start;

  if (closed) {
    verts.pop();
  }

  return { verts, edges, closed } satisfies Chain<C>;
}

/**
 * Repeat a closed chain's first vertex at its end, so every chain is a path — §5's cut point.
 *
 * A cyclic band is not a band, and §5 offers two ways out: bordering, or promoting a vertex
 * to a node and letting the Schur complement carry the wrap-around. This is the second. The
 * cut vertex then appears twice, once at each end, with a block at each, and the two are
 * reconciled through the junction machinery exactly as two chains meeting at a node are.
 *
 * `closed` survives the cut, because "these two ends are the same vertex" stays true and is
 * how a caller knows the extra vertex is not a real one. Open chains are returned unchanged
 * and uncopied.
 */
export function cutOpen<C extends Curve>(chain: Chain<C>): Chain<C> {
  if (!chain.closed) {
    return chain;
  }

  return { verts: [...chain.verts, chain.verts[0]], edges: chain.edges, closed: true };
}

/**
 * Decompose a mesh into nodes and chains — `docs/plans/spower-solver.md` §5.
 *
 * Nodes are the vertices of valence other than 2; chains are the maximal paths between
 * them. Every edge lands in exactly one chain, and valence-0 vertices land in none. Chains
 * are emitted node-first, so an all-valence-2 component — which has no node to start from —
 * comes out last, as a closed chain cut at an arbitrary vertex.
 *
 * A cycle hanging off a node is reported closed too, since it does return to its start. A
 * caller that wants a linear system out of it has to cut it somewhere regardless, and
 * `verts[0]` is as good a place as any.
 *
 * This is the *topological* decomposition of §5's table, the static one. Splitting further
 * at authored corners to get solve blocks is the solver's job, not this function's.
 */
export function chains<C extends Curve>(mesh: SolvableMesh<C>): Chain<C>[] {
  const out: Chain<C>[] = [];
  const seen = new Set<SolvableEdge<C>>();

  for (const v of mesh.verts) {
    if (v.edges.length === 2) {
      continue;
    }

    for (const e of v.edges as SolvableEdge<C>[]) {
      if (!seen.has(e)) {
        out.push(follow(v, e, seen));
      }
    }
  }

  for (const e of mesh.edges) {
    if (!seen.has(e)) {
      out.push(follow(e.v1, e, seen));
    }
  }

  return out;
}
