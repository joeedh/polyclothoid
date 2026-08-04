/**
 * What couples one chain to another — the interface of §5's substructuring, and Phase 8 of
 * `docs/plans/spower-solver.md` §12.
 *
 * `topology.ts` decomposes a mesh into nodes and chains. That decomposition is *topological*
 * and static; this one is the solve-side reading of it, and it answers a narrower question:
 * which unknowns does more than one chain touch? Everything else stays inside a chain, where
 * the band is, and the answer here is what has to be eliminated onto instead.
 *
 * §5 names two things that break naive decoupling, and both are handled here rather than
 * anywhere they could be forgotten:
 *
 * 1. **A level-1 pairing shares no block but still emits a G1 row.** A junction authored
 *    G1-only couples two chains through a single row touching both, so {@link GammaSlot} has
 *    a `"row"` case alongside its `"entry"` case and the interface is the union of the two.
 * 2. **`Rᵥ` is shared** by every edge incident on a node. That one needs no machinery,
 *    because §3 pins `Rᵥ` to chord lengths and chord lengths are solver *input* —
 *    `referenceLengths` computes it once for the whole mesh before any chain is built. It is
 *    a coupling in the model and a constant in the solve.
 *
 * The interface is deliberately *not* "every node DOF". A node group read by a single
 * edge-end is a private block like any other and stays in its chain's band; only a group with
 * two or more ends in it, or a row spanning two ends, is promoted. That is what keeps the
 * common case free: §3's default at valence ≥ 3 is fully split, so a plain junction produces
 * no interface at all and the chains fall apart into exactly the independent banded solves
 * Phase 2 already ran.
 */
import { type Curve } from "./curve.js";
import { type SolvableEdge, type SolvableVertex } from "./mesh_types.js";
import { curvatureChannel, vertexPartitions } from "./pairing.js";
import { type Chain } from "./topology.js";

/** Where one chain meets a junction: the outermost edge, and which of the chain's two ends. */
export interface ChainEnd {
  /** Index into the chain list this end was derived from. */
  chain: number;

  /** `0` for `verts[0]`, `1` for the last vertex. */
  end: 0 | 1;

  vertex: SolvableVertex;
  edge: SolvableEdge;

  /**
   * Whether this end's chain runs against its group's canonical direction through the vertex.
   *
   * A vertex block is written along the direction the chain travels *through* the vertex, and
   * two chains meeting at a node need not travel the same way. Only {@link interfaceSlots}
   * decides this, because "the same way" is a statement about a whole group and not about one
   * end; {@link chainEnds} leaves it `false`.
   */
  reversed: boolean;
}

/**
 * One interface unknown.
 *
 * `"entry"` is a block entry that several edge-ends read out of the same scalar — order `n`
 * of §3's `n`th partition, restricted to groups with more than one member. `"row"` is the
 * multiplier of a G1 row at a junction, one per edge of a spanning tree of a tangent group,
 * which is `m − 1` rows for a group of `m` and not `m(m−1)/2`.
 */
export type GammaSlot =
  | { kind: "entry"; vertex: SolvableVertex; order: number; ends: ChainEnd[] }
  | { kind: "row"; vertex: SolvableVertex; a: ChainEnd; b: ChainEnd };

/**
 * Chains that have to be solved together, and the unknowns that make it so.
 *
 * A component with an empty {@link gamma} is a single chain and nothing else — the Schur
 * complement is `0 × 0` and the solve is the banded one it always was.
 */
export interface Component<C extends Curve = Curve> {
  chains: Chain<C>[];

  /** Positions of {@link chains} in the list they came from, which is the diagnostics' chain index. */
  indices: number[];

  /** Interface unknowns, numbered component-locally. {@link ChainEnd.chain} indexes {@link chains}. */
  gamma: GammaSlot[];
}

function find(parent: Int32Array, i: number): number {
  while (parent[i] !== i) {
    parent[i] = parent[parent[i]];
    i = parent[i];
  }

  return i;
}

/** Every chain end, keyed by the vertex it lands on, in chain order. */
export function chainEnds<C extends Curve>(all: Chain<C>[]) {
  const out = new Map<SolvableVertex, ChainEnd[]>();

  const add = (end: ChainEnd) => {
    const list = out.get(end.vertex);

    if (list) {
      list.push(end);
    } else {
      out.set(end.vertex, [end]);
    }
  };

  for (let c = 0; c < all.length; c++) {
    const { verts, edges } = all[c];

    add({ chain: c, end: 0, vertex: verts[0], edge: edges[0] as SolvableEdge, reversed: false });
    add({
      chain   : c,
      end     : 1,
      vertex  : verts[verts.length - 1],
      edge    : edges[edges.length - 1] as SolvableEdge,
      reversed: false,
    });
  }

  return out;
}

/** The members of each group of a partition, as index lists, in first-appearance order. */
function groupsOf(labels: Int32Array) {
  const out = new Map<number, number[]>();

  for (let i = 0; i < labels.length; i++) {
    const list = out.get(labels[i]);

    if (list) {
      list.push(i);
    } else {
      out.set(labels[i], [i]);
    }
  }

  return [...out.values()];
}

/**
 * Every interface unknown among `all`, with {@link ChainEnd.chain} indexing `all`.
 *
 * `all` must be paths — run every chain through `cutOpen` first, or a closed one will
 * report the same vertex at both ends without a block at the second and the wrap-around will
 * be silently dropped. That cut is exactly why this takes a chain *list* rather than a mesh:
 * the promotion of a cut vertex to a junction is a decision about the solve, and this
 * function is where it takes effect rather than where it is made.
 *
 * `p` is needed because the partitions are: a level of `3` shares two entries at `p = 1` and
 * three at `p = 2`, so how many slots a junction contributes is a function of the order being
 * solved at. Continuation (§9) rebuilds this per rung for that reason — but see
 * {@link components} for why the *grouping* need not be rebuilt with it.
 */
export function interfaceSlots<C extends Curve>(all: Chain<C>[], p: number, channel = curvatureChannel) {
  const slots: GammaSlot[] = [];

  for (const [vertex, list] of chainEnds(all)) {
    if (list.length < 2) {
      continue;
    }

    const byEdge = new Map<SolvableEdge, ChainEnd>();

    for (const end of list) {
      byEdge.set(end.edge, end);
    }

    // An incident edge always has an end here — chains stop at nodes and are cut at cut
    // points — but a group is filtered rather than asserted, since a caller may pass a subset.
    const members = (group: number[]) =>
      group.map((i) => byEdge.get(vertex.edges[i])).filter((e): e is ChainEnd => e !== undefined);

    const parts = vertexPartitions(vertex, p, channel);
    const tangents = groupsOf(parts[0]).map(members);

    /*
      Orientation, along the same spanning tree the G1 rows use. Two ends agree on the direction
      through the node exactly when one arrives and the other leaves, so equal `end` is a flip —
      and a group of three or more cannot satisfy that pairwise, which is why it has to be read
      off a tree rather than off every pair. Every other partition subdivides the coarsest one,
      so a tree over that one covers them all; width's coarsest is order 0, not the tangents.
    */
    for (const group of groupsOf(parts[channel.coarsest])) {
      const ends = members(group);

      for (let k = 1; k < ends.length; k++) {
        ends[k].reversed = ends[k - 1].reversed !== (ends[k].end === ends[k - 1].end);
      }
    }

    for (let n = 0; n <= p; n++) {
      for (const group of groupsOf(parts[n + 1])) {
        const shared = members(group);

        if (shared.length < 2) {
          continue;
        }

        slots.push({ kind: "entry", vertex, order: n, ends: shared });
      }
    }

    for (const tangent of tangents) {
      for (let k = 1; k < tangent.length; k++) {
        slots.push({ kind: "row", vertex, a: tangent[k - 1], b: tangent[k] });
      }
    }
  }

  return slots;
}

/** The chains a slot touches. */
function touched(slot: GammaSlot) {
  return slot.kind === "entry" ? slot.ends.map((e) => e.chain) : [slot.a.chain, slot.b.chain];
}

/**
 * Group chains into independently solvable components, with the interface between them.
 *
 * The grouping does not depend on `p`, even though the slot *list* does, so a continuation
 * ladder may call this once at the top order and re-derive only {@link interfaceSlots} per
 * rung. Tangent grouping needs level 1 and order-`0` sharing needs level 2, both of which
 * every `p ≥ 0` can represent; lowering `p` only drops the *higher* shared orders, which
 * never groups two chains that were not already grouped by order `0`.
 */
export function components<C extends Curve>(all: Chain<C>[], p: number, channel = curvatureChannel): Component<C>[] {
  const slots = interfaceSlots(all, p, channel);
  const parent = new Int32Array(all.length);

  for (let i = 0; i < all.length; i++) {
    parent[i] = i;
  }

  for (const slot of slots) {
    const chains = touched(slot);

    for (let k = 1; k < chains.length; k++) {
      parent[find(parent, chains[0])] = find(parent, chains[k]);
    }
  }

  const grouped = new Map<number, number[]>();

  for (let c = 0; c < all.length; c++) {
    const root = find(parent, c);
    const list = grouped.get(root);

    if (list) {
      list.push(c);
    } else {
      grouped.set(root, [c]);
    }
  }

  const out: Component<C>[] = [];

  for (const indices of grouped.values()) {
    const local = new Map(indices.map((c, i) => [c, i]));
    const rebase = (end: ChainEnd) => ({ ...end, chain: local.get(end.chain) ?? -1 });

    const gamma = slots
      .filter((slot) => local.has(slot.kind === "entry" ? slot.ends[0].chain : slot.a.chain))
      .map((slot) =>
        slot.kind === "entry"
          ? { ...slot, ends: slot.ends.map(rebase) }
          : { ...slot, a: rebase(slot.a), b: rebase(slot.b) }
      );

    out.push({ chains: indices.map((c) => all[c]), indices, gamma });
  }

  return out;
}
