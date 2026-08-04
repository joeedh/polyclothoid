/**
 * Continuity levels, and the per-order partitions they induce on a vertex's edge-ends —
 * `docs/plans/spower-solver.md` §3.
 *
 * The authored primitive is a {@link Pairing}: two edge-ends at a vertex and a level.
 *
 * | Level | Meaning | Structure |
 * |---|---|---|
 * | 0 | position only — a corner | no G1 row, no shared block entry |
 * | 1 | G1 | G1 row, no shared entry |
 * | 2 | G2 | G1 row + shared entry `0` (`κ`) |
 * | `k` | G`k` | G1 row + shared entries `0 … k−2` |
 *
 * A single partition of the edge-ends cannot express that, because "broken tangent" and
 * "broken curvature" are different amounts of breaking and a vertex may carry both at once.
 * What can is **one partition per order, each a coarsening of the next** — `p + 2` of them,
 * indexed `n = −1 … p`, where `n = −1` says which ends are tangent-coupled (the G1 rows) and
 * `n ≥ 0` says which ends share block entry `n`.
 *
 * The `n = −1` partition is the one that is easy to omit, and its omission is silent: without
 * it a level-1 pairing groups nothing, so the G1 row at a G1-only joint is simply never
 * emitted and the joint comes apart with no error anywhere.
 *
 * Transitivity is why these are union-finds rather than pairwise lookups. Level-3 A–B and
 * level-2 B–C put A and C in one group at order `0` without anyone having authored A–C, which
 * is the correct reading — they share the same `κ` through B — and is not visible pairwise.
 * At valence 2 every group has at most two members and the union-find is a formality; it is
 * written for the junctions Phase 8 adds.
 */
import { type Pairing, type SolvableEdge, type SolvableVertex } from "./mesh_types.js";

/** The highest level order `p` can express: a G1 row plus all `p + 1` block entries shared. */
export function maxLevel(p: number) {
  return p + 2;
}

/**
 * The level a pairing must reach to group two ends in partition `n`.
 *
 * `n = −1` is the tangent partition and needs only a G1 pairing; sharing block entry `n`
 * needs level `n + 2`.
 */
export function levelFor(n: number) {
  return n < 0 ? 1 : n + 2;
}

/**
 * The default level between two edge-ends at a vertex nobody has authored a pairing for.
 *
 * Valence ≤ 2 defaults to the ceiling: the pairing is unambiguous, so there is nothing to
 * guess and no reason to deliver less than was asked for. Valence ≥ 3 defaults to `0` —
 * fully split — because no canonical pairing exists there and any other default would be the
 * solver inferring structure, which §3 forbids. Fully split is the honest answer at a
 * junction, not a degraded one.
 */
export function defaultLevel(v: SolvableVertex, p: number) {
  return v.edges.length <= 2 ? maxLevel(p) : 0;
}

/** The authored pairing of `a` and `b` at `v`, in either order, if there is one. */
export function findPairing(v: SolvableVertex, a: SolvableEdge, b: SolvableEdge): Pairing | undefined {
  for (const pairing of v.pairings ?? []) {
    if ((pairing.a === a && pairing.b === b) || (pairing.a === b && pairing.b === a)) {
      return pairing;
    }
  }

  return undefined;
}

/**
 * The level in force between two edge-ends: authored if it is, defaulted if it is not, and
 * clamped to what order `p` can represent.
 *
 * Clamping down is not a break — a level of `5` at `p = 1` asks for derivatives that do not
 * exist in the state space, so there is nothing to report and nothing was lost.
 */
export function pairingLevel(v: SolvableVertex, a: SolvableEdge, b: SolvableEdge, p: number) {
  const authored = findPairing(v, a, b);
  const level = authored ? authored.level : defaultLevel(v, p);

  return Math.max(0, Math.min(level, maxLevel(p)));
}

function find(parent: Int32Array, i: number): number {
  while (parent[i] !== i) {
    parent[i] = parent[parent[i]];
    i = parent[i];
  }

  return i;
}

/**
 * The `p + 2` partitions of `v.edges`, indexed `[n + 1]` for `n = −1 … p`.
 *
 * Each entry labels every incident edge-end with the representative of its group, so two ends
 * are grouped exactly when their labels match. Each partition is a coarsening of the one
 * below it, since {@link levelFor} is increasing.
 */
export function vertexPartitions(v: SolvableVertex, p: number) {
  const m = v.edges.length;
  const out: Int32Array[] = [];

  for (let n = -1; n <= p; n++) {
    const parent = new Int32Array(m);

    for (let i = 0; i < m; i++) {
      parent[i] = i;
    }

    const needed = levelFor(n);

    for (let i = 0; i < m; i++) {
      for (let j = i + 1; j < m; j++) {
        if (pairingLevel(v, v.edges[i], v.edges[j], p) >= needed) {
          parent[find(parent, i)] = find(parent, j);
        }
      }
    }

    for (let i = 0; i < m; i++) {
      parent[i] = find(parent, i);
    }

    out.push(parent);
  }

  return out;
}

/** What two edge-ends at a vertex actually hold in common. See {@link sharing}. */
export interface Sharing {
  /** Leading block entries the two ends read out of the same DOF: `0 … entries − 1`. */
  entries: number;

  /** Whether a G1 row couples them. */
  tangent: boolean;
}

/**
 * Read {@link vertexPartitions} for one pair of ends.
 *
 * `entries` counts agreeing orders rather than deriving them from a level, because at a
 * junction the group at order `n` may include ends nobody paired directly — see the module
 * comment. The coarsening property makes the agreeing orders a prefix, so a count is enough.
 */
export function sharing(v: SolvableVertex, a: SolvableEdge, b: SolvableEdge, p: number): Sharing {
  const parts = vertexPartitions(v, p);

  const ia = v.edges.indexOf(a);
  const ib = v.edges.indexOf(b);

  if (ia < 0 || ib < 0) {
    return { entries: 0, tangent: false };
  }

  let entries = 0;

  for (let n = 0; n <= p; n++) {
    if (parts[n + 1][ia] === parts[n + 1][ib]) {
      entries++;
    }
  }

  return { entries, tangent: parts[0][ia] === parts[0][ib] };
}
