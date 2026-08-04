/**
 * Width as a second scalar profile on the same chains — `docs/plans/spower-solver.md` §10 and
 * Phase 9 of §12.
 *
 * Everything structural is shared with `κ`. Width lives in the same s-power basis, on the same
 * vertex-owned blocks, split at joints by the same authored levels, eliminated across junctions
 * by the same Schur complement. What is *not* shared is the three things §10 says are genuinely
 * different, and they are the whole of this file:
 *
 * 1. **There are no constraint rows.** `κ`'s G1 rows exist because sharing the curvature value
 *    at a joint does not make the tangents agree; two edges can read the same `κ` and still meet
 *    at a kink, so the tangent is a solved condition. A width joint has no such gap — the two
 *    ends reading the same order-0 entry *is* the width being continuous. So the system is
 *    symmetric positive definite rather than a KKT saddle, `interfaceSlots` produces only
 *    `"entry"` slots, and there is no Newton loop: it is one quadratic program.
 * 2. **The levels shift down by one.** Level `k` on the width channel shares entries `0 … k−1`,
 *    against `κ`'s `0 … k−2`, for the same reason — there is no row to buy at level 1. See
 *    `ChannelRules` in `pairing.ts`.
 * 3. **There is a data term and there are bounds.** Width is *given*, per sample, by pressure
 *    input, so the objective is `Σ ωᵢ(w(uᵢ) − pᵢ)² + λ·Σ_e E_e` rather than pure energy; and an
 *    unconstrained fit overshoots between samples, where a negative width turns the offset
 *    inside out. The bound is `w ≥ w_min` on the order-0 entries, run as an active set.
 *
 * ## What the bounds actually bound
 *
 * Order-0 entries are `w` at a vertex, so clamping them is clamping the width *at the vertices*
 * and nowhere else. Between two vertices the profile is a degree-`2p+1` polynomial and nothing
 * stops it dipping below `w_min` in the interior — the true constraint is semi-infinite and this
 * one is its collocation at the joints. Rather than pretend otherwise, {@link WidthSolver}
 * samples each edge afterwards and *reports* interior undershoot as
 * {@link WidthSolverReport.undershoot}. Adding the offending point as a new bound and re-solving
 * is the standard cutting-plane fix and would fit the active set unchanged; it is not done here
 * because it costs a factorization per cut and §10 does not ask for it.
 *
 * That the bound survives the interface at all is a small piece of luck worth naming: a width
 * does not change sign when you walk the other way, so its orientation sign is `(−1)ⁿ` and
 * order 0 is even. Every end mirroring a shared order-0 entry mirrors it with `+1`, so
 * `γ ≥ w_min` and `w ≥ w_min` are the same inequality. The odd entries do flip, and nothing is
 * bounded there.
 *
 * ## Ordering
 *
 * §10's one-directional coupling: `w`'s domain is arclength, arclength is `KSCALE`, and `KSCALE`
 * is a `κ` output. Run `SPowerSolver` first and this second. There is no outer loop — width does
 * not feed back into the geometry — but running this on curves that have not been solved reads
 * whatever `L_e` they happen to hold.
 */
import { BandedSymmetric, DenseLU, KKTFactorization, defaultKKTSolveOptions } from "../math/index.js";
import { type KKTSolveOptions } from "../math/banded.js";
import { type EdgeFrame, widthDOF } from "./blocks.js";
import { type ChainEnd, type GammaSlot, components } from "./junctions.js";
import { type SolvableEdge, type SolvableMesh, type SolvableVertex } from "./mesh_types.js";
import { sharing, widthChannel } from "./pairing.js";
import { SPOWER_ORDER, SPowerClothoid } from "./spower_clothoid.js";
import { evalSPower } from "./spower.js";
import { referenceLengths } from "./spower_solver.js";
import { type Chain, chains, cutOpen } from "./topology.js";

/**
 * One pressure reading: a width the fit should pass near, at a point on an edge.
 *
 * `u` is the edge's *own* fraction of arclength, `0` at `v1` and `1` at `v2`, so a client that
 * knows a world arclength `s` along the edge passes `s / edge.curve.length`. Chain orientation
 * is not the client's business and the frame folds it in.
 *
 * `width` is whatever quantity the stroker will offset by twice — a full stroke width, on the
 * convention that the offset curves sit at `±width/2`. {@link WidthSolverReport.cusps} reads it
 * that way and nothing else in here cares.
 */
export interface WidthSample {
  edge: SolvableEdge;
  u: number;
  width: number;

  /** Relative confidence, defaulting to `1`. Enters the normal equations as `ω`. */
  weight?: number;
}

export interface WidthSolverOptions {
  /** The s-power order `p`, which need not be `κ`'s. */
  order: number;

  /**
   * Edit-locality parameter, as §6 — the foundation term in `widthEnergy`.
   *
   * Far smaller than `κ`'s default, because the two are not the same knob wearing one name. For
   * curvature the foundation pulls towards `κ = 0`, which is a straight line and a perfectly
   * reasonable thing to be pulled towards. For width it pulls towards `w = 0`, which is a
   * degenerate stroke. It is here for the decay length it buys the smoothing operator and for
   * definiteness on a chain that carries no samples at all, not as a prior worth having.
   */
  alpha: number;

  /**
   * Weight on the smoothing energy against the data term, in units of the mean chord.
   *
   * The trade is dimensional and cannot be made otherwise: the data term is a squared width and
   * `∫(dw/ds)² ds` is a width squared per length, so a bare number cannot stand between them.
   * The mean chord of the mesh supplies the missing length, which makes this dimensionless and
   * scale-invariant — the same number gives the same fit on a drawing zoomed by any factor.
   */
  smoothing: number;

  /** `w_min` — the lower bound the active set enforces at every vertex. */
  minimum: number;

  /** Active-set iterations before giving up and reporting what it had. */
  iterations: number;

  /** Bound violations below this are not worth a factorization. Relative to the mean sampled width. */
  tolerance: number;

  /** Points per edge at which the cusp and undershoot checks look. */
  probes: number;

  kkt: Partial<KKTSolveOptions>;
}

export const defaultWidthSolverOptions: WidthSolverOptions = {
  order     : SPOWER_ORDER,
  alpha     : 1e-2,
  smoothing : 0.05,
  minimum   : 0.0,
  iterations: 16,
  tolerance : 1e-9,
  probes    : 17,
  kkt       : {},
};

/**
 * A point where the offset at half-width cusps — §10's second difference, reported and not acted
 * on.
 *
 * The offset of a curve at distance `h` is regular only where `h·|κ| < 1`; at `h = 1/κ` its
 * tangent vanishes and past it the outline turns inside out. Two individually valid profiles can
 * produce this, which is exactly why it is checked after both are solved. What to *do* about it
 * — clamp the width, subdivide, let the outline self-intersect and rely on a nonzero winding
 * fill — is a stroker policy question, so this says where and how badly and stops.
 */
export interface WidthCusp {
  edge: SolvableEdge;
  u: number;

  /** World arclength along the edge, i.e. `u · L_e`. */
  s: number;

  width: number;
  curvature: number;

  /** `|κ|·width/2`. At least `1` by construction, and how far past tells you how bad. */
  product: number;
}

/** A point where the fitted profile dips under `minimum` between two vertices — see the module doc. */
export interface WidthUndershoot {
  edge: SolvableEdge;
  u: number;
  width: number;
}

export interface WidthSolverReport {
  chains: number;
  components: number;

  /** Total unknowns across every component, interface included. */
  unknowns: number;

  /** Samples that landed on an edge the solver owns. */
  samples: number;

  /** Active-set iterations, summed over components. */
  iterations: number;

  /** Order-0 entries sitting on their bound when the solve finished. */
  active: number;

  /** False if any component's factorization failed. The rest of the report still describes what ran. */
  ok: boolean;

  /** Worst `‖Kz − b‖∞` over the components. */
  residual: number;

  /** Weighted RMS of `w(uᵢ) − pᵢ`, or `0` with no samples. */
  rms: number;

  cusps: WidthCusp[];
  undershoot: WidthUndershoot[];
}

interface ChainSample {
  edge: number;
  u: number;
  width: number;
  weight: number;
}

/**
 * One chain's width unknowns and its assembled normal equations.
 *
 * The layout is `ChainSystem`'s with the multipliers taken out: `[b₀][b₁]…[b_m]`, each interior
 * block split into the entries the two ends share and a private tail apiece. It stays banded for
 * the same reason — an edge touches two adjacent blocks and nothing further.
 */
export class WidthChain {
  p: number;
  block: number;

  frames: EdgeFrame[] = [];
  curves: SPowerClothoid[] = [];

  blockAt: Int32Array;
  arrivingAt: Int32Array;
  leavingAt: Int32Array;

  /** Entries the two ends at `verts[i]` read out of the same scalar. */
  shared: Int32Array;

  n = 0;
  z: Float64Array;

  kkt: BandedSymmetric;

  /** The data term's right-hand side, `Σ ω·p·r`, in local DOF order. */
  rhs: Float64Array;

  /** Interface index of each local unknown, or `-1` where it is private. Set by {@link WidthComponent}. */
  localToGlobal: Int32Array;

  /** `±1` relating each mirrored unknown to the interface value it mirrors. */
  localSign: Float64Array;

  samples: ChainSample[] = [];

  constructor(
    public chain: Chain,
    refs: Map<SolvableVertex, number>,
    public options: WidthSolverOptions,

    /** The length that makes {@link WidthSolverOptions.smoothing} dimensionless. */
    public reference: number
  ) {
    const { verts, edges } = chain;

    this.p = options.order;
    this.block = widthDOF.blockLength(this.p);

    for (let i = 0; i < edges.length; i++) {
      const e = edges[i] as SolvableEdge;
      const curve = e.curve as SPowerClothoid;

      curve.update(e);
      this.curves.push(curve);

      const chord = verts[i].vectorDistance(verts[i + 1]);

      this.frames.push({
        chord,
        arclength: curve.length > 0.0 ? curve.length : chord,
        rEarlier : refs.get(verts[i]) ?? 1.0,
        rLater   : refs.get(verts[i + 1]) ?? 1.0,
        forward  : e.v1 === verts[i],
      });
    }

    this.blockAt = new Int32Array(verts.length);
    this.arrivingAt = new Int32Array(verts.length).fill(-1);
    this.leavingAt = new Int32Array(verts.length).fill(-1);
    this.shared = new Int32Array(verts.length).fill(this.block);

    for (let i = 1; i < edges.length; i++) {
      const { entries } = sharing(verts[i], edges[i - 1], edges[i], this.p, widthChannel);

      this.shared[i] = Math.min(entries, this.block);
    }

    this.layout();

    this.z = new Float64Array(this.n);
    this.rhs = new Float64Array(this.n);
    this.localToGlobal = new Int32Array(this.n).fill(-1);
    this.localSign = new Float64Array(this.n).fill(1.0);
    this.kkt = new BandedSymmetric(this.n, this.bandwidth());
  }

  layout() {
    const m = this.chain.edges.length;
    let at = 0;

    for (let i = 0; i < this.chain.verts.length; i++) {
      const interior = i > 0 && i < m;

      this.blockAt[i] = at;
      at += this.shared[i];

      if (interior) {
        const rest = this.block - this.shared[i];

        this.arrivingAt[i] = at;
        at += rest;

        this.leavingAt[i] = at;
        at += rest;
      }
    }

    this.n = at;
  }

  /** Widest offset any edge's DOF span, measured off the layout rather than predicted from it. */
  bandwidth() {
    const d = 2 * this.block;
    let w = 0;

    for (let i = 0; i < this.chain.edges.length; i++) {
      let lo = Infinity;
      let hi = -Infinity;

      for (let k = 0; k < d; k++) {
        const at = this.at(i, k);

        lo = Math.min(lo, at);
        hi = Math.max(hi, at);
      }

      w = Math.max(w, hi - lo);
    }

    return Math.min(w, Math.max(0, this.n - 1));
  }

  /** Global index of entry `n` of `verts[i]`, as read by the end arriving at it or leaving it. */
  slot(i: number, n: number, arriving: boolean) {
    if (n < this.shared[i]) {
      return this.blockAt[i] + n;
    }

    return (arriving ? this.arrivingAt[i] : this.leavingAt[i]) + n - this.shared[i];
  }

  /** Global index of local DOF `k` of edge `i`, whose DOF are its two blocks in chain order. */
  at(i: number, k: number) {
    return k < this.block ? this.slot(i, k, false) : this.slot(i + 1, k - this.block, true);
  }

  endVertex(end: 0 | 1) {
    return end === 1 ? this.frames.length : 0;
  }

  endSlot(end: 0 | 1, n: number) {
    return this.slot(this.endVertex(end), n, end === 1);
  }

  /** Local indices holding a width *value*, i.e. an order-0 entry — the ones the bound applies to. */
  bounded(into: number[]) {
    for (let i = 0; i < this.chain.verts.length; i++) {
      const leaving = this.slot(i, 0, false);
      const arriving = this.slot(i, 0, true);

      into.push(leaving);

      if (arriving !== leaving) {
        into.push(arriving);
      }
    }

    return into;
  }

  /** Scatter `λ·Σ M_eᵀE_eM_e` and the sample normal equations into the band and `rhs`. */
  assemble() {
    const { p, kkt, rhs, frames } = this;
    const d = 2 * this.block;
    const lambda = this.options.smoothing * this.reference;

    kkt.zero();
    rhs.fill(0.0);

    for (let i = 0; i < frames.length; i++) {
      const h = widthDOF.hessian(p, frames[i], this.options.alpha);

      // Local and global index order agree, so the local lower triangle is the global one.
      for (let r = 0; r < d; r++) {
        for (let c = 0; c <= r; c++) {
          kkt.add(this.at(i, r), this.at(i, c), lambda * h[r * d + c]);
        }
      }
    }

    for (const sample of this.samples) {
      const row = widthDOF.pullback(p, frames[sample.edge], widthDOF.valueWeights(p, sample.u));
      const w = sample.weight;

      for (let r = 0; r < d; r++) {
        rhs[this.at(sample.edge, r)] += w * sample.width * row[r];

        for (let c = 0; c <= r; c++) {
          kkt.add(this.at(sample.edge, r), this.at(sample.edge, c), w * row[r] * row[c]);
        }
      }
    }
  }

  /** Push the current DOF through the width transform into each curve's `ws`. */
  write() {
    const { p, frames, curves } = this;
    const count = widthDOF.coefficientLength(p);
    const scratch = new Float64Array(2 * this.block);

    for (let i = 0; i < curves.length; i++) {
      for (let k = 0; k < 2 * this.block; k++) {
        scratch[k] = this.z[this.at(i, k)];
      }

      if (curves[i].ws.length !== count) {
        curves[i].ws = new Float64Array(count);
      }

      widthDOF.coefficients(p, frames[i], scratch, curves[i].ws);
    }
  }
}

/**
 * Chains that share width unknowns, eliminated onto the interface between them — §5's
 * substructuring with no multipliers in it.
 *
 * `S = D − Σ Bᵢᵀ Aᵢ⁻¹ Bᵢ`, the same complement `ComponentSystem` forms, over a symmetric
 * positive definite `Aᵢ` instead of a quasi-definite KKT block. An empty interface makes this
 * one banded solve and nothing else.
 */
export class WidthComponent {
  /** Position of each chain-local unknown among its chain's private ones, `-1` for interface. */
  interior: Int32Array[] = [];

  offsets: number[] = [];
  counts: number[] = [];

  gammaAt = 0;
  n = 0;

  /** Fixed-at-bound flags and values, in global indices — the active set. */
  clamped: Uint8Array;
  bound: Float64Array;

  /** Global indices carrying a width value, which is where {@link clamped} may be set. */
  bounds: number[] = [];

  blocks: BandedSymmetric[] = [];
  coupling: Float64Array[] = [];
  d: Float64Array;
  schur: Float64Array;
  ys: Float64Array[] = [];

  b: Float64Array;

  /**
   * The data term's right-hand side before {@link build} folds the active set into it.
   *
   * Kept because the two are not interchangeable where it matters most: `b` has had the clamped
   * columns subtracted out of the free rows and the clamped rows overwritten by their bounds,
   * and the multiplier of a clamped unknown is `(Kz − b)` against the *unclamped* right-hand
   * side. Reading it off `b` instead compares the gradient to `w_min` and cycles.
   */
  b0: Float64Array;
  z: Float64Array;

  private factors: KKTFactorization[] = [];
  private solution: Float64Array[] = [];
  private column: Float64Array[] = [];
  private reduced: Float64Array;
  private work: Float64Array;
  private masked: Float64Array;
  private correction: Float64Array;
  private lu?: DenseLU;

  constructor(
    public systems: WidthChain[],
    public gamma: GammaSlot[]
  ) {
    for (const s of systems) {
      s.localToGlobal.fill(-1);
      s.localSign.fill(1.0);
    }

    for (let g = 0; g < gamma.length; g++) {
      const slot = gamma[g];

      if (slot.kind !== "entry") {
        continue;
      }

      for (const end of slot.ends) {
        const s = systems[end.chain];
        const at = s.endSlot(end.end, slot.order);

        s.localToGlobal[at] = g;
        s.localSign[at] = WidthComponent.orientation(end, slot.order);
      }
    }

    const m = gamma.length;
    let at = 0;

    for (const s of systems) {
      const map = new Int32Array(s.n).fill(-1);
      let count = 0;

      for (let i = 0; i < s.n; i++) {
        if (s.localToGlobal[i] < 0) {
          map[i] = count++;
        }
      }

      this.interior.push(map);
      this.offsets.push(at);
      this.counts.push(count);

      this.blocks.push(new BandedSymmetric(count, Math.min(s.kkt.bandwidth, Math.max(0, count - 1))));
      this.coupling.push(new Float64Array(count * m));
      this.ys.push(new Float64Array(count * m));
      this.solution.push(new Float64Array(count));
      this.column.push(new Float64Array(count));

      at += count;
    }

    this.gammaAt = at;
    this.n = at + m;

    for (let c = 0; c < systems.length; c++) {
      for (const local of systems[c].bounded([])) {
        this.bounds.push(this.globalOf(c, local));
      }
    }

    this.bounds = [...new Set(this.bounds)];

    this.d = new Float64Array(m * m);
    this.schur = new Float64Array(m * m);
    this.reduced = new Float64Array(m);
    this.b = new Float64Array(this.n);
    this.b0 = new Float64Array(this.n);
    this.z = new Float64Array(this.n);
    this.work = new Float64Array(this.n);
    this.masked = new Float64Array(this.n);
    this.correction = new Float64Array(this.n);
    this.clamped = new Uint8Array(this.n);
    this.bound = new Float64Array(this.n);
  }

  /**
   * The sign relating entry `order` at `end` to the interface value it mirrors.
   *
   * `(−1)ⁿ`, where `κ`'s is `(−1)ⁿ⁺¹`: reversing the direction of travel flips `d/ds` for both,
   * and flips the field for `κ` alone. So width's *even* entries — the value among them — come
   * through a reversed end unchanged, and its odd ones flip.
   */
  static orientation(end: ChainEnd, order: number) {
    return end.reversed && (order & 1) === 1 ? -1.0 : 1.0;
  }

  /** Global index of chain `c`'s local unknown `i`. */
  globalOf(c: number, i: number) {
    const g = this.systems[c].localToGlobal[i];

    return g >= 0 ? this.gammaAt + g : this.offsets[c] + this.interior[c][i];
  }

  /** Distribute the global solution back into the chains, mirroring the interface as it goes. */
  take() {
    for (let c = 0; c < this.systems.length; c++) {
      const s = this.systems[c];

      for (let i = 0; i < s.n; i++) {
        s.z[i] = s.localSign[i] * this.z[this.globalOf(c, i)];
      }
    }
  }

  /**
   * Split every chain's band into `Aᵢ`, `Bᵢ` and `D`, and form `b`, with the active set folded in.
   *
   * A clamped unknown leaves the system entirely: its row and column become the identity, and
   * every entry that used to couple it to a free unknown moves to that unknown's right-hand
   * side. Doing it here rather than after assembly is what keeps `Aᵢ` banded and positive
   * definite — a fixed variable is not a constraint row, it is one fewer variable.
   */
  build() {
    const m = this.gamma.length;

    this.d.fill(0.0);
    this.b.fill(0.0);
    this.b0.fill(0.0);

    for (const a of this.blocks) {
      a.zero();
    }

    for (const b of this.coupling) {
      b.fill(0.0);
    }

    for (let c = 0; c < this.systems.length; c++) {
      const s = this.systems[c];

      s.assemble();

      for (let i = 0; i < s.n; i++) {
        this.b[this.globalOf(c, i)] += s.localSign[i] * s.rhs[i];
      }
    }

    this.b0.set(this.b);

    for (let c = 0; c < this.systems.length; c++) {
      const s = this.systems[c];
      const bw = s.kkt.bandwidth;

      for (let j = 0; j < s.n; j++) {
        const hi = Math.min(s.n - 1, j + bw);

        for (let i = j; i <= hi; i++) {
          const v = s.kkt.get(i, j);

          if (v !== 0.0) {
            this.scatter(c, i, j, s.localSign[i] * s.localSign[j] * v);
          }
        }
      }
    }

    for (const g of this.bounds) {
      if (!this.clamped[g]) {
        continue;
      }

      this.b[g] = this.bound[g];

      if (g >= this.gammaAt) {
        this.d[(g - this.gammaAt) * m + (g - this.gammaAt)] = 1.0;
      } else {
        const c = this.chainOf(g);

        this.blocks[c].set(g - this.offsets[c], g - this.offsets[c], 1.0);
      }
    }
  }

  /** Which chain a private global index belongs to. */
  private chainOf(g: number) {
    let c = 0;

    while (c + 1 < this.offsets.length && this.offsets[c + 1] <= g) {
      c++;
    }

    return c;
  }

  /** One symmetric matrix entry, routed by where its two unknowns live and whether they are clamped. */
  private scatter(c: number, i: number, j: number, v: number) {
    const m = this.gamma.length;
    const s = this.systems[c];

    const gi = this.globalOf(c, i);
    const gj = this.globalOf(c, j);

    if (this.clamped[gi] || this.clamped[gj]) {
      /* Diagonal and clamped-to-clamped are carried by the identity rows `build` writes. */
      if (!this.clamped[gj]) {
        this.b[gj] -= v * this.bound[gi];
      } else if (!this.clamped[gi]) {
        this.b[gi] -= v * this.bound[gj];
      }

      return;
    }

    const li = s.localToGlobal[i];
    const lj = s.localToGlobal[j];
    const map = this.interior[c];

    if (li < 0 && lj < 0) {
      this.blocks[c].add(map[i], map[j], v);
    } else if (li < 0) {
      this.coupling[c][map[i] * m + lj] += v;
    } else if (lj < 0) {
      this.coupling[c][map[j] * m + li] += v;
    } else {
      this.d[li * m + lj] += v;

      if (i !== j) {
        this.d[lj * m + li] += v;
      }
    }
  }

  /**
   * `out = K x` for the *unclamped* global matrix, walked off the chain bands.
   *
   * The active set needs this and {@link build}'s `A`/`B`/`D` cannot give it: those describe the
   * reduced problem, in which a clamped unknown has no row left to read a multiplier off.
   */
  applyUnfixed(x: Float64Array, out: Float64Array) {
    out.fill(0.0);

    for (let c = 0; c < this.systems.length; c++) {
      const s = this.systems[c];
      const bw = s.kkt.bandwidth;

      for (let j = 0; j < s.n; j++) {
        const hi = Math.min(s.n - 1, j + bw);
        const gj = this.globalOf(c, j);

        for (let i = j; i <= hi; i++) {
          const v = s.kkt.get(i, j);

          if (v === 0.0) {
            continue;
          }

          const gi = this.globalOf(c, i);
          const both = s.localSign[i] * s.localSign[j] * v;

          out[gi] += both * x[gj];

          if (i !== j) {
            out[gj] += both * x[gi];
          }
        }
      }
    }

    return out;
  }

  /** Factor every `Aᵢ`, then form and factor `S`. False if a factorization underflowed. */
  factor() {
    const m = this.gamma.length;

    this.factors.length = 0;

    for (let c = 0; c < this.systems.length; c++) {
      const f = new KKTFactorization(this.blocks[c], [], this.systems[c].options.kkt);

      if (!f.ok) {
        return false;
      }

      this.factors.push(f);
    }

    if (m === 0) {
      return true;
    }

    this.schur.set(this.d);

    for (let c = 0; c < this.systems.length; c++) {
      const k = this.counts[c];
      const b = this.coupling[c];
      const y = this.ys[c];
      const col = this.column[c];
      const sol = this.solution[c];

      for (let g = 0; g < m; g++) {
        for (let i = 0; i < k; i++) {
          col[i] = b[i * m + g];
        }

        this.factors[c].solve(col, sol);

        for (let i = 0; i < k; i++) {
          y[i * m + g] = sol[i];
        }
      }

      for (let r = 0; r < m; r++) {
        for (let cc = 0; cc < m; cc++) {
          let sum = 0.0;

          for (let i = 0; i < k; i++) {
            sum += b[i * m + r] * y[i * m + cc];
          }

          this.schur[r * m + cc] -= sum;
        }
      }
    }

    this.lu = new DenseLU(this.schur, m);

    return this.lu.ok;
  }

  /** One pass of the substructured solve against the factorization {@link factor} left behind. */
  solveOnce(b: Float64Array, out: Float64Array) {
    const m = this.gamma.length;

    for (let c = 0; c < this.systems.length; c++) {
      const off = this.offsets[c];
      const k = this.counts[c];

      this.factors[c].solve(b.subarray(off, off + k), this.solution[c]);
      out.set(this.solution[c], off);
    }

    if (m === 0) {
      return;
    }

    const rs = this.reduced;

    for (let g = 0; g < m; g++) {
      rs[g] = b[this.gammaAt + g];
    }

    for (let c = 0; c < this.systems.length; c++) {
      const off = this.offsets[c];
      const k = this.counts[c];
      const bc = this.coupling[c];

      for (let i = 0; i < k; i++) {
        const xi = out[off + i];

        for (let g = 0; g < m; g++) {
          rs[g] -= bc[i * m + g] * xi;
        }
      }
    }

    this.lu?.solveInPlace(rs);

    for (let g = 0; g < m; g++) {
      out[this.gammaAt + g] = rs[g];
    }

    for (let c = 0; c < this.systems.length; c++) {
      const off = this.offsets[c];
      const k = this.counts[c];
      const y = this.ys[c];

      for (let i = 0; i < k; i++) {
        let sum = 0.0;

        for (let g = 0; g < m; g++) {
          sum += y[i * m + g] * rs[g];
        }

        out[off + i] -= sum;
      }
    }
  }

  /**
   * Solve the reduced system into {@link z}, refining against the unclamped operator.
   *
   * The refinement measures `b − K z` with the clamped unknowns held at their bounds and their
   * own rows zeroed, which is the residual of the problem actually being solved rather than of
   * the one the bands describe.
   */
  solveReduced(refinement: number) {
    if (!this.factor()) {
      this.z.fill(0.0);

      return Infinity;
    }

    this.solveOnce(this.b, this.z);

    for (const g of this.bounds) {
      if (this.clamped[g]) {
        this.z[g] = this.bound[g];
      }
    }

    let residual = Infinity;

    for (let pass = 0; pass <= refinement; pass++) {
      residual = this.residual(this.work);

      if (pass === refinement) {
        break;
      }

      this.solveOnce(this.work, this.correction);

      for (let i = 0; i < this.n; i++) {
        this.z[i] += this.correction[i];
      }
    }

    return residual;
  }

  /** `out = b − Kz` over the free unknowns, and `0` on the clamped ones. */
  residual(out: Float64Array) {
    for (let i = 0; i < this.n; i++) {
      this.masked[i] = this.clamped[i] ? 0.0 : this.z[i];
    }

    this.applyUnfixed(this.masked, out);

    let worst = 0.0;

    for (let i = 0; i < this.n; i++) {
      out[i] = this.clamped[i] ? 0.0 : this.b[i] - out[i];
      worst = Math.max(worst, Math.abs(out[i]));
    }

    return worst;
  }

  /**
   * Run the bounded QP: solve, clamp what went under, release what wants to come off — §10's
   * "active-set on bounds only, warm-startable, each iteration the same banded solve".
   *
   * The multiplier of a clamped unknown is `(Kz − b)` on its own row of the *unclamped* matrix,
   * and a negative one means the objective wants to move it up off its floor, which is allowed.
   * Adding every violator at once rather than one per pass is the usual gradient-projection
   * shortcut; the release step is one-at-a-time-safe because it only ever shrinks the set.
   */
  run(minimum: number, tolerance: number, iterations: number, refinement: number) {
    let residual = Infinity;
    let steps = 0;
    let ok = true;

    for (const g of this.bounds) {
      this.bound[g] = minimum;
    }

    for (let pass = 0; pass < Math.max(1, iterations); pass++) {
      steps++;

      this.build();
      residual = this.solveReduced(refinement);

      if (!Number.isFinite(residual)) {
        ok = false;

        break;
      }

      let added = 0;

      for (const g of this.bounds) {
        if (!this.clamped[g] && this.z[g] < minimum - tolerance) {
          this.clamped[g] = 1;
          added++;
        }
      }

      if (added > 0) {
        continue;
      }

      const gradient = this.applyUnfixed(this.z, this.work);
      let released = 0;

      for (const g of this.bounds) {
        if (this.clamped[g] && gradient[g] - this.b0[g] < -tolerance) {
          this.clamped[g] = 0;
          released++;
        }
      }

      if (released === 0) {
        break;
      }
    }

    this.take();

    let active = 0;

    for (const g of this.bounds) {
      active += this.clamped[g];
    }

    return { residual, steps, active, ok };
  }
}

/**
 * §10's width solve: `κ` first, then this.
 *
 * Not a {@link CurveSolver} — it does not move the curve, and its report answers different
 * questions. Construct it with the samples the fit should honour and call {@link solve} once;
 * there is no outer iteration to drive.
 */
export class WidthSolver {
  options: WidthSolverOptions;

  constructor(
    public mesh: SolvableMesh,
    public samples: WidthSample[] = [],
    options: Partial<WidthSolverOptions> = {}
  ) {
    this.options = { ...defaultWidthSolverOptions, ...options };
  }

  solve(): WidthSolverReport {
    const opts = this.options;
    const paths = chains(this.mesh).map(cutOpen);
    const refs = referenceLengths(this.mesh);
    const kkt = { ...defaultKKTSolveOptions, ...opts.kkt };

    const report: WidthSolverReport = {
      chains    : paths.length,
      components: 0,
      unknowns  : 0,
      samples   : 0,
      iterations: 0,
      active    : 0,
      ok        : true,
      residual  : 0.0,
      rms       : 0.0,
      cusps     : [],
      undershoot: [],
    };

    let chords = 0.0;
    let count = 0;

    for (const path of paths) {
      for (let i = 0; i < path.edges.length; i++) {
        chords += path.verts[i].vectorDistance(path.verts[i + 1]);
        count++;
      }
    }

    const reference = count === 0 ? 1.0 : chords / count;

    for (const component of components(paths, opts.order, widthChannel)) {
      const systems = component.chains.map((chain) => new WidthChain(chain, refs, opts, reference));
      const owner = new Map<SolvableEdge, { system: WidthChain; edge: number }>();

      for (const system of systems) {
        for (let i = 0; i < system.chain.edges.length; i++) {
          owner.set(system.chain.edges[i] as SolvableEdge, { system, edge: i });
        }
      }

      for (const sample of this.samples) {
        const at = owner.get(sample.edge);

        if (at === undefined) {
          continue;
        }

        const u = Math.min(Math.max(sample.u, 0.0), 1.0);

        at.system.samples.push({
          edge: at.edge,
          u,
          width : sample.width,
          weight: sample.weight ?? 1.0,
        });

        report.samples++;
      }

      const system = new WidthComponent(systems, component.gamma);
      const run = system.run(opts.minimum, opts.tolerance, opts.iterations, kkt.refinement);

      for (const chain of systems) {
        chain.write();
      }

      report.components++;
      report.unknowns += system.n;
      report.iterations += run.steps;
      report.active += run.active;
      report.ok &&= run.ok;
      report.residual = Math.max(report.residual, Number.isFinite(run.residual) ? run.residual : Infinity);
    }

    this.measure(report);

    return report;
  }

  /** The three post-hoc checks: the fit's RMS, §10's cusps, and interior undershoot. */
  private measure(report: WidthSolverReport) {
    const { minimum, probes } = this.options;

    let weighted = 0.0;
    let total = 0.0;

    for (const sample of this.samples) {
      const curve = sample.edge.curve as SPowerClothoid;

      if (!(curve instanceof SPowerClothoid) || curve.ws.length === 0) {
        continue;
      }

      const w = sample.weight ?? 1.0;
      const gap = evalSPower(curve.ws, curve.ws.length, sample.u) - sample.width;

      weighted += w * gap * gap;
      total += w;
    }

    report.rms = total > 0.0 ? Math.sqrt(weighted / total) : 0.0;

    for (const edge of this.mesh.edges) {
      const curve = edge.curve as SPowerClothoid;

      if (!(curve instanceof SPowerClothoid) || curve.ws.length === 0) {
        continue;
      }

      let cusp: WidthCusp | undefined;
      let dip: WidthUndershoot | undefined;

      for (let i = 0; i < probes; i++) {
        const u = i / (probes - 1);
        const width = evalSPower(curve.ws, curve.ws.length, u);
        const k = curve.curvature(u * curve.length);
        const product = Math.abs(k) * Math.abs(width) * 0.5;

        if (product >= 1.0 && (cusp === undefined || product > cusp.product)) {
          cusp = { edge, u, s: u * curve.length, width, curvature: k, product };
        }

        if (width < minimum && (dip === undefined || width < dip.width)) {
          dip = { edge, u, width };
        }
      }

      if (cusp !== undefined) {
        report.cusps.push(cusp);
      }

      if (dip !== undefined) {
        report.undershoot.push(dip);
      }
    }
  }
}
