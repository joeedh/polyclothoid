import { IDGen, Vector3, Vector4, applyVec3Mixin, listRemove, type VecLike, type Vec3Mixin } from "../math/index.js";
import {
  Clothoid,
  ClothoidSolver,
  type Curve,
  type CurveSolver,
  type PairingChannel,
  type SolvableMesh,
  type SolverTuning,
  emptySolveReport,
} from "../curve/index.js";
import * as nstructjs from "nstructjs";

export const MeshTypes = {
  VERTEX  : 1,
  EDGE    : 2,
  HANDLE  : 4,
  LOOP    : 8,
  LOOPLIST: 16,
  FACE    : 32,
};

export const MeshFlags = {
  SELECT: 1,
  HIDE  : 2,
};

export const RecalcFlags = {
  SOLVE: 1,
};

/** Curve types are pluggable; see {@link Mesh.switchSplineType}. */
export type CurveConstructor = new (v1?: Vertex, v2?: Vertex) => Curve;
export type CurveSolverConstructor = new (mesh: SolvableMesh, tuning?: SolverTuning) => CurveSolver;

export class Element {
  static STRUCT = nstructjs.inlineRegister(
    this,
    `
mesh.Element {
  type     : int;
  flag     : int;
  index    : int;
  eid      : int;
}
`
  );

  flag = 0;
  index = 0;
  eid = -1;

  /** Slot in the owning {@link ElementArray}, or -1 once removed. */
  _index = -1;

  constructor(public type: number) {}
}

/**
 * An authored continuity pairing between two of a vertex's edge-ends — `curve/pairing.ts`.
 *
 * Document state, not solve state: it is the ceiling the solver works under, so it is
 * serialized and a load/solve/save round-trip must not change it.
 */
export class VertexPairing {
  static STRUCT = nstructjs.inlineRegister(
    this,
    `
mesh.VertexPairing {
  a       : int | this.a.eid;
  b       : int | this.b.eid;
  level   : int;
  channel : string;
}
`
  );

  a!: Edge;
  b!: Edge;
  level = 0;

  /** Which scalar profile this level speaks about — `curve/mesh_types.ts`. */
  channel: PairingChannel = "curvature";

  constructor(a?: Edge, b?: Edge, level = 0, channel: PairingChannel = "curvature") {
    if (a && b) {
      this.a = a;
      this.b = b;
    }

    this.level = level;
    this.channel = channel;
  }
}

export class Vertex extends Element {
  static override STRUCT = nstructjs.inlineRegister(
    this,
    `
mesh.Vertex {
  0           : float;
  1           : float;
  2           : float;
  edges       : array(e, int) | e.eid;
  pairings    : array(mesh.VertexPairing);
}
`
  );

  edges: Edge[] = [];

  /** Empty means unspecified, which `curve/pairing.ts` reads as the valence default. */
  pairings: VertexPairing[] = [];

  constructor(co?: VecLike) {
    super(MeshTypes.VERTEX);

    this.initVector3();

    if (co !== undefined) {
      this.load(co);
    }
  }

  /** The authored pairing of `a` and `b` here on one channel, in either order. */
  pairing(a: Edge, b: Edge, channel: PairingChannel = "curvature") {
    return this.pairings.find((p) => p.channel === channel && ((p.a === a && p.b === b) || (p.a === b && p.b === a)));
  }

  /**
   * Author a continuity level between two edge-ends. This is where corner detection lands.
   *
   * A level at or above the ceiling is stored as authored and clamped when it is read, so
   * raising `p` later restores the continuity that was asked for rather than the one that
   * happened to fit.
   *
   * The channels are independent: a level-0 curvature pairing is a corner in the curve and
   * says nothing about whether the width steps there.
   */
  setPairing(a: Edge, b: Edge, level: number, channel: PairingChannel = "curvature") {
    const existing = this.pairing(a, b, channel);

    if (existing) {
      existing.level = level;
    } else {
      this.pairings.push(new VertexPairing(a, b, level, channel));
    }

    return this;
  }

  /** Only meaningful on 2-valence vertices, which is what the curve solvers walk. */
  otherEdge(e: Edge) {
    if (this.edges.length !== 2) {
      throw new Error("otherEdge only works on 2-valence vertices");
    }

    if (e === this.edges[0]) {
      return this.edges[1];
    }

    if (e === this.edges[1]) {
      return this.edges[0];
    }

    return undefined;
  }
}

export interface Vertex extends Vec3Mixin {}
applyVec3Mixin(Vertex);

/** A bezier-style control point owned by an {@link Edge}. Display only. */
export class Handle extends Element {
  static override STRUCT = nstructjs.inlineRegister(
    this,
    `
mesh.Handle {
  0           : float;
  1           : float;
  2           : float;
  owner       : int | this.owner ? this.owner.eid : -1;
}
`
  );

  owner?: Edge;

  constructor(co?: VecLike) {
    super(MeshTypes.HANDLE);

    this.initVector3();

    if (co !== undefined) {
      this.load(co);
    }
  }
}

export interface Handle extends Vec3Mixin {}
applyVec3Mixin(Handle);

/** Shrinks the domain away from the endpoints, where curvature is least reliable. */
const DERIVATIVE_EPS = 0.00005;

export class Edge extends Element {
  static override STRUCT = nstructjs.inlineRegister(
    this,
    `
mesh.Edge {
  v1           : int | this.v1.eid;
  v2           : int | this.v2.eid;
  h1           : int | this.h1.eid;
  h2           : int | this.h2.eid;
  l            : int | this.l ? this.l.eid : -1;
  curve        : abstract(Curve);
}
`
  );

  v1!: Vertex;
  v2!: Vertex;
  h1!: Handle;
  h2!: Handle;
  curve!: Curve;

  /** Head of the radial cycle of loops using this edge. */
  l?: Loop;

  constructor() {
    super(MeshTypes.EDGE);
  }

  get length() {
    return this.curve.length;
  }

  /** Walks the radial cycle. Guarded, because a corrupt cycle would never terminate. */
  get loops() {
    const start = this.l;

    return (function* () {
      if (!start) {
        return;
      }

      let l = start;

      for (let i = 0; i < 100; i++) {
        yield l;

        l = l.radial_next!;

        if (l === start) {
          return;
        }
      }

      console.warn("Infinite loop detected in Edge.loops");
    })();
  }

  update() {
    this.curve.update(this);

    return this;
  }

  evaluate(s: number) {
    return this.curve.evaluate(s);
  }

  derivative(s: number) {
    return this.curve.derivative(s * (1.0 - DERIVATIVE_EPS * 2.0) + DERIVATIVE_EPS);
  }

  derivative2(s: number) {
    return this.curve.derivative2(s * (1.0 - DERIVATIVE_EPS * 2.0) + DERIVATIVE_EPS);
  }

  normal(s: number) {
    return this.curve.normal(s);
  }

  curvature(s: number) {
    return this.curve.curvature(s);
  }

  /** The handle at `v`'s end of the edge. */
  handle(v: Vertex) {
    return v === this.v1 ? this.h1 : this.h2;
  }

  /** Inverse of {@link handle}. */
  vertex(h: Handle) {
    return h === this.h1 ? this.v1 : this.v2;
  }

  otherVertex(v: Vertex) {
    if (v === this.v1) {
      return this.v2;
    }

    if (v === this.v2) {
      return this.v1;
    }

    throw new Error("vertex " + v.eid + " not in edge");
  }
}

export class Loop extends Element {
  static override STRUCT = nstructjs.inlineRegister(
    this,
    `
mesh.Loop {
  v : int | this.v.eid;
  e : int | this.e.eid;
}
`
  );

  v!: Vertex;
  e!: Edge;
  f!: Face;
  list!: LoopList;

  radial_next?: Loop;
  radial_prev?: Loop;
  next!: Loop;
  prev!: Loop;

  constructor() {
    super(MeshTypes.LOOP);
  }
}

/** One boundary cycle of a {@link Face} — the outer contour or a hole. */
export class LoopList extends Element {
  static override STRUCT = nstructjs.inlineRegister(
    this,
    `
mesh.LoopList {
  _loops : iter(int) | this._save_loops();
}
`
  );

  length = 0;
  l?: Loop;

  /** Populated only while serializing; holds loop eids. */
  _loops?: number[];

  constructor() {
    super(MeshTypes.LOOPLIST);
  }

  get verts() {
    return this.iterVerts();
  }

  *iterVerts() {
    for (const l of this) {
      yield l.v;
    }
  }

  *[Symbol.iterator]() {
    const start = this.l;

    if (!start) {
      return;
    }

    let l = start;

    for (let i = 0; i < 10000; i++) {
      yield l;

      l = l.next;

      if (l === start) {
        return;
      }
    }

    console.warn("Infinite loop detected in LoopList");
  }

  _save_loops() {
    return [...this].map((l) => l.eid);
  }
}

export class Face extends Element {
  static override STRUCT = nstructjs.inlineRegister(
    this,
    `
mesh.Face {
  lists     : iter(list, int) | list.eid;
  fillColor : vec4;
  blur      : float;
}
`
  );

  lists: LoopList[] = [];
  blur = 0.0;
  center = new Vector3();
  fillColor = new Vector4([0.5, 0.5, 0.5, 1]);

  constructor() {
    super(MeshTypes.FACE);
  }

  get loops() {
    return this.iterLoops();
  }

  *iterLoops() {
    for (const list of this.lists) {
      yield* list;
    }
  }

  get verts() {
    return this.iterVerts();
  }

  *iterVerts() {
    for (const l of this.loops) {
      yield l.v;
    }
  }

  calcCenter() {
    this.center.zero();

    let tot = 0;

    for (const l of this.loops) {
      this.center.add(l.v);
      tot++;
    }

    if (tot) {
      this.center.mulScalar(1.0 / tot);
    }

    return this.center;
  }
}

export class ElementSet<T extends Element> extends Set<T> {
  constructor(public type: number) {
    super();
  }

  get editable() {
    return this.iterVisible();
  }

  *iterVisible() {
    for (const item of this) {
      if (!(item.flag & MeshFlags.HIDE)) {
        yield item;
      }
    }
  }

  get length() {
    return this.size;
  }

  remove(item: T) {
    this.delete(item);
  }
}

/**
 * A sparse, eid-stable list of one element type.
 *
 * Removal punches a hole rather than compacting, so `_index` stays valid for every other
 * element. `length` is the live count, not `list.length`.
 */
export class ElementArray<T extends Element = Element> {
  static STRUCT = nstructjs.inlineRegister(
    this,
    `
mesh.ElementArray {
  this        : iter(abstract(mesh.Element));
  highlight   : int | this.highlight ? this.highlight.eid : -1;
  active      : int | this.active ? this.active.eid : -1;
  type        : int;
}
`
  );

  list: (T | undefined)[] = [];
  length = 0;
  freelist: number[] = [];

  selected: ElementSet<T>;
  active?: T;
  highlight?: T;

  constructor(public type: number) {
    this.selected = new ElementSet<T>(type);
  }

  get visible() {
    return this.iterVisible();
  }

  *iterVisible() {
    for (const item of this) {
      if (!(item.flag & MeshFlags.HIDE)) {
        yield item;
      }
    }
  }

  get editable() {
    return this.visible;
  }

  *[Symbol.iterator]() {
    for (const item of this.list) {
      if (item !== undefined) {
        yield item;
      }
    }
  }

  concat(b: Iterable<T>) {
    return [...this, ...b];
  }

  push(v: T) {
    v._index = this.list.length;
    this.list.push(v);
    this.length++;

    if (v.flag & MeshFlags.SELECT) {
      this.selected.add(v);
    }

    return this;
  }

  remove(v: T) {
    this.selected.delete(v);

    if (this.active === v) {
      this.active = undefined;
    }

    if (this.highlight === v) {
      this.highlight = undefined;
    }

    if (v._index < 0 || this.list[v._index] !== v) {
      throw new Error("element not in array");
    }

    this.freelist.push(v._index);
    this.list[v._index] = undefined;

    v._index = -1;
    this.length--;

    return this;
  }

  setSelect(v: T, state: boolean) {
    if (state) {
      v.flag |= MeshFlags.SELECT;
      this.selected.add(v);
    } else {
      v.flag &= ~MeshFlags.SELECT;
      this.selected.delete(v);
    }

    return this;
  }

  selectNone() {
    for (const e of this) {
      this.setSelect(e, false);
    }
  }

  selectAll() {
    for (const e of this) {
      this.setSelect(e, true);
    }
  }

  loadSTRUCT(reader: (obj: this) => void) {
    reader(this);

    for (const elem of this) {
      if (elem.flag & MeshFlags.SELECT) {
        this.selected.add(elem);
      }
    }
  }
}

/**
 * A half-edge mesh whose edges carry pluggable {@link Curve} geometry.
 *
 * The face/loop machinery is inherited from the original editor and is not exercised by
 * the stroker path, which only ever builds vertex-edge chains.
 */
export class Mesh implements SolvableMesh {
  static STRUCT = nstructjs.inlineRegister(
    this,
    `
mesh.Mesh {
  elists : array(mesh.ElementArray) | this.getElists();
  eidgen : IDGen;
}
`
  );

  eidgen = new IDGen();
  eidMap = new Map<number, Element>();

  recalc = RecalcFlags.SOLVE;

  /** What the last {@link solve} reported. Not serialized — it describes a run, not a document. */
  report = emptySolveReport();

  CurveCls: CurveConstructor = Clothoid;
  SolverCls: CurveSolverConstructor = ClothoidSolver;

  /** The live {@link SolverCls} instance, kept for its cross-solve state — see {@link solve}. */
  solver?: CurveSolver;

  /**
   * Handed to {@link SolverCls} when it is built. Change it through {@link setSolverTuning}.
   *
   * Not serialized: it describes how this session wants curves fitted, not the document.
   */
  solverTuning: SolverTuning = {};

  elists = new Map<number, ElementArray>();

  verts!: ElementArray<Vertex>;
  handles!: ElementArray<Handle>;
  edges!: ElementArray<Edge>;
  loops!: ElementArray<Loop>;
  lists!: ElementArray<LoopList>;
  faces!: ElementArray<Face>;

  constructor() {
    this.makeElists();
  }

  get elements() {
    return this.eidMap.values();
  }

  get hasHighlight() {
    for (const elist of this.elists.values()) {
      if (elist.highlight) {
        return true;
      }
    }

    return false;
  }

  getElists() {
    return [...this.elists.values()];
  }

  makeElists() {
    this.elists.clear();

    for (const type of Object.values(MeshTypes)) {
      this.elists.set(type, new ElementArray(type));
    }

    this.addElistAliases();
  }

  addElistAliases() {
    this.verts = this.elists.get(MeshTypes.VERTEX) as ElementArray<Vertex>;
    this.handles = this.elists.get(MeshTypes.HANDLE) as ElementArray<Handle>;
    this.edges = this.elists.get(MeshTypes.EDGE) as ElementArray<Edge>;
    this.loops = this.elists.get(MeshTypes.LOOP) as ElementArray<Loop>;
    this.lists = this.elists.get(MeshTypes.LOOPLIST) as ElementArray<LoopList>;
    this.faces = this.elists.get(MeshTypes.FACE) as ElementArray<Face>;
  }

  ensureSolve() {
    if (this.recalc & RecalcFlags.SOLVE) {
      this.solve();
    }
  }

  regenSolve() {
    this.recalc |= RecalcFlags.SOLVE;

    return this;
  }

  /**
   * Run {@link SolverCls} over this mesh and hand back its §8 report.
   *
   * Also kept on {@link report}, because {@link ensureSolve} is the path most callers take
   * and it has nowhere to return a value to.
   *
   * The solver instance is kept rather than rebuilt, because §8's hysteresis is state that
   * has to outlive one solve to mean anything — a joint that lost a continuity level has to
   * remember it in order to be held to the stricter threshold when it asks for it back. It is
   * dropped whenever {@link SolverCls} changes, since the levels a solver withheld say nothing
   * about a different curve type.
   */
  solve() {
    this.recalc &= ~RecalcFlags.SOLVE;

    if (!this.solver || !(this.solver instanceof this.SolverCls)) {
      this.solver = new this.SolverCls(this, this.solverTuning);
    }

    this.report = this.solver.solve();

    return this.report;
  }

  /**
   * Reconfigure {@link SolverCls}, and re-solve if that actually changed anything.
   *
   * A solver reads its options once, when it is constructed, so a change means a new instance
   * and the loss of §8's hysteresis along with it. That is the right trade when the bound
   * really moved and the wrong one on every unrelated UI event, which is what the comparison
   * is for — a control panel can call this on every change without paying for it.
   */
  setSolverTuning(tuning: SolverTuning) {
    const before = this.solverTuning as Record<string, unknown>;
    const after = tuning as Record<string, unknown>;
    const keys = new Set([...Object.keys(before), ...Object.keys(after)]);

    if ([...keys].every((k) => before[k] === after[k])) {
      return this;
    }

    this.solverTuning = { ...tuning };
    this.solver = undefined;

    return this.regenSolve();
  }

  /** Replace every edge's curve with a fresh instance of `CurveCls` and re-solve. */
  switchSplineType(CurveCls: CurveConstructor, SolverCls: CurveSolverConstructor) {
    this.CurveCls = CurveCls;
    this.SolverCls = SolverCls;
    this.solver = undefined;

    for (const e of this.edges) {
      e.curve = new CurveCls(e.v1, e.v2);
      e.curve.init(e);
    }

    return this.regenSolve();
  }

  _element_init(e: Element) {
    e.eid = this.eidgen.next();
    this.eidMap.set(e.eid, e);
  }

  getList(type: number) {
    return this.elists.get(type)!;
  }

  setSelect(e: Element, state: boolean) {
    this.getList(e.type).setSelect(e, state);
  }

  selectNone() {
    for (const elist of this.elists.values()) {
      elist.selectNone();
    }
  }

  selectAll() {
    for (const elist of this.elists.values()) {
      elist.selectAll();
    }
  }

  setActive(elem?: Element) {
    if (!elem) {
      for (const elist of this.elists.values()) {
        elist.active = undefined;
      }
    } else {
      this.getList(elem.type).active = elem;
    }

    return this;
  }

  /** Returns whether the highlight actually changed, so callers can skip a redraw. */
  setHighlight(elem?: Element) {
    if (!elem) {
      return this.clearHighlight();
    }

    const elist = this.getList(elem.type);
    const changed = elist.highlight !== elem;

    elist.highlight = elem;

    return changed;
  }

  clearHighlight() {
    const existed = this.hasHighlight;

    for (const elist of this.elists.values()) {
      elist.highlight = undefined;
    }

    return existed;
  }

  makeVertex(co?: VecLike) {
    const v = new Vertex(co);

    this._element_init(v);
    this.verts.push(v);

    return v;
  }

  makeHandle(co?: VecLike) {
    const h = new Handle(co);

    this._element_init(h);
    this.handles.push(h);

    return h;
  }

  getEdge(v1: Vertex, v2: Vertex) {
    for (const e of v1.edges) {
      if (e.otherVertex(v1) === v2) {
        return e;
      }
    }

    return undefined;
  }

  makeEdge(v1: Vertex, v2: Vertex) {
    const e = new Edge();

    e.v1 = v1;
    e.v2 = v2;

    e.h1 = this.makeHandle(v1);
    e.h1.interp(v2, 1.0 / 3.0);
    e.h1.owner = e;

    e.h2 = this.makeHandle(v1);
    e.h2.interp(v2, 2.0 / 3.0);
    e.h2.owner = e;

    e.curve = new this.CurveCls(v1, v2);
    e.curve.init(e);

    v1.edges.push(e);
    v2.edges.push(e);

    this._element_init(e);
    this.edges.push(e);

    return e;
  }

  reverseEdge(e: Edge) {
    [e.v1, e.v2] = [e.v2, e.v1];
    [e.h1, e.h2] = [e.h2, e.h1];
  }

  killVertex(v: Vertex) {
    if (v.eid === -1) {
      console.warn("Vertex", v.eid, "already freed");
      return;
    }

    while (v.edges.length > 0) {
      this.killEdge(v.edges[0]);
    }

    this.eidMap.delete(v.eid);
    this.verts.remove(v);

    v.eid = -1;
  }

  killEdge(e: Edge) {
    if (e.eid === -1) {
      console.warn("Edge", e.eid, "already freed");
      return;
    }

    while (e.l) {
      this.killFace(e.l.f);
    }

    this.edges.remove(e);
    this.eidMap.delete(e.eid);

    this.eidMap.delete(e.h1.eid);
    this.handles.remove(e.h1);

    this.eidMap.delete(e.h2.eid);
    this.handles.remove(e.h2);

    e.eid = -1;

    listRemove(e.v1.edges, e);
    listRemove(e.v2.edges, e);

    for (const v of [e.v1, e.v2]) {
      v.pairings = v.pairings.filter((pairing) => pairing.a !== e && pairing.b !== e);
    }
  }

  _killLoop(l: Loop) {
    this.eidMap.delete(l.eid);
    this.loops.remove(l);

    l.eid = -1;
  }

  _killList(list: LoopList) {
    this.eidMap.delete(list.eid);
    this.lists.remove(list);

    list.eid = -1;
  }

  killFace(f: Face) {
    for (const list of f.lists) {
      for (const l of list) {
        this.radialLoopRemove(l.e, l);
        this._killLoop(l);
      }

      this._killList(list);
    }

    this.eidMap.delete(f.eid);
    this.faces.remove(f);

    f.eid = -1;
  }

  radialLoopInsert(e: Edge, l: Loop) {
    if (!e.l) {
      e.l = l;
      l.radial_next = l.radial_prev = l;

      return;
    }

    l.radial_prev = e.l;
    l.radial_next = e.l.radial_next;

    e.l.radial_next!.radial_prev = l;
    e.l.radial_next = l;
  }

  radialLoopRemove(e: Edge, l: Loop) {
    if (e.l === l) {
      e.l = e.l.radial_next;
    }

    if (e.l === l) {
      e.l = undefined;
      return;
    }

    l.radial_next!.radial_prev = l.radial_prev;
    l.radial_prev!.radial_next = l.radial_next;
  }

  addLoopList(f: Face, vs: Vertex[]) {
    const list = new LoopList();

    this._element_init(list);
    this.lists.push(list);

    let firstl: Loop | undefined;
    let lastl: Loop | undefined;

    for (let i = 0; i < vs.length; i++) {
      const v1 = vs[i];
      const v2 = vs[(i + 1) % vs.length];

      const e = this.getEdge(v1, v2) ?? this.makeEdge(v1, v2);
      const l = new Loop();

      this._element_init(l);
      this.loops.push(l);

      l.v = v1;
      l.e = e;
      l.f = f;
      l.list = list;

      this.radialLoopInsert(e, l);

      if (!firstl) {
        firstl = l;
      } else {
        lastl!.next = l;
        l.prev = lastl!;
      }

      lastl = l;
      list.length++;
    }

    if (firstl && lastl) {
      firstl.prev = lastl;
      lastl.next = firstl;
      list.l = firstl;
    }

    f.lists.push(list);

    return list;
  }

  makeFace(vs: Vertex[]) {
    const f = new Face();

    this._element_init(f);
    this.faces.push(f);

    this.addLoopList(f, vs);

    return f;
  }

  unlinkFace(f: Face) {
    for (const list of f.lists) {
      for (const l of list) {
        this.radialLoopRemove(l.e, l);
      }
    }
  }

  linkFace(f: Face, forceRelink = true) {
    for (const list of f.lists) {
      for (const l of list) {
        if (forceRelink || !l.e) {
          l.e = this.getEdge(l.v, l.next.v) ?? this.makeEdge(l.v, l.next.v);
        }

        this.radialLoopInsert(l.e, l);
      }
    }
  }

  /**
   * Insert a vertex at `t` along `e`, splitting it into `[newEdge, newVertex]`.
   *
   * Loops on both sides of the edge are duplicated so every face touching it keeps a
   * closed cycle.
   */
  splitEdge(e: Edge, t = 0.5): [Edge, Vertex] {
    const nv = this.makeVertex(e.v1).interp(e.v2, t);
    const ne = this.makeEdge(nv, e.v2);

    listRemove(e.v2.edges, e);
    e.v2 = nv;
    nv.edges.push(e);

    ne.h1.load(nv).interp(ne.v2, 1.0 / 3.0);
    ne.h2.load(nv).interp(ne.v2, 2.0 / 3.0);

    e.h1.load(e.v1).interp(nv, 1.0 / 3.0);
    e.h2.load(e.v1).interp(nv, 2.0 / 3.0);

    e.update();
    ne.update();

    if (e.flag & MeshFlags.SELECT) {
      this.edges.setSelect(ne, true);
      this.verts.setSelect(nv, true);
    }

    if (e.l) {
      for (const l of [...e.loops]) {
        const l2 = new Loop();

        this._element_init(l2);
        this.loops.push(l2);

        l2.f = l.f;
        l2.list = l.list;
        l2.v = nv;

        if (l.e === e) {
          l2.e = ne;

          this.radialLoopInsert(ne, l2);
        } else {
          /* This loop runs against the edge, so the new half is the original edge. */
          this.radialLoopRemove(e, l);

          l.e = ne;
          l2.e = e;

          this.radialLoopInsert(ne, l);
          this.radialLoopInsert(e, l2);
        }

        l2.prev = l;
        l2.next = l.next;
        l.next.prev = l2;
        l.next = l2;

        l.list.length++;
      }
    }

    return [ne, nv];
  }

  /** Merge the two edges either side of a 2-valence vertex, removing it. */
  dissolveVertex(v: Vertex) {
    if (v.edges.length !== 2) {
      throw new Error("can't dissolve vertex without exactly two edges");
    }

    const loops = new Set<Loop>();
    const faces = new Set<Face>();

    for (const e of v.edges) {
      for (const l of e.loops) {
        const lv = l.v === v ? l : l.next;

        loops.add(lv);
        faces.add(lv.f);
      }
    }

    for (const f of faces) {
      this.unlinkFace(f);
    }

    for (const l of loops) {
      l.prev.next = l.next;
      l.next.prev = l.prev;

      if (l === l.list.l) {
        l.list.l = l.list.l.next;
      }

      if (l === l.list.l) {
        /* The cycle collapsed to nothing; the whole boundary goes with it. */
        listRemove(l.f.lists, l.list);
        this._killList(l.list);

        if (l.f.lists.length === 0) {
          faces.delete(l.f);
          this.killFace(l.f);
          continue;
        }
      } else {
        l.list.length--;
      }

      this._killLoop(l);
    }

    const [e1, e2] = v.edges;
    const v1 = e1.otherVertex(v);
    const v2 = e2.otherVertex(v);

    const flag = (e1.flag | e2.flag) & ~MeshFlags.HIDE;

    this.killVertex(v);

    const e3 = this.makeEdge(v1, v2);
    e3.flag |= flag;

    if (flag & MeshFlags.SELECT) {
      this.edges.setSelect(e3, true);
    }

    for (const f of faces) {
      this.linkFace(f, true);
    }

    return e3;
  }

  selectFlush(selmode: number) {
    if (selmode & MeshTypes.VERTEX) {
      this.edges.selectNone();

      const active = this.edges.active;
      const setActive = !active || !((active.v1.flag | active.v2.flag) & MeshFlags.SELECT);

      for (const e of this.edges) {
        if (!(e.v1.flag & MeshFlags.SELECT) || !(e.v2.flag & MeshFlags.SELECT)) {
          continue;
        }

        this.edges.setSelect(e, true);
        this.handles.setSelect(e.h1, true);
        this.handles.setSelect(e.h2, true);

        if (setActive) {
          this.edges.active = e;
        }
      }

      for (const f of this.faces) {
        let ok = true;

        for (const l of f.loops) {
          if (!(l.e.flag & MeshFlags.SELECT)) {
            ok = false;
            break;
          }
        }

        if (ok) {
          this.faces.setSelect(f, true);
        }
      }
    } else if (selmode & MeshTypes.EDGE) {
      this.verts.selectNone();

      for (const v of this.verts) {
        if (v.edges.some((e) => e.flag & MeshFlags.SELECT)) {
          this.verts.setSelect(v, true);
        }
      }
    }
  }

  /**
   * Rebuild every object reference from the eids nstructjs wrote out.
   *
   * Fields typed `int` in the STRUCT blocks come back as eids, so each pass below swaps
   * an eid for the element it names. The radial cycles are rebuilt last because they
   * depend on loops, faces and lists all being resolved first.
   */
  loadSTRUCT(reader: (obj: this) => void) {
    reader(this);

    const elists = this.elists as unknown as ElementArray[];

    this.elists = new Map();
    for (const elist of elists) {
      this.elists.set(elist.type, elist);
    }

    this.addElistAliases();

    const eidMap = (this.eidMap = new Map<number, Element>());

    for (const elist of this.getElists()) {
      for (const elem of elist) {
        eidMap.set(elem.eid, elem);
      }
    }

    const get = <T extends Element>(eid: unknown) => eidMap.get(eid as number) as T | undefined;

    for (const v of this.verts) {
      const eids = v.edges as unknown as number[];

      v.edges = eids.map((eid) => get<Edge>(eid)!).filter(Boolean);

      for (const pairing of v.pairings) {
        pairing.a = get<Edge>(pairing.a)!;
        pairing.b = get<Edge>(pairing.b)!;

        /* Files written before §10's width channel have no field here at all. */
        pairing.channel = pairing.channel || "curvature";
      }

      /* A pairing naming an edge that did not survive is not a corner, it is a dangling ref. */
      v.pairings = v.pairings.filter((pairing) => pairing.a && pairing.b);
    }

    for (const h of this.handles) {
      h.owner = get<Edge>(h.owner);
    }

    /* Edge.l is stashed and reapplied after the radial cycles are rebuilt below. */
    const edgeLoops = new Map<Edge, Loop | undefined>();

    for (const e of this.edges) {
      e.v1 = get<Vertex>(e.v1)!;
      e.v2 = get<Vertex>(e.v2)!;
      e.h1 = get<Handle>(e.h1)!;
      e.h2 = get<Handle>(e.h2)!;

      edgeLoops.set(e, get<Loop>(e.l));
      e.l = undefined;
    }

    for (const l of this.loops) {
      l.v = get<Vertex>(l.v)!;
      l.e = get<Edge>(l.e)!;
    }

    for (const list of this.lists) {
      const loops = (list._loops ?? []).map((eid) => get<Loop>(eid)!).filter(Boolean);

      list._loops = undefined;
      list.l = loops[0];

      for (let i = 0; i < loops.length; i++) {
        const l1 = loops[i];
        const l2 = loops[(i + 1) % loops.length];

        l1.next = l2;
        l2.prev = l1;
      }
    }

    for (const f of this.faces) {
      f.lists = f.lists.map((eid) => get<LoopList>(eid)!).filter(Boolean);

      for (const list of f.lists) {
        list.length = 0;

        for (const l of list) {
          l.f = f;
          l.list = list;

          this.radialLoopInsert(l.e, l);
          list.length++;
        }
      }
    }

    for (const [e, l] of edgeLoops) {
      e.l = l;
    }

    for (const elist of this.getElists()) {
      elist.active = get(elist.active);
      elist.highlight = get(elist.highlight);
    }

    for (const e of this.edges) {
      if (e.curve) {
        e.curve.update(e);
      } else {
        e.curve = new this.CurveCls(e.v1, e.v2);
        e.curve.init(e);
      }
    }

    this.regenSolve();
  }
}
