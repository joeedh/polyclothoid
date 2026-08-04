# Architecture

Design document. For the underlying math see `docs/research/clothoids.md`; for port
progress see `docs/plans/typescript-port.md`.

## Shape of the thing

polyclothoid is a **library** that turns a sparse sequence of input points into a smooth,
evenly-spaced sequence of brush dabs. The interesting part is the middle: instead of
fitting cubic beziers, it fits curves whose *curvature profile* is the free variable, then
solves a small constraint system so curvature agrees across segment joints.

```
input points  ->  Mesh (verts + edges)  ->  Solver  ->  Curve.evaluate(s)  ->  dabs
                        ^                     |
                        +---- curvature ------+
                             profiles (ks)
```

The property that motivates the design: because `s` is arclength, the unit tangent is
exactly `(cos θ, sin θ)`. No normalization, no arclength reparameterization table. Evenly
spaced dabs are then just a fixed step in `s`, which is what `Stroker` exploits.

## Layers

Strictly one-directional. Nothing lower may import from something higher.

```
  demo/          canvas harness — DOM lives here and nowhere else
    |
  stroke.ts      Stroker: input points -> evenly spaced dabs
    |
  mesh/          Mesh, Vertex, Edge, Handle, Face; owns the curve per edge
    |
  curve/         Curve base + Clothoid / CubicBezier / BSpline, and their solvers
    |
  math/          Vector2/3/4, cachering, IDGen, Constraint/Solver, binomial, fract
```

| Layer | Files | Lines | Owns |
| --- | --- | --- | --- |
| `src/math` | `vector.ts`, `vec3_mixin.ts`, `util.ts`, `solver.ts` | ~760 | The slice of path.ux the stroker actually used, vendored. |
| `src/curve` | `curve.ts`, `mesh_types.ts`, `clothoid.ts`, `bezier.ts`, `bspline.ts` | ~1090 | Curve primitives and their per-type solvers. Knows no topology beyond "an edge has two vertices". |
| `src/mesh` | `mesh.ts`, `colors.ts` | ~990 | Half-edge topology, selection and flag state, serialization, curve-type dispatch. |
| `src/stroke.ts` | — | ~130 | The only stateful, streaming thing in the library. |
| `src/demo` | `main.ts`, `render.ts`, `ui.ts`, `dabs.ts` | ~570 | Canvas, events, controls. Excluded from the library bundle. |

`math/` has no project imports at all. `demo/` sits outside the stack and consumes the
public barrel exactly as an external caller would, which is the cheapest available check
that the barrel is actually sufficient.

## The two seams

Most of the port's structural work went into two boundaries. Both are load-bearing.

### 1. `curve/mesh_types.ts` — mesh depends on curve, never the reverse

Curve solvers need mesh-shaped things: `ClothoidSolver` walks every edge, finds 2-valence
vertices, and asks each for its neighbour so it can detect corners. That is a genuine
dependency on topology, and importing `Mesh` into `curve/` to express it would make the
cycle real.

Instead `curve/mesh_types.ts` declares what the solvers need **structurally**:

```ts
interface SolvableVertex extends CurveVertex { edges: SolvableEdge[]; otherEdge(e): ... }
interface SolvableEdge<C extends Curve = Curve> { v1; v2; curve: C; length; evaluate(s); ... }
interface SolvableMesh<C extends Curve = Curve> { verts: Iterable<...>; edges: Iterable<...> }
interface CurveSolver { solve(): void }
```

`Mesh` declares `implements SolvableMesh` and satisfies it by construction. The arrow
points one way, mesh → curve, and a solver can be exercised against a hand-built object in
a test without pulling in the mesh at all.

Solver constructors take the **unparameterized** `SolvableMesh` and cast at the point of
use (`e.curve as Clothoid`). Narrowing to `SolvableMesh<Clothoid>` reads better but makes
`Mesh` — whose edges hold `Curve` — fail to fit, since the curve type is a runtime choice.

### 2. `Canvas2DLike` — curves can draw without the library knowing about the DOM

Curves keep a `draw()` hook, because seeing a clothoid's control polygon is how you debug
a bad fit. It takes a structural type declared in `curve/curve.ts`, not
`CanvasRenderingContext2D`:

```ts
export type CanvasPaint = string | object;

export interface Canvas2DLike {
  fillStyle: CanvasPaint; strokeStyle: CanvasPaint; lineWidth: number;
  beginPath(); moveTo(x, y); lineTo(x, y); rect(x, y, w, h); arc(x, y, r, a0, a1);
  fill(); stroke();
}
```

`CanvasPaint` is deliberately wide. TypeScript checks mutable properties covariantly, so
declaring `fillStyle: string` would make a real 2D context — whose `fillStyle` is
`string | CanvasGradient | CanvasPattern` — *fail* to satisfy the interface. Widening the
declared type is what lets the demo pass a real context in and a test pass a stub in.

## Why the curve type is pluggable

`Mesh.switchSplineType(CurveCls, SolverCls)` swaps every edge's curve implementation. This
exists because the project is a comparison harness as much as a library: the whole point is
judging clothoid output against bezier and b-spline output on identical input. Keep the
three implementations interface-compatible even when that costs a little elegance.

```ts
export type CurveConstructor = new (v1?: Vertex, v2?: Vertex) => Curve;
export type CurveSolverConstructor = new (mesh: SolvableMesh) => CurveSolver;

class Mesh {
  CurveCls: CurveConstructor = Clothoid;
  SolverCls: CurveSolverConstructor = ClothoidSolver;
}
```

`makeEdge` news a `CurveCls`, `solve()` news a `SolverCls` and runs it. The original picked
the type with a module-level re-export alias in `spline.js`, so it was fixed at import time
and a saved mesh could only ever round-trip clothoids. It is now a runtime choice.

`Stroker` exposes the same knob through `StrokerOptions`, applied to the throwaway mesh it
fits each segment with:

```ts
new Stroker(cb, { CurveCls: CubicBezier, SolverCls: BezierSolver });
```

The demo runs two strokers over one input stream rather than re-solving one mesh twice.

The shared curve interface is deliberately small:

```ts
evaluate(s)     // point at arclength s
derivative(s)   // unit tangent — analytic for clothoids, since s is arclength
derivative2(s)
curvature(s)
get length()
update(edge)
```

## Solve scheduling

Solving is deferred and flag-driven, not eager.

| Call | Effect |
| --- | --- |
| `regenSolve()` | sets `RecalcFlags.SOLVE`; cheap, call it freely from edits |
| `ensureSolve()` | solves only if the flag is set — what a renderer calls each frame |
| `solve()` | clears the flag, runs the solver unconditionally, returns and keeps its `SolveReport` as `Mesh.report` |

The solver instance is cached on the mesh and reused across solves, rebuilt only when
`SolverCls` changes. That is not just an allocation saving: a solver may carry state between
solves — `SPowerSolver` remembers which joints it broke, so hysteresis needs the same
instance to see the next solve.

Individual curves are lazy the same way: `Curve.update()` marks dirty and recomputation
happens on the next query.

## The public API surface

`src/index.ts` re-exports all four barrels flat:

```ts
export * from "./math/index.js";
export * from "./curve/index.js";
export * from "./mesh/index.js";
export { Stroker, type StrokeCallback, type StrokerOptions } from "./stroke.js";
```

| Group | Exports |
| --- | --- |
| Entry point | `Stroker`, `StrokeCallback`, `StrokerOptions` |
| Mesh | `Mesh`, `Element`, `Vertex`, `Handle`, `Edge`, `Loop`, `LoopList`, `Face`, `ElementArray`, `ElementSet`, `MeshTypes`, `MeshFlags`, `RecalcFlags`, `CurveConstructor`, `CurveSolverConstructor`, `ElemColors`, `getElemColor` |
| Curves | `Curve`, `Clothoid`, `ClothoidSolver`, `CubicBezier`, `BezierSolver`, `BSpline`, `BSplinePoint`, `BSplineSolver` |
| Curve support | `Canvas2DLike`, `CanvasPaint`, `CurveVertex`, `CurveEdge`, `SolvableMesh`, `SolvableEdge`, `SolvableVertex`, `CurveSolver`, `CurvatureProfile`, `activeProfile`, `setCurvatureProfile`, `piecewiseLinear`, `circleArc`, `bernsteinCurvature`, the `K*` parameter-slot indices, `ClothoidSolverOptions` |
| Math | `Vector2/3/4`, `VecLike`, `Vec3Mixin`, `applyVec3Mixin`, `CacheRing`, `IDGen`, `Constraint`, `Solver`, `fract`, `clamp`, `binomial`, `listRemove`, `time_ms` |

Most consumers want `Stroker` and nothing else. The rest is exported because the demo needs
it, and because a curve solver is a reasonable thing to want on its own.

## The one runtime dependency

`nstructjs`, for serialization. Every persistent class carries a `STRUCT` and is registered
with it. This is not a UI concern and not something to vendor — it is a real serialization
format with its own `.d.ts`.

```ts
static override STRUCT = nstructjs.inlineRegister(this, `mesh.Vertex { ... }`);
```

`inlineRegister` supersedes the deprecated `inherit()`: it registers at class-definition
time and walks inheritance itself, so the separate `register()` calls and the repeated
parent fields both go away.

`Edge.curve` serializes as `abstract(Curve)`, which is why `Curve` registers an empty
struct despite carrying no state — it makes the base a known type so a saved mesh
round-trips whichever curve type it was built with.

Everything else the library once took from path.ux now lives in `src/math`, at a fraction
of the size, because only a thin slice was ever used.

## Numerical conventions

These are easy to violate by accident and expensive to debug.

- **`s` is arclength, not a normalized parameter.** `Curve.evaluate(0.5)` means "half a
  unit along", not "halfway". Normalized parameters appear only inside a curve
  implementation, scaled by `KARCSCALE`.
- **Curvature samples live in canonical space.** A `ks[i]` value is only comparable across
  segments after dividing by that segment's `KSCALE`. The solver does this on every
  cross-segment comparison; new constraints must too.
- **Cacherings return borrowed vectors.** `evaluate`, `derivative`, and friends return a
  vector from a fixed-size ring buffer. It is valid until roughly 64-128 further calls, so
  callers that need to keep a value must copy it. This is a deliberate allocation
  optimization for the stroker's inner loop, kept from the original.
- **A constraint returning `0.0` is meaningful.** It short-circuits twice: in
  `Constraint.evaluate`, before any gradient is finite-differenced, and again in
  `Solver.solveStep`, before any parameter is touched. `curv_c` returns `0.0`
  unconditionally and does its work as a **direct projection** — an under-relaxed
  assignment, `fac = 0.5` — interleaved Gauss-Seidel style with `tan_c`'s gradient descent,
  at one call per iteration instead of `2 × paramCount`. Do not "fix" it into returning a
  residual.

## Deliberately not included

- Undo, tool operators, and a datapath/property system. These came from path.ux and served
  the old editor; a stroker library has no business with them.
- Global prototype patching. path.ux defined `Math.fract` and `Array.prototype.remove`
  globally. Both are now explicit imports — see `CLAUDE.md`.
- A second serialization path. The old `toJSON`/`loadJSON` pair duplicated the `STRUCT`
  blocks and was buggier than them; nstructjs is the one we kept.

The vertex mixin exists for the same anti-monkeypatching reason: `Vertex` and `Handle`
extend `Element` and so cannot extend `Vector3`, and the original copied `Array.prototype`
methods onto them with a `Reflect.ownKeys` walk. `applyVec3Mixin` plus declaration merging
replaces that, which is also why `eslint.config.js` disables
`no-unsafe-declaration-merging` and `no-empty-object-type` for `src/mesh/mesh.ts` alone.

## Build

`tools/build.mjs` emits two bundles with esbuild; `tsgo` typechecks and emits declarations
but never emits JS.

| Output | Entry | Notes |
| --- | --- | --- |
| `dist/polyclothoid.js` | `src/index.ts` | `platform: "neutral"`, no externals. Needs explicit `mainFields`/`conditions` because neutral ignores package.json `main`, and nstructjs ships only that field. |
| `dist/demo.js` | `src/demo/main.ts` | `platform: "browser"`; what `index.html` loads. |
| `dist/types/**.d.ts` | — | `pnpm types`, via `tsgo --emitDeclarationOnly`. |
