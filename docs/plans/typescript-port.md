# Master plan — TypeScript port and de-path.ux-ing

Status legend: `[ ]` not started · `[~]` in progress · `[x]` done · `[-]` dropped

Living document. Update as work lands, not at the end.

---

## Goal

Turn `scripts/stroker` from a JS module embedded in a path.ux application into a
standalone, dependency-light TypeScript **library**, with a minimal vanilla-TS canvas
harness for visual verification.

## Decisions taken

| Decision | Choice | Rationale |
| --- | --- | --- |
| path.ux submodule | **Delete** | The library used only `Vector2/3/4`, `cachering`, `IDGen`, `Constraint`/`Solver`, `binomial`, `Math.fract`. All vendored in well under 500 lines. A UI toolkit is the wrong dependency for a stroker library. |
| Demo app (`scripts/core`, ~2500 lines) | **Rebuild minimal, vanilla** | Confirmed with the user. Drops `AppState`, `ToolOp`, undo stack, datapath `Context`, property templates, screen/sarea layout. Keeps canvas draw, vertex dragging, toggles, stroke preview. |
| `nstructjs` | **Keep as an npm dep** | Not a UI library — it is the serialization layer behind every `STRUCT` block in `mesh.ts` and the curve classes. Ships its own `.d.ts`. Dropping it would mean deleting save/load, which is real non-UI functionality. This is the library's only runtime dependency. |
| Typechecker | `tsgo` (`@typescript/native-preview`) | Per project direction. |
| Bundler / server | esbuild + esbuild's own HTTP server | Per project direction. Replaces rollup and the hand-rolled `serv.js`. |
| Formatter | `@pathtx/prettier` | Fork providing `alignPropertyValues: "group"`; matches existing house style. Config lifted from path.ux's own `.prettierrc`. |
| Port method | Two-pass | Pass 1 annotates from reasoning alone with the typechecker never invoked; pass 2 drives errors to zero. Keeps pass 1 honest about intent rather than letting the checker drive the design. |

## Phases

### 1. Toolchain scaffold — `[x]`

- [x] `package.json` on pnpm, ESM, scripts for build/serve/typecheck/format/lint
- [x] `tsconfig.json` — strict, `verbatimModuleSyntax`, `noEmit` (esbuild emits)
- [x] `tools/build.mjs`, `tools/serve.mjs` (esbuild)
- [x] `.prettierrc` from path.ux house style
- [x] `eslint.config.js` — flat config, typescript-eslint, `no-inferrable-types`,
      and `no-warning-comments` set to fail on any surviving `CLAUDENOTE:`
- [x] Verify all four tools actually run

### 2. Docs and conventions — `[x]`

- [x] `CLAUDE.md` — toolchain, TS style, comment rules, doc layout
- [x] `docs/plans/typescript-port.md` (this file)
- [x] `docs/architecture.md` — module boundaries and the public API surface
- [x] `docs/research/clothoids.md` — pre-existing, still accurate

### 3. Vendored math layer — `[x]`

Replaces the entire path.ux surface. Target `src/math/`. ~700 lines total, versus a UI
toolkit submodule.

- [x] `vector.ts` — `Vector2`/`Vector3`/`Vector4`, hand-written and unrolled per
      component. The original generated these with `eval` at module load.
- [x] `vec3_mixin.ts` — vector behaviour for `Vertex`/`Handle`, which extend `Element`
      and so cannot extend `Vector3`. Replaces a `Reflect.ownKeys` prototype walk that
      copied `Array.prototype` methods onto non-array objects. Typed via declaration
      merging, so `this`-returning methods stay correctly typed on the subclass.
- [x] `util.ts` — `CacheRing` (generic, holds storage in a field rather than
      `extends Array`), `IDGen`, `fract`, `clamp`, `time_ms`, `listRemove`, `binomial`
- [x] `solver.ts` — `Constraint`/`Solver`, generic over the constraint param type
- [x] `index.ts` barrel

**Do not break:** a constraint returning `0.0` short-circuits both the gradient
finite-differencing (`Constraint.evaluate`) and the parameter update (`Solver.solveStep`).
`curv_c` in the clothoid solver depends on this to act as a direct projection interleaved
with gradient descent. It is a designed feature, not an accident.

### 4. Library port, pass 1 — annotate by reasoning — `[x]`

Typechecker **must not be run** during this phase.

- [x] `curve/curve.ts` — `Curve` base, `Canvas2DLike`, `CurveVertex`/`CurveEdge`
- [x] `curve/mesh_types.ts` — structural `Solvable*` interfaces, breaking the mesh/curve cycle
- [x] `curve/clothoid.ts` (613 lines — the core; exempt from the 3-line comment cap)
- [x] `curve/bezier.ts` (272 lines)
- [x] `curve/bspline.ts` (507 lines)
- [x] `curve/spline.ts` — folded into `curve/index.ts`; the type is now chosen at runtime
      via `Mesh.switchSplineType` rather than by a module-level re-export alias
- [x] `mesh/mesh.ts` (1731 lines — the big one)
- [x] `mesh/colors.ts` — `ElemColors`/`getElemColor` split out of `mesh.ts`
- [x] `stroke.ts` (182 lines)
- [x] `index.ts` — public API surface

### 5. Library port, pass 2 — drive tsgo to zero — `[x]`

- [x] `pnpm typecheck` clean
- [x] Every `CLAUDENOTE:` from pass 1 resolved or escalated to the user — none survived

Four clusters of errors, all structural rather than incidental:

| Cluster | Fix |
| --- | --- |
| `nstructjs` has no default export (5 files) | `import * as nstructjs` |
| `Index signature for type 'number' is missing` on every mixed-in vertex (~30) | added `[index: number]: number` to `Vec3Mixin`, which is what makes a mixed-in class an `ArrayLike` |
| `nstructjs.inherit()` is deprecated; `static STRUCT` needs `override` (9 classes) | `static override STRUCT = nstructjs.inlineRegister(this, ...)`, which registers inline and walks inheritance itself, so the separate `register()` calls and the repeated parent fields both go away |
| Solver constructors were narrowed to `SolvableMesh<Clothoid>` etc., so `Mesh` did not fit | widened to `SolvableMesh`, with the cast at the point of use |

### 6. Demo harness — `[x]`

- [x] Canvas renderer: mesh, curves, curvature combs, faces, vertices
- [x] Vertex dragging / selection / highlight, plus split, dissolve and delete keys
- [x] Brush stroke preview driven by `Stroker`, clothoid and bezier side by side
- [x] Plain HTML controls replacing the path.ux property panel, persisted to localStorage
- [x] `index.html` pointing at the esbuild bundle

`src/demo/` is four files, ~600 lines including the control panel: `main.ts` (boot and
events), `render.ts`, `ui.ts`, `dabs.ts`. It is excluded from the library bundle.

To compare the two fits, `Stroker` gained an optional `StrokerOptions` second argument
that overrides the curve type of its throwaway mesh. The demo runs two strokers over the
same input rather than re-solving one mesh twice, which is what the original did.

### 7. Cleanup and gates — `[x]`

- [x] path.ux submodule removed: `deinit`, unstaged, `.gitmodules` deleted
- [x] Old JS sources deleted — `scripts/core`, `scripts/stroker`, `rollup.config.js`,
      `serv.js`, `build.sh`, `git_pull.sh`
- [x] `pnpm format`
- [x] `pnpm lint` clean — which by construction means zero `CLAUDENOTE:` left
- [x] `pnpm build` and `pnpm serve` verified end to end
- [x] `Readme.Md` rewritten

Lint needed four real fixes (`prefer-const`, a type-only import, two empty methods that
now say why they are empty, and five `this` aliases replaced by named generator methods)
plus one documented per-file override: `no-unsafe-declaration-merging` and
`no-empty-object-type` are off for `src/mesh/mesh.ts`, because the `Vertex`/`Handle`
vector mixin is typed by declaration merging and there is no alternative spelling.

---

## Known defects to carry across the port

From `docs/research/clothoids.md` §8. The port should not silently preserve these, and
should not silently fix them either — each needs a deliberate call.

| Defect | Where | Disposition |
| --- | --- | --- |
| `mesh.order` is undefined; `ks[NaN] = 0` silently no-ops the `v2` branch of corner clearing | `clothoid.js:543`, 599-611 | **Fix during port** — it is live and silent |
| `changeOrder` reads `temp[i]` where it means `temp[i1]` | `clothoid.js:405` | Fix; dormant code behind `if (0)` |
| `changeOrder` builds a `Float32Array` view over a `Float64Array` buffer, and assigns to `e._ks` not `e.curve._ks` | `clothoid.js:413` | Fix; same dormant block |
| `bstep` sums basis functions without weighting by `ks[i]` | `clothoid.js:45-52` | Fix or delete with the Bernstein experiment |
| Dead `dx`/`dy` pair overwritten two lines later | `clothoid.js:193-197` | Delete |
| Unreachable code after early `return` in `dstep`/`istep`/`istep2` | `clothoid.js` | Delete, preserving the live branch |
| `binomial` calls undefined `bin()` | path.ux `curve1d_bspline.js:159` | **Fixed** while vendoring — `src/math/util.ts` |
| `Solver.solveStepSimple` reads `con.wlst[j]` with `j` out of scope; throws a TDZ `ReferenceError` if ever reached. Unreachable only because `simple` defaults to `false` | path.ux `solver.js:171` | **Fixed** while vendoring — `wlst[i]`, `src/math/solver.ts` |
| `Clothoid.derivative2` finite-differenced `derivative`, feeding noise into `curvature` | `clothoid.js` | **Fixed** — analytic; `dth/ds` is `k(s)`, so it is `k` times the tangent rotated a quarter turn |
| `BSpline.update()` sets `this.flag`, which is never initialized, instead of `this.regen` — so `update()` never actually dirtied anything | `bspline.js:242` | **Fixed** |
| `BSpline.loadStruct` is lowercase (nstructjs looks for `loadSTRUCT`) and ORs in the undefined `BSplineFlags.FULL`, making `regen` NaN | `bspline.js:473` | **Fixed** — renamed, uses `BSplineRecalc.FULL` |
| `BSpline.STRUCT` declares `degere : int` | `bspline.js:491` | **Fixed** — `degree` |
| `dbasis` has three stacked implementations after the live `return` | `bspline.js:66-101` | Deleted; the reduce derivation is preserved as a doc comment |
| `Face.verts` generator closes over `this` instead of `this2`, so it throws | `mesh.js:630` | **Fixed** — delegates to `Face.loops` |
| `ElementArray.loadJSON` assigns `this.list.length = []` | `mesh.js:811` | Moot — the JSON path is dropped, see below |
| `Mesh.splitEdge` loads `ne.h1` twice; `ne.h2` is never placed | `mesh.js:1409-1410` | **Fixed** |
| `Mesh.splitEdge` computes `vector`/`h1`/`h2` and never uses them | `mesh.js:1400-1403` | Deleted |
| `Mesh.loadSTRUCT` ends with `switchSplineType(...)`, which discards every curve it just deserialized and re-ran `afterSTRUCT` on | `mesh.js:1720` | **Fixed** — ends with `regenSolve()` |
| `Mesh.makeEdge` places `h1` at 1/2 and `h2` at 2/3 | `mesh.js:1121` | **Fixed** — 1/3 and 2/3 |
| `Mesh.regen_render` calls `window.redraw_all()`; `Stroker.onInput` does too | `mesh.js:1613`, `stroke.js:181` | Deleted — the library does not reach into the host page |

Experiments currently parked behind `if (0)` are **not** defects — preserve them as
documented, switchable options rather than deleting them. They are the research record.

## Deliberate scope changes made during pass 1

| Change | Rationale |
| --- | --- |
| Dropped the `toJSON`/`loadJSON` pair on every mesh element | A second, unreferenced serialization path duplicating the `STRUCT` blocks, and buggier than them. nstructjs is the serialization layer we kept. |
| `Edge.curve` serializes as `abstract(Curve)` rather than a hardcoded struct name | The original interpolated `${Curve.structName}` from the `spline.js` alias, so a saved mesh could only ever round-trip clothoids. `Curve` is now registered as an empty base struct. |
| `LoopList[Symbol.iterator]` uses a generator, not the 1024-deep pooled `LoopListIter` | The generator was already written in the file, below an unreachable `return`. Same iteration order, no pool depth limit, no `stack.cur` leak. |
| `Edge.normal` delegates to `Curve.normal` | The original returned an unnormalized normal scaled by 0.01, in the opposite handedness from `Curve.normal`. Scaling for display is the renderer's job. |
| `Stroker` keeps a position ring and a count instead of `mpos1..6` / `v1..v6` fields | Same behaviour, and the "have we seen N samples yet" test becomes readable. |
| Curvature profiles are named exports selected by `setCurvatureProfile` | Replaces the module-level swap that had to be done by editing which function `bstep`/`istep` pointed at. |
