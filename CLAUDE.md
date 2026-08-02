# polyclothoid

A polynomial clothoid-based spline solver optimized for brush strokes. Curves are defined
by a sampled **curvature profile** `k(s)` and recovered by integration, not by control
points. See `docs/research/clothoids.md` for the math and `docs/architecture.md` for the
module boundaries and public API surface.

This is a **library first**. Everything under `src/` except `src/demo` must stay free of
DOM and UI concerns; the canvas demo is a harness for eyeballing curve quality, not a
product. Curves do expose a `draw()` hook, but it takes a structural `Canvas2DLike`
declared by the library rather than the DOM's `CanvasRenderingContext2D`.

## Toolchain

Package manager is **pnpm**. Do not use npm or yarn — the lockfile is pnpm's.

| Task | Command |
| --- | --- |
| Typecheck | `pnpm typecheck` (native `tsgo`, not `tsc`) |
| Bundle | `pnpm build` (esbuild) |
| Dev server | `pnpm serve` (esbuild's own HTTP server, port 8080) |
| Pages bundle | `pnpm build:site` (minified build, then assemble `site/`) |
| Format | `pnpm format` (`@pathtx/prettier`, a fork — not upstream prettier) |
| Lint | `pnpm lint` (typescript-eslint, flat config) |
| Declarations | `pnpm types` (`tsgo --emitDeclarationOnly` to `dist/types`) |
| All gates | `pnpm check` |

`@pathtx/prettier` supplies the `alignPropertyValues: "group"` option that produces the
aligned-colon object style used throughout this codebase. Upstream prettier will reformat
those blocks incorrectly — do not swap it out.

## Deployment

The demo is published to GitHub Pages at <https://joeedh.github.io/polyclothoid/> by
`.github/workflows/pages.yml`, on every push to `master` plus `workflow_dispatch`. The
build job gates on `pnpm typecheck`, runs `pnpm build:site`, and uploads `site/`; a second
job deploys it. `site/` is generated and gitignored — never commit it.

Three things about that setup are load-bearing:

1. **`site/` mirrors the repo root rather than flattening.** `index.html` loads
   `./dist/demo.js` by relative path, so `tools/site.mjs` copies the bundle to
   `site/dist/demo.js`. That relative path is also what makes the `/polyclothoid/` Pages
   subpath work without a base href. Do not hoist the bundle to `site/`.

2. **`tools/site.mjs` writes a `.nojekyll` marker.** Without it Pages runs the output
   through Jekyll, which drops paths beginning with an underscore.

3. **`pnpm/action-setup` must run before `actions/setup-node`.** `setup-node`'s
   `cache: pnpm` shells out to pnpm to locate the store, so the reverse order fails the
   job outright.

The Pages source is configured as **GitHub Actions** in repo settings, not a branch.
`actions/configure-pages` fails with a bare `Not Found` if that ever gets switched back.

## TypeScript style

**Prefer implicit type annotation over explicit typings.** Let inference do the work;
annotate only where it genuinely cannot infer — public API boundaries, empty container
initializers, and cases where the inferred type is wrong or uselessly wide. Redundant
annotations are a lint error (`@typescript-eslint/no-inferrable-types`).

`verbatimModuleSyntax` is on, so type-only imports need the `type` keyword. Use the inline
form: `import { type Curve, Mesh } from "./mesh.js"`.

Import paths carry the `.js` extension even in `.ts` sources — that is what
`moduleResolution: "bundler"` plus ESM output expects.

## Invariants

Four things that look like bugs and are not. Breaking any of them fails quietly.

1. **A constraint returning `0.0` is meaningful.** It short-circuits twice: in
   `Constraint.evaluate`, before any gradient is finite-differenced, and again in
   `Solver.solveStep`, before any parameter is touched. `curv_c` in the clothoid solver
   returns `0.0` unconditionally on purpose and does its work as a direct projection
   interleaved with `tan_c`'s gradient descent. Do not "fix" it into returning a residual.

2. **`s` is arclength, not a normalized parameter.** `evaluate(0.5)` is half a unit along
   the curve, not halfway along it. Normalized parameters exist only inside a curve
   implementation.

3. **Vectors from `evaluate`/`derivative`/`normal` are borrowed** from a `CacheRing` and
   stay valid for a bounded number of further calls. Copy anything you keep.

4. **Never patch built-in prototypes.** See History below.

`curve/` must not import `mesh/`. Solvers reach mesh-shaped things through the structural
`Solvable*` interfaces in `curve/mesh_types.ts`; the dependency runs one way, mesh → curve.

## Serialization

`nstructjs` is the only runtime dependency. Register structs inline — `nstructjs.inherit()`
is deprecated:

```ts
static override STRUCT = nstructjs.inlineRegister(this, `mesh.Vertex { ... }`);
```

`inlineRegister` walks inheritance itself, so there is no separate `register()` call and
child structs need not repeat parent fields. It has no default export, so import it as
`import * as nstructjs from "nstructjs"`.

## Code comments

Two hard rules, both enforced in review:

1. **In-progress comments start with `CLAUDENOTE:`.** Anything provisional — an assumption
   being carried, a type that needs revisiting, a question for a later pass — gets that
   prefix so it can be found and stripped mechanically. `pnpm lint` fails while any
   `CLAUDENOTE:` remains, so a branch cannot merge with them in place.

2. **Non-doc comments are capped at 3 lines.** If an explanation needs more room, it
   belongs in a doc comment (`/** ... */`) or in `docs/`. The sole exception is the core
   clothoid solver, `src/curve/clothoid.ts`, where the math genuinely needs prose and
   long-form derivations are welcome.

Doc comments (`/** ... */`) are exempt from the line cap everywhere.

## Documentation layout

| Kind | Location |
| --- | --- |
| Research, derivations, algorithm notes | `docs/research/` |
| Design and architecture | `docs/` (root) |
| Plans and task tracking | `docs/plans/` |

Put new documents in the right directory the first time. `docs/plans/typescript-port.md`
tracked the TypeScript port and is now complete — keep it as the record of why things are
shaped the way they are, and start a new plan doc for new work rather than reopening it.
Update plan docs as work lands rather than at the end.

`docs/research/spower.md` is background reading, not a description of anything in `src/`:
Sánchez-Reyes' s-power series, a two-point Taylor analogue whose truncations are Hermite
interpolants, and the 2003 result that approximates clothoids with them. Nothing is
implemented — treat it as a survey of a possible direction, and note that its `s` is a
symmetric parameter `(1−u)u`, not this codebase's arclength `s`.

## Layout

```
src/
  index.ts        public API surface — re-exports all four barrels flat
  math/           vendored vector/util/solver layer (no external deps)
  curve/          curve primitives: clothoid, bezier, bspline, and their solvers
  mesh/           the spline mesh the solver runs over
  stroke.ts       the brush stroker itself
  demo/           canvas harness, vanilla TS, no UI library
tools/            esbuild build and serve scripts, plus site.mjs (Pages assembly)
docs/             see above
.github/workflows/  pages.yml — build and deploy the demo
```

The curve type is a runtime choice, not an import-time one: `Mesh` holds `CurveCls` and
`SolverCls`, `Mesh.switchSplineType()` swaps them, and `Stroker` takes the same pair via
`StrokerOptions`. Keep the clothoid, bezier and b-spline implementations
interface-compatible — comparing them on identical input is half the point of the project.

## History

The project previously depended on **path.ux** (the author's UI toolkit) as a git
submodule, and the stroker imported its vector math, `cachering`, `Constraint`/`Solver`,
and `nstructjs` re-export. That submodule has been removed: the library used only a small
slice of it, now vendored into `src/math`. Do not reintroduce the dependency.

The original JS sources under `scripts/`, along with the rollup config and the hand-rolled
`serv.js`, were deleted once the port landed. Git history has them if a behaviour needs
checking against the original.

path.ux also monkeypatched globals — `Math.fract`, `Array.prototype.remove`. The port
replaces these with explicit imports. Never patch built-in prototypes in this codebase.
