# polyclothoid

A polynomial clothoid-based spline solver optimized for brush strokes, with a small
vanilla-canvas demo for looking at what it produces.

The curves are defined by a sampled curvature profile `k(s)` rather than by control
points. Geometry is recovered by integrating `θ(s) = ∫k ds`, then `x = ∫cos θ ds`,
`y = ∫sin θ ds`. Because `s` is arclength, the unit tangent is exactly `(cos θ, sin θ)`,
which is what makes even dab spacing cheap. See `docs/research/clothoids.md`.

WIP.

General ideas:

- Try out shortest path decomposition paper?
  - Would be fairly easy to do with a bspline curvature function
  - No need to have separate arc/clothoid/line primitives in that case, you'd just use a
    degree 1 (linear) bspline.

## Setup

Requires [pnpm](https://pnpm.io) and Node 20+.

```
git clone https://github.com/joeedh/polyclothoid.git
cd polyclothoid
pnpm install
```

## Running the demo

```
pnpm serve
```

esbuild watches `src/`, rebuilds `dist/demo.js`, and serves the page on
`http://127.0.0.1:8080/` (set `PORT` to change it).

The demo has two modes, switched from the panel or with `Space`:

- **Edit mesh** — drag vertices, shift-click to extend the selection. `A` selects all or
  none, `E` splits selected edges, `D` dissolves 2-valence vertices, `X` deletes.
- **Paint stroke** — drag to lay down dabs. Every input position is fitted twice, once
  with a clothoid and once with a cubic bezier, and both dab trails are drawn so their
  spacing can be compared directly. `C` clears them.

## Building

```
pnpm build          # dist/polyclothoid.js (library) and dist/demo.js
pnpm build --minify
pnpm types          # dist/types/**.d.ts
```

## Checks

```
pnpm check          # typecheck + lint + format:check
pnpm format
```

## Using the library

```js
import { Stroker } from "polyclothoid";

const stroker = new Stroker((x, y, dx, dy) => paintDab(x, y));

canvas.addEventListener("pointermove", (e) => {
  stroker.onInput(e.offsetX, e.offsetY, brushRadius, brushSpacing);
});
```

`Stroker` takes a second argument to swap the curve type it fits with — see
`StrokerOptions`. The lower-level pieces (`Mesh`, `Clothoid`, `ClothoidSolver`,
`CubicBezier`, `BSpline`) are exported too.

## Layout

| Path | What |
| --- | --- |
| `src/math` | vectors, `CacheRing`, `IDGen`, the Gauss-Seidel `Solver` |
| `src/curve` | `Curve` base plus the clothoid, bezier and b-spline implementations |
| `src/mesh` | half-edge mesh whose edges carry pluggable curve geometry |
| `src/stroke.ts` | `Stroker` — pointer positions in, evenly spaced dabs out |
| `src/demo` | the canvas harness; not part of the published library |
| `docs/` | design notes, `docs/research/` findings, `docs/plans/` work plans |
