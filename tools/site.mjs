import { cp, mkdir, rm, writeFile } from "node:fs/promises";

/**
 * Assembles the static site published to GitHub Pages.
 *
 * `index.html` loads `./dist/demo.js`, so the layout inside `site/` mirrors the repo
 * root rather than flattening. Run `pnpm build` first — this only copies.
 */
const outdir = "site";

await rm(outdir, { recursive: true, force: true });
await mkdir(`${outdir}/dist`, { recursive: true });

await cp("index.html", `${outdir}/index.html`);
await cp("dist/demo.js", `${outdir}/dist/demo.js`);
await cp("dist/demo.js.map", `${outdir}/dist/demo.js.map`);

/** Keeps Pages from running the output through Jekyll. */
await writeFile(`${outdir}/.nojekyll`, "");

console.log(`assembled -> ${outdir}/`);
