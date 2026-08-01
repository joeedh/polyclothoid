import { build } from "esbuild";
import { rm } from "node:fs/promises";

const outdir = "dist";
const minify = process.argv.includes("--minify");

await rm(outdir, { recursive: true, force: true });

/**
 * The library bundle: no externals, since the library owns all of its deps.
 *
 * The "neutral" platform ignores package.json `main` unless told otherwise, and nstructjs
 * ships only that field, so resolution has to be spelled out.
 */
await build({
  entryPoints: ["src/index.ts"],
  outfile    : `${outdir}/polyclothoid.js`,
  bundle     : true,
  format     : "esm",
  target     : "es2022",
  platform   : "neutral",
  mainFields : ["module", "main"],
  conditions : ["import", "module", "default"],
  sourcemap  : true,
  minify,
});

/** The demo bundle, loaded by index.html. */
await build({
  entryPoints: ["src/demo/main.ts"],
  outfile    : `${outdir}/demo.js`,
  bundle     : true,
  format     : "esm",
  target     : "es2022",
  platform   : "browser",
  sourcemap  : true,
  minify,
});

console.log(`built -> ${outdir}/ ${minify ? "(minified)" : ""}`);
