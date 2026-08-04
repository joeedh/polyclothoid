import { build } from "esbuild";
import { rm, readdir } from "node:fs/promises";
import { spawn } from "node:child_process";
import path from "node:path";

/*
  Node's own type stripping cannot run these sources directly: they import with a `.js`
  extension that only a bundler resolves back to `.ts`. So bundle each test to
  `.test-build/` and hand the directory to `node --test`.
*/
const outdir = ".test-build";
const srcdir = "tests";

async function walk(dir) {
  const out = [];

  for (const ent of await readdir(dir, { withFileTypes: true })) {
    const full = path.join(dir, ent.name);

    if (ent.isDirectory()) {
      out.push(...(await walk(full)));
    } else if (ent.name.endsWith(".test.ts")) {
      out.push(full);
    }
  }

  return out;
}

const entryPoints = await walk(srcdir);

if (entryPoints.length === 0) {
  console.error(`no *.test.ts under ${srcdir}/`);
  process.exit(1);
}

await rm(outdir, { recursive: true, force: true });

await build({
  entryPoints,
  outdir,
  outbase   : srcdir,
  bundle    : true,
  format    : "esm",
  target    : "es2022",
  platform  : "node",
  mainFields: ["module", "main"],
  conditions: ["import", "module", "default"],
  sourcemap : "inline",
});

const args = ["--test", "--enable-source-maps"];

// Everything after `--` on our own command line is forwarded (e.g. `--test-name-pattern`).
const sep = process.argv.indexOf("--");
if (sep !== -1) {
  args.push(...process.argv.slice(sep + 1));
}

// A directory argument is resolved as a module; the glob form is what enumerates files.
args.push(`${outdir}/**/*.js`);

const child = spawn(process.execPath, args, { stdio: "inherit" });
child.on("exit", (code) => process.exit(code ?? 1));
