import { context } from "esbuild";
import { spawn } from "node:child_process";
import { mkdir, writeFile } from "node:fs/promises";
import { findpath } from "nw";

/** Unset by default: esbuild then finds a free port itself, so `pnpm serve` can coexist. */
const port = process.env.PORT ? Number(process.env.PORT) : undefined;
const cdpPort = Number(process.env.CDP_PORT ?? 9222);

const ctx = await context({
  entryPoints: ["src/demo/main.ts"],
  outfile    : "dist/demo.js",
  bundle     : true,
  format     : "esm",
  target     : "es2022",
  platform   : "browser",
  sourcemap  : true,
});

await ctx.watch();

const { port: served } = await ctx.serve({
  servedir: ".",
  fallback: "index.html",
  ...(port !== undefined && { port }),
});

const url = `http://127.0.0.1:${served}/`;

/**
 * NW.js wants an app directory with its own manifest, and the served port is only known
 * at runtime, so the manifest is generated under dist/ rather than checked in.
 */
const appDir = "dist/nwjs";
await mkdir(appDir, { recursive: true });
await writeFile(
  `${appDir}/package.json`,
  JSON.stringify(
    {
      name  : "polyclothoid-demo",
      main  : url,
      window: { title: "polyclothoid", width: 1280, height: 800 },
    },
    null,
    2
  )
);

const binary = await findpath();

console.log(`serving ${url}`);
console.log(`CDP endpoint http://127.0.0.1:${cdpPort}/`);

const child = spawn(binary, [`--remote-debugging-port=${cdpPort}`, appDir], {
  stdio: "inherit",
});

child.on("close", async (code) => {
  await ctx.dispose();
  process.exit(code ?? 0);
});
