import { context } from "esbuild";

const port = Number(process.env.PORT ?? 8080);

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

const { hosts, port: served } = await ctx.serve({
  servedir: ".",
  fallback: "index.html",
  port,
});

for (const host of hosts) {
  console.log(`serving http://${host}:${served}/`);
}
