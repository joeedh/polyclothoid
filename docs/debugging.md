# Debugging the demo

The demo runs in two harnesses: a plain browser against the dev server, and an NW.js
desktop shell. Both serve the same `index.html` + `dist/demo.js` bundle with sourcemaps
and esbuild watch, so edits to `src/` rebuild on save and a reload picks them up.

## Browser

```
pnpm serve
```

Starts esbuild's dev server on port 8080 (`PORT` overrides) and watches `src/demo`. Open
<http://localhost:8080/> and use the browser's own devtools.

## NW.js shell

```
pnpm nwjs
```

`tools/nwjs.mjs` does three things:

1. Starts its own esbuild watch + serve, on whatever free port esbuild finds (it scans
   upward from 8000), so it never contends with a `pnpm serve` already running on 8080.
   Set `PORT` to pin it. The served URL is printed at launch.
2. Generates an NW.js app manifest under `dist/nwjs/` whose `main` points at the served
   URL — NW.js needs a manifest directory of its own, and the port is only known at
   runtime, so it is not checked in.
3. Launches the NW.js **SDK** binary (installed by the `nw` dev dependency, pinned to an
   `-sdk` version) with `--remote-debugging-port=9222` (`CDP_PORT` overrides).

Closing the NW.js window shuts the dev server down with it.

Because it is the SDK flavor, DevTools are available in-window: press `F12` (or
right-click → Inspect).

## Connecting over CDP

The shell always starts with the Chrome DevTools Protocol listening on
`http://127.0.0.1:9222/`. Ways in:

- **List targets:** `curl http://127.0.0.1:9222/json/list` — each target carries a
  `webSocketDebuggerUrl` for raw CDP clients.
- **Remote DevTools from a browser:** open `http://127.0.0.1:9222/` in Chrome and click
  the page target, or add `localhost:9222` under `chrome://inspect` → *Configure…* and
  attach from there.
- **Puppeteer / Playwright:** attach to the running shell rather than launching a
  browser:

  ```js
  // puppeteer
  const browser = await puppeteer.connect({ browserURL: "http://127.0.0.1:9222" });

  // playwright
  const browser = await chromium.connectOverCDP("http://127.0.0.1:9222");
  ```

The demo page is the target whose URL matches the served URL printed at launch — NW.js
also exposes its own background page as a target; ignore it.

Note that CDP here debugs the *page* (DOM, console, JS breakpoints in the demo bundle).
There is no Node context to attach to — the shell loads the app over HTTP and grants it
no Node access.
