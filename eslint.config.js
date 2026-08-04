import js from "@eslint/js";
import tseslint from "typescript-eslint";

export default tseslint.config([
  {
    ignores: ["dist/**", "node_modules/**", ".test-build/**"],
  },

  js.configs.recommended,
  tseslint.configs.recommended,
  tseslint.configs.stylistic,

  {
    files: ["src/**/*.ts", "tests/**/*.ts"],
    languageOptions: {
      parserOptions: {
        projectService: true,
        tsconfigRootDir: import.meta.dirname,
      },
    },
    rules: {
      // House style: let inference do the work, annotate only where it cannot.
      "@typescript-eslint/no-inferrable-types": "error",
      "@typescript-eslint/consistent-type-imports": ["error", { fixStyle: "inline-type-imports" }],
      "@typescript-eslint/no-unused-vars": ["error", { argsIgnorePattern: "^_", varsIgnorePattern: "^_" }],

      // In-progress markers must not survive into a finished branch.
      "no-warning-comments": ["error", { terms: ["CLAUDENOTE:"], location: "start" }],
    },
  },

  {
    // Vertex/Handle extend Element, so they get vector behaviour from a mixin instead of
    // a base class. Declaration merging with an empty interface is how that mixin is
    // typed; both rules below fire on exactly that pattern and have no alternative here.
    files: ["src/mesh/mesh.ts"],
    rules: {
      "@typescript-eslint/no-unsafe-declaration-merging": "off",
      "@typescript-eslint/no-empty-object-type": "off",
    },
  },

  {
    files: ["tools/**/*.mjs"],
    languageOptions: {
      globals: { process: "readonly", console: "readonly" },
    },
  },
]);
