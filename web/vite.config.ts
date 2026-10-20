import { defineConfig } from "vitest/config";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  build: {
    outDir: "../src/bespoke_bio/demo/static",
    emptyOutDir: true,
  },
  server: {
    host: "127.0.0.1",
    proxy: {
      // Preserve the browser's Host/Origin pair for the gateway's same-origin check.
      "/api": { target: "http://127.0.0.1:8787", changeOrigin: false },
    },
  },
  test: { include: ["src/**/*.test.ts"] },
});
