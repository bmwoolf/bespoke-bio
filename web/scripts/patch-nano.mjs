import { readFile, writeFile } from "node:fs/promises";

// Compatibility fixes for Nano 0.0.15's chain IDs, light background, and camera.
const bundle = new URL(
  "../node_modules/@juliocesar-io/nano-protein-viewer-react/dist/index.esm.js",
  import.meta.url,
);

const replacements = [
  [
    `              const P = I.unit, B = P?.label_asym_id ?? P?.chainGroupId ?? "A", U = x.chainColors && (x.chainColors[B] ?? x.chainColors[String(B)]) || 5164484;
              return typeof U == "string" ? parseInt(String(U).replace("#", ""), 16) : U;`,
    `              const P = I.unit, B = P?.model?.atomicHierarchy, U = B?.chainAtomSegments?.index?.[I.element], z = Number.isInteger(U) ? B?.chains?.auth_asym_id?.value?.(U) : void 0, L = z ?? P?.label_asym_id ?? P?.chainGroupId ?? "A", F = x.chainColors && (x.chainColors[L] ?? x.chainColors[String(L)]) || 5164484;
              return typeof F == "string" ? parseInt(String(F).replace("#", ""), 16) : F;`,
  ],
  [
    `            const w = D?.label_asym_id ?? D?.chainGroupId ?? "A";
            w && C.add(String(w));`,
    `            const w = D?.model?.atomicHierarchy, I = D?.elements?.[0], P = Number.isInteger(I) ? w?.chainAtomSegments?.index?.[I] : void 0, B = Number.isInteger(P) ? w?.chains?.auth_asym_id?.value?.(P) : void 0, U = B ?? D?.label_asym_id ?? D?.chainGroupId ?? "A";
            U && C.add(String(U));`,
  ],
  [
    `  async function f() {
    e && e.plugin.managers.camera.reset();
  }`,
    `  async function f() {
    if (!e) return;
    const g = e.plugin, x = g.canvas3d?.scene?.boundingSphereVisible;
    g.canvas3d?.setProps({ camera: { helper: { axes: { name: "off", params: {} } } } });
    x ? g.managers.camera.focusSphere(x, { extraRadius: 0, durationMs: 0 }) : g.managers.camera.reset();
  }`,
  ],
  [
    `          const C = g === "dark" ? { r: 0, g: 0, b: 0 } : { r: 255, g: 255, b: 255 };
          x.canvas3d.setProps({
            renderer: {
              backgroundColor: C
            }
          });`,
    `          const C = g === "dark" ? 0 : 16777215;
          x.canvas3d.setProps({
            renderer: {
              backgroundColor: C
            },
            camera: {
              helper: { axes: { name: "off", params: {} } }
            }
          });`,
  ],
];

let source = await readFile(bundle, "utf8");
let changed = false;

for (const [original, replacement] of replacements) {
  if (source.includes(replacement)) continue;
  if (!source.includes(original)) {
    throw new Error(
      "Nano viewer 0.0.15 changed; review scripts/patch-nano.mjs before installing.",
    );
  }
  source = source.replace(original, replacement);
  changed = true;
}

if (changed) await writeFile(bundle, source);
