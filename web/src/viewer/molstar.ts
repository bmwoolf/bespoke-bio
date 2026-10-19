import { PluginContext } from "molstar/lib/mol-plugin/context";
import { DefaultPluginSpec } from "molstar/lib/mol-plugin/spec";
import { PluginConfig } from "molstar/lib/mol-plugin/config";
import { Color } from "molstar/lib/mol-util/color";
import { createStructureRepresentationParams } from "molstar/lib/mol-plugin-state/helpers/structure-representation-params";
import { StateTransforms } from "molstar/lib/mol-plugin-state/transforms";
import {
  Structure,
  StructureProperties,
} from "molstar/lib/mol-model/structure";
import type { StateObjectSelector } from "molstar/lib/mol-state";
import type { PluginStateObject } from "molstar/lib/mol-plugin-state/objects";

export type ViewStyle = "surface" | "cartoon" | "atoms";
export type ViewSettings = {
  style: ViewStyle;
  metallic: boolean;
  sequenceGradient: boolean;
  ambientOcclusion: boolean;
  shadows: boolean;
  outlines: boolean;
  depthOfField: boolean;
};
export type ViewerInfo = {
  atoms: number;
  residues: number;
  chains: number;
  webgl: string;
  surfacePath: "GPU" | "CPU" | null;
  colors: { label: string; color: string }[];
};
export const defaultSettings: ViewSettings = {
  style: "surface",
  metallic: true,
  sequenceGradient: false,
  ambientOcclusion: true,
  shadows: true,
  outlines: true,
  depthOfField: true,
};

/** One plugin owns one canvas. Serialize representation changes and release GPU resources on unmount. */
export class MolecularScene {
  private pending: Promise<void> = Promise.resolve();
  private disposed = false;
  private representation:
    | StateObjectSelector<PluginStateObject.Molecule.Structure.Representation3D>
    | undefined;

  constructor(
    private plugin: PluginContext,
    private structure: StateObjectSelector<PluginStateObject.Molecule.Structure>,
  ) {}

  setStyle(settings: ViewSettings): Promise<ViewerInfo | undefined> {
    const operation = this.pending.then(async () => {
      if (this.disposed) return;
      const structure = this.structure.cell?.obj?.data;
      if (!structure) throw new Error("The structure is no longer available");
      const material = {
        metalness: settings.metallic ? 0.8 : 0,
        roughness: settings.metallic ? 0.22 : 0.65,
        bumpiness: 0,
      };
      const palette = {
        name: "colors" as const,
        params: {
          list: {
            kind: "set" as const,
            colors: [
              Color(0xb9cbdc),
              Color(0xd0a46a),
              Color(0x77b7ab),
              Color(0xab99cc),
            ],
          },
        },
      };
      const color = settings.sequenceGradient
        ? {
            color: "sequence-id" as const,
            colorParams: {
              list: {
                kind: "interpolate" as const,
                colors: [Color(0x7185e8), Color(0x55c2b8), Color(0xf0b35a)],
              },
            },
          }
        : { color: "chain-id" as const, colorParams: { palette } };
      this.plugin.canvas3d?.setProps({
        postprocessing: {
          occlusion: settings.ambientOcclusion
            ? {
                name: "on",
                params: {
                  samples: 32,
                  multiScale: { name: "off", params: {} },
                  radius: 5,
                  bias: 0.8,
                  blurKernelSize: 15,
                  blurDepthBias: 0.5,
                  resolutionScale: 1,
                  color: Color(0x080c10),
                  transparentThreshold: 0.4,
                },
              }
            : { name: "off", params: {} },
          shadow: settings.shadows
            ? { name: "on", params: { steps: 4, maxDistance: 3, tolerance: 1 } }
            : { name: "off", params: {} },
          outline: settings.outlines
            ? {
                name: "on",
                params: {
                  scale: 1,
                  threshold: 0.75,
                  color: Color(0x080c10),
                  includeTransparent: true,
                },
              }
            : { name: "off", params: {} },
          dof: settings.depthOfField
            ? {
                name: "on",
                params: {
                  blurSize: 9,
                  blurSpread: 0.5,
                  inFocus: 0,
                  PPM: 60,
                  center: "camera-target",
                  mode: "sphere",
                },
              }
            : { name: "off", params: {} },
        },
      });
      const params =
        settings.style === "surface"
          ? createStructureRepresentationParams(this.plugin, structure, {
              type: "gaussian-surface",
              typeParams: {
                material,
                tryUseGpu: true,
                quality: "custom",
                resolution: 0.8,
              },
              ...color,
            })
          : settings.style === "cartoon"
            ? createStructureRepresentationParams(this.plugin, structure, {
                type: "cartoon",
                typeParams: { material },
                ...color,
              })
            : createStructureRepresentationParams(this.plugin, structure, {
                type: "ball-and-stick",
                typeParams: { material },
                ...color,
              });
      const update = this.plugin.build();
      if (this.representation) {
        await update.to(this.representation).update(params).commit();
      } else {
        this.representation = await update
          .to(this.structure)
          .apply(
            StateTransforms.Representation.StructureRepresentation3D,
            params,
          )
          .commit();
      }
      if (this.disposed) return;
      const objects =
        this.representation.cell?.obj?.data.repr.renderObjects ?? [];
      if (!objects.length)
        throw new Error("Mol* did not produce a molecular representation");
      const legend =
        this.representation.cell?.obj?.data.repr.theme.color.legend;
      const chainIds = new Set<string>();
      Structure.eachAtomicHierarchyElement(structure, {
        chain: (location) =>
          chainIds.add(StructureProperties.chain.auth_asym_id(location)),
      });
      return {
        colors:
          !settings.sequenceGradient && legend?.kind === "table-legend"
            ? legend.table.map(([label, color]) => ({
                label,
                color: `#${color.toString(16).padStart(6, "0")}`,
              }))
            : [],
        atoms: structure.elementCount,
        residues: structure.polymerResidueCount,
        chains: chainIds.size,
        webgl: this.plugin.canvas3d?.webgl.isWebGL2 ? "WebGL 2" : "WebGL 1",
        surfacePath:
          settings.style === "surface"
            ? objects.some((object) => object.type === "texture-mesh")
              ? ("GPU" as const)
              : ("CPU" as const)
            : null,
      };
    });
    this.pending = operation.then(
      () => undefined,
      () => undefined,
    );
    return operation;
  }

  resize() {
    if (!this.disposed) {
      this.plugin.handleResize();
      // Sync the camera viewport before computing its new aspect-dependent distance.
      // requestResize alone runs after the camera-reset queue in Mol*'s render loop.
      this.plugin.canvas3d?.handleResize();
      // Keep the full molecule framed when the viewport changes aspect ratio.
      this.plugin.canvas3d?.requestCameraReset({ durationMs: 0 });
    }
  }
  resetCamera() {
    if (!this.disposed)
      this.plugin.canvas3d?.requestCameraReset({
        durationMs: window.matchMedia("(prefers-reduced-motion: reduce)")
          .matches
          ? 0
          : 250,
      });
  }
  dispose() {
    if (this.disposed) return;
    this.disposed = true;
    void this.pending.then(() => this.plugin.dispose());
  }
}

export async function loadScene(
  canvas: HTMLCanvasElement,
  parent: HTMLDivElement,
  url: string,
  signal: AbortSignal,
): Promise<MolecularScene> {
  const plugin = new PluginContext({
    ...DefaultPluginSpec(),
    config: [[PluginConfig.VolumeStreaming.Enabled, false]],
  });
  try {
    await plugin.init();
    signal.throwIfAborted();
    if (!(await plugin.initViewerAsync(canvas, parent)))
      throw new Error(
        "WebGL is unavailable. Download the structure to inspect it in another viewer.",
      );
    signal.throwIfAborted();
    plugin.canvas3d?.setProps({
      renderer: { backgroundColor: Color(0x141b22), ambientIntensity: 0.6 },
      camera: {
        manualReset: true,
        helper: { axes: { name: "off", params: {} } },
      },
    });
    const response = await fetch(url, { signal });
    if (!response.ok)
      throw new Error(`Structure download failed (${response.status})`);
    const contents = await response.text();
    signal.throwIfAborted();
    const data = await plugin.builders.data.rawData({
      data: contents,
      label: "Design artifact",
    });
    const trajectory = await plugin.builders.structure.parseTrajectory(
      data,
      url.endsWith(".pdb") ? "pdb" : "mmcif",
    );
    if (trajectory.cell?.obj?.data.frameCount !== 1)
      throw new Error(
        "This viewer requires one coordinate model per artifact; download multi-model files for inspection.",
      );
    const model = await plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.builders.structure.createStructure(model, {
      name: "model",
      params: {},
    });
    if (!structure.cell?.obj?.data.elementCount)
      throw new Error("The artifact contains no atomic coordinates");
    signal.throwIfAborted();
    return new MolecularScene(plugin, structure);
  } catch (error) {
    plugin.dispose();
    throw error;
  }
}
