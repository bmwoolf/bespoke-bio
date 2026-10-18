import { PluginContext } from "molstar/lib/mol-plugin/context";
import { DefaultPluginSpec } from "molstar/lib/mol-plugin/spec";
import { PluginConfig } from "molstar/lib/mol-plugin/config";
import { Color } from "molstar/lib/mol-util/color";
import { createStructureRepresentationParams } from "molstar/lib/mol-plugin-state/helpers/structure-representation-params";
import { StateTransforms } from "molstar/lib/mol-plugin-state/transforms";
import type { StateObjectSelector } from "molstar/lib/mol-state";
import type { PluginStateObject } from "molstar/lib/mol-plugin-state/objects";

export type ViewStyle = "surface" | "cartoon" | "atoms";
export type ViewSettings = {
  style: ViewStyle;
  metallic: boolean;
  colorByChain: boolean;
};
export type ViewerInfo = {
  atoms: number;
  residues: number;
  chains: number;
  webgl: string;
  surfacePath: "GPU" | "CPU" | null;
};
export const defaultSettings: ViewSettings = {
  style: "surface",
  metallic: true,
  colorByChain: true,
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
      const color = settings.colorByChain
        ? { color: "chain-id" as const, colorParams: { palette } }
        : {
            color: "uniform" as const,
            colorParams: { value: Color(0xc4ced8) },
          };
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
      return {
        atoms: structure.elementCount,
        residues: structure.polymerResidueCount,
        chains: structure.models.reduce(
          (sum, model) => sum + model.atomicHierarchy.chains._rowCount,
          0,
        ),
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
    if (!this.disposed) this.plugin.handleResize();
  }
  resetCamera() {
    if (!this.disposed)
      this.plugin.canvas3d?.requestCameraReset({ durationMs: 250 });
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
      renderer: { backgroundColor: Color(0x101b27), ambientIntensity: 0.6 },
      camera: { manualReset: true },
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
