import { useEffect, useRef, useState } from "react";
import type { StructureSelection } from "../api";
import { ErrorMessage, Loading, formatCount } from "../components/Status";
import { defaultSettings, loadScene } from "./molstar";
import type { MolecularScene, ViewerInfo, ViewSettings } from "./molstar";

export default function MolecularViewer({
  selection,
}: {
  selection: StructureSelection;
}) {
  const framed = useRef(false);
  const parent = useRef<HTMLDivElement>(null);
  const canvas = useRef<HTMLCanvasElement>(null);
  const [scene, setScene] = useState<MolecularScene | null>(null);
  const [info, setInfo] = useState<ViewerInfo | null>(null);
  const [error, setError] = useState<Error | null>(null);
  const [busy, setBusy] = useState(true);
  const [settings, setSettings] = useState<ViewSettings>(defaultSettings);

  useEffect(() => {
    const abort = new AbortController();
    let loaded: MolecularScene | undefined;
    let resize: ResizeObserver | undefined;
    const container = parent.current;
    const element = canvas.current;
    if (!container || !element) return;
    void loadScene(element, container, selection.url, abort.signal)
      .then((value) => {
        if (abort.signal.aborted) {
          value.dispose();
          return;
        }
        loaded = value;
        resize = new ResizeObserver(() => value.resize());
        resize.observe(container);
        setScene(value);
      })
      .catch((reason: unknown) => {
        if (!abort.signal.aborted) {
          setError(
            reason instanceof Error
              ? reason
              : new Error("Could not load structure"),
          );
          setBusy(false);
        }
      });
    return () => {
      abort.abort();
      resize?.disconnect();
      loaded?.dispose();
    };
  }, [selection.url]);

  useEffect(() => {
    if (!scene) return;
    let current = true;
    setBusy(true);
    setError(null);
    void scene
      .setStyle(settings)
      .then((value) => {
        if (!current || !value) return;
        if (!framed.current) {
          scene.resetCamera();
          framed.current = true;
        }
        setInfo(value);
        setBusy(false);
      })
      .catch((reason: unknown) => {
        if (current) {
          setError(
            reason instanceof Error
              ? reason
              : new Error("Could not build representation"),
          );
          setBusy(false);
        }
      });
    return () => {
      current = false;
    };
  }, [scene, settings]);

  return (
    <section className="molecule" aria-label="Molecular structure">
      <div className="viewer-heading">
        <h3>{selection.title}</h3>
        <a href={selection.url}>Download structure</a>
      </div>
      <div className="viewer-controls">
        <div className="segmented" aria-label="Representation">
          {(
            [
              ["surface", "Surface"],
              ["cartoon", "Ribbon"],
              ["atoms", "Atoms"],
            ] as const
          ).map(([style, label]) => (
            <button
              className="secondary"
              aria-pressed={settings.style === style}
              key={style}
              onClick={() => setSettings({ ...settings, style })}
            >
              {label}
            </button>
          ))}
        </div>
        <label className="check">
          <input
            type="checkbox"
            checked={settings.metallic}
            onChange={(e) =>
              setSettings({ ...settings, metallic: e.target.checked })
            }
          />
          Metallic finish
        </label>
        <label className="check">
          <input
            type="checkbox"
            checked={settings.colorByChain}
            onChange={(e) =>
              setSettings({ ...settings, colorByChain: e.target.checked })
            }
          />
          Color by chain
        </label>
        <button
          className="secondary"
          disabled={!scene || busy}
          onClick={() => scene?.resetCamera()}
        >
          Reset view
        </button>
      </div>
      <ErrorMessage error={error} />
      <div className="molstar-viewport" ref={parent} aria-busy={busy}>
        <canvas ref={canvas} aria-label="Interactive molecular structure" />
        {busy && (
          <div className="viewer-loading">
            <Loading>Building molecular representation…</Loading>
          </div>
        )}
      </div>
      <p className="viewer-info" role="status" data-testid="viewer-info">
        {info &&
          `${formatCount(info.atoms)} atoms · ${formatCount(info.residues)} residues · ${info.chains} ${info.chains === 1 ? "chain" : "chains"} · ${info.webgl}${info.surfacePath ? ` · ${info.surfacePath} surface computation` : ""}`}
      </p>
      <p className="muted">
        Drag to rotate · Scroll to zoom. The surface is calculated from atomic
        coordinates; the metallic finish changes lighting only.
      </p>
    </section>
  );
}
