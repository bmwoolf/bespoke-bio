import { useEffect, useRef, useState } from "react";
import type { StructureSelection } from "../api";
import { ErrorMessage, Loading, formatCount } from "../components/Status";
import { Popover } from "../components/Popover";
import { defaultSettings, loadScene } from "./molstar";
import type { MolecularScene, ViewerInfo, ViewSettings } from "./molstar";

export default function MolecularViewer({
  selection,
}: {
  selection: StructureSelection;
}) {
  const framed = useRef(false);
  const frame = useRef<HTMLElement>(null);
  const [fullscreen, setFullscreen] = useState(false);
  const [fullscreenError, setFullscreenError] = useState<Error | null>(null);
  useEffect(() => {
    const changed = () =>
      setFullscreen(document.fullscreenElement === frame.current);
    document.addEventListener("fullscreenchange", changed);
    return () => document.removeEventListener("fullscreenchange", changed);
  }, []);
  async function toggleFullscreen() {
    try {
      setFullscreenError(null);
      if (document.fullscreenElement === frame.current)
        await document.exitFullscreen();
      else await frame.current?.requestFullscreen();
    } catch {
      setFullscreenError(
        new Error("Fullscreen is unavailable in this browser."),
      );
    }
  }
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
    <section className="molecule" aria-label="Molecular structure" ref={frame}>
      <ErrorMessage error={error} />
      <ErrorMessage error={fullscreenError} />
      {error && (
        <a className="fallback-download" href={selection.url}>
          Download structure
        </a>
      )}
      <div className="molstar-viewport" ref={parent} aria-busy={busy}>
        <canvas ref={canvas} aria-label={selection.title} />
        <div className="viewer-topline">
          <span className="viewer-hint">Drag to rotate · Scroll to zoom</span>
          <Popover className="appearance-popover" label="Appearance" icon="◐">
            <label className="check">
              <input
                type="checkbox"
                checked={settings.metallic}
                onChange={(event) =>
                  setSettings({ ...settings, metallic: event.target.checked })
                }
              />
              Metallic finish
            </label>
            <label className="check">
              <input
                type="checkbox"
                checked={settings.colorByChain}
                onChange={(event) =>
                  setSettings({
                    ...settings,
                    colorByChain: event.target.checked,
                  })
                }
              />
              Color by chain
            </label>
            <p className="viewer-info" role="status" data-testid="viewer-info">
              {info &&
                `${formatCount(info.atoms)} atoms · ${formatCount(info.residues)} residues · ${info.chains} ${info.chains === 1 ? "chain" : "chains"} · ${info.webgl}${info.surfacePath ? ` · ${info.surfacePath} surface computation` : ""}`}
            </p>
            <p className="muted">
              Surface from atomic coordinates. Metallic finish changes lighting
              only.
            </p>
          </Popover>
        </div>
        {busy && (
          <div className="viewer-loading">
            <Loading>Loading structure…</Loading>
          </div>
        )}
        <div className="viewer-bottomline">
          <div className="chain-legend" aria-label="Chain colors">
            {info?.colors.map((entry) => (
              <span key={entry.label}>
                <i style={{ background: entry.color }} />
                {entry.label}
              </span>
            ))}
          </div>
          <div className="segmented" role="group" aria-label="Representation">
            {(
              [
                ["surface", "Surface"],
                ["cartoon", "Ribbon"],
                ["atoms", "Atoms"],
              ] as const
            ).map(([style, label]) => (
              <button
                aria-pressed={settings.style === style}
                key={style}
                onClick={() => setSettings({ ...settings, style })}
              >
                {label}
              </button>
            ))}
          </div>
          <div className="camera-controls">
            <button
              className="icon-button"
              aria-label="Reset view"
              title="Reset view"
              disabled={!scene || busy}
              onClick={() => scene?.resetCamera()}
            >
              <svg
                viewBox="0 0 24 24"
                fill="none"
                stroke="currentColor"
                strokeWidth="1.6"
                aria-hidden="true"
              >
                <path d="M4 10a8 8 0 1 1 1 7M4 4v6h6" />
              </svg>
            </button>
            {document.fullscreenEnabled && (
              <button
                className="icon-button"
                aria-label={fullscreen ? "Exit fullscreen" : "Fullscreen"}
                title={fullscreen ? "Exit fullscreen" : "Fullscreen"}
                onClick={() => void toggleFullscreen()}
              >
                <svg
                  viewBox="0 0 24 24"
                  fill="none"
                  stroke="currentColor"
                  strokeWidth="1.6"
                  aria-hidden="true"
                >
                  <path d="M4 9V4h5M15 4h5v5M20 15v5h-5M9 20H4v-5" />
                </svg>
              </button>
            )}
          </div>
        </div>
      </div>
    </section>
  );
}
