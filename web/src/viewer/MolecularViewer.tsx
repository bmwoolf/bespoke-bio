import { useEffect, useRef, useState } from "react";
import type { MolstarViewerHandle } from "@juliocesar-io/nano-protein-viewer-react";
import type { StructureSelection } from "../api";
import { ErrorMessage, Loading, formatCount } from "../components/Status";
import { Popover } from "../components/Popover";
import { defaultSettings, loadScene } from "./molstar";
import type { MolecularScene, ViewerInfo, ViewSettings } from "./molstar";

type AppearanceToggle = Exclude<keyof ViewSettings, "style">;
const appearanceToggles: readonly [AppearanceToggle, string][] = [
  ["metallic", "Metallic finish"],
  ["sequenceGradient", "Sequence gradient (N → C)"],
  ["ambientOcclusion", "Ambient occlusion"],
  ["shadows", "Shadows"],
  ["outlines", "Outlines"],
  ["depthOfField", "Depth of field"],
];

function viewerHref(mode: "styled" | "nano") {
  const url = new URL(window.location.href);
  if (mode === "nano") url.searchParams.set("viewer", mode);
  else url.searchParams.delete("viewer");
  return `${url.pathname}${url.search}${url.hash}`;
}

function ViewerMode({ active }: { active: "styled" | "nano" }) {
  return (
    <nav className="viewer-mode" aria-label="Viewer">
      {(["styled", "nano"] as const).map((mode) => (
        <a
          aria-current={active === mode ? "page" : undefined}
          href={viewerHref(mode)}
          key={mode}
        >
          {mode === "styled" ? "Styled" : "Nano"}
        </a>
      ))}
    </nav>
  );
}

function NanoViewer({ selection }: { selection: StructureSelection }) {
  const parent = useRef<HTMLDivElement>(null);
  const viewer = useRef<MolstarViewerHandle | null>(null);
  const [busy, setBusy] = useState(true);
  const [error, setError] = useState<Error | null>(null);
  const [plddtAvailable, setPlddtAvailable] = useState(false);
  const [color, setColor] = useState<
    "chain" | "secondary" | "rainbow" | "element" | "residue" | "plddt"
  >("chain");
  const [illustrative, setIllustrative] = useState(true);
  const [surface, setSurface] = useState(true);
  const pathname = new URL(selection.url, window.location.href).pathname;

  useEffect(() => {
    const abort = new AbortController();
    let current = true;
    let loaded: MolstarViewerHandle | undefined;
    void Promise.all([
      import("@juliocesar-io/nano-protein-viewer-react"),
      import("molstar/lib/mol-plugin-ui/skin/dark.scss"),
      fetch(selection.url, { signal: abort.signal }),
    ])
      .then(async ([module, , response]) => {
        if (!response.ok)
          throw new Error(`Structure download failed (${response.status})`);
        const contents = await response.text();
        abort.signal.throwIfAborted();
        if (!parent.current) return;
        loaded = module.createMolstarViewer();
        await loaded.mount(parent.current);
        await loaded.loadStructureText(
          contents,
          pathname.endsWith(".pdb") ? "pdb" : "mmcif",
          false,
          true,
        );
        loaded.setBackgroundColor("dark");
        await loaded.updateColorTheme("chain");
        await loaded.applyIllustrativeStyle(true);
        await loaded.applySurface(true, { opacity: 8, inherit: true });
        if (!current) return;
        viewer.current = loaded;
        setPlddtAvailable(await loaded.isPLDDTAvailable());
        setBusy(false);
      })
      .catch((reason: unknown) => {
        if (!abort.signal.aborted && current) {
          setError(
            reason instanceof Error
              ? reason
              : new Error("Could not load Nano viewer"),
          );
          setBusy(false);
        }
      });
    return () => {
      current = false;
      abort.abort();
      viewer.current = null;
      if (loaded) void loaded.clear();
    };
  }, [pathname, selection.url]);

  function reportFailure(reason: unknown) {
    setError(
      reason instanceof Error
        ? reason
        : new Error("Could not update Nano viewer"),
    );
  }

  function changeColor(next: typeof color) {
    setColor(next);
    const active = viewer.current;
    if (!active) return;
    const update =
      next === "plddt"
        ? active.applyPLDDTTheme(true)
        : active
            .applyPLDDTTheme(false)
            .then(() => active.updateColorTheme(next));
    void update.catch(reportFailure);
  }

  return (
    <section className="molecule" aria-label="Molecular structure">
      <ErrorMessage error={error} />
      {error && (
        <a className="fallback-download" href={selection.url}>
          Download structure
        </a>
      )}
      <div className="molstar-viewport nano-viewport" aria-busy={busy}>
        <div className="nano-canvas" ref={parent} />
        <div className="viewer-topline nano-topline">
          <span className="viewer-hint">Drag to rotate · Scroll to zoom</span>
          <ViewerMode active="nano" />
        </div>
        <div className="nano-controls" aria-label="Nano appearance">
          <strong>Nano controls</strong>
          <label>
            Color
            <select
              value={color}
              onChange={(event) =>
                changeColor(event.target.value as typeof color)
              }
              disabled={busy}
            >
              <option value="chain">Chain</option>
              <option value="secondary">Secondary structure</option>
              <option value="rainbow">Sequence rainbow</option>
              <option value="element">Element</option>
              <option value="residue">Residue</option>
              {plddtAvailable && <option value="plddt">pLDDT</option>}
            </select>
          </label>
          <label className="check">
            <input
              type="checkbox"
              checked={illustrative}
              disabled={busy}
              onChange={(event) => {
                setIllustrative(event.target.checked);
                void viewer.current
                  ?.applyIllustrativeStyle(event.target.checked)
                  .catch(reportFailure);
              }}
            />
            Illustrative
          </label>
          <label className="check">
            <input
              type="checkbox"
              checked={surface}
              disabled={busy}
              onChange={(event) => {
                setSurface(event.target.checked);
                void viewer.current
                  ?.applySurface(event.target.checked, {
                    opacity: 8,
                    inherit: true,
                  })
                  .catch(reportFailure);
              }}
            />
            Molecular surface
          </label>
          <button
            disabled={busy}
            onClick={() => void viewer.current?.resetView()}
          >
            Reset view
          </button>
        </div>
        {busy && (
          <div className="viewer-loading">
            <Loading>Loading Nano viewer…</Loading>
          </div>
        )}
      </div>
    </section>
  );
}

function StyledViewer({ selection }: { selection: StructureSelection }) {
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
          <ViewerMode active="styled" />
          <Popover className="appearance-popover" label="Appearance" icon="◐">
            {appearanceToggles.map(([key, label]) => (
              <label className="check" key={key}>
                <input
                  type="checkbox"
                  checked={settings[key]}
                  onChange={(event) =>
                    setSettings({ ...settings, [key]: event.target.checked })
                  }
                />
                {label}
              </label>
            ))}
            <p className="viewer-info" role="status" data-testid="viewer-info">
              {info &&
                `${formatCount(info.atoms)} atoms · ${formatCount(info.residues)} residues · ${info.chains} ${info.chains === 1 ? "chain" : "chains"} · ${info.webgl}${info.surfacePath ? ` · ${info.surfacePath} surface computation` : ""}`}
            </p>
            <p className="muted">
              Visual effects and sequence colors do not change the structure or
              report prediction confidence.
            </p>
          </Popover>
        </div>
        {busy && (
          <div className="viewer-loading">
            <Loading>Loading structure…</Loading>
          </div>
        )}
        <div className="viewer-bottomline">
          <div
            className="chain-legend"
            aria-label={
              settings.sequenceGradient ? "Sequence colors" : "Chain colors"
            }
          >
            {settings.sequenceGradient ? (
              <span>
                <i className="sequence-gradient" />
                Sequence · N → C
              </span>
            ) : (
              info?.colors.map((entry) => (
                <span
                  key={
                    selection.chainRoles?.target.includes(entry.label)
                      ? `Target · ${entry.label}`
                      : selection.chainRoles?.binder.includes(entry.label)
                        ? `Key · ${entry.label}`
                        : entry.label
                  }
                >
                  <i style={{ background: entry.color }} />
                  {selection.chainRoles?.target.includes(entry.label)
                    ? `Target · ${entry.label}`
                    : selection.chainRoles?.binder.includes(entry.label)
                      ? `Key · ${entry.label}`
                      : entry.label}
                </span>
              ))
            )}
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

export default function MolecularViewer({
  selection,
}: {
  selection: StructureSelection;
}) {
  return new URLSearchParams(window.location.search).get("viewer") ===
    "nano" ? (
    <NanoViewer selection={selection} />
  ) : (
    <StyledViewer selection={selection} />
  );
}
