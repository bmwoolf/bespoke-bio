import { useEffect, useRef, useState } from "react";
import type { MolstarViewerHandle } from "@juliocesar-io/nano-protein-viewer-react";
import type { StructureSelection } from "../api";
import { ErrorMessage, Loading } from "../components/Status";

function NanoViewer({ selection }: { selection: StructureSelection }) {
  const parent = useRef<HTMLDivElement>(null);
  const [busy, setBusy] = useState(true);
  const [error, setError] = useState<Error | null>(null);
  const pathname = new URL(selection.url, window.location.href).pathname;
  const targetChains = selection.chainRoles?.target.join(",") ?? "";
  const keyChains = selection.chainRoles?.binder.join(",") ?? "";

  useEffect(() => {
    let current = true;
    let loaded: MolstarViewerHandle | undefined;
    setBusy(true);
    setError(null);
    void Promise.all([
      import("@juliocesar-io/nano-protein-viewer-react"),
      import("molstar/lib/mol-plugin-ui/skin/dark.scss"),
      fetch(selection.url),
    ])
      .then(async ([module, , response]) => {
        if (!response.ok)
          throw new Error(`Structure download failed (${response.status})`);
        const contents = await response.text();
        if (!current) return;
        if (!parent.current) return;
        loaded = module.createMolstarViewer();
        await loaded.mount(parent.current);
        await loaded.loadStructureText(
          contents,
          pathname.endsWith(".pdb") ? "pdb" : "mmcif",
          false,
          true,
        );
        loaded.setBackgroundColor("light");
        const chains = await loaded.listChains();
        const fallbackColor = selection.title.toLowerCase().startsWith("key")
          ? "#F47A1F"
          : "#35B99A";
        const targetIds = new Set(targetChains.split(",").filter(Boolean));
        const keyIds = new Set(keyChains.split(",").filter(Boolean));
        const hasChainRoles = targetIds.size > 0 || keyIds.size > 0;
        const chainColors = Object.fromEntries(
          chains.map((chain) => [
            chain,
            hasChainRoles
              ? targetIds.has(chain)
                ? "#35B99A"
                : keyIds.has(chain)
                  ? "#F47A1F"
                  : "#78909C"
              : fallbackColor,
          ]),
        );
        await loaded.updateColorTheme("chain", { chainColors });
        await loaded.applyIllustrativeStyle(false);
        await loaded.applySurface(true, { opacity: 70, inherit: true });
        await loaded.resetView();
        if (!current) return;
        setBusy(false);
      })
      .catch((reason: unknown) => {
        if (current) {
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
      if (loaded) void loaded.clear().catch(() => undefined);
    };
  }, [keyChains, pathname, selection.title, selection.url, targetChains]);

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
        {busy && (
          <div className="viewer-loading">
            <Loading>Loading Nano viewer…</Loading>
          </div>
        )}
      </div>
    </section>
  );
}

export default function MolecularViewer({
  selection,
}: {
  selection: StructureSelection;
}) {
  return <NanoViewer selection={selection} />;
}
