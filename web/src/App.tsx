import { useState } from "react";
import { useQuery, useQueryClient } from "@tanstack/react-query";
import { api } from "./api";
import type { Job } from "./api";
import { TargetCatalog } from "./components/TargetCatalog";
import { ExampleWorkspace, TargetWorkspace } from "./components/Workspace";
import { Dialog } from "./components/Dialog";
import { ErrorMessage, Loading, jobLabel } from "./components/Status";

type Selection =
  | { kind: "example" }
  | { kind: "target"; id: string; jobId?: string };

export function App() {
  const cache = useQueryClient();
  const [selection, setSelection] = useState<Selection>({ kind: "example" });
  const [dialog, setDialog] = useState<"proteins" | "history" | null>(null);
  const catalog = useQuery({
    queryKey: ["catalog"],
    queryFn: ({ signal }) => api.catalog(signal),
  });
  const jobs = useQuery({
    queryKey: ["jobs"],
    queryFn: ({ signal }) => api.jobs(signal),
  });
  const meta = catalog.data;
  const selectedId = selection.kind === "target" ? selection.id : null;
  function selectTarget(id: string) {
    setSelection({ kind: "target", id });
    setDialog(null);
  }
  function selectJob(job: Job) {
    cache.setQueryData(["job", job.job_id], job);
    setSelection(
      job.target_id
        ? { kind: "target", id: job.target_id, jobId: job.job_id }
        : { kind: "example" },
    );
    setDialog(null);
  }
  return (
    <div className="app-shell">
      <header className="app-header">
        <a className="brand" href="/">
          bespoke<span>bio</span>
          <span className="brand-dot" />
        </a>
        <nav aria-label="Workspace">
          <button
            className="quiet mobile-browse"
            onClick={() => setDialog("proteins")}
          >
            Proteins
          </button>
          <button className="quiet" onClick={() => setDialog("history")}>
            History
          </button>
          <button
            className="quiet"
            aria-pressed={selection.kind === "example"}
            onClick={() => setSelection({ kind: "example" })}
          >
            Example
          </button>
        </nav>
      </header>
      <div className="app-body">
        <aside className="desktop-catalog">
          <TargetCatalog
            selectedId={selectedId}
            onSelect={selectTarget}
            catalog={meta}
          />
        </aside>
        <main className="workspace" id="workspace">
          <ErrorMessage error={catalog.error} />
          {catalog.isPending && <Loading>Loading workspace…</Loading>}
          {meta &&
            (selection.kind === "example" ? (
              <ExampleWorkspace
                configured={meta.service_configured}
                saved={jobs.data?.find((job) => job.target_id === null)}
                historyReady={!jobs.isPending}
              />
            ) : (
              <TargetWorkspace
                key={`${selection.id}:${selection.jobId ?? "new"}`}
                targetId={selection.id}
                initialJobId={selection.jobId}
                catalog={meta}
              />
            ))}
        </main>
      </div>
      <Dialog
        open={dialog === "proteins"}
        title="Proteins"
        onClose={() => setDialog(null)}
      >
        {dialog === "proteins" && (
          <TargetCatalog
            selectedId={selectedId}
            onSelect={selectTarget}
            catalog={meta}
            compact
          />
        )}
      </Dialog>
      <Dialog
        open={dialog === "history"}
        title="History"
        onClose={() => setDialog(null)}
      >
        <ErrorMessage error={jobs.error} />
        {jobs.isPending && <Loading />}
        {jobs.data?.length === 0 && (
          <p className="muted">
            Your prepared targets and saved results will appear here.
          </p>
        )}
        <div className="history-list">
          {jobs.data?.map((job) => (
            <button
              className="history-item"
              key={job.job_id}
              onClick={() => selectJob(job)}
            >
              <strong>
                {job.target_id
                  ? `Protein · ${job.target_id.slice(-8)}`
                  : "PD-L1 · Saved example"}
              </strong>
              <span>{jobLabel(job.status)}</span>
              <small>
                {job.target_id ?? "Previously generated candidates"}
              </small>
            </button>
          ))}
        </div>
      </Dialog>
    </div>
  );
}
