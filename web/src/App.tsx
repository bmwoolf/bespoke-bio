import { useEffect, useState } from "react";
import { useQuery, useQueryClient } from "@tanstack/react-query";
import { api, isActive, isComplete } from "./api";
import type { Job } from "./api";
import { TargetCatalog } from "./components/TargetCatalog";
import { ExampleWorkspace, TargetWorkspace } from "./components/Workspace";
import { Dialog } from "./components/Dialog";
import { ErrorMessage, Loading, jobLabel } from "./components/Status";

type Selection =
  | { kind: "example" }
  | { kind: "target"; id: string; jobId?: string };

function readSelection(): Selection {
  const params = new URLSearchParams(location.search);
  const id = params.get("target");
  return id
    ? { kind: "target", id, jobId: params.get("job") ?? undefined }
    : { kind: "example" };
}

export function App() {
  const cache = useQueryClient();
  const [selection, setSelection] = useState<Selection>(readSelection);
  function navigate(next: Selection) {
    const url = new URL(location.href);
    url.search = "";
    if (next.kind === "target") {
      url.searchParams.set("target", next.id);
      if (next.jobId) url.searchParams.set("job", next.jobId);
    }
    if (url.href !== location.href) history.pushState(null, "", url);
    setSelection(next);
  }
  useEffect(() => {
    const back = () => setSelection(readSelection());
    addEventListener("popstate", back);
    return () => removeEventListener("popstate", back);
  }, []);
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
  const savedJob =
    jobs.data?.find(
      (job) => job.target_id === selectedId && isComplete(job.status),
    ) ??
    jobs.data?.find(
      (job) =>
        job.target_id === selectedId &&
        (isActive(job.status) ||
          job.status === "prepared_pending_gpu_execution"),
    );
  const selectedJobId =
    selection.kind === "target"
      ? (selection.jobId ?? savedJob?.job_id)
      : undefined;
  const modeledTargets = new Set(
    jobs.data
      ?.filter((job) => job.target_id && isComplete(job.status))
      .map((job) => job.target_id as string),
  );
  function selectTarget(id: string) {
    navigate({ kind: "target", id });
    setDialog(null);
  }
  function selectJob(job: Job) {
    cache.setQueryData(["job", job.job_id], job);
    navigate(
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
            onClick={() => navigate({ kind: "example" })}
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
            modeledTargets={modeledTargets}
          />
        </aside>
        <main className="workspace" id="workspace">
          <ErrorMessage error={catalog.error} />
          {catalog.isPending && <Loading>Loading workspace…</Loading>}
          {meta &&
            (selection.kind === "example" ? (
              <ExampleWorkspace
                configured={meta.service_configured}
                readOnly={meta.read_only}
                saved={jobs.data?.find((job) => job.target_id === null)}
                historyReady={!jobs.isPending}
              />
            ) : (
              <TargetWorkspace
                key={`${selection.id}:${selectedJobId ?? "new"}`}
                targetId={selection.id}
                initialJobId={selectedJobId}
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
            modeledTargets={modeledTargets}
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
