import { useEffect, useState } from "react";
import { useQuery, useQueryClient } from "@tanstack/react-query";
import { api, appPath, isActive, isComplete, isOperator } from "./api";
import type { Job } from "./api";
import { DesignCatalog } from "./components/DesignCatalog";
import { DesignWorkspace } from "./components/DesignWorkspace";
import { TargetCatalog } from "./components/TargetCatalog";
import { ExampleWorkspace, TargetWorkspace } from "./components/Workspace";
import { Dialog } from "./components/Dialog";
import { ErrorMessage, Loading, jobLabel } from "./components/Status";

type Selection =
  | { kind: "example" }
  | { kind: "design"; jobId?: string }
  | { kind: "target"; id: string; jobId?: string };

function readSelection(): Selection {
  const params = new URLSearchParams(location.search);
  const id = params.get("target");
  if (params.has("design"))
    return { kind: "design", jobId: params.get("design") || undefined };
  if (params.has("example") || isOperator())
    return id
      ? { kind: "target", id, jobId: params.get("job") ?? undefined }
      : { kind: "example" };
  return id
    ? { kind: "target", id, jobId: params.get("job") ?? undefined }
    : { kind: "design" };
}

export function App() {
  const cache = useQueryClient();
  const [selection, setSelection] = useState<Selection>(readSelection);
  const [collection, setCollection] = useState<"designs" | "proteins">(() =>
    readSelection().kind === "target" || isOperator() ? "proteins" : "designs",
  );
  function navigate(next: Selection) {
    const url = new URL(location.href);
    url.search = "";
    if (next.kind === "design" && next.jobId)
      url.searchParams.set("design", next.jobId);
    if (next.kind === "example") url.searchParams.set("example", "1");
    if (next.kind === "target") {
      url.searchParams.set("target", next.id);
      if (next.jobId) url.searchParams.set("job", next.jobId);
    }
    if (url.href !== location.href) history.pushState(null, "", url);
    setSelection(next);
    if (next.kind !== "example")
      setCollection(next.kind === "target" ? "proteins" : "designs");
  }
  useEffect(() => {
    const back = () => {
      const next = readSelection();
      setSelection(next);
      setCollection(
        next.kind === "target" || isOperator() ? "proteins" : "designs",
      );
    };
    addEventListener("popstate", back);
    return () => removeEventListener("popstate", back);
  }, []);
  const [dialog, setDialog] = useState<"proteins" | "history" | null>(null);
  const catalog = useQuery({
    queryKey: ["catalog"],
    queryFn: ({ signal }) => api.catalog(signal),
    refetchInterval: isOperator() ? 30_000 : false,
  });
  const jobs = useQuery({
    queryKey: ["jobs"],
    queryFn: ({ signal }) => api.jobs(signal),
  });
  const firstDesigns = useQuery({
    queryKey: ["designs", "", 0],
    queryFn: ({ signal }) => api.designs("", 0, signal),
    enabled: !isOperator(),
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
  function selectDesign(jobId: string) {
    navigate({ kind: "design", jobId });
    setDialog(null);
  }
  function selectJob(job: Job) {
    if (meta?.read_only) {
      selectDesign(job.job_id);
      return;
    }
    cache.setQueryData(["job", job.job_id], job);
    navigate(
      job.target_id
        ? { kind: "target", id: job.target_id, jobId: job.job_id }
        : { kind: "example" },
    );
    setDialog(null);
  }
  function browse(compact = false) {
    return (
      <div className="catalog-browser">
        {!isOperator() && (
          <div
            className="catalog-tabs"
            role="group"
            aria-label="Catalog collection"
          >
            <button
              aria-pressed={collection === "designs"}
              onClick={() => setCollection("designs")}
            >
              Saved keys
            </button>
            <button
              aria-pressed={collection === "proteins"}
              onClick={() => setCollection("proteins")}
            >
              Proteins
            </button>
          </div>
        )}
        {collection === "designs" ? (
          <DesignCatalog
            selectedId={
              selection.kind === "design"
                ? (selection.jobId ??
                  firstDesigns.data?.records[0]?.job_id ??
                  null)
                : null
            }
            onSelect={selectDesign}
          />
        ) : (
          <TargetCatalog
            selectedId={selectedId}
            onSelect={selectTarget}
            catalog={meta}
            modeledTargets={modeledTargets}
            compact={compact}
          />
        )}
      </div>
    );
  }
  return (
    <div className="app-shell">
      <header className="app-header">
        <a className="brand" href={appPath("/")}>
          bespoke<span>bio</span>
          <span className="brand-dot" />
        </a>
        <nav aria-label="Workspace">
          {isOperator() && <span className="muted">Operator</span>}
          <button
            className="quiet mobile-browse"
            onClick={() => setDialog("proteins")}
          >
            Browse
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
        <aside className="desktop-catalog">{browse()}</aside>
        <main className="workspace" id="workspace">
          <ErrorMessage error={catalog.error} />
          {catalog.isPending && <Loading>Loading workspace…</Loading>}
          {meta &&
            (selection.kind === "design" ? (
              <DesignWorkspace
                key={selection.jobId ?? "first"}
                jobId={selection.jobId}
              />
            ) : selection.kind === "example" ? (
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
        title="Browse"
        onClose={() => setDialog(null)}
      >
        {dialog === "proteins" && browse(true)}
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
