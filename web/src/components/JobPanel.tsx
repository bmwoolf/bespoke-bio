import { lazy, Suspense, useState } from "react";
import { useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import { api, isActive, isComplete } from "../api";
import type { ProteinStructure, Target } from "../api";
import { ProteinView } from "./ProteinView";
import { Popover } from "./Popover";
import { ViewerBoundary } from "./ViewerBoundary";
import { CandidateList, CandidateDetails } from "./CandidateList";
import { TargetDetails } from "./TargetDetails";
import { ErrorMessage, Loading, jobLabel } from "./Status";

const MolecularViewer = lazy(() => import("../viewer/MolecularViewer"));

export function JobPanel({
  jobId,
  allowRun,
  target,
  proteinStructure,
}: {
  jobId: string;
  allowRun: boolean;
  target?: Target;
  proteinStructure?: ProteinStructure | null;
}) {
  const cache = useQueryClient();
  const [candidateId, setCandidateId] = useState<string | null>(null);
  const [structureKind, setStructureKind] = useState<
    "complex" | "binder" | "protein"
  >("complex");
  const status = useQuery({
    queryKey: ["job", jobId],
    queryFn: ({ signal }) => api.job(jobId, signal),
    staleTime: 0,
    refetchInterval: (query) =>
      isActive(query.state.data?.status ?? "") ? 2500 : false,
  });
  const results = useQuery({
    queryKey: ["results", jobId],
    queryFn: ({ signal }) => api.results(jobId, signal),
    enabled: isComplete(status.data?.status ?? ""),
  });
  const action = useMutation({
    mutationFn: (name: "run" | "cancel") => api[name](jobId),
    onSuccess: (job) => {
      cache.setQueryData(["job", jobId], job);
      void cache.invalidateQueries({ queryKey: ["jobs"] });
    },
  });
  const job = status.data;
  // A shared URL must never display another target's job or result under this heading.
  if (
    (job && job.target_id !== (target?.target_id ?? null)) ||
    (results.data && results.data.target_id !== (target?.target_id ?? null))
  ) {
    return (
      <div className="error" role="alert">
        This result belongs to a different protein. Choose it from History.
      </div>
    );
  }
  const candidate =
    results.data?.candidates.find((value) => value.id === candidateId) ??
    results.data?.candidates[0];
  const structure =
    candidate && structureKind !== "protein"
      ? {
          url:
            structureKind === "complex"
              ? candidate.complex_structure_url
              : candidate.binder_structure_url,
          title: `${structureKind === "complex" ? "Complex" : "Binder"} · ${candidate.id}`,
        }
      : null;
  const ready = job?.status === "prepared_pending_gpu_execution";
  const canCancel =
    job && ["queued", "running", "cancel_requested"].includes(job.status);
  return (
    <>
      <ErrorMessage error={status.error} />
      <ErrorMessage error={results.error} />
      <ErrorMessage error={action.error} />
      {status.error && (
        <button className="quiet" onClick={() => void status.refetch()}>
          Refresh status
        </button>
      )}
      {results.error && (
        <button className="quiet" onClick={() => void results.refetch()}>
          Reload results
        </button>
      )}
      {!job && status.isPending && <Loading>Loading job…</Loading>}
      {job && (
        <>
          {job.error && (
            <div className="error" role="alert">
              {job.error}
            </div>
          )}
          <div className="result-toolbar">
            {candidate ? (
              <label className="structure-picker">
                <span className="sr-only">Structure</span>
                <select
                  value={structureKind}
                  onChange={(event) =>
                    setStructureKind(
                      event.target.value === "binder"
                        ? "binder"
                        : event.target.value === "protein"
                          ? "protein"
                          : "complex",
                    )
                  }
                >
                  <option value="complex">Complex</option>
                  <option value="binder">Key</option>
                  {proteinStructure && <option value="protein">Protein</option>}
                </select>
              </label>
            ) : (
              <span className="status-label" data-testid="job-status">
                {results.data && !results.data.candidates.length
                  ? "No keys returned"
                  : jobLabel(job.status)}
              </span>
            )}
            <div className="result-actions">
              {ready && (
                <>
                  <span className="muted">
                    {!allowRun
                      ? "Generation unavailable in this demo"
                      : "Generation uses paid compute"}
                  </span>
                  <button
                    disabled={
                      !allowRun || !!job.dispatch_attempted || action.isPending
                    }
                    onClick={() => {
                      if (
                        window.confirm(
                          "Generate keys for this protein using paid GPU compute?",
                        )
                      )
                        action.mutate("run");
                    }}
                  >
                    Generate keys
                  </button>
                </>
              )}
              {canCancel && (
                <button
                  className="quiet"
                  disabled={
                    job.status === "cancel_requested" || action.isPending
                  }
                  onClick={() => action.mutate("cancel")}
                >
                  Cancel
                </button>
              )}
              {results.data && (
                <Popover className="download-popover" label="Download">
                  <a href={results.data.download_url}>All results (.zip)</a>
                  {candidate && (
                    <>
                      <a href={candidate.complex_structure_url}>
                        Complex structure
                      </a>
                      <a href={candidate.binder_structure_url}>
                        Binder structure
                      </a>
                    </>
                  )}
                </Popover>
              )}
            </div>
          </div>
          {isActive(job.status) && (
            <Loading>
              {job.progress?.stage?.replaceAll("_", " ") ||
                jobLabel(job.status)}
            </Loading>
          )}
          {structure ? (
            <ViewerBoundary key={structure.url} downloadUrl={structure.url}>
              <Suspense
                fallback={
                  <div className="empty-structure">
                    <Loading>Loading molecular viewer…</Loading>
                  </div>
                }
              >
                <MolecularViewer key={structure.url} selection={structure} />
              </Suspense>
            </ViewerBoundary>
          ) : proteinStructure ? (
            <ProteinView
              structure={proteinStructure}
              title={target?.gene_symbol || "Target"}
            />
          ) : (
            <div className="empty-structure">
              <div>
                <span className="empty-mark" aria-hidden="true">
                  ◇
                </span>
                <h2>
                  {isActive(job.status)
                    ? jobLabel(job.status)
                    : ready
                      ? "Protein ready"
                      : results.data
                        ? "No candidates returned"
                        : "No structure available"}
                </h2>
                <p>
                  {results.isFetching
                    ? "Loading candidates…"
                    : ready
                      ? "Sequence validated. Generate binders to create structures."
                      : job.progress?.stage
                        ? job.progress.stage.replaceAll("_", " ")
                        : "Structures appear here when results are available."}
                </p>
              </div>
            </div>
          )}
          {candidate && results.data && (
            <CandidateList
              candidates={results.data.candidates}
              selectedId={candidate.id}
              onSelect={setCandidateId}
            />
          )}
          <div className="workspace-footer">
            <details className="workspace-details">
              <summary>Details</summary>
              {candidate && <CandidateDetails candidate={candidate} />}
              {target && <TargetDetails target={target} />}
              <details>
                <summary>Job record</summary>
                <pre>{JSON.stringify(job, null, 2)}</pre>
              </details>
            </details>
            <span className="research-note">Research prediction</span>
          </div>
        </>
      )}
    </>
  );
}
