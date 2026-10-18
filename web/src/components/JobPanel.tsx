import { lazy, Suspense, useState } from "react";
import { useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import { api, isActive, isComplete } from "../api";
import type { StructureSelection } from "../api";
import { ViewerBoundary } from "./ViewerBoundary";
import { CandidateList } from "./CandidateList";
import { ErrorMessage, Loading } from "./Status";

const MolecularViewer = lazy(() => import("../viewer/MolecularViewer"));

export function JobPanel({
  jobId,
  allowRun,
}: {
  jobId: string;
  allowRun: boolean;
}) {
  const cache = useQueryClient();
  const [structure, setStructure] = useState<StructureSelection | null>(null);
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
  return (
    <>
      <ErrorMessage error={status.error} />
      <ErrorMessage error={results.error} />
      <ErrorMessage error={action.error} />
      {!job && status.isPending && <Loading>Loading job…</Loading>}
      {job && (
        <>
          <div className="job-line">
            <strong data-testid="job-status">
              {job.status.replaceAll("_", " ")}
            </strong>
            <span className="tag">
              {job.mode === "actual_cached_replay"
                ? "CACHED REPLAY · NO NEW INFERENCE"
                : job.execution_performed
                  ? "MODEL EXECUTION REPORTED"
                  : "NO MODEL EXECUTION"}
            </span>
          </div>
          <p className="muted job-target">
            {job.target_id
              ? `Target: ${job.target_id}`
              : "Public PD-L1 example — separate from the selected donor protein."}
          </p>
          {job.progress?.stage && <p>Current stage: {job.progress.stage}</p>}
          {job.error && (
            <div className="error" role="alert">
              {job.error}
            </div>
          )}
          <div className="actions">
            <button
              className="secondary"
              disabled={status.isFetching}
              onClick={() => {
                void status.refetch();
              }}
            >
              Refresh status
            </button>
            <button
              disabled={
                !allowRun ||
                job.status !== "prepared_pending_gpu_execution" ||
                !!job.dispatch_attempted ||
                action.isPending
              }
              onClick={() => {
                if (
                  window.confirm(
                    "Start a paid GPU design job for this prepared target?",
                  )
                )
                  action.mutate("run");
              }}
            >
              Start GPU job
            </button>
            <button
              className="secondary"
              disabled={
                !["queued", "running", "cancel_requested"].includes(
                  job.status,
                ) || action.isPending
              }
              onClick={() => action.mutate("cancel")}
            >
              Cancel job
            </button>
            {results.data && (
              <a href={results.data.download_url}>Download complete ZIP</a>
            )}
          </div>
          {results.isFetching && <Loading>Loading candidates…</Loading>}
          {results.data ? (
            <CandidateList
              candidates={results.data.candidates}
              onView={setStructure}
            />
          ) : (
            !results.isFetching && (
              <p className="empty">No candidate results for this job yet.</p>
            )
          )}
          {structure && (
            <ViewerBoundary key={structure.url} downloadUrl={structure.url}>
              <Suspense fallback={<Loading>Loading molecular viewer…</Loading>}>
                <MolecularViewer key={structure.url} selection={structure} />
              </Suspense>
            </ViewerBoundary>
          )}
        </>
      )}
    </>
  );
}
