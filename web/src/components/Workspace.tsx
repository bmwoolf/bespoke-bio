import { useState } from "react";
import { useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import { api } from "../api";
import type { Catalog, Job } from "../api";
import { Popover } from "./Popover";
import { JobPanel } from "./JobPanel";
import { TargetDetails } from "./TargetDetails";
import { ErrorMessage, Loading, formatCount } from "./Status";
import { ProteinView } from "./ProteinView";

export function EmptyStructure({
  title = "No structure yet",
  children,
}: {
  title?: string;
  children?: React.ReactNode;
}) {
  return (
    <div className="empty-structure">
      <div>
        <span className="empty-mark" aria-hidden="true">
          ◇
        </span>
        <h2>{title}</h2>
        {children}
      </div>
    </div>
  );
}

export function ExampleWorkspace({
  configured,
  readOnly,
  saved,
  historyReady,
}: {
  configured: boolean;
  readOnly: boolean;
  saved?: Job;
  historyReady: boolean;
}) {
  // Cached replay is the only automatic POST. It cannot start inference, and the gateway
  // reuses its receipt. Disable all automatic retries/refetches, including Strict Mode remounts.
  const cache = useQueryClient();
  const example = useQuery({
    queryKey: ["saved-example"],
    queryFn: async () => {
      if (saved) return saved;
      const job = await api.demo();
      void cache.invalidateQueries({ queryKey: ["jobs"] });
      return job;
    },
    enabled: historyReady && (!!saved || (configured && !readOnly)),
    retry: false,
    retryOnMount: false,
    staleTime: Infinity,
    refetchOnMount: false,
    refetchOnWindowFocus: false,
    refetchOnReconnect: false,
  });
  return (
    <>
      <div className="workspace-heading">
        <div>
          <div className="eyebrow">Saved example</div>
          <h1>PD-L1</h1>
          <p>Previously generated binder candidates</p>
        </div>
        <span className="subtle-tag">Example data</span>
      </div>
      <ErrorMessage error={example.error} />
      {historyReady && !saved && (!configured || readOnly) ? (
        <EmptyStructure title="Choose a protein">
          <p>
            The catalog is available. No saved example has been published here
            yet.
          </p>
        </EmptyStructure>
      ) : example.data ? (
        <JobPanel
          key={example.data.job_id}
          jobId={example.data.job_id}
          allowRun={false}
        />
      ) : !example.error ? (
        <EmptyStructure title="Opening example">
          <Loading>Loading saved structures…</Loading>
        </EmptyStructure>
      ) : null}
    </>
  );
}

export function TargetWorkspace({
  targetId,
  initialJobId,
  catalog,
}: {
  targetId: string;
  initialJobId?: string;
  catalog: Catalog;
}) {
  const cache = useQueryClient();
  const [jobId, setJobId] = useState(initialJobId);
  const [length, setLength] = useState("75");
  const [count, setCount] = useState("2");
  const target = useQuery({
    queryKey: ["target", targetId],
    queryFn: ({ signal }) => api.target(targetId, signal),
  });
  const structure = useQuery({
    queryKey: ["protein-structure", targetId],
    queryFn: ({ signal }) => api.structure(targetId, signal),
    staleTime: Infinity,
  });
  const prepare = useMutation({
    mutationFn: async (selection: Parameters<typeof api.prepare>[0]) => {
      const job = await api.prepare(selection);
      return catalog.gpu_execution_enabled &&
        job.status === "prepared_pending_gpu_execution" &&
        !job.dispatch_attempted
        ? api.run(job.job_id)
        : job;
    },
    onSuccess: (job) => {
      cache.setQueryData(["job", job.job_id], job);
      void cache.invalidateQueries({ queryKey: ["jobs"] });
    },
    onError: () => void cache.invalidateQueries({ queryKey: ["jobs"] }),
  });
  const data = target.data;
  return (
    <>
      <ErrorMessage error={target.error} />
      <ErrorMessage error={prepare.error} />
      <ErrorMessage error={structure.error} />
      {data &&
        structure.data &&
        structure.data.sequence_sha256 !== data.sequence_sha256 && (
          <div className="error" role="alert">
            The saved structure belongs to a different sequence.
          </div>
        )}
      {!data && target.isPending && <Loading>Loading protein…</Loading>}
      {data && (
        <>
          <div className="workspace-heading">
            <div>
              <div className="eyebrow">{data.donor}</div>
              <h1>{data.gene_symbol || data.gene_id}</h1>
              <p>
                {data.transcript_id} ·{" "}
                {formatCount(data.amino_acid_sequence.length)} aa ·{" "}
                <span title={data.target_id}>{data.target_id.slice(-8)}</span>
              </p>
            </div>
            {!jobId && !catalog.read_only && (
              <form
                className="prepare-actions"
                onInvalidCapture={(event) => {
                  const settings = event.currentTarget.querySelector("details");
                  if (settings) settings.open = true;
                }}
                onSubmit={(event) => {
                  event.preventDefault();
                  if (
                    catalog.gpu_execution_enabled &&
                    !window.confirm(
                      "Generate keys for this protein using paid GPU compute?",
                    )
                  )
                    return;
                  // Per-call success only runs while this observer is mounted; a late response
                  // cannot replace a different selected protein or the example.
                  prepare.mutate(
                    {
                      target_id: targetId,
                      binder_length: Number(length),
                      num_designs: Number(count),
                    },
                    { onSuccess: (job) => setJobId(job.job_id) },
                  );
                }}
              >
                <Popover className="settings-popover" label="Settings">
                  <label>
                    Binder length
                    <input
                      type="number"
                      min={25}
                      max={300}
                      required
                      value={length}
                      disabled={prepare.isPending}
                      onChange={(event) => setLength(event.target.value)}
                    />
                  </label>
                  <label>
                    Candidates
                    <input
                      type="number"
                      min={1}
                      max={32}
                      required
                      value={count}
                      disabled={prepare.isPending}
                      onChange={(event) => setCount(event.target.value)}
                    />
                  </label>
                  <p className="muted">
                    Small demo defaults; adjust before preparing.
                  </p>
                </Popover>
                <button
                  disabled={prepare.isPending || !catalog.service_configured}
                >
                  {prepare.isPending
                    ? "Submitting…"
                    : catalog.gpu_execution_enabled
                      ? "Generate keys"
                      : "Prepare target"}
                </button>
              </form>
            )}
          </div>
          {!catalog.service_configured && !catalog.read_only && (
            <p className="notice">
              Connect a design service to prepare this protein.
            </p>
          )}
          {jobId ? (
            <JobPanel
              key={jobId}
              jobId={jobId}
              allowRun={catalog.gpu_execution_enabled}
              target={data}
              proteinStructure={
                structure.data?.sequence_sha256 === data.sequence_sha256
                  ? structure.data
                  : null
              }
            />
          ) : (
            <>
              {structure.data?.sequence_sha256 === data.sequence_sha256 ? (
                <ProteinView
                  structure={structure.data}
                  title={data.gene_symbol || data.gene_id}
                />
              ) : structure.isPending ? (
                <Loading>Loading protein structure…</Loading>
              ) : (
                <EmptyStructure title="Structure pending">
                  <p>
                    The complete sequence is available. Its structure has not
                    been published yet.
                    <br />
                    {catalog.read_only
                      ? "New structures must be generated by an operator."
                      : "Prepare it, then generate binders to create structures."}
                  </p>
                  <a className="quiet" href="/">
                    View the saved 3D example
                  </a>
                </EmptyStructure>
              )}
              <details className="workspace-details">
                <summary>Details</summary>
                <TargetDetails target={data} />
              </details>
            </>
          )}
        </>
      )}
    </>
  );
}
