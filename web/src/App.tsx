import { useState } from "react";
import { useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import { api } from "./api";
import type { Job } from "./api";
import { TargetCatalog } from "./components/TargetCatalog";
import { TargetDetails } from "./components/TargetDetails";
import { JobPanel } from "./components/JobPanel";
import { ErrorMessage, Loading, formatCount } from "./components/Status";

export function App() {
  const cache = useQueryClient();
  const [targetId, setTargetId] = useState<string | null>(null);
  const [jobId, setJobId] = useState<string | null>(null);
  const catalog = useQuery({
    queryKey: ["catalog"],
    queryFn: ({ signal }) => api.catalog(signal),
  });
  const jobs = useQuery({
    queryKey: ["jobs"],
    queryFn: ({ signal }) => api.jobs(signal),
  });
  function selectJob(job: Job) {
    cache.setQueryData(["job", job.job_id], job);
    setJobId(job.job_id);
    void cache.invalidateQueries({ queryKey: ["jobs"] });
  }
  const demo = useMutation({ mutationFn: api.demo, onSuccess: selectJob });
  const meta = catalog.data;
  return (
    <>
      <header>
        <a className="brand" href="/">
          BESPOKE <span>BIO</span>
        </a>
        <span className="tag">RESEARCH DEMO</span>
      </header>
      <main>
        <div className="intro">
          <div>
            <p className="eyebrow">FROM DONOR SEQUENCE TO MOLECULAR DESIGN</p>
            <h1>Protein → candidate</h1>
            <p>
              Explore the complete protein catalog. Prepare a target. Inspect
              candidate designs.
            </p>
          </div>
          <div className="mode">
            {meta?.gpu_execution_enabled
              ? "GPU inference is available through an explicit start action."
              : "Integration mode · GPU inference disabled"}
          </div>
        </div>
        <ErrorMessage error={catalog.error} />
        {catalog.isPending && <Loading>Loading catalog…</Loading>}
        {meta && (
          <>
            {!meta.service_configured && (
              <p className="error" role="alert">
                Catalog ready. Configure the design service to prepare targets
                or inspect cached results.
              </p>
            )}
            <section className="stats" aria-label="Catalog totals">
              <div>
                <strong>{formatCount(meta.target_records)}</strong>
                <span>biological records</span>
              </div>
              <div>
                <strong>{formatCount(meta.unique_sequences)}</strong>
                <span>unique sequences</span>
              </div>
              <div>
                <strong>{formatCount(meta.nonreference_records)}</strong>
                <span>non-reference records</span>
              </div>
              <div>
                <strong>Intact</strong>
                <span>original sequences & provenance</span>
              </div>
            </section>
            <div className="workspace">
              <TargetCatalog selectedId={targetId} onSelect={setTargetId} />
              <div className="detail-column">
                <section className="panel">
                  <div className="panel-title">
                    <span className="step">02</span>
                    <h2>Prepare target</h2>
                  </div>
                  {targetId ? (
                    <TargetDetails
                      key={targetId}
                      targetId={targetId}
                      onPrepared={selectJob}
                    />
                  ) : (
                    <p className="empty">
                      Select a protein to inspect its exact sequence and donor
                      evidence.
                    </p>
                  )}
                </section>
                <section className="panel">
                  <div className="panel-title">
                    <span className="step">03</span>
                    <h2>Inspect results</h2>
                    <button
                      className="secondary"
                      disabled={demo.isPending || !meta.service_configured}
                      onClick={() => demo.mutate()}
                    >
                      {demo.isPending
                        ? "Loading example…"
                        : "Load cached PD-L1 example"}
                    </button>
                  </div>
                  <p className="muted">
                    Candidates are computational predictions. Filter passes and
                    confidence scores do not establish measured binding.
                  </p>
                  <ErrorMessage error={demo.error} />
                  {jobId ? (
                    <JobPanel
                      key={jobId}
                      jobId={jobId}
                      allowRun={meta.gpu_execution_enabled}
                    />
                  ) : (
                    <p className="empty">
                      Prepared jobs appear here. The cached example shows saved
                      model outputs without running inference.
                    </p>
                  )}
                </section>
                <section className="panel history">
                  <h2>Local job history</h2>
                  <ErrorMessage error={jobs.error} />
                  <div className="history-list">
                    {jobs.data?.length === 0 && (
                      <p className="muted">No jobs yet.</p>
                    )}
                    {jobs.data?.map((job) => (
                      <button
                        className="secondary"
                        key={job.job_id}
                        onClick={() => selectJob(job)}
                      >
                        {job.target_id || "Cached PD-L1 example"} ·{" "}
                        {job.status.replaceAll("_", " ")}
                      </button>
                    ))}
                  </div>
                </section>
              </div>
            </div>
          </>
        )}
        <footer>
          <span>
            Full Stage 1 output stays on the server. Only selected records are
            sent to the configured design service.
          </span>
          <a href="https://molstar.org/" target="_blank" rel="noreferrer">
            Molecular viewer: Mol*
          </a>
        </footer>
      </main>
    </>
  );
}
