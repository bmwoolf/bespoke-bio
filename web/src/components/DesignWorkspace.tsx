import { useQuery } from "@tanstack/react-query";
import { api } from "../api";
import { JobPanel } from "./JobPanel";
import { EmptyStructure } from "./Workspace";
import { ErrorMessage, Loading } from "./Status";

export function DesignWorkspace({ jobId }: { jobId?: string }) {
  const first = useQuery({
    queryKey: ["designs", "", 0],
    queryFn: ({ signal }) => api.designs("", 0, signal),
    enabled: !jobId,
  });
  const selectedId = jobId ?? first.data?.records[0]?.job_id;
  const design = useQuery({
    queryKey: ["design", selectedId],
    queryFn: ({ signal }) => api.design(selectedId!, signal),
    enabled: !!selectedId,
  });
  const data = design.data;
  return (
    <>
      <ErrorMessage error={first.error || design.error} />
      {data ? (
        <>
          <div className="workspace-heading">
            <div>
              <div className="eyebrow">{data.source}</div>
              <h1>{data.name}</h1>
              <p>
                {data.description}
                {data.length ? ` · ${data.length} aa` : ""}
              </p>
            </div>
            <span className="subtle-tag">
              {data.candidate_count} candidate keys
            </span>
          </div>
          <JobPanel
            key={data.job_id}
            jobId={data.job_id}
            expectedTargetId={data.target_id}
            allowRun={false}
            initialView="complex"
            proteinStructure={data.protein_structure}
          />
          {Object.keys(data.provenance).length > 0 && (
            <details className="workspace-details">
              <summary>Target provenance</summary>
              <dl>
                {[
                  ["Target state", data.provenance.state],
                  ["Sequence SHA-256", data.provenance.sequence_sha256],
                  ["Structure context", data.provenance.target_view],
                ]
                  .filter(([, value]) => typeof value === "string")
                  .map(([label, value]) => (
                    <div className="metadata-row" key={String(label)}>
                      <dt>{String(label)}</dt>
                      <dd>{String(value)}</dd>
                    </div>
                  ))}
              </dl>
              {typeof data.provenance.amino_acid_sequence === "string" && (
                <>
                  <h3>Target sequence · {data.length} aa</h3>
                  <pre>{data.provenance.amino_acid_sequence}</pre>
                </>
              )}
              {data.provenance.complete_study_url ===
                "https://github.com/bmwoolf/bespoke-bio/releases/tag/protein-results-2026-10-18" && (
                <p>
                  <a
                    href={data.provenance.complete_study_url}
                    target="_blank"
                    rel="noreferrer"
                  >
                    Complete study, controls and source archives ↗
                  </a>
                </p>
              )}
              <details>
                <summary>Complete source record</summary>
                <pre>{JSON.stringify(data.provenance, null, 2)}</pre>
              </details>
            </details>
          )}
        </>
      ) : first.data?.total === 0 && !jobId ? (
        <EmptyStructure title="No batch results published yet">
          <p>The Proteins tab contains the existing folded catalog.</p>
        </EmptyStructure>
      ) : !design.error && !first.error ? (
        <Loading>Loading designed target…</Loading>
      ) : null}
    </>
  );
}
