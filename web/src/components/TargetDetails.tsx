import { useState } from "react";
import { useMutation, useQuery } from "@tanstack/react-query";
import { api } from "../api";
import type { Job } from "../api";
import { ErrorMessage, Loading, formatCount } from "./Status";

export function TargetDetails({
  targetId,
  onPrepared,
}: {
  targetId: string;
  onPrepared: (job: Job) => void;
}) {
  const { data: target, error } = useQuery({
    queryKey: ["target", targetId],
    queryFn: ({ signal }) => api.target(targetId, signal),
  });
  const [length, setLength] = useState("75");
  const [count, setCount] = useState("2");
  const prepare = useMutation({
    mutationFn: api.prepare,
    onSuccess: onPrepared,
  });
  return (
    <>
      <ErrorMessage error={error} />
      <ErrorMessage error={prepare.error} />
      {!target && !error && <Loading>Loading complete sequence…</Loading>}
      {target && (
        <>
          <div className="job-line">
            <h3>{target.gene_symbol || target.gene_id}</h3>
            <span className="tag">
              {target.sequence_sha256 === target.reference_sequence_sha256
                ? "REFERENCE-MATCHING"
                : "NON-REFERENCE"}
            </span>
          </div>
          <dl>
            {Object.entries({
              Donor: target.donor,
              Transcript: target.transcript_id,
              "Target ID": target.target_id,
              "Phase scenario": target.phase_scenario,
              Length: `${formatCount(target.amino_acid_sequence.length)} amino acids`,
              Expression: "Not measured",
              "Sequence SHA-256": target.sequence_sha256,
            }).map(([key, value]) => (
              <div className="metadata-row" key={key}>
                <dt>{key}</dt>
                <dd>{value}</dd>
              </div>
            ))}
          </dl>
          <p className="field-label">Complete amino-acid sequence</p>
          <pre>{target.amino_acid_sequence}</pre>
          <details>
            <summary>Original target JSON</summary>
            <pre>{JSON.stringify(target, null, 2)}</pre>
          </details>
          <form
            className="settings"
            onSubmit={(event) => {
              event.preventDefault();
              prepare.mutate({
                target_id: targetId,
                binder_length: Number(length),
                num_designs: Number(count),
              });
            }}
          >
            <label>
              Binder length
              <input
                type="number"
                min={25}
                max={300}
                required
                value={length}
                onChange={(e) => setLength(e.target.value)}
              />
            </label>
            <label>
              Requested designs
              <input
                type="number"
                min={1}
                max={32}
                required
                value={count}
                onChange={(e) => setCount(e.target.value)}
              />
            </label>
            <button disabled={prepare.isPending}>
              {prepare.isPending ? "Preparing…" : "Prepare selected target"}
            </button>
          </form>
          <p className="muted">
            Preparation validates inputs. It does not start GPU inference.
          </p>
        </>
      )}
    </>
  );
}
