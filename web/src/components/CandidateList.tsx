import type { Candidate } from "../api";

const metricLabels: Record<string, string> = {
  design_ptm: "Key folding confidence (pTM)",
  design_to_target_iptm: "Key → target interface (ipTM)",
  esmfold2_ipsae_min: "Interface confidence (minimum ipSAE)",
  filter_rmsd: "Native filter RMSD",
  filter_rmsd_design: "Key filter RMSD",
  "designfolding-filter_rmsd": "Independent key fold RMSD",
  quality_score: "Within-batch ranking score",
  N_target_clearance_angstrom: "N-terminus clearance",
  C_target_clearance_angstrom: "C-terminus clearance",
};

export function filterLabel(candidate: Candidate) {
  return candidate.native_filter_passed === true
    ? "Passed"
    : candidate.native_filter_passed === false
      ? "Failed"
      : "Not assessed";
}
export function CandidateList({
  candidates,
  selectedId,
  onSelect,
}: {
  candidates: Candidate[];
  selectedId: string;
  onSelect: (id: string) => void;
}) {
  return (
    <div className="candidate-strip" aria-label="Candidates">
      <span className="strip-label">
        Candidates <span>{candidates.length}</span>
      </span>
      <div className="candidate-options">
        {candidates.map((candidate, index) => (
          <button
            className="candidate-option"
            key={candidate.id}
            aria-pressed={candidate.id === selectedId}
            onClick={() => onSelect(candidate.id)}
            title={candidate.id}
          >
            <span>{String(index + 1).padStart(2, "0")}</span>
            <span
              className={
                candidate.native_filter_passed === true
                  ? "pass"
                  : candidate.native_filter_passed === false
                    ? "fail"
                    : "muted"
              }
            >
              {filterLabel(candidate)}
            </span>
          </button>
        ))}
      </div>
      <span className="strip-note">Native filters</span>
    </div>
  );
}
export function CandidateDetails({ candidate }: { candidate: Candidate }) {
  return (
    <div className="candidate-details">
      <h3>Candidate</h3>
      <p className="identifier">{candidate.id}</p>
      <p>
        Original native filters: {filterLabel(candidate).toLowerCase()}
        {candidate.failed_filters.length
          ? ` · ${candidate.failed_filters.join(", ")}`
          : ""}
      </p>
      <dl>
        {Object.entries(candidate.confidence_metrics).map(([name, value]) => (
          <div className="metadata-row" key={name}>
            <dt title={name}>
              {metricLabels[name] ?? name.replaceAll("_", " ")}
            </dt>
            <dd>
              {value === null ? "Not available" : value.toFixed(3)}{" "}
              {candidate.confidence_units?.[name]}
            </dd>
          </div>
        ))}
      </dl>
      {candidate.comparison_scores.length > 0 && (
        <>
          <h3>Fixed-key comparisons</h3>
          <p className="muted">
            The same key was scored against each listed target state. Choose a
            context and seed in Prediction to view its pose.
          </p>
          <div className="comparison-table">
            <table>
              <thead>
                <tr>
                  <th>Target state</th>
                  <th>Seed</th>
                  <th>Interface ipTM</th>
                  <th>Key pTM</th>
                </tr>
              </thead>
              <tbody>
                {candidate.comparison_scores.map((row) => (
                  <tr key={row.id}>
                    <td>{row.scored_on}</td>
                    <td>{row.seed}</td>
                    <td>
                      {row.metrics.design_to_target_iptm?.toFixed(3) ??
                        "Not available"}
                    </td>
                    <td>
                      {row.metrics.design_ptm?.toFixed(3) ?? "Not available"}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <details>
            <summary>All comparison scores</summary>
            <pre>
              {JSON.stringify(
                candidate.comparison_scores.map(
                  ({ id, seed, scored_on, metrics }) => ({
                    id,
                    seed,
                    scored_on,
                    metrics,
                  }),
                ),
                null,
                2,
              )}
            </pre>
          </details>
        </>
      )}
      <h3>Binder sequence · {candidate.sequence.length} aa</h3>
      <pre data-testid="binder-sequence">{candidate.sequence}</pre>
      <p className="muted">
        Confidence and filter passes are model evidence, not measured binding.
      </p>
    </div>
  );
}
