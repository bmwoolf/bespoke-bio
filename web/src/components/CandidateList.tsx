import type { Candidate } from "../api";

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
            <dt>{name.replaceAll("_", " ")}</dt>
            <dd>
              {value === null ? "Not available" : value.toFixed(3)}{" "}
              {candidate.confidence_units?.[name]}
            </dd>
          </div>
        ))}
      </dl>
      <h3>Binder sequence · {candidate.sequence.length} aa</h3>
      <pre data-testid="binder-sequence">{candidate.sequence}</pre>
      <p className="muted">
        Confidence and filter passes are model evidence, not measured binding.
      </p>
    </div>
  );
}
