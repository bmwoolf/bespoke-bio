import type { Candidate, StructureSelection } from "../api";

export function CandidateList({
  candidates,
  onView,
}: {
  candidates: Candidate[];
  onView: (value: StructureSelection) => void;
}) {
  if (!candidates.length)
    return (
      <p className="empty">The job completed with no exported candidates.</p>
    );
  return (
    <div className="candidates">
      {candidates.map((candidate) => (
        <article className="candidate" key={candidate.id}>
          <h3>{candidate.id}</h3>
          <p
            className={
              candidate.native_filter_passed === true ? "pass" : "fail"
            }
          >
            Native computational filters:{" "}
            {candidate.native_filter_passed === true
              ? "passed"
              : candidate.native_filter_passed === false
                ? "failed"
                : "unknown"}
          </p>
          <div className="metrics">
            {Object.entries(candidate.confidence_metrics).map(
              ([name, value]) => (
                <div key={name}>
                  {name.replaceAll("_", " ")}
                  <strong>
                    {value === null ? "Not available" : value.toFixed(3)}
                  </strong>
                  <small>{candidate.confidence_units?.[name]}</small>
                </div>
              ),
            )}
          </div>
          <pre>{candidate.sequence}</pre>
          {!!candidate.failed_filters.length && (
            <p className="muted">
              Failed filters: {candidate.failed_filters.join(", ")}
            </p>
          )}
          <button
            className="secondary"
            onClick={() =>
              onView({
                url: candidate.complex_structure_url,
                title: `Complex · ${candidate.id}`,
              })
            }
          >
            View complex
          </button>
          <button
            className="secondary"
            onClick={() =>
              onView({
                url: candidate.binder_structure_url,
                title: `Binder · ${candidate.id}`,
              })
            }
          >
            View binder
          </button>
        </article>
      ))}
    </div>
  );
}
