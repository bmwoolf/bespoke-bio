import type { Target } from "../api";
import { formatCount } from "./Status";

export function TargetDetails({ target }: { target: Target }) {
  return (
    <div className="target-details">
      <h3>Protein record</h3>
      <dl>
        {Object.entries({
          Donor: target.donor,
          Transcript: target.transcript_id,
          "Target ID": target.target_id,
          "Phase scenario": target.phase_scenario,
          Length: `${formatCount(target.amino_acid_sequence.length)} amino acids`,
          Sequence:
            target.sequence_sha256 === target.reference_sequence_sha256
              ? "Matches reference"
              : "Differs from reference",
          Expression: "Not measured",
          "Sequence SHA-256": target.sequence_sha256,
        }).map(([key, value]) => (
          <div className="metadata-row" key={key}>
            <dt>{key}</dt>
            <dd>{value}</dd>
          </div>
        ))}
      </dl>
      <h3>Amino-acid sequence</h3>
      <pre>{target.amino_acid_sequence}</pre>
      <details>
        <summary>Original record</summary>
        <pre>{JSON.stringify(target, null, 2)}</pre>
      </details>
    </div>
  );
}
