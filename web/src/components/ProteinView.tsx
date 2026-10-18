import { lazy, Suspense } from "react";
import type { ProteinStructure } from "../api";
import { Loading } from "./Status";
import { ViewerBoundary } from "./ViewerBoundary";

const MolecularViewer = lazy(() => import("../viewer/MolecularViewer"));

export function ProteinView({
  structure,
  title,
}: {
  structure: ProteinStructure;
  title: string;
}) {
  return (
    <>
      <ViewerBoundary
        key={structure.structure_url}
        downloadUrl={structure.structure_url}
      >
        <Suspense fallback={<Loading>Loading molecular viewer…</Loading>}>
          <MolecularViewer
            key={structure.structure_url}
            selection={{
              url: structure.structure_url,
              title: `Protein · ${title}`,
            }}
          />
        </Suspense>
      </ViewerBoundary>
      <div className="workspace-footer">
        <span className="muted">
          {structure.source}
          {structure.mean_plddt !== null &&
            ` · mean pLDDT ${structure.mean_plddt.toFixed(1)}/100`}
        </span>
        {!!structure.geometry_screen_flags && (
          <span className="research-note">Geometry review needed</span>
        )}
        <details className="workspace-details">
          <summary>Structure data</summary>
          <p>
            <a href={structure.structure_url}>Coordinates</a> ·{" "}
            <a href={structure.confidence_url}>Confidence data</a> ·{" "}
            <a href={structure.provenance_url}>Provenance</a>
          </p>
          <p className="muted">
            pLDDT estimates local confidence, not binding or biological
            activity. Full confidence data includes uncertainty between regions.
            {!!structure.geometry_screen_flags &&
              ` The peptide geometry screen flagged ${structure.geometry_screen_flags} bond distance(s); inspect before design.`}
          </p>
        </details>
        <span className="research-note">Research prediction</span>
      </div>
    </>
  );
}
