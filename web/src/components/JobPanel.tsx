import { lazy, Suspense, useState } from "react";
import { useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import { api, isActive, isComplete } from "../api";
import type { ProteinStructure, Target } from "../api";
import { ProteinView } from "./ProteinView";
import { Popover } from "./Popover";
import { ViewerBoundary } from "./ViewerBoundary";
import { CandidateList, CandidateDetails } from "./CandidateList";
import { TargetDetails } from "./TargetDetails";
import { ErrorMessage, Loading, jobLabel } from "./Status";

const MolecularViewer = lazy(() => import("../viewer/MolecularViewer"));

export function JobPanel({
  jobId,
  allowRun,
  target,
  proteinStructure,
  expectedTargetId,
  initialView = "complex",
}: {
  jobId: string;
  allowRun: boolean;
  target?: Target;
  proteinStructure?: ProteinStructure | null;
  expectedTargetId?: string | null;
  initialView?: "complex" | "binder" | "protein";
}) {
  const cache = useQueryClient();
  const [candidateId, setCandidateId] = useState<string | null>(null);
  const [evaluation, setEvaluation] = useState("evaluated");
  const [seed, setSeed] = useState(11);
  const [structureKind, setStructureKind] = useState<
    "complex" | "binder" | "protein"
  >(initialView);
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
  // A shared URL must never display another target's job or result under this heading.
  if (
    (job &&
      job.target_id !== (expectedTargetId ?? target?.target_id ?? null)) ||
    (results.data &&
      results.data.target_id !==
        (expectedTargetId ?? target?.target_id ?? null))
  ) {
    return (
      <div className="error" role="alert">
        This result belongs to a different protein. Choose it from History.
      </div>
    );
  }
  const candidate =
    results.data?.candidates.find((value) => value.id === candidateId) ??
    results.data?.candidates[0];
  const prediction =
    evaluation === "evaluated"
      ? (candidate?.predictions.find((p) => p.seed === seed) ??
        candidate?.predictions[0])
      : undefined;
  const comparisons =
    candidate?.comparison_scores.filter(
      (p) => evaluation === `comparison:${p.scored_on}`,
    ) ?? [];
  const comparison = comparisons.find((p) => p.seed === seed) ?? comparisons[0];
  const model = prediction ?? comparison ?? candidate;
  const context = prediction
    ? "evaluated"
    : comparison
      ? `comparison:${comparison.scored_on}`
      : "native";
  const seeds = prediction ? candidate!.predictions : comparisons;
  const structure =
    candidate &&
    model &&
    (structureKind !== "protein" || model.target_structure_url)
      ? {
          url:
            structureKind === "protein"
              ? model.target_structure_url!
              : structureKind === "complex"
                ? model.complex_structure_url
                : model.binder_structure_url,
          title: `${structureKind === "complex" ? "Complex" : structureKind === "protein" ? "Target" : "Key"} · ${candidate.id}`,
          chainRoles:
            structureKind === "complex" ? model.complex_chains : undefined,
        }
      : null;
  const ready = job?.status === "prepared_pending_gpu_execution";
  const canCancel =
    job && ["queued", "running", "cancel_requested"].includes(job.status);
  return (
    <>
      <ErrorMessage error={status.error} />
      <ErrorMessage error={results.error} />
      <ErrorMessage error={action.error} />
      {status.error && (
        <button className="quiet" onClick={() => void status.refetch()}>
          Refresh status
        </button>
      )}
      {results.error && (
        <button className="quiet" onClick={() => void results.refetch()}>
          Reload results
        </button>
      )}
      {!job && status.isPending && <Loading>Loading job…</Loading>}
      {job && (
        <>
          {job.error && (
            <div className="error" role="alert">
              {job.error}
            </div>
          )}
          <div className="result-toolbar">
            {candidate ? (
              <div
                className="structure-tabs"
                role="group"
                aria-label="Structure"
              >
                {(
                  [
                    ["protein", "Target"],
                    ["binder", "Key"],
                    ["complex", "Complex"],
                  ] as const
                ).map(([value, label]) => (
                  <button
                    key={value}
                    aria-pressed={structureKind === value}
                    disabled={
                      value === "protein" &&
                      !proteinStructure &&
                      !model?.target_structure_url
                    }
                    onClick={() => setStructureKind(value)}
                  >
                    {label}
                  </button>
                ))}
              </div>
            ) : (
              <span className="status-label" data-testid="job-status">
                {results.data && !results.data.candidates.length
                  ? "No keys returned"
                  : jobLabel(job.status)}
              </span>
            )}
            {candidate && (
              <Popover className="prediction-popover" label="Prediction">
                {candidate.predictions.length > 0 ||
                candidate.comparison_scores.length > 0 ? (
                  <>
                    <label>
                      Prediction context
                      <select
                        value={context}
                        onChange={(event) => setEvaluation(event.target.value)}
                      >
                        {candidate.predictions.length > 0 && (
                          <option value="evaluated">Standard refold</option>
                        )}
                        <option value="native">Original design</option>
                        {[
                          ...new Set(
                            candidate.comparison_scores.map((p) => p.scored_on),
                          ),
                        ].map((state) => (
                          <option key={state} value={`comparison:${state}`}>
                            Fixed key · {state}
                          </option>
                        ))}
                      </select>
                    </label>
                    {(prediction || comparison) && (
                      <label>
                        Prediction seed
                        <select
                          value={(prediction ?? comparison)!.seed}
                          onChange={(event) =>
                            setSeed(Number(event.target.value))
                          }
                        >
                          {seeds.map((p) => (
                            <option key={p.id} value={p.seed}>
                              Seed {p.seed}
                            </option>
                          ))}
                        </select>
                      </label>
                    )}
                  </>
                ) : (
                  <p>Original design · No standard refold in this export.</p>
                )}
                <p className="prediction-note">
                  {prediction
                    ? `${prediction.context}. Target and Key show independent folds; Complex shows their separately predicted binding pose.`
                    : comparison
                      ? `Fixed key scored against ${comparison.scored_on}, seed ${comparison.seed}. Target and Key are extracted from this complex. Scores across prediction protocols are not directly comparable.`
                      : "Complex shows the original predicted pose. Target and Key are displayed separately."}
                </p>
              </Popover>
            )}
            <div className="result-actions">
              {ready && (
                <>
                  <span className="muted">
                    {!allowRun
                      ? "Generation unavailable in this demo"
                      : "Generation uses paid compute"}
                  </span>
                  <button
                    disabled={
                      !allowRun || !!job.dispatch_attempted || action.isPending
                    }
                    onClick={() => {
                      if (
                        window.confirm(
                          "Generate keys for this protein using paid GPU compute?",
                        )
                      )
                        action.mutate("run");
                    }}
                  >
                    Generate keys
                  </button>
                </>
              )}
              {canCancel && (
                <button
                  className="quiet"
                  disabled={
                    job.status === "cancel_requested" || action.isPending
                  }
                  onClick={() => action.mutate("cancel")}
                >
                  Cancel
                </button>
              )}
              {results.data && (
                <Popover className="download-popover" label="Download">
                  <a href={results.data.download_url}>All results (.zip)</a>
                  {candidate && (
                    <>
                      <a href={model?.complex_structure_url}>
                        Complex structure
                      </a>
                      <a href={model?.binder_structure_url}>Binder structure</a>
                      {model?.target_structure_url && (
                        <a href={model.target_structure_url}>
                          Target structure
                        </a>
                      )}
                      {prediction && (
                        <a href={prediction.evidence_url}>
                          Refolding evidence (.json)
                        </a>
                      )}
                    </>
                  )}
                </Popover>
              )}
            </div>
          </div>
          {isActive(job.status) && (
            <Loading>
              {job.progress?.stage?.replaceAll("_", " ") ||
                jobLabel(job.status)}
            </Loading>
          )}
          {structure ? (
            <ViewerBoundary key={structure.url} downloadUrl={structure.url}>
              <Suspense
                fallback={
                  <div className="empty-structure">
                    <Loading>Loading molecular viewer…</Loading>
                  </div>
                }
              >
                <MolecularViewer key={structure.url} selection={structure} />
              </Suspense>
            </ViewerBoundary>
          ) : proteinStructure ? (
            <ProteinView
              structure={proteinStructure}
              title={target?.gene_symbol || "Target"}
            />
          ) : (
            <div className="empty-structure">
              <div>
                <span className="empty-mark" aria-hidden="true">
                  ◇
                </span>
                <h2>
                  {isActive(job.status)
                    ? jobLabel(job.status)
                    : ready
                      ? "Protein ready"
                      : results.data
                        ? "No candidates returned"
                        : "No structure available"}
                </h2>
                <p>
                  {results.isFetching
                    ? "Loading candidates…"
                    : ready
                      ? "Sequence validated. Generate binders to create structures."
                      : job.progress?.stage
                        ? job.progress.stage.replaceAll("_", " ")
                        : "Structures appear here when results are available."}
                </p>
              </div>
            </div>
          )}
          {candidate && results.data && (
            <CandidateList
              candidates={results.data.candidates}
              selectedId={candidate.id}
              onSelect={setCandidateId}
            />
          )}
          {prediction && (
            <div
              className="prediction-metrics"
              aria-label="Computational confidence"
            >
              <span>Computational confidence</span>
              {structureKind === "complex" ? (
                <>
                  <span>
                    Target{" "}
                    <strong>
                      {prediction.confidence_metrics.target_plddt?.toFixed(1)}
                    </strong>{" "}
                    pLDDT
                  </span>
                  <span>
                    Key{" "}
                    <strong>
                      {prediction.confidence_metrics.binder_plddt?.toFixed(1)}
                    </strong>{" "}
                    pLDDT
                  </span>
                  <span>
                    Interface{" "}
                    <strong>
                      {prediction.confidence_metrics.target_to_binder_iptm?.toFixed(
                        3,
                      )}{" "}
                      /{" "}
                      {prediction.confidence_metrics.binder_to_target_iptm?.toFixed(
                        3,
                      )}
                    </strong>{" "}
                    ipTM
                  </span>
                  <span>
                    Close overlaps{" "}
                    <strong>
                      {prediction.confidence_metrics.close_atomic_overlaps}
                    </strong>
                  </span>
                </>
              ) : (
                <span>
                  {structureKind === "binder" ? "Key" : "Target"}{" "}
                  <strong>
                    {(structureKind === "binder"
                      ? prediction.binder_alone_plddt
                      : prediction.target_alone_plddt
                    ).toFixed(1)}
                  </strong>{" "}
                  pLDDT
                </span>
              )}
            </div>
          )}
          {candidate && !prediction && (
            <div
              className="prediction-metrics"
              aria-label="Computational confidence"
            >
              <span>
                {comparison
                  ? `Fixed key · ${comparison.scored_on} · seed ${comparison.seed}`
                  : "Original design confidence"}
              </span>
              <span title="Predicted interface confidence, not binding probability">
                Interface{" "}
                <strong>
                  {(
                    comparison?.metrics ?? candidate.confidence_metrics
                  ).design_to_target_iptm?.toFixed(3) ?? "Not available"}
                </strong>{" "}
                ipTM
              </span>
              {!comparison && (
                <span title="Minimum interface predicted aligned-error score">
                  <strong>
                    {candidate.confidence_metrics.esmfold2_ipsae_min?.toFixed(
                      3,
                    ) ?? "Not available"}
                  </strong>{" "}
                  ipSAE
                </span>
              )}
              <span
                title={
                  comparison
                    ? "Predicted global confidence of the designed key"
                    : "Original pipeline's filter RMSD"
                }
              >
                {comparison ? "Key pTM" : "Filter RMSD"}{" "}
                <strong>
                  {(comparison
                    ? comparison.metrics.design_ptm
                    : candidate.confidence_metrics.filter_rmsd
                  )?.toFixed(3) ?? "Not available"}
                </strong>
                {comparison ? "" : " Å"}
              </span>
            </div>
          )}
          {candidate?.assembly_assessment && (
            <p className="assembly-note">
              Assembly check:{" "}
              {{
                intact_tetramer_incompatible_rigid_pose:
                  "the predicted key overlaps the intact tetramer",
                substantial_overlap_all_placements:
                  "the predicted key overlaps the native assembly in every tested placement",
              }[candidate.assembly_assessment] ??
                candidate.assembly_assessment.replaceAll("_", " ")}
              . Native filter status is unchanged.
            </p>
          )}
          <div className="workspace-footer">
            <details className="workspace-details">
              <summary>Details</summary>
              {candidate && <CandidateDetails candidate={candidate} />}
              {target && <TargetDetails target={target} />}
              <details>
                <summary>Job record</summary>
                <pre>{JSON.stringify(job, null, 2)}</pre>
              </details>
            </details>
            <span className="research-note">Research prediction</span>
          </div>
        </>
      )}
    </>
  );
}
