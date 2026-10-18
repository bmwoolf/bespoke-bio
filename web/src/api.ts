import { z } from "zod";

const count = z.number().int().nonnegative();
const digest = z.string().regex(/^[a-f0-9]{64}$/);
const jobId = z.string().regex(/^[a-f0-9]{32}$/);
const artifactUrl = z
  .string()
  .regex(/^\/api\/jobs\/[a-f0-9]{32}\/(artifacts\/[^?#]+|download)$/);
export const catalogSchema = z.object({
  target_records: count,
  unique_sequences: count,
  nonreference_records: count,
  service_configured: z.boolean(),
  gpu_execution_enabled: z.boolean(),
});
export const targetSchema = z
  .object({
    target_id: z.string(),
    donor: z.string(),
    gene_id: z.string(),
    gene_symbol: z.string(),
    transcript_id: z.string(),
    phase_scenario: z.string(),
    amino_acid_sequence: z.string().min(1),
    sequence_sha256: digest,
    reference_sequence_sha256: digest,
    expression: z.null(),
  })
  .passthrough();
export const pageSchema = z.object({
  total: count,
  offset: count,
  limit: count,
  records: z.array(
    z.object({
      target_id: z.string(),
      gene: z.string(),
      transcript: z.string(),
      donor: z.string(),
      changed: z.union([z.literal(0), z.literal(1)]),
      length: count,
    }),
  ),
});
export const jobSchema = z.object({
  job_id: jobId,
  target_id: z.string().nullable(),
  status: z.string(),
  mode: z.string().nullable(),
  execution_performed: z.boolean().nullable(),
  dispatch_attempted: z.boolean().nullable(),
  progress: z.object({ stage: z.string().optional() }).passthrough().nullable(),
  error: z.string().nullable(),
});
const metrics = z.record(z.string(), z.number().finite().nullable());
export const candidateSchema = z.object({
  id: z.string(),
  sequence: z.string(),
  native_filter_passed: z.boolean().nullable(),
  failed_filters: z.array(z.string()),
  confidence_metrics: metrics,
  confidence_units: z.record(z.string(), z.string()).optional(),
  binder_structure_url: artifactUrl,
  complex_structure_url: artifactUrl,
});
export const resultsSchema = z.object({
  run_id: jobId,
  target_id: z.string().nullable(),
  candidates: z.array(candidateSchema),
  download_url: artifactUrl,
});
export type Catalog = z.infer<typeof catalogSchema>;
export type Target = z.infer<typeof targetSchema>;
export type TargetPage = z.infer<typeof pageSchema>;
export type Job = z.infer<typeof jobSchema>;
export type Candidate = z.infer<typeof candidateSchema>;
export type Results = z.infer<typeof resultsSchema>;
export type Selection = {
  target_id: string;
  binder_length: number;
  num_designs: number;
};
export type StructureSelection = { url: string; title: string };

export async function request<T>(
  path: string,
  schema: z.ZodType<T>,
  options: RequestInit = {},
): Promise<T> {
  const response = await fetch(path, options);
  const body: unknown = await response.json();
  if (!response.ok) {
    const error = z.object({ detail: z.string() }).safeParse(body);
    throw new Error(
      error.success ? error.data.detail : `Request failed (${response.status})`,
    );
  }
  const result = schema.safeParse(body);
  if (!result.success)
    throw new Error(
      `The service returned an unexpected response for ${path.split("?")[0]}`,
    );
  return result.data;
}
const post = <T>(path: string, schema: z.ZodType<T>, body: unknown = {}) =>
  request(path, schema, {
    method: "POST",
    headers: { "Content-Type": "application/json", "X-Bespoke-Request": "1" },
    body: JSON.stringify(body),
  });
export const api = {
  catalog: (signal?: AbortSignal) =>
    request("/api/catalog", catalogSchema, { signal }),
  targets: (
    q: string,
    changed: boolean,
    offset: number,
    signal?: AbortSignal,
  ) =>
    request(
      "/api/targets?" +
        new URLSearchParams({
          q,
          changed: String(changed),
          offset: String(offset),
          limit: "25",
        }),
      pageSchema,
      { signal },
    ),
  target: (id: string, signal?: AbortSignal) =>
    request("/api/targets/" + encodeURIComponent(id), targetSchema, { signal }),
  jobs: (signal?: AbortSignal) =>
    request("/api/jobs", z.array(jobSchema), { signal }),
  job: (id: string, signal?: AbortSignal) =>
    request(`/api/jobs/${id}`, jobSchema, { signal }),
  results: (id: string, signal?: AbortSignal) =>
    request(`/api/jobs/${id}/results`, resultsSchema, { signal }),
  prepare: (selection: Selection) => post("/api/prepare", jobSchema, selection),
  demo: () => post("/api/demo", jobSchema),
  run: (id: string) => post(`/api/jobs/${id}/run`, jobSchema),
  cancel: (id: string) => post(`/api/jobs/${id}/cancel`, jobSchema),
};
export const isComplete = (status: string) =>
  ["completed", "cached_demo_completed"].includes(status);
export const isActive = (status: string) =>
  ["queued", "running", "cancel_requested", "dispatch_requested"].includes(
    status,
  );
