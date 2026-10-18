import { afterEach, describe, expect, it, vi } from "vitest";
import {
  api,
  candidateSchema,
  isActive,
  isComplete,
  resultsSchema,
  targetSchema,
} from "./api";

const target = {
  target_id: "synthetic-1",
  donor: "synthetic",
  gene_id: "G1",
  gene_symbol: "DEMO",
  transcript_id: "T1",
  phase_scenario: "explicit-1",
  amino_acid_sequence: "ME",
  sequence_sha256: "a".repeat(64),
  reference_sequence_sha256: "b".repeat(64),
  expression: null,
};
const job = {
  job_id: "c".repeat(32),
  target_id: target.target_id,
  status: "prepared_pending_gpu_execution",
  mode: "prepare",
  execution_performed: false,
  progress: null,
  error: null,
  dispatch_attempted: null,
};
const candidate = {
  id: "binder-1",
  sequence: "ME",
  native_filter_passed: false,
  failed_filters: ["rmsd"],
  confidence_metrics: { iptm: 0.2, unavailable: null },
  binder_structure_url: `/api/jobs/${job.job_id}/artifacts/binder.cif`,
  complex_structure_url: `/api/jobs/${job.job_id}/artifacts/complex.cif`,
};

afterEach(() => vi.unstubAllGlobals());
describe("typed handoff boundary", () => {
  it("keeps operator requests and artifact links inside the protected prefix", async () => {
    vi.stubGlobal("location", { pathname: "/operator/" });
    const fetch = vi.fn().mockResolvedValue(new Response(JSON.stringify(job)));
    vi.stubGlobal("fetch", fetch);
    await api.job(job.job_id);
    expect(fetch).toHaveBeenCalledWith(`/operator/api/jobs/${job.job_id}`, {
      signal: undefined,
    });
    const parsed = candidateSchema.parse(candidate);
    expect(parsed.complex_structure_url).toBe(
      "/operator" + candidate.complex_structure_url,
    );
    expect(parsed.binder_structure_url).toBe(
      "/operator" + candidate.binder_structure_url,
    );
    expect(
      candidateSchema.safeParse({
        ...candidate,
        complex_structure_url: "//untrusted.example/model.cif",
      }).success,
    ).toBe(false);
  });
  it("retains unmeasured expression and additional provenance verbatim", () => {
    const original = { ...target, evidence: { source: "synthetic fixture" } };
    expect(targetSchema.parse(original)).toEqual(original);
    expect(
      targetSchema.safeParse({ ...target, sequence_sha256: "broken" }).success,
    ).toBe(false);
  });
  it("retains failed candidates and missing confidence without inventing scores", () => {
    expect(candidateSchema.parse(candidate)).toEqual({
      ...candidate,
      predictions: [],
    });
    expect(
      candidateSchema.safeParse({
        ...candidate,
        confidence_metrics: { iptm: NaN },
      }).success,
    ).toBe(false);
  });
  it("rejects external artifact links and malformed results", () => {
    expect(
      candidateSchema.safeParse({
        ...candidate,
        complex_structure_url: "https://elsewhere.example/file.cif",
      }).success,
    ).toBe(false);
    expect(
      resultsSchema.safeParse({ run_id: job.job_id, candidates: [] }).success,
    ).toBe(false);
  });
  it("prepares through an explicit mutation with no run request or browser token", async () => {
    const fetch = vi
      .fn()
      .mockResolvedValue(new Response(JSON.stringify(job), { status: 201 }));
    vi.stubGlobal("fetch", fetch);
    const selection = {
      target_id: target.target_id,
      binder_length: 75,
      num_designs: 2,
    };
    expect(await api.prepare(selection)).toEqual(job);
    expect(fetch).toHaveBeenCalledExactlyOnceWith("/api/prepare", {
      method: "POST",
      headers: { "Content-Type": "application/json", "X-Bespoke-Request": "1" },
      body: JSON.stringify(selection),
    });
  });
  it("propagates ambiguous submissions without retrying", async () => {
    const fetch = vi
      .fn()
      .mockRejectedValue(new TypeError("Network interrupted"));
    vi.stubGlobal("fetch", fetch);
    await expect(
      api.prepare({
        target_id: target.target_id,
        binder_length: 75,
        num_designs: 2,
      }),
    ).rejects.toThrow("Network interrupted");
    expect(fetch).toHaveBeenCalledTimes(1);
  });
  it("rejects bad successful responses instead of rendering made-up defaults", async () => {
    vi.stubGlobal(
      "fetch",
      vi.fn().mockResolvedValue(new Response('{"target_records":1}')),
    );
    await expect(api.catalog()).rejects.toThrow("unexpected response");
  });
  it("keeps unknown submission distinct from running and cached results", () => {
    expect(isActive("submission_unknown")).toBe(false);
    expect(isComplete("prepared_pending_gpu_execution")).toBe(false);
    expect(isComplete("cached_demo_completed")).toBe(true);
    expect(isActive("running")).toBe(true);
  });
});
