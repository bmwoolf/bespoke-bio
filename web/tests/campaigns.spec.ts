import { expect, test } from "@playwright/test";

test("read-only release preserves all campaigns, variants, predictions and pagination", async ({
  page,
  request,
}, info) => {
  test.skip(
    !process.env.BESPOKE_DEMO_READ_ONLY,
    "Opt-in complete read-only release",
  );
  const meta = await (await request.get("/api/catalog")).json();
  expect(meta).toMatchObject({
    read_only: true,
    operator_access: false,
    service_configured: false,
    gpu_execution_enabled: false,
  });
  for (const path of [
    "/operator",
    "/operator/",
    "/operator/api/jobs",
    "/operator/openapi.json",
  ]) {
    expect((await request.get(path)).status()).toBe(410);
  }
  expect((await request.post("/api/prepare", { data: {} })).status()).toBe(403);
  const first = await (await request.get("/api/designs?limit=25")).json();
  const second = await (
    await request.get("/api/designs?limit=25&offset=25")
  ).json();
  expect(first).toMatchObject({
    total: 35,
    candidate_total: 204,
    structured_candidates: 204,
  });
  const records = [...first.records, ...second.records];
  expect(new Set(records.map((r: { job_id: string }) => r.job_id)).size).toBe(
    35,
  );
  const gdf = records.find(
    (r: { name: string; description: string }) =>
      r.name === "GDF15" && r.description.startsWith("Reference"),
  );
  const results = await (
    await request.get(`/api/jobs/${gdf.job_id}/results`)
  ).json();
  const key = results.candidates.find(
    (c: { comparison_scores: unknown[] }) => c.comparison_scores.length > 0,
  );
  const writes: string[] = [],
    errors: string[] = [];
  page.on("request", (r) => {
    if (r.method() !== "GET") writes.push(r.url());
  });
  page.on("pageerror", (e) => errors.push(e.message));
  await page.goto(`/?design=${gdf.job_id}`);
  await expect(
    page.getByRole("heading", { name: "GDF15", exact: true }),
  ).toBeVisible();
  await expect(page.locator(".molstar-viewport")).toHaveAttribute(
    "aria-busy",
    "false",
  );
  await page.locator(`.candidate-option[title="${key.id}"]`).click();
  await page.getByLabel("Prediction", { exact: true }).click();
  await page
    .getByRole("combobox", { name: "Prediction context" })
    .selectOption("comparison:H202D");
  await page
    .getByRole("combobox", { name: "Prediction seed" })
    .selectOption("22");
  await expect(page.locator(".molstar-viewport")).toHaveAttribute(
    "aria-busy",
    "false",
  );
  const selected = key.comparison_scores.find(
    (c: { scored_on: string; seed: number }) =>
      c.scored_on === "H202D" && c.seed === 22,
  );
  await expect(
    page.getByLabel("Computational confidence", { exact: true }),
  ).toContainText(selected.metrics.design_to_target_iptm.toFixed(3));
  await expect(page.getByTestId("viewer-info")).toContainText(
    "187 residues · 2 chains",
  );
  await expect(page.getByLabel("Chain colors")).toContainText("Target · A");
  await expect(page.getByLabel("Chain colors")).toContainText("Key · B");
  await page.getByLabel("Prediction", { exact: true }).click();
  await page.screenshot({
    path: info.outputPath("campaign-comparison-desktop.png"),
  });
  await page
    .getByRole("group", { name: "Structure", exact: true })
    .getByRole("button", { name: "Key", exact: true })
    .click();
  await expect(page.getByTestId("viewer-info")).toContainText(
    "75 residues · 1 chain",
  );
  await page.getByRole("button", { name: "Next", exact: true }).click();
  await expect(page.locator(".pagination")).toContainText("26–35");
  const pf4 = records.find((r: { name: string }) => r.name === "PF4");
  const pf4Result = await (
    await request.get(`/api/jobs/${pf4.job_id}/results`)
  ).json();
  const flagged = pf4Result.candidates.find(
    (c: { assembly_assessment?: string }) => c.assembly_assessment,
  );
  expect(flagged.native_filter_passed).toBe(true);
  await page.goto(`/?design=${pf4.job_id}`);
  await page.locator(`.candidate-option[title="${flagged.id}"]`).click();
  await expect(page.locator(".assembly-note")).toContainText(
    "overlaps the intact tetramer",
  );
  await expect(page.locator(".molstar-viewport")).toHaveAttribute(
    "aria-busy",
    "false",
  );
  await page.setViewportSize({ width: 375, height: 812 });
  await page.screenshot({ path: info.outputPath("campaign-mobile.png") });
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= innerWidth,
    ),
  ).toBe(true);
  await expect(
    page.getByRole("button", { name: "Generate keys", exact: true }),
  ).toHaveCount(0);
  await expect(
    page.getByRole("link", { name: "Generate keys", exact: true }),
  ).toHaveCount(0);
  expect(writes).toEqual([]);
  expect(errors).toEqual([]);
});
