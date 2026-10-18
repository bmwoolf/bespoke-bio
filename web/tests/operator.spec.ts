import { readFileSync } from "node:fs";
import { expect, test } from "@playwright/test";

const passwordFile = process.env.BESPOKE_DEMO_OPERATOR_PASSWORD_FILE;
const targetId = "protein_e055e4bd687a0c5bd5ad8afc76a18d42";

test("operator access protects every route while public browsing stays read-only", async ({
  request,
}) => {
  test.skip(!passwordFile, "Opt-in deployed operator check");
  for (const path of [
    "/operator/",
    "/operator/api/catalog",
    "/operator/api/jobs",
    "/operator/openapi.json",
  ]) {
    expect((await request.get(path)).status()).toBe(401);
  }
  expect(
    (await request.post("/operator/api/prepare", { data: {} })).status(),
  ).toBe(401);
  expect((await request.post("/api/prepare", { data: {} })).status()).toBe(403);
  expect((await request.get("/api/catalog")).ok()).toBe(true);
});

test("operator button preserves the selected protein and dispatches only after confirmation", async ({
  browser,
}, testInfo) => {
  test.skip(!passwordFile, "Opt-in deployed operator check");
  const context = await browser.newContext({
    baseURL: testInfo.project.use.baseURL,
    httpCredentials: {
      username: "operator",
      password: readFileSync(passwordFile!, "utf8").trim(),
    },
    viewport: { width: 1440, height: 1000 },
  });
  const page = await context.newPage();
  const errors: string[] = [];
  const writes: string[] = [];
  page.on("pageerror", (error) => errors.push(error.message));
  // Never start a real GPU in browser tests. Both mutations and status are intercepted.
  const id = "f".repeat(32);
  const job = {
    job_id: id,
    target_id: targetId,
    status: "running",
    mode: "test",
    execution_performed: false,
    dispatch_attempted: true,
    progress: { stage: "folding" },
    error: null,
  };
  await page.route("**/operator/api/jobs", (route) =>
    route.fulfill({ json: [] }),
  );
  await page.route("**/operator/api/prepare", (route) => {
    writes.push("prepare");
    expect(route.request().postDataJSON().target_id).toBe(targetId);
    return route.fulfill({
      json: {
        ...job,
        status: "prepared_pending_gpu_execution",
        dispatch_attempted: false,
      },
    });
  });
  await page.route("**/operator/api/jobs/*/run", (route) => {
    writes.push("run");
    return route.fulfill({ json: job });
  });
  await page.route(`**/operator/api/jobs/${id}`, (route) =>
    route.fulfill({ json: job }),
  );
  await page.goto(`/?target=${targetId}`);
  const link = page.getByRole("link", { name: "Generate keys", exact: true });
  await expect(link).toHaveAttribute("href", `/operator/?target=${targetId}`);
  await link.click();
  await expect(page).toHaveURL(new RegExp(`/operator/\\?target=${targetId}`));
  await expect(
    page.getByRole("heading", { name: "A1BG", exact: true }),
  ).toBeVisible();
  await expect(page.locator(".molstar-viewport")).toHaveAttribute(
    "aria-busy",
    "false",
  );
  const generate = page.getByRole("button", {
    name: "Generate keys",
    exact: true,
  });
  await expect(generate).toBeEnabled();
  expect(writes).toEqual([]);
  await page.screenshot({
    path: testInfo.outputPath("operator-ready-desktop.png"),
  });
  page.once("dialog", (dialog) => dialog.dismiss());
  await generate.click();
  expect(writes).toEqual([]);
  await page.setViewportSize({ width: 320, height: 700 });
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= innerWidth,
    ),
  ).toBe(true);
  await page.screenshot({
    path: testInfo.outputPath("operator-ready-mobile.png"),
  });
  page.once("dialog", (dialog) => dialog.accept());
  await generate.click();
  await expect(page.getByTestId("job-status")).toHaveText("Generating");
  expect(writes).toEqual(["prepare", "run"]);
  await expect(page.getByTestId("viewer-info")).toContainText(
    "495 residues · 1 chain",
  );
  expect(errors).toEqual([]);
  await context.close();
});
