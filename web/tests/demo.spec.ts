import { expect, test } from "@playwright/test";
import type { Page } from "@playwright/test";

const publicDeployment = process.env.BESPOKE_DEMO_PUBLIC === "1";

async function hidePreparedHistory(page: Page) {
  await page.route("**/api/jobs", async (route) => {
    const response = await route.fetch();
    const jobs = await response.json();
    await route.fulfill({
      json: jobs.filter(
        (job: { target_id: string | null }) => job.target_id === null,
      ),
    });
  });
}

// Opt-in checks against the running HG002 gateway. Preparation and cached replay only.
async function ready(page: Page) {
  await expect(page.locator(".molstar-viewport")).toHaveAttribute(
    "aria-busy",
    "false",
  );
}
async function searchTarget(page: Page) {
  await page.getByLabel("Search gene or transcript").fill("A1BG");
  await page.getByLabel("Changed sequences", { exact: true }).check();
  await expect(
    page.getByText("8 matching records", { exact: true }),
  ).toBeVisible();
  await page
    .getByRole("button")
    .filter({ hasText: "HG002 · 45bf5095" })
    .click();
  await expect(
    page.getByRole("heading", { name: "A1BG", exact: true }),
  ).toBeVisible();
}

test("public catalog preserves selection across refresh and never submits jobs", async ({
  page,
}, testInfo) => {
  test.skip(!publicDeployment, "Public deployment check");
  const writes: string[] = [];
  page.on("request", (request) => {
    if (!["GET", "HEAD", "OPTIONS"].includes(request.method()))
      writes.push(request.url());
  });
  await page.goto("/");
  await ready(page);
  await searchTarget(page);
  await expect(
    page.getByRole("heading", {
      name: "No 3D structure for this protein",
      exact: true,
    }),
  ).toBeVisible();
  await expect(
    page.getByRole("button", { name: "Prepare target", exact: true }),
  ).toHaveCount(0);
  await expect(page).toHaveURL(
    /target=protein_dcd5494921477cd9218a2a6745bf5095/,
  );
  await page.reload();
  await expect(
    page.getByRole("heading", { name: "A1BG", exact: true }),
  ).toBeVisible();
  await expect(page.locator(".molecule")).toHaveCount(0);
  await page.screenshot({
    path: testInfo.outputPath("sequence-only-desktop.png"),
  });
  await page.setViewportSize({ width: 390, height: 844 });
  await page.getByRole("button", { name: "Proteins", exact: true }).click();
  const drawer = page.getByRole("dialog", { name: "Proteins", exact: true });
  await expect(drawer).toBeVisible();
  await page.screenshot({
    path: testInfo.outputPath("public-mobile-catalog.png"),
  });
  await page.keyboard.press("Escape");
  await page.setViewportSize({ width: 320, height: 700 });
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= innerWidth,
    ),
  ).toBe(true);
  await page.screenshot({
    path: testInfo.outputPath("sequence-only-mobile.png"),
  });
  await page.getByRole("button", { name: "Example", exact: true }).click();
  await ready(page);
  await page.goBack();
  await expect(
    page.getByRole("heading", { name: "A1BG", exact: true }),
  ).toBeVisible();
  await page.goForward();
  await ready(page);
  expect(writes).toEqual([]);
});

test("existing jobs reopen from the sidebar, and cross-target URLs never display a model", async ({
  page,
}) => {
  const targetId = "protein_dcd5494921477cd9218a2a6745bf5095";
  const jobId = "d".repeat(32);
  const prepared = {
    job_id: jobId,
    target_id: targetId,
    status: "prepared_pending_gpu_execution",
    mode: "prepare",
    execution_performed: false,
    dispatch_attempted: null,
    progress: null,
    error: null,
  };
  // Synthetic status only: never attach example coordinates to a donor protein.
  await page.route("**/api/jobs", async (route) => {
    const response = await route.fetch();
    const jobs = await response.json();
    await route.fulfill({
      json: [
        ...jobs.filter(
          (job: { target_id: string | null }) => job.target_id !== targetId,
        ),
        prepared,
      ],
    });
  });
  await page.route(`**/api/jobs/${jobId}`, (route) =>
    route.fulfill({ json: prepared }),
  );
  await page.goto("/");
  await searchTarget(page);
  await expect(page.getByTestId("job-status")).toHaveText("Ready to generate");
  await expect(
    page.getByRole("button", { name: "Prepare target", exact: true }),
  ).toHaveCount(0);
  await page.reload();
  await expect(page.getByTestId("job-status")).toHaveText("Ready to generate");
  const savedJobs = await (await page.request.get("/api/jobs")).json();
  const example = savedJobs.find(
    (job: { target_id: string | null }) => job.target_id === null,
  );
  await page.goto(`/?target=${targetId}&job=${example.job_id}`);
  await expect(page.getByRole("alert")).toContainText("different protein");
  await expect(page.locator(".molecule")).toHaveCount(0);
});

test("public page with no published example never attempts automatic replay", async ({
  page,
}) => {
  test.skip(!publicDeployment, "Public deployment check");
  const writes: string[] = [];
  page.on("request", (request) => {
    if (request.method() === "POST") writes.push(request.url());
  });
  await page.route("**/api/jobs", (route) => route.fulfill({ json: [] }));
  await page.goto("/");
  await expect(
    page.getByText("No saved example has been published here yet.", {
      exact: false,
    }),
  ).toBeVisible();
  expect(writes).toEqual([]);
});

test("saved example opens directly; representations, candidates and downloads stay together", async ({
  page,
}, testInfo) => {
  const errors: string[] = [];
  const unexpected: string[] = [];
  page.on("pageerror", (error) => errors.push(error.message));
  page.on("request", (request) => {
    if (
      request.url().endsWith("/run") ||
      new URL(request.url()).origin !==
        new URL(testInfo.project.use.baseURL!).origin
    )
      unexpected.push(request.url());
  });
  await page.setViewportSize({ width: 1440, height: 900 });
  await page.goto("/");
  await expect(
    page.getByRole("heading", { name: "PD-L1", exact: true }),
  ).toBeVisible();
  await ready(page);
  await expect(page.getByTestId("viewer-info")).toContainText(
    "1,496 atoms · 190 residues · 2 chains",
  );
  await expect(page.getByTestId("viewer-info")).toContainText(
    "GPU surface computation",
  );
  await expect(page.locator(".candidate-option")).toHaveCount(2);
  await expect(
    page.getByRole("button", { name: "01 Passed", exact: true }),
  ).toHaveAttribute("aria-pressed", "true");
  expect(
    await page.evaluate(
      () => document.documentElement.scrollHeight <= innerHeight,
    ),
  ).toBe(true);
  await page.screenshot({ path: testInfo.outputPath("desktop-surface.png") });
  await page.getByRole("button", { name: "Ribbon", exact: true }).click();
  await ready(page);
  await expect(
    page.getByRole("button", { name: "Ribbon", exact: true }),
  ).toHaveAttribute("aria-pressed", "true");
  await page.screenshot({ path: testInfo.outputPath("desktop-ribbon.png") });
  await page.getByRole("button", { name: "Atoms", exact: true }).click();
  await ready(page);
  await page.getByLabel("Appearance", { exact: true }).click();
  await page.getByLabel("Metallic finish").uncheck();
  await page.getByLabel("Color by chain").uncheck();
  await ready(page);
  await expect(page.locator(".chain-legend span")).toHaveCount(0);
  await page.getByLabel("Appearance", { exact: true }).click();
  await page.getByText("Details", { exact: true }).click();
  const binderLength = (await page.getByTestId("binder-sequence").innerText())
    .length;
  await page.getByText("Details", { exact: true }).click();
  await page
    .getByRole("combobox", { name: "Structure", exact: true })
    .selectOption("binder");
  await ready(page);
  await expect(page.getByTestId("viewer-info")).toContainText(
    `${binderLength} residues · 1 chain`,
  );
  await page
    .getByRole("combobox", { name: "Structure", exact: true })
    .selectOption("complex");
  await ready(page);
  await page.getByRole("button", { name: "02 Failed", exact: true }).click();
  await ready(page);
  await expect(
    page.getByRole("button", { name: "02 Failed", exact: true }),
  ).toHaveAttribute("aria-pressed", "true");
  await page.getByText("Details", { exact: true }).click();
  await expect(page.locator(".candidate-details")).toContainText(
    "Computational filters: failed",
  );
  await page.getByText("Details", { exact: true }).click();
  await page.getByRole("button", { name: "01 Passed", exact: true }).click();
  await ready(page);
  await page.getByRole("button", { name: "Fullscreen", exact: true }).click();
  await expect(
    page.getByRole("button", { name: "Exit fullscreen", exact: true }),
  ).toBeVisible();
  await page.getByRole("button", { name: "Reset view", exact: true }).click();
  await page
    .getByRole("button", { name: "Exit fullscreen", exact: true })
    .click();
  await page.getByText("Download", { exact: true }).click();
  const download = page.waitForEvent("download");
  await page
    .getByRole("link", { name: "All results (.zip)", exact: true })
    .click();
  expect(await (await download).failure()).toBeNull();
  await expect(page.locator(".download-popover")).not.toHaveAttribute(
    "open",
    "",
  );
  await page.setViewportSize({ width: 390, height: 844 });
  await ready(page);
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= innerWidth,
    ),
  ).toBe(true);
  await page.screenshot({ path: testInfo.outputPath("mobile-surface.png") });
  expect(errors).toEqual([]);
  expect(unexpected).toEqual([]);
});

test("donor selection removes example; preparation, details and mobile history preserve identity", async ({
  page,
}, testInfo) => {
  test.skip(
    publicDeployment,
    "Preparation belongs to the private operator gateway",
  );
  await hidePreparedHistory(page);
  const runs: string[] = [];
  page.on("request", (request) => {
    if (request.url().endsWith("/run")) runs.push(request.url());
  });
  await page.goto("/");
  await ready(page);
  await searchTarget(page);
  await expect(page.locator(".molecule")).toHaveCount(0);
  await expect(
    page.getByText("No 3D structure for this protein", { exact: true }),
  ).toBeVisible();
  await expect(
    page.getByRole("heading", { name: "PD-L1", exact: true }),
  ).toHaveCount(0);
  await page.getByText("Details", { exact: true }).click();
  await expect(page.getByText("Not measured", { exact: true })).toBeVisible();
  await expect(page.locator(".target-details")).toContainText(
    "protein_dcd5494921477cd9218a2a6745bf5095",
  );
  await page.getByText("Details", { exact: true }).click();
  const preparation = page.waitForRequest((request) =>
    request.url().endsWith("/api/prepare"),
  );
  await page
    .getByRole("button", { name: "Prepare target", exact: true })
    .click();
  expect((await preparation).postDataJSON()).toEqual({
    target_id: "protein_dcd5494921477cd9218a2a6745bf5095",
    binder_length: 75,
    num_designs: 2,
  });
  await expect(page.getByTestId("job-status")).toHaveText("Ready to generate");
  await expect(
    page.getByRole("button", { name: "Generate binders", exact: true }),
  ).toBeDisabled();
  await page.screenshot({ path: testInfo.outputPath("prepared-target.png") });
  await page.getByRole("button", { name: "Example", exact: true }).click();
  await ready(page);
  await page.setViewportSize({ width: 390, height: 844 });
  await page.getByRole("button", { name: "Proteins", exact: true }).click();
  const drawer = page.getByRole("dialog", { name: "Proteins", exact: true });
  await expect(drawer).toBeVisible();
  await drawer.getByLabel("Search gene or transcript").fill("A1BG");
  await drawer.getByLabel("Changed sequences", { exact: true }).check();
  await expect(
    drawer.getByText("8 matching records", { exact: true }),
  ).toBeVisible();
  await page.screenshot({ path: testInfo.outputPath("mobile-proteins.png") });
  await drawer
    .getByRole("button")
    .filter({ hasText: "HG002 · 45bf5095" })
    .click();
  await expect(drawer).not.toBeVisible();
  await expect(
    page.getByRole("heading", { name: "A1BG", exact: true }),
  ).toBeVisible();
  await expect(page.locator(".molecule")).toHaveCount(0);
  await page.setViewportSize({ width: 320, height: 700 });
  await page.getByLabel("Settings", { exact: true }).click();
  const settingsBox = await page
    .locator(".settings-popover .popover-content")
    .boundingBox();
  expect(settingsBox).not.toBeNull();
  expect(settingsBox!.x).toBeGreaterThanOrEqual(0);
  expect(settingsBox!.x + settingsBox!.width).toBeLessThanOrEqual(320);
  await page.screenshot({ path: testInfo.outputPath("mobile-settings.png") });
  await page.getByLabel("Binder length", { exact: true }).fill("20");
  await page.keyboard.press("Escape");
  await page
    .getByRole("button", { name: "Prepare target", exact: true })
    .click();
  await expect(page.getByLabel("Binder length", { exact: true })).toBeVisible();
  await expect(page.getByLabel("Binder length", { exact: true })).toBeFocused();
  await page.getByLabel("Binder length", { exact: true }).fill("75");
  await page.keyboard.press("Escape");
  await page.setViewportSize({ width: 390, height: 844 });
  await page.getByRole("button", { name: "History", exact: true }).click();
  const history = page.getByRole("dialog", { name: "History", exact: true });
  await history
    .getByRole("button")
    .filter({ hasText: "PD-L1 · Saved example" })
    .click();
  await expect(history).not.toBeVisible();
  await expect(
    page.getByRole("heading", { name: "PD-L1", exact: true }),
  ).toBeVisible();
  await ready(page);
  await page.getByRole("button", { name: "History", exact: true }).click();
  await page.keyboard.press("Escape");
  await expect(history).not.toBeVisible();
  await expect(
    page.getByRole("button", { name: "History", exact: true }),
  ).toBeFocused();
  await page.setViewportSize({ width: 320, height: 700 });
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= innerWidth,
    ),
  ).toBe(true);
  await page.screenshot({ path: testInfo.outputPath("small-mobile.png") });
  expect(runs).toEqual([]);
});

test("a late preparation response cannot replace a newly selected example", async ({
  page,
}) => {
  test.skip(
    publicDeployment,
    "Preparation belongs to the private operator gateway",
  );
  await hidePreparedHistory(page);
  await page.goto("/");
  await ready(page);
  await searchTarget(page);
  let release!: () => void;
  const pending = new Promise<void>((resolve) => {
    release = resolve;
  });
  await page.route("**/api/prepare", async (route) => {
    const response = await route.fetch();
    await pending;
    await route.fulfill({ response });
  });
  const request = page.waitForRequest("**/api/prepare");
  await page
    .getByRole("button", { name: "Prepare target", exact: true })
    .click();
  await request;
  await page.getByRole("button", { name: "Example", exact: true }).click();
  const response = page.waitForResponse("**/api/prepare");
  release();
  await response;
  await ready(page);
  await expect(
    page.getByRole("heading", { name: "PD-L1", exact: true }),
  ).toBeVisible();
  await expect(
    page.getByRole("button", { name: "Prepare target", exact: true }),
  ).toHaveCount(0);
});

test("surface computation falls back to CPU without float buffers", async ({
  page,
}) => {
  await page.addInitScript(() => {
    const original = WebGL2RenderingContext.prototype.getExtension;
    Object.defineProperty(WebGL2RenderingContext.prototype, "getExtension", {
      value: function (this: WebGL2RenderingContext, name: string) {
        if (
          name === "EXT_color_buffer_float" ||
          name === "WEBGL_color_buffer_float"
        )
          return null;
        return Reflect.apply(original, this, [name]);
      },
    });
  });
  await page.goto("/");
  await ready(page);
  await expect(page.getByTestId("viewer-info")).toContainText("1,496 atoms");
  await expect(page.getByTestId("viewer-info")).toContainText(
    "CPU surface computation",
  );
});

test("without WebGL the original structure remains downloadable", async ({
  page,
}) => {
  await page.addInitScript(() => {
    const original = HTMLCanvasElement.prototype.getContext;
    Object.defineProperty(HTMLCanvasElement.prototype, "getContext", {
      value: function (
        this: HTMLCanvasElement,
        name: string,
        options?: unknown,
      ) {
        if (name.includes("webgl")) return null;
        return Reflect.apply(original, this, [name, options]);
      },
    });
  });
  await page.goto("/");
  await expect(page.getByRole("alert")).toContainText("WebGL is unavailable");
  await expect(
    page.getByRole("link", { name: "Download structure", exact: true }),
  ).toBeVisible();
  await expect(page.locator(".candidate-option")).toHaveCount(2);
});

test("a failed viewer chunk leaves candidates and downloads usable", async ({
  page,
}) => {
  await page.route(/MolecularViewer(?:-[^/]+\.js|\.tsx)(?:\?|$)/, (route) =>
    route.abort(),
  );
  await page.goto("/");
  await expect(page.getByRole("alert")).toContainText(
    "molecular viewer could not load",
  );
  await expect(
    page.getByRole("link", { name: "download the original structure" }),
  ).toBeVisible();
  await page.getByText("Download", { exact: true }).click();
  await expect(
    page.getByRole("link", { name: "All results (.zip)", exact: true }),
  ).toBeVisible();
  await expect(page.locator(".candidate-option")).toHaveCount(2);
});

test("empty results do not invent a candidate or substitute a structure", async ({
  page,
}) => {
  await page.route("**/api/jobs/*/results", async (route) => {
    const response = await route.fetch();
    const body = await response.json();
    await route.fulfill({ json: { ...body, candidates: [] } });
  });
  await page.goto("/");
  await expect(
    page.getByRole("heading", { name: "No candidates returned", exact: true }),
  ).toBeVisible();
  await expect(page.locator(".molecule")).toHaveCount(0);
  await expect(page.locator(".candidate-option")).toHaveCount(0);
  await page.getByText("Download", { exact: true }).click();
  await expect(
    page.getByRole("link", { name: "All results (.zip)", exact: true }),
  ).toBeVisible();
});

test("automatic cached replay never retries a failed submission", async ({
  page,
}) => {
  test.skip(publicDeployment, "Public snapshots never submit cached replays");
  let submissions = 0;
  await page.route("**/api/jobs", (route) => route.fulfill({ json: [] }));
  await page.route("**/api/demo", (route) => {
    submissions += 1;
    return route.fulfill({
      status: 503,
      json: { detail: "Saved example unavailable" },
    });
  });
  await page.goto("/");
  await expect(page.getByRole("alert")).toHaveText("Saved example unavailable");
  await searchTarget(page);
  await page.getByRole("button", { name: "Example", exact: true }).click();
  await expect(page.getByRole("alert")).toHaveText("Saved example unavailable");
  expect(submissions).toBe(1);
});
