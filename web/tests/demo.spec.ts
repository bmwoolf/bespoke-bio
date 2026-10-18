import { expect, test } from "@playwright/test";

// Opt-in integration check against an already running gateway. Creates only preparation/cached jobs.
test("full catalog, exact preparation and Mol* representations", async ({
  page,
}, testInfo) => {
  const errors: string[] = [];
  const runRequests: string[] = [];
  const externalRequests: string[] = [];
  page.on("pageerror", (error) => errors.push(error.message));
  page.on("request", (request) => {
    if (request.method() === "POST" && request.url().endsWith("/run"))
      runRequests.push(request.url());
    if (
      new URL(request.url()).origin !==
      new URL(testInfo.project.use.baseURL!).origin
    )
      externalRequests.push(request.url());
  });
  await page.goto("/");
  await expect(
    page.getByRole("region", { name: "Catalog totals" }),
  ).toContainText("383,350");
  await page.getByLabel("Gene or transcript prefix").fill("A1BG");
  await page.getByLabel("Non-reference proteins only").check();
  await page.getByRole("button", { name: "Search", exact: true }).click();
  await expect(
    page.getByText("1–8 of 8 matching records", { exact: true }),
  ).toBeVisible();
  await page.getByRole("button").filter({ hasText: "Record 45bf5095" }).click();
  await expect(page.getByText("Not measured", { exact: true })).toBeVisible();
  await page.getByRole("button", { name: "Prepare selected target" }).click();
  await expect(page.getByTestId("job-status")).toHaveText(
    "prepared pending gpu execution",
  );
  await expect(
    page.getByRole("button", { name: "Start GPU job" }),
  ).toBeDisabled();
  await page.getByRole("button", { name: "Load cached PD-L1 example" }).click();
  await expect(page.getByTestId("job-status")).toHaveText(
    "cached demo completed",
  );
  await expect(page.locator(".candidate")).toHaveCount(2);
  await expect(
    page.getByText("Native computational filters: failed", { exact: true }),
  ).toBeVisible();
  const first = page.locator(".candidate").first();
  const binderLength = (await first.locator("pre").innerText()).trim().length;
  await first
    .getByRole("button", { name: "View complex", exact: true })
    .click();
  await expect(page.getByTestId("viewer-info")).toContainText("1,496 atoms");
  await expect(page.getByTestId("viewer-info")).toContainText(
    "GPU surface computation",
  );
  await expect(page.locator(".molstar-viewport")).toHaveAttribute(
    "aria-busy",
    "false",
  );
  await page
    .getByRole("region", { name: "Molecular structure" })
    .scrollIntoViewIfNeeded();
  await page.screenshot({ path: testInfo.outputPath("metallic-surface.png") });
  await page.getByRole("button", { name: "Ribbon", exact: true }).click();
  await expect(
    page.getByRole("button", { name: "Ribbon", exact: true }),
  ).toHaveAttribute("aria-pressed", "true");
  await expect(page.locator(".molstar-viewport")).toHaveAttribute(
    "aria-busy",
    "false",
  );
  await page.screenshot({ path: testInfo.outputPath("ribbon.png") });
  await page.getByRole("button", { name: "Atoms", exact: true }).click();
  await expect(page.locator(".molstar-viewport")).toHaveAttribute(
    "aria-busy",
    "false",
  );
  await page.getByLabel("Metallic finish").uncheck();
  await page.getByLabel("Color by chain").uncheck();
  await expect(page.locator(".molstar-viewport")).toHaveAttribute(
    "aria-busy",
    "false",
  );
  await first.getByRole("button", { name: "View binder", exact: true }).click();
  await expect(page.getByTestId("viewer-info")).toContainText(
    `${binderLength} residues · 1 chain`,
  );
  // Replace and unmount the viewer to exercise lifecycle and stale asynchronous loads.
  await first
    .getByRole("button", { name: "View complex", exact: true })
    .click();
  await expect(page.getByTestId("viewer-info")).toContainText(
    "190 residues · 2 chains",
  );
  const download = page.waitForEvent("download");
  await page.getByRole("link", { name: "Download complete ZIP" }).click();
  expect(await (await download).failure()).toBeNull();
  await page.setViewportSize({ width: 390, height: 844 });
  await expect(page.locator(".molstar-viewport")).toHaveAttribute(
    "aria-busy",
    "false",
  );
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= innerWidth,
    ),
  ).toBe(true);
  await page
    .getByRole("region", { name: "Molecular structure" })
    .scrollIntoViewIfNeeded();
  await page.screenshot({ path: testInfo.outputPath("mobile-surface.png") });
  expect(errors).toEqual([]);
  expect(runRequests).toEqual([]);
  expect(externalRequests).toEqual([]);
});

test("surface computation falls back to CPU when float buffers are unavailable", async ({
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
  await page.getByRole("button", { name: "Load cached PD-L1 example" }).click();
  await expect(page.locator(".candidate")).toHaveCount(2);
  await page
    .locator(".candidate")
    .first()
    .getByRole("button", { name: "View complex", exact: true })
    .click();
  await expect(page.getByTestId("viewer-info")).toContainText("1,496 atoms");
  await expect(page.getByTestId("viewer-info")).toContainText(
    "CPU surface computation",
  );
  await expect(page.locator(".molstar-viewport")).toHaveAttribute(
    "aria-busy",
    "false",
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
  await page.getByRole("button", { name: "Load cached PD-L1 example" }).click();
  await expect(page.locator(".candidate")).toHaveCount(2);
  await page
    .locator(".candidate")
    .first()
    .getByRole("button", { name: "View complex", exact: true })
    .click();
  await expect(page.getByRole("alert")).toContainText("WebGL is unavailable");
  await expect(
    page.getByRole("link", { name: "Download structure", exact: true }),
  ).toBeVisible();
});

test("a failed viewer chunk leaves candidates and downloads usable", async ({
  page,
}) => {
  await page.route("**/MolecularViewer-*.js", (route) => route.abort());
  await page.goto("/");
  await page.getByRole("button", { name: "Load cached PD-L1 example" }).click();
  await expect(page.locator(".candidate")).toHaveCount(2);
  await page
    .locator(".candidate")
    .first()
    .getByRole("button", { name: "View complex", exact: true })
    .click();
  await expect(page.getByRole("alert")).toContainText(
    "molecular viewer could not load",
  );
  await expect(
    page.getByRole("link", { name: "download the original structure" }),
  ).toBeVisible();
  await expect(
    page.getByRole("link", { name: "Download complete ZIP" }),
  ).toBeVisible();
  await expect(page.locator(".candidate")).toHaveCount(2);
});
