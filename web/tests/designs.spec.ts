import { expect, test } from "@playwright/test";
import type { Page } from "@playwright/test";

async function ready(page: Page, _residues: string) {
  await expect(page.locator(".molstar-viewport")).toHaveAttribute(
    "aria-busy",
    "false",
  );
}

test("saved keys open first and every view stays associated with its target", async ({
  page,
}, info) => {
  const writes: string[] = [],
    errors: string[] = [];
  page.on("request", (request) => {
    if (request.method() === "POST") writes.push(request.url());
  });
  page.on("pageerror", (error) => errors.push(error.message));
  await page.goto("/");
  await expect(
    page.getByRole("button", { name: "Saved keys", exact: true }),
  ).toHaveAttribute("aria-pressed", "true");
  await expect(
    page.getByRole("heading", { name: "FGF2", exact: true }),
  ).toBeVisible();
  await ready(page, "199 residues · 2 chains");
  await expect(page.locator(".candidate-option")).toHaveCount(32);
  await page.screenshot({
    path: info.outputPath("saved-complex-default-desktop.png"),
  });
  const structure = page.getByRole("group", {
    name: "Structure",
    exact: true,
  });
  await structure.getByRole("button", { name: "Complex", exact: true }).click();
  await ready(page, "199 residues · 2 chains");
  await page.getByLabel("Prediction", { exact: true }).click();
  await expect(
    page.getByRole("combobox", { name: "Prediction context" }),
  ).toHaveValue("evaluated");
  const scores = page.getByLabel("Computational confidence", { exact: true });
  await expect(scores).toContainText("97.1");
  await expect(scores).toContainText("94.1");
  await expect(scores).toContainText("0.868 / 0.899");
  await page
    .getByRole("combobox", { name: "Prediction seed" })
    .selectOption("22");
  await ready(page, "199 residues · 2 chains");
  await page
    .getByRole("combobox", { name: "Prediction seed" })
    .selectOption("33");
  await ready(page, "199 residues · 2 chains");
  await page
    .getByRole("combobox", { name: "Prediction context" })
    .selectOption("native");
  await ready(page, "199 residues · 2 chains");
  await expect(scores).toContainText("Original design confidence");
  await page
    .getByRole("combobox", { name: "Prediction context" })
    .selectOption("evaluated");
  await page
    .getByRole("combobox", { name: "Prediction seed" })
    .selectOption("11");
  await ready(page, "199 residues · 2 chains");
  await page.getByLabel("Prediction", { exact: true }).click();
  await expect(page.getByRole("group", { name: "Complex chains" })).toHaveCount(
    0,
  );
  await page.screenshot({ path: info.outputPath("saved-complex-desktop.png") });
  await structure.getByRole("button", { name: "Target", exact: true }).click();
  await ready(page, "124 residues · 1 chain");
  await page.locator(".candidate-option").last().click();
  await structure.getByRole("button", { name: "Key", exact: true }).click();
  await ready(page, "75 residues · 1 chain");
  await page
    .getByRole("button")
    .filter({ hasText: /^IL1RN/ })
    .first()
    .click();
  await expect(
    page.getByRole("heading", { name: "IL1RN", exact: true }),
  ).toBeVisible();
  await ready(page, "220 residues · 2 chains");
  await page.reload();
  await expect(
    page.getByRole("heading", { name: "IL1RN", exact: true }),
  ).toBeVisible();
  await ready(page, "220 residues · 2 chains");
  await structure.getByRole("button", { name: "Complex", exact: true }).click();
  await ready(page, "220 residues · 2 chains");
  await page.getByRole("button", { name: "Proteins", exact: true }).click();
  await page
    .getByRole("button")
    .filter({ hasText: "HG002 · 45bf5095" })
    .click();
  await expect(
    page.getByRole("heading", { name: "A1BG", exact: true }),
  ).toBeVisible();
  await expect(page.locator(".candidate-option")).toHaveCount(2);
  await structure.getByRole("button", { name: "Key", exact: true }).click();
  await ready(page, "75 residues · 1 chain");
  await structure.getByRole("button", { name: "Target", exact: true }).click();
  await ready(page, "495 residues · 1 chain");
  await page.screenshot({ path: info.outputPath("a1bg-keys-desktop.png") });
  await page.setViewportSize({ width: 320, height: 700 });
  await page.getByRole("button", { name: "Browse", exact: true }).click();
  const drawer = page.getByRole("dialog", { name: "Browse", exact: true });
  await drawer.getByRole("button", { name: "Saved keys", exact: true }).click();
  await page.screenshot({
    path: info.outputPath("saved-keys-mobile-sidebar.png"),
  });
  await drawer.getByRole("button").filter({ hasText: /^FGF2/ }).first().click();
  await ready(page, "199 residues · 2 chains");
  await page.screenshot({ path: info.outputPath("saved-complex-mobile.png") });
  await page.getByLabel("Prediction", { exact: true }).click();
  await expect(
    page.getByRole("combobox", { name: "Prediction seed" }),
  ).toHaveValue("11");
  const box = await page
    .locator(".prediction-popover .popover-content")
    .boundingBox();
  expect(box).not.toBeNull();
  expect(box!.x).toBeGreaterThanOrEqual(0);
  expect(box!.x + box!.width).toBeLessThanOrEqual(320);
  await page.screenshot({
    path: info.outputPath("prediction-menu-mobile.png"),
  });
  await page.getByLabel("Prediction", { exact: true }).click();
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= innerWidth,
    ),
  ).toBe(true);
  expect(writes).toEqual([]);
  expect(errors).toEqual([]);
});
