import { expect, test } from "@playwright/test";
import type { Page } from "@playwright/test";

async function fillsScreen(page: Page) {
  await expect
    .poll(() =>
      page.locator(".molstar-viewport canvas").evaluate((canvas) => {
        const box = canvas.getBoundingClientRect();
        return {
          left: Math.round(box.left),
          top: Math.round(box.top),
          width: Math.round(box.width - innerWidth),
          height: Math.round(box.height - innerHeight),
        };
      }),
    )
    .toEqual({ left: 0, top: 0, width: 0, height: 0 });
}

for (const [name, path] of [
  ["saved complex", "/"],
  ["donor protein", "/?target=protein_dcd5494921477cd9218a2a6745bf5095"],
] as const) {
  test(`${name} fills and centers its canvas in fullscreen`, async ({
    page,
  }, testInfo) => {
    await page.setViewportSize({ width: 1440, height: 900 });
    await page.goto(path);
    const viewport = page.locator(".molstar-viewport");
    await expect(viewport).toHaveAttribute("aria-busy", "false");
    const embedded = await viewport.boundingBox();
    expect(embedded).not.toBeNull();
    await page.getByRole("button", { name: "Fullscreen", exact: true }).click();
    const exit = page.getByRole("button", {
      name: "Exit fullscreen",
      exact: true,
    });
    await expect(exit).toBeVisible();
    await fillsScreen(page);
    await page.screenshot({ path: testInfo.outputPath("fullscreen.png") });
    await page.getByRole("button", { name: "Ribbon", exact: true }).click();
    await expect(viewport).toHaveAttribute("aria-busy", "false");
    await fillsScreen(page);
    await page.setViewportSize({ width: 390, height: 844 });
    await fillsScreen(page);
    await expect(exit).toBeInViewport();
    await page.screenshot({
      path: testInfo.outputPath("fullscreen-mobile.png"),
    });
    await page.setViewportSize({ width: 1440, height: 900 });
    await exit.click();
    await expect
      .poll(async () => Math.round((await viewport.boundingBox())!.height))
      .toBe(Math.round(embedded!.height));
    await expect(
      page.getByRole("button", { name: "Fullscreen", exact: true }),
    ).toBeInViewport();
  });
}
