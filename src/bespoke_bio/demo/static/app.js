"use strict";
const $ = (id) => document.getElementById(id);
let selection = null,
  job = null,
  offset = 0,
  total = 0,
  allowRun = false,
  timer = null,
  viewer = null;
const format = (value) => Number(value).toLocaleString();
function error(message) {
  $("alert").textContent = message;
  $("alert").hidden = false;
}
async function api(path, body) {
  const response = await fetch(
    path,
    body === undefined
      ? {}
      : {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
            "X-Bespoke-Request": "1",
          },
          body: JSON.stringify(body),
        },
  );
  const value = await response.json();
  if (!response.ok)
    throw Error(
      typeof value.detail === "string"
        ? value.detail
        : JSON.stringify(value.detail),
    );
  return value;
}
async function action(task) {
  $("alert").hidden = true;
  try {
    await task();
  } catch (e) {
    error(e.message);
  }
}
function text(tag, value, className) {
  const el = document.createElement(tag);
  el.textContent = value;
  if (className) el.className = className;
  return el;
}
async function search() {
  const params = new URLSearchParams({
    q: $("query").value,
    changed: $("nonreference").checked,
    offset,
    limit: 25,
  });
  const result = await api("/api/targets?" + params);
  total = result.total;
  $("page-info").textContent = total
    ? `${format(offset + 1)}–${format(Math.min(offset + 25, total))} of ${format(total)} matching records`
    : "No matching records.";
  $("targets").replaceChildren();
  for (const row of result.records) {
    const button = text(
      "button",
      "",
      "target" + (selection?.target_id === row.target_id ? " selected" : ""),
    );
    button.append(
      text("strong", row.gene || "Unnamed gene"),
      text("small", row.transcript),
      text(
        "small",
        `${row.donor} · ${format(row.length)} aa · ${row.changed ? "non-reference" : "reference-matching"}`,
        "allele",
      ),
    );
    button.setAttribute(
      "aria-label",
      `Select ${row.gene} ${row.transcript} ${row.target_id}`,
    );
    button.onclick = () => action(() => select(row.target_id));
    $("targets").append(button);
  }
  $("previous").disabled = offset === 0;
  $("next").disabled = offset + 25 >= total;
}
async function select(id) {
  selection = await api("/api/targets/" + encodeURIComponent(id));
  $("empty-target").hidden = true;
  $("target-detail").hidden = false;
  $("gene-name").textContent = selection.gene_symbol || selection.gene_id;
  $("selected-badge").textContent =
    selection.sequence_sha256 === selection.reference_sequence_sha256
      ? "REFERENCE-MATCHING"
      : "NON-REFERENCE";
  $("metadata").replaceChildren();
  for (const [name, value] of [
    ["Donor", selection.donor],
    ["Transcript", selection.transcript_id],
    ["Phase scenario", selection.phase_scenario],
    ["Length", `${format(selection.amino_acid_sequence.length)} amino acids`],
    [
      "Expression",
      selection.expression === null
        ? "Not measured"
        : JSON.stringify(selection.expression),
    ],
    ["Sequence SHA-256", selection.sequence_sha256],
  ]) {
    $("metadata").append(text("dt", name), text("dd", value));
  }
  $("sequence").textContent = selection.amino_acid_sequence;
  $("target-json").textContent = JSON.stringify(selection, null, 2);
  await search();
}
async function history() {
  const rows = await api("/api/jobs");
  $("history").replaceChildren();
  if (!rows.length) $("history").textContent = "No jobs yet.";
  for (const row of rows) {
    const b = text(
      "button",
      `${row.target_id || "Cached PD-L1 example"} · ${row.status}`,
      "secondary",
    );
    b.onclick = () => action(() => showJob(row));
    $("history").append(b);
  }
}
async function showJob(record) {
  clearTimeout(timer);
  job = record;
  $("job-detail").hidden = false;
  $("structure-panel").hidden = true;
  $("candidates").replaceChildren();
  $("download").hidden = true;
  $("results-empty").hidden = false;
  $("results-empty").textContent = "No candidate results for this job yet.";
  await refresh();
  await history();
}
async function refresh() {
  clearTimeout(timer);
  if (!job) return;
  const requested = job.job_id;
  const current = await api(`/api/jobs/${requested}`);
  if (job.job_id !== requested) return;
  job = current;
  $("job-status").textContent = job.status.replaceAll("_", " ");
  $("job-mode").textContent =
    job.mode === "actual_cached_replay"
      ? "CACHED REPLAY · NO NEW INFERENCE"
      : job.execution_performed
        ? "MODEL EXECUTION REPORTED"
        : "NO MODEL EXECUTION";
  $("job-target").textContent = job.target_id
    ? `Target: ${job.target_id}`
    : "Public PD-L1 example — separate from the selected donor protein.";
  $("progress").textContent = job.progress?.stage
    ? `Current stage: ${job.progress.stage}`
    : "";
  $("run").disabled =
    !allowRun ||
    job.dispatch_attempted ||
    job.status !== "prepared_pending_gpu_execution";
  $("cancel").disabled = !["queued", "running", "cancel_requested"].includes(
    job.status,
  );
  if (["completed", "cached_demo_completed"].includes(job.status)) {
    const results = await api(`/api/jobs/${requested}/results`);
    if (job.job_id !== requested) return;
    renderResults(results);
  } else if (
    ["queued", "running", "cancel_requested", "dispatch_requested"].includes(
      job.status,
    )
  ) {
    timer = setTimeout(() => action(refresh), 2500);
  }
}
function renderResults(result) {
  $("results-empty").hidden = result.candidates.length > 0;
  $("results-empty").textContent =
    "The job completed with no exported candidates.";
  $("download").hidden = false;
  $("download").href = result.download_url;
  $("candidates").replaceChildren();
  for (const c of result.candidates) {
    const card = text("article", "", "candidate");
    card.append(
      text("h3", c.id),
      text(
        "p",
        c.native_filter_passed === true
          ? "Native computational filters: passed"
          : c.native_filter_passed === false
            ? "Native computational filters: failed"
            : "Native filter outcome: unknown",
        c.native_filter_passed ? "pass" : "fail",
      ),
    );
    const metrics = text("div", "", "metrics");
    for (const [key, value] of Object.entries(c.confidence_metrics || {})) {
      const item = text("div", key.replaceAll("_", " "));
      item.append(
        text(
          "strong",
          value === null ? "Not available" : Number(value).toFixed(3),
        ),
      );
      metrics.append(item);
    }
    card.append(metrics, text("pre", c.sequence));
    if (c.failed_filters?.length)
      card.append(
        text("p", "Failed filters: " + c.failed_filters.join(", "), "muted"),
      );
    for (const [label, path] of [
      ["View complex", c.complex_structure_url],
      ["View binder", c.binder_structure_url],
    ]) {
      const b = text("button", label, "secondary");
      b.onclick = () => action(() => structure(path, label + " · " + c.id));
      card.append(b);
    }
    $("candidates").append(card);
  }
}
async function structure(path, title) {
  const response = await fetch(path);
  if (!response.ok) throw Error("Structure download failed");
  const contents = await response.text();
  $("structure-panel").hidden = false;
  $("structure-title").textContent = title;
  $("structure-download").href = path;
  if (!window.$3Dmol) {
    $("viewer-note").textContent =
      "The viewer library could not load. The structure download is available.";
    return;
  }
  if (!viewer)
    viewer = $3Dmol.createViewer($("viewer"), { backgroundColor: "#edf1ed" });
  viewer.clear();
  viewer.addModel(contents, path.endsWith(".pdb") ? "pdb" : "cif");
  if (!viewer.getModel().selectedAtoms({ atom: "CA" }).length)
    throw Error(
      "The viewer could not resolve the protein backbone. Download the original structure for inspection.",
    );
  viewer.setStyle({}, { cartoon: { colorscheme: "chain" } });
  viewer.zoomTo();
  viewer.resize();
  viewer.render();
  $("viewer-note").textContent =
    "Predicted coordinates · Drag to rotate · Scroll to zoom · Colors distinguish chains";
}
$("search").onsubmit = (e) => {
  e.preventDefault();
  offset = 0;
  action(search);
};
$("nonreference").onchange = () => {
  offset = 0;
  action(search);
};
$("previous").onclick = () => {
  offset = Math.max(0, offset - 25);
  action(search);
};
$("next").onclick = () => {
  offset += 25;
  action(search);
};
$("prepare").onclick = () =>
  action(async () => {
    if (!selection) return;
    $("prepare").disabled = true;
    try {
      await showJob(
        await api("/api/prepare", {
          target_id: selection.target_id,
          binder_length: Number($("binder-length").value),
          num_designs: Number($("design-count").value),
        }),
      );
    } finally {
      $("prepare").disabled = false;
    }
  });
$("cached").onclick = () =>
  action(async () => {
    $("cached").disabled = true;
    try {
      await showJob(await api("/api/demo", {}));
    } finally {
      $("cached").disabled = false;
    }
  });
$("refresh").onclick = () => action(refresh);
$("run").onclick = () =>
  action(async () => {
    if (
      !job ||
      !allowRun ||
      !confirm("Start a paid GPU design job for this prepared target?")
    )
      return;
    $("run").disabled = true;
    await api(`/api/jobs/${job.job_id}/run`, {});
    await refresh();
  });
$("cancel").onclick = () =>
  action(async () => {
    await api(`/api/jobs/${job.job_id}/cancel`, {});
    await refresh();
  });
action(async () => {
  const meta = await api("/api/catalog");
  allowRun = meta.gpu_execution_enabled;
  $("records").textContent = format(meta.target_records);
  $("unique").textContent = format(meta.unique_sequences);
  $("changed").textContent = format(meta.nonreference_records);
  $("mode").textContent = allowRun
    ? "GPU execution is available through an explicit start action."
    : "Integration mode · GPU execution disabled";
  if (!meta.service_configured)
    error(
      "Catalog is ready. Configure the design service to prepare targets or inspect cached results.",
    );
  await search();
  await history();
});
