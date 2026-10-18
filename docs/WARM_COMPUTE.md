# Reusable BoltzGen compute

The executable entrypoint is `examples/boltzgen_smoke/bootstrap.sh warmup` on a Runpod pod. `bootstrap.sh` without an argument runs the existing two-design, no-epitope PD-L1 pipeline after setup.

Use a **100 GB network volume** mounted at `/workspace`, the official `runpod-torch-v280` template, an 80 GB A100 or H100, a **CUDA 13.0 minimum host driver**, and a 100 GB local container disk. The Boltz environment inherits Torch 2.8/CUDA12.8 from this exact template; the isolated ESM environment needs a CUDA13-compatible driver. Do not change template underneath a saved runtime snapshot.

The live 2026-10-17 catalog lists A100 SXM80 at $1.79/hour and H100 SXM80 at $3.99/hour. US-MO-1 lists both. US-MD-1 has A100 stock but **rejects network-volume creation**, so cannot host this cache. Availability changes; successful creation with the CUDA floor, not catalog GPU stock alone, confirms placement.

[Official standard network storage pricing](https://www.runpod.io/pricing) is $0.07/GB/month under 1 TB: 100 GB is approximately $7/month, or $0.0096/hour. It continues billing after pod deletion. Keep this volume for the hackathon; delete disposable GPU pods after outputs are downloaded.

## What persists

- `/workspace/cache/boltzgen`: BoltzGen, Boltz2, SolubleMPNN, chemical components.
- `/workspace/cache/huggingface`: pinned ESMFold2 and ESMC-6B snapshots, approximately 27 GB on first download.
- `/workspace/cache/runtimes/<revision>-torch280-py312/runtime.tar.zst`: versioned, SHA256-checked local-runtime snapshot. Warm restoration hashes the compressed stream during extraction in one network read and refuses execution if the digest differs.
- `/workspace/cache/warm_manifest.json`: actual downloaded file paths and sizes; explicitly does not claim inference success.
- Job inputs, logs and outputs in separate `/workspace/jobs/<job-id>` directories.

Package installation happens on local `/opt` and `/tmp`, avoiding many small writes to the network filesystem. The first successful CUDA runtime check saves both environments as a compressed tar archive. Later pods restore it to exactly the original absolute paths, including any uv-managed Python installation. Runtime restoration is restricted to trusted account-local cache artifacts. No account credentials are stored in that archive.

Four OpenMP/MKL threads and one OpenBLAS thread avoid the severe CPU oversubscription seen in the initial test. Kernels stay disabled, matching the completed reference run.

## Checks and use

```bash
cd /workspace/boltzgen_smoke
bash bootstrap.sh warmup
# Input check, real CUDA allocations, and complete required model downloads.
# No new scientific sampling is done by warmup.
bash bootstrap.sh
# Full two-design reference pipeline; retains all seven stages.
```

`warmup` checks both environments' imports and GPU allocations. A cache manifest proves files exist and have nonzero size; the completed earlier PD-L1 run establishes this source revision's inference behavior. A new worker is only demonstrated end to end after a real job succeeds. Cold runtime and model cache preparation succeeded on H100 pod d6f79xbzrghsag in US-MO-1. The compressed runtime is about 6.7 GB; the populated cache occupies approximately 37 GB. Native imports and CUDA allocation checks passed in both environments. The original uncompressed snapshot hash was intentionally canceled and superseded by the verified compressed snapshot. A fresh A100 pod restored the checksum-verified snapshot, passed both CUDA probes and native input checks, and reused every model weight in **76 seconds** (2026-10-18 01:52:59–01:54:15 UTC). Cold preparation took **566 seconds**, including the retired uncompressed archive hash. This measured warm setup was approximately 7.4 times faster. Probe logs and manifests are in `outputs/warm_compute/restore/`; the probe pod was deleted and deletion confirmed. This is runtime readiness timing, not end-to-end binder design latency.

Do not share one output directory across concurrent jobs. Keep one GPU worker serial for the demo; queue incoming upstream target packages rather than loading multiple design models concurrently.

Both test GPU pods were deleted after results were downloaded, and both deletions returned 404 on read-back. The 100 GB cache volume `sjltm1vrwh` remains in US-MO-1. Estimated GPU charges for cold setup, live API inference and the independent warm restore were approximately $1.14 combined, plus small storage charges; finalized billing was not queried.
