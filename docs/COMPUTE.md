# Current compute operation

Local CPU work includes input validation, chain/residue mapping, upstream package import,
HTTP/client integration, saved-result export, geometry checks and optional bounded
Foldseek retrieval. Fresh BoltzGen design and prediction use a CUDA GPU.

The tested worker uses the official `runpod-torch-v280` image, 100 GB local container disk,
a 100 GB network volume at `/workspace`, and an 80 GB A100 or H100 with a CUDA 13-capable
host driver. The main environment uses Torch 2.8; upstream ESMFold2 uses its isolated
Torch 2.11/CUDA13 runtime. Model source is pinned, standard Torch kernels are used, and
OMP/MKL thread counts are capped at four.

Python environments execute on local disk. Model files and a compressed, checksummed
runtime snapshot persist on the network volume. A fresh A100 restored the runtime,
passed both CUDA checks and reused the weights in 76 seconds. This is startup readiness,
not model inference time. See [the operational instructions and measurements](WARM_COMPUTE.md).

Both test GPU pods were deleted after results were downloaded. The 100 GB cache remains
in US-MO-1 and costs approximately $7/month until deleted. Fresh inference requires
provisioning and configuring another worker; the API does not silently create one.
Keep a pod-scoped spending guard, collect partial outputs on failure, and verify cleanup.

The current API binds to loopback. Feed/API deployment and ownership of worker lifecycle
remain team integration decisions. Credentials, dependency caches and live run directories
are excluded from the Git repository and handoff archive. Public benchmark results are
retained under `examples/boltzgen_results/` and `examples/boltzgen_api_results/`.

Older Mosaic/RFD3 setup notes are [historical](COMPUTE_HISTORY.md).
