#!/bin/bash
# Local orchestrator: loop iterations (Python → R → metrics) then visualise.
# Reads ../config.yaml directly. No Slurm.

set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(cd "${SCRIPTS_DIR}/.." && pwd)"
CONFIG_YAML="${BENCH_DIR}/config.yaml"

mkdir -p "${BENCH_DIR}/logs" "${BENCH_DIR}/checkpoints" "${BENCH_DIR}/tmp" "${BENCH_DIR}/results"

# Resolve settings from config.yaml
eval "$(
  CONFIG_YAML="${CONFIG_YAML}" /data-master/pure-workspace/hmami/conda_envs/atlas/bin/python - <<'PY'
import os
import yaml

with open(os.environ["CONFIG_YAML"], encoding="utf-8") as fh:
    cfg = yaml.safe_load(fh)

methods = cfg["methods"]
order = [
    "pca", "seurat_v3", "seurat_v4", "scvi", "scanvi", "scpoli",
    "bbknn", "harmony", "combat", "scanorama", "css_pearson", "css_spearman",
]
enabled = [m for m in order if methods.get(m)]
r_methods = {"seurat_v3", "seurat_v4", "css_pearson", "css_spearman"}
gpu_methods = {"scvi", "scanvi", "scpoli"}

print(f'N_ITERS={int(cfg["benchmark"]["n_iters"])}')
print(f'PY="{cfg["paths"]["python"]}"')
print(f'RSCRIPT="{cfg["paths"].get("rscript", "Rscript")}"')
print(f'ANY_R={"1" if any(methods.get(m) for m in r_methods) else "0"}')
print(f'ANY_GPU={"1" if any(methods.get(m) for m in gpu_methods) else "0"}')
print(f'BENCH_METHODS="{",".join(enabled)}"')
print(f'ENABLED_LIST="{" ".join(enabled)}"')
PY
)"

if [[ -z "${BENCH_METHODS}" ]]; then
  echo "FATAL: no methods enabled in ${CONFIG_YAML}"
  exit 1
fi

echo "============================================================"
echo "HKOCA Level_2 benchmark (local orchestrator)"
echo "config:  ${CONFIG_YAML}"
echo "iters:   ${N_ITERS}"
echo "methods: ${ENABLED_LIST}"
echo "GPU methods enabled: ${ANY_GPU}"
echo "R methods enabled:   ${ANY_R}"
echo "============================================================"

export BENCH_METHODS
export BENCH_DIR
export CONFIG_YAML
export PYTHONUNBUFFERED=1
export TMPDIR="${BENCH_DIR}/tmp"
export MPLCONFIGDIR="${BENCH_DIR}/tmp/mpl"
export NUMBA_CACHE_DIR="${BENCH_DIR}/tmp/numba"
export JAX_PLATFORMS=cpu
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"

if [[ -f /data-master/pure-workspace/hmami/miniconda3/etc/profile.d/conda.sh ]]; then
  # shellcheck disable=SC1091
  source /data-master/pure-workspace/hmami/miniconda3/etc/profile.d/conda.sh
  conda activate atlas 2>/dev/null || true
fi
export LD_PRELOAD="${CONDA_PREFIX:-/data-master/pure-workspace/hmami/conda_envs/atlas}/lib/libstdc++.so.6"
export LD_LIBRARY_PATH="${CONDA_PREFIX:-/data-master/pure-workspace/hmami/conda_envs/atlas}/lib:${LD_LIBRARY_PATH:-}"

cd "${SCRIPTS_DIR}"

for ITER_ID in $(seq 1 "${N_ITERS}"); do
  export BENCHMARK_ITER_ID="${ITER_ID}"
  CHECKPOINT_DIR="${BENCH_DIR}/checkpoints/subsample_iter_${ITER_ID}"
  LOG="${BENCH_DIR}/logs/iter_${ITER_ID}.log"
  mkdir -p "${CHECKPOINT_DIR}"

  echo
  echo "---------- ITERATION ${ITER_ID}/${N_ITERS} ----------"
  {
    echo "Level_2 BENCHMARK — ITERATION ${ITER_ID}/${N_ITERS}"
    echo "host=$(hostname)  time=$(date -Is)"
    echo "methods=${BENCH_METHODS}"

    echo "[STEP 1] Python integrations"
    "${PY}" -u "${SCRIPTS_DIR}/run_integrations.py"

    if [[ "${ANY_R}" == "1" ]]; then
      echo "[STEP 2] R integrations (enabled subset)"
      "${RSCRIPT}" "${SCRIPTS_DIR}/run_r_integrations.R"
    else
      echo "[STEP 2] SKIP R — no R methods enabled"
    fi

    echo "[STEP 3] scIB metrics"
    "${PY}" -u "${SCRIPTS_DIR}/run_metrics.py"

    echo "SUCCESS iter ${ITER_ID}  $(date -Is)"
  } 2>&1 | tee "${LOG}"
done

echo
echo "[STEP 4] Aggregate + visualise"
"${PY}" -u "${SCRIPTS_DIR}/visualisation_results.py"

echo
echo "DONE $(date -Is)"
echo "Results: ${BENCH_DIR}/results/"
ls -lah "${BENCH_DIR}/results/"*.csv "${BENCH_DIR}/results/"*.png 2>/dev/null || true
