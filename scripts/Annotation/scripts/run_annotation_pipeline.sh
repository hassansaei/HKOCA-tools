#!/usr/bin/env bash
# Optional RDS→h5ad → concat → Leiden cluster → Snapseed v5
#
#   bash scripts/run_annotation.sh
#   bash scripts/run_annotation.sh --from concat
#   bash scripts/run_annotation.sh --only annotate
#   bash scripts/run_annotation.sh --dry-run
#
# RDS→h5ad is optional: if filtered singlet .h5ad files already exist
# (e.g. written by QC), conversion is skipped. If only .rds are present,
# they are converted.

set -euo pipefail

PIPE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS="${PIPE_ROOT}/scripts"
CONFIG="${PIPE_ROOT}/config"
MARKERS="${PIPE_ROOT}/markers/snapseed_markers_v5.yaml"
LOGDIR="${PIPE_ROOT}/logs"
mkdir -p "${LOGDIR}"

PY="${PY:-/data-master/pure-workspace/hmami/conda_envs/atlas/bin/python}"
RSCRIPT="${RSCRIPT:-$(command -v Rscript || true)}"
[[ -n "${RSCRIPT}" ]] || RSCRIPT="/data-master/pure-workspace/hmami/conda_envs/atlas/bin/Rscript"

H5AD_GLOB="*_harmonized_singlets_filtered.h5ad"

FROM="rds"
ONLY=""
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --from) FROM="${2:-}"; shift 2 ;;
    --only) ONLY="${2:-}"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

ts() { date -Is; }
log() { echo "[$(ts)] $*"; }

yaml_get() {
  "${PY}" - "$1" "$2" <<'PY'
import sys
from pathlib import Path
try:
    import yaml
    cfg = yaml.safe_load(Path(sys.argv[1]).read_text())
    cur = cfg
    for part in sys.argv[2].split("."):
        cur = cur[part]
    print(cur)
except Exception:
    key = sys.argv[2]
    for line in Path(sys.argv[1]).read_text().splitlines():
        s = line.strip()
        if s.startswith("#") or ":" not in s:
            continue
        k, v = s.split(":", 1)
        if k.strip() == key:
            print(v.strip().strip('"').strip("'"))
            raise SystemExit(0)
    raise SystemExit(f"key not found: {key}")
PY
}

run_cmd() {
  local title="$1"; shift
  log "START ${title}"
  log "CMD   $*"
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    log "dry-run skip"
    return 0
  fi
  local logfile="${LOGDIR}/${title}_$(date +%Y%m%d_%H%M%S).log"
  "$@" 2>&1 | tee "${logfile}"
  local st=${PIPESTATUS[0]}
  log "DONE  ${title} exit=${st} log=${logfile}"
  return "${st}"
}

should_run() {
  local step="$1"
  if [[ -n "${ONLY}" ]]; then
    [[ "${ONLY}" == "${step}" ]]; return
  fi
  case "${FROM}" in
    rds) return 0 ;;
    concat) [[ "${step}" != "rds" ]] ;;
    cluster) [[ "${step}" == "cluster" || "${step}" == "annotate" ]] ;;
    annotate) [[ "${step}" == "annotate" ]] ;;
    *) echo "Bad --from ${FROM} (rds|concat|cluster|annotate)" >&2; return 1 ;;
  esac
}

count_files() {
  local dir="$1" pattern="$2"
  if [[ ! -d "${dir}" ]]; then
    echo 0
    return
  fi
  find "${dir}" -maxdepth 1 -type f -name "${pattern}" 2>/dev/null | wc -l | tr -d ' '
}

sync_h5ad_to_concat() {
  local src="$1" dst="$2"
  [[ "${src}" == "${dst}" ]] && return 0
  [[ ! -d "${src}" ]] && return 0
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    log "dry-run would sync ${H5AD_GLOB} from ${src} → ${dst}"
    return 0
  fi
  mkdir -p "${dst}"
  local n
  n="$(find "${src}" -maxdepth 1 -type f -name "${H5AD_GLOB}" -exec cp -n {} "${dst}/" \; -print | wc -l | tr -d ' ')"
  log "synced ${n} h5ad file(s) → ${dst}"
}

RDS_CFG="${CONFIG}/rds_to_h5ad_parameters.yaml"
CONCAT_CFG="${CONFIG}/concatenation_parameters.yaml"
ANN_CFG="${CONFIG}/snapseed_annotation_parameters.yaml"

INPUT_RDS_DIR="$(yaml_get "${RDS_CFG}" input_rds_dir)"
OUTPUT_H5AD_DIR="$(yaml_get "${RDS_CFG}" output_h5ad_dir)"
CONCAT_INPUT_DIR="$(yaml_get "${CONCAT_CFG}" input_dir)"
CONCAT_OUT="$(yaml_get "${CONCAT_CFG}" output_h5ad)"
ANN_OUTDIR="$(yaml_get "${ANN_CFG}" output_dir)"
ANN_SCORE="$(yaml_get "${ANN_CFG}" score_flag_threshold)"
ANN_RES="$(yaml_get "${ANN_CFG}" resolution)"
CLUSTERED="$(yaml_get "${ANN_CFG}" clustered_h5ad)"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-8}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-8}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-8}"
export PYTHONUNBUFFERED=1
export LD_LIBRARY_PATH="/data-master/pure-workspace/hmami/conda_envs/atlas/lib:${LD_LIBRARY_PATH:-}"

log "annotation root=${PIPE_ROOT} from=${FROM} only=${ONLY:-all} dry_run=${DRY_RUN}"

if should_run rds; then
  N_H5_CONCAT="$(count_files "${CONCAT_INPUT_DIR}" "${H5AD_GLOB}")"
  N_H5_OUT="$(count_files "${OUTPUT_H5AD_DIR}" "${H5AD_GLOB}")"
  N_RDS="$(count_files "${INPUT_RDS_DIR}" "*.rds")"

  log "detect inputs: h5ad(concat_dir)=${N_H5_CONCAT} h5ad(convert_out)=${N_H5_OUT} rds=${N_RDS}"

  if [[ "${N_H5_CONCAT}" -gt 0 || "${N_H5_OUT}" -gt 0 ]]; then
    log "h5ad already present (e.g. from QC) — skip RDS→h5ad conversion"
    if [[ "${N_H5_CONCAT}" -eq 0 && "${N_H5_OUT}" -gt 0 ]]; then
      sync_h5ad_to_concat "${OUTPUT_H5AD_DIR}" "${CONCAT_INPUT_DIR}"
    fi
  elif [[ "${N_RDS}" -gt 0 ]]; then
    log "RDS detected, no filtered h5ad yet — running conversion"
    if [[ ! -x "${RSCRIPT}" && ! -f "${RSCRIPT}" ]]; then
      log "ERROR: Rscript not found (${RSCRIPT})"
      exit 1
    fi
    TMP_R="${LOGDIR}/rds_to_h5ad_configured_$$.R"
    cat > "${TMP_R}" <<REOF
suppressPackageStartupMessages({ library(Seurat); library(SeuratDisk) })
INPUT_DIR  <- "${INPUT_RDS_DIR}"
OUTPUT_DIR <- "${OUTPUT_H5AD_DIR}"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
rds_files <- list.files(INPUT_DIR, pattern = "\\\\.rds\$", full.names = TRUE)
if (length(rds_files) == 0) stop("No .rds files in: ", INPUT_DIR)
for (rds_path in rds_files) {
  base_name <- tools::file_path_sans_ext(basename(rds_path))
  h5seurat_path <- file.path(OUTPUT_DIR, paste0(base_name, ".h5Seurat"))
  h5ad_path <- file.path(OUTPUT_DIR, paste0(base_name, ".h5ad"))
  message("Processing: ", basename(rds_path))
  obj <- readRDS(rds_path)
  if (!inherits(obj, "Seurat")) next
  if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- as(obj[["RNA"]], "Assay")
  SaveH5Seurat(obj, filename = h5seurat_path, overwrite = TRUE)
  rm(obj); gc()
  Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)
  produced <- sub("\\\\.h5Seurat\$", ".h5ad", h5seurat_path)
  if (file.exists(produced) && normalizePath(produced) != normalizePath(h5ad_path, mustWork = FALSE))
    file.rename(produced, h5ad_path)
}
REOF
    run_cmd "1_rds_to_h5ad" "${RSCRIPT}" "${TMP_R}"
    rm -f "${TMP_R}"
    sync_h5ad_to_concat "${OUTPUT_H5AD_DIR}" "${CONCAT_INPUT_DIR}"
  else
    log "ERROR: no ${H5AD_GLOB} in concat/convert dirs and no .rds in ${INPUT_RDS_DIR}"
    exit 1
  fi
fi

if should_run concat; then
  run_cmd "2_concatenate" "${PY}" -u "${SCRIPTS}/concatenate_filtered_singlets.py"
fi

if should_run cluster; then
  mkdir -p "${ANN_OUTDIR}"
  run_cmd "3_cluster" \
    "${PY}" -u "${SCRIPTS}/annotate_snapseed_v5.py" \
      --input "${CONCAT_OUT}" \
      --markers "${MARKERS}" \
      --outdir "${ANN_OUTDIR}" \
      --clustered-out "${CLUSTERED}" \
      --resolution "${ANN_RES}" \
      --cluster-only
fi

if should_run annotate; then
  mkdir -p "${ANN_OUTDIR}"
  run_cmd "4_snapseed" \
    "${PY}" -u "${SCRIPTS}/annotate_snapseed_v5.py" \
      --input "${CLUSTERED}" \
      --markers "${MARKERS}" \
      --outdir "${ANN_OUTDIR}" \
      --score-threshold "${ANN_SCORE}"
fi

log "finished"
