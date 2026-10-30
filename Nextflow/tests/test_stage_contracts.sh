#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

fail=0
ok()  { echo "OK  $*"; }
bad() { echo "FAIL $*"; fail=1; }

echo "== Required layout =="
for path in \
  main.nf \
  nextflow.config \
  workflows/hkoca.nf \
  modules/cellbender.nf \
  modules/harmonization.nf \
  modules/qc.nf \
  modules/annotation.nf \
  modules/integration.nf \
  modules/projection.nf \
  modules/prepare_meta.nf \
  modules/prepare_harmonize_csv.py \
  conda/environment_cellbender.yaml \
  conda/environment_harmonize.yaml \
  conda/environment_qc.yaml \
  conda/environment_integration.yaml \
  conda/environment_projection.yaml \
  docs/AUDIT.md \
  docs/DEPENDENCY_MAP.md \
  docs/nextflow.md \
  use_cases/GSM9712682/README.md
do
  if [[ -e "$path" ]]; then ok "$path"; else bad "missing $path"; fi
done

echo "== prepare_harmonize_csv.py =="
TMP=$(mktemp -d)
cat > "$TMP/sample_info.csv" <<'CSV'
data_path,sample_id,study,source,diff_protocol,sc_protocol,sequencing,genome_build,Age,type,sample_dir,run_cellbender
/tmp/raw/sampleA/raw_feature_bc_matrix.h5,sampleA,Demo,Lab,organoid,10x,v3,GRCh38,30,organoid,/tmp/raw/sampleA,True
CSV
echo -e "sampleA\t/tmp/cb/sampleA_filtered_filtered.h5" > "$TMP/cb_map.tsv"
python3 modules/prepare_harmonize_csv.py \
  --sample-info "$TMP/sample_info.csv" \
  --cellbender-map "$TMP/cb_map.tsv" \
  --output "$TMP/harmonize_metadata.csv" \
  --working-dir "$TMP"
if grep -q '/tmp/cb/sampleA_filtered_filtered.h5' "$TMP/harmonize_metadata.csv"; then
  ok "CellBender path rewrite"
else
  bad "CellBender path rewrite"
fi
rm -rf "$TMP"

echo "== Nextflow config (optional) =="
if command -v nextflow >/dev/null 2>&1; then
  if nextflow config -profile conda >/dev/null 2>&1; then
    ok "nextflow config -profile conda"
  else
    echo "WARN nextflow config failed (often missing Java 17+); skip"
  fi
else
  echo "WARN nextflow not on PATH; skip"
fi

if [[ "$fail" -ne 0 ]]; then
  echo "Contract checks FAILED"
  exit 1
fi
echo "Contract checks PASSED"
