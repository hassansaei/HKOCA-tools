#!/usr/bin/env bash
# Run this once after: conda env create -f conda/environment_harmonize.yaml
# It installs R packages that are not available on conda-forge.
set -euo pipefail

ENV_NAME="${1:-hkoca_harmonize}"

RSCRIPT=$(conda run -n "$ENV_NAME" which Rscript 2>/dev/null || echo "")
if [[ -z "$RSCRIPT" ]]; then
    echo "ERROR: Rscript not found in conda env '$ENV_NAME'." >&2
    exit 1
fi

echo "Installing kBET from GitHub into env '$ENV_NAME' ..."
conda run -n "$ENV_NAME" "$RSCRIPT" -e \
    "remotes::install_github('theislab/kBET', upgrade = 'never', quiet = FALSE)"

echo "kBET installation complete."
