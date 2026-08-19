#!/usr/bin/env bash
# Run once after: conda env create -f conda/environment_harmonize.yaml
#             and: conda env create -f conda/environment_integration.yaml
#
# Installs kBET from GitHub into both R environments.
# kBET is not on conda-forge so it cannot be declared in the yaml directly.
#
# Usage:
#   bash conda/post_install_harmonize.sh
#   bash conda/post_install_harmonize.sh hkoca_harmonize hkoca_integration
#
# Inside an Apptainer/Singularity container where conda run is unavailable,
# run manually using the env's own Rscript:
#   /opt/conda/envs/hkoca_harmonize/bin/Rscript  -e "install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_github('theislab/kBET', upgrade='never')"
#   /opt/conda/envs/hkoca_integration/bin/Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_github('theislab/kBET', upgrade='never')"
set -euo pipefail

ENVS=("${@:-hkoca_harmonize hkoca_integration}")

_install_kbet() {
    local env_name="$1"
    local rscript
    rscript=$(conda run -n "$env_name" which Rscript 2>/dev/null || true)
    if [[ -z "$rscript" ]]; then
        echo "WARNING: Rscript not found in conda env '$env_name' — skipping." >&2
        return 0
    fi
    echo "[$env_name] Installing remotes if missing ..."
    conda run -n "$env_name" "$rscript" -e \
        "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='https://cloud.r-project.org')"
    echo "[$env_name] Installing kBET from GitHub ..."
    conda run -n "$env_name" "$rscript" -e \
        "remotes::install_github('theislab/kBET', upgrade='never', quiet=FALSE)"
    echo "[$env_name] kBET installed."
}

for env in "${ENVS[@]}"; do
    _install_kbet "$env"
done
