#!/usr/bin/env bash
# Activate the selected HKOCA conda env, then run hkoca (or an explicit command).
set -euo pipefail

readonly ALLOWED_ENVS="hkoca_harmonize hkoca_cellbender sc_qc_pipeline"
HKOCA_ENV="${HKOCA_ENV:-hkoca_harmonize}"

case " ${ALLOWED_ENVS} " in
    *" ${HKOCA_ENV} "*) ;;
    *)
        echo "error: unknown HKOCA_ENV='${HKOCA_ENV}' (expected one of: ${ALLOWED_ENVS})" >&2
        exit 2
        ;;
esac

# shellcheck source=/dev/null
source /opt/conda/etc/profile.d/conda.sh
conda activate "${HKOCA_ENV}"

# Default: treat argv as hkoca subcommands/flags.
# Use `--` to run an arbitrary command inside the env, e.g.:
#   docker run --rm hkoca:local -- bash
#   docker run --rm hkoca:local -- Rscript -e 'packageVersion("Seurat")'
if [[ "${1:-}" == "--" ]]; then
    shift
    exec "$@"
fi

if [[ "${1:-}" == "hkoca" ]]; then
    shift
fi

exec hkoca "$@"
