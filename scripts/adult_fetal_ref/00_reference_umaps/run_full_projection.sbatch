#!/bin/bash
#SBATCH --job-name=ova_wknn_map
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --time=24:00:00
#SBATCH --output=/data-master/pure-workspace/labss/hmami/atlas_shortcut_new_start/organoid_vs_primary_atlas/logs/slurm_mapping_full_%j.out
#SBATCH --error=/data-master/pure-workspace/labss/hmami/atlas_shortcut_new_start/organoid_vs_primary_atlas/logs/slurm_mapping_full_%j.err

set -euo pipefail
ROOT=/data-master/pure-workspace/labss/hmami/atlas_shortcut_new_start/organoid_vs_primary_atlas
PY=/data-master/pure-workspace/hmami/conda_envs/atlas/bin/python
mkdir -p "$ROOT/logs"
export MPLCONFIGDIR=/tmp/mpl_ova_map_${SLURM_JOB_ID:-$$}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export NUMBA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

cd "$ROOT"
WHICH=${1:-both}
echo "=== START full projection which=$WHICH $(date) job=${SLURM_JOB_ID:-NA} ==="
$PY -u scripts/00_reference_umaps/compute_mapping.py --which "$WHICH"
echo "=== DONE $(date) ==="
