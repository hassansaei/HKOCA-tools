#!/bin/bash
#SBATCH -p gpu
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=1

apptainer exec \
    --nv \
    --bind /data-master/workspace/labss/akiselev/dir/HKOCA-tools:/out \
    /data-master/workspace/labss/akiselev/dir/HKOCA-tools/containers/hkoca-tools.sif \
    python /out/scripts/integration_scvi.py
