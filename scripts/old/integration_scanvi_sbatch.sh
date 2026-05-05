#!/bin/bash
#SBATCH -p gpu
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=1

apptainer exec \
    --bind /data-master/workspace/labss/akiselev/dir/HKOCA-tools:/out \
    --nv \
    /data-master/workspace/labss/akiselev/dir/HKOCA-tools/containers/hkoca-tools.sif \
    python /out/scripts/integration_scanvi.py
