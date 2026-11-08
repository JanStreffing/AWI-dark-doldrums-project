#!/bin/bash
#SBATCH --job-name=pm_tco95
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/perfect_model_tco95_%j.log
set -euo pipefail
# sklearn 1.1.1 for the v5 solar pickle; its shapely needs the environment's own libstdc++
P=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z
export PATH=$P/bin:$PATH LD_LIBRARY_PATH=$P/lib OPENBLAS_NUM_THREADS=1
cd /work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/scripts/tco95
for run in ctl1950g hi1950g; do
    N_WORKERS=60 python3 -u perfect_model_series.py $run 1952 2099
done
echo "=== $(date) DONE ==="
