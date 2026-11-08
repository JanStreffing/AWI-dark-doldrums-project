#!/bin/bash
#SBATCH --job-name=tco319_pipe
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/tco319_pipeline_%j.log
set -euo pipefail
# sklearn 1.1.1 for the solar pickles; its shapely needs the environment's own libstdc++
P=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z
export PATH=$P/bin:$PATH LD_LIBRARY_PATH=$P/lib OPENBLAS_NUM_THREADS=1
R=/work/ab0246/a270092/software/dark-doldrums-project
echo "=== v7 vs v6 vs v5, Germany, held-out folds"
N_WORKERS=16 python3 -u $R/solar/scripts/evaluation/evaluate_v6_germany.py
echo "=== 31 km Germany series"
for run in ctl1950d ssp585d; do
    N_WORKERS=20 python3 -u $R/dunkelflaute/scripts/tco319/compute_germany_tco319.py $run
done
echo "=== events, 100 / 31 / 9 km"
python3 -u $R/dunkelflaute/scripts/tco319/tco319_events.py
echo "=== $(date) DONE ==="
