#!/bin/bash
# Daily region 2t for the 9 km and 31 km runs (interannual transfer test).
#SBATCH --job-name=interannual_ours
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/interannual_extract_ours_%j.log
set -uo pipefail
P=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z
export PATH=$P/bin:$PATH LD_LIBRARY_PATH=$P/lib OPENBLAS_NUM_THREADS=1
R=/work/ab0246/a270092/software/dark-doldrums-project
echo "=== $(date) tco319"
N_WORKERS=20 python3 -u $R/dunkelflaute/scripts/regimes/interannual_extract_ours.py tco319 || echo FAILED tco319
echo "=== $(date) dart9km"
N_WORKERS=48 python3 -u $R/dunkelflaute/scripts/regimes/interannual_extract_ours.py dart9km || echo FAILED dart9km
echo "=== $(date) DONE"
