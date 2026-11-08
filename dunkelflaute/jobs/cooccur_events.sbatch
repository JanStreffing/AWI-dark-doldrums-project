#!/bin/bash
# Europe-wide co-occurrence of low-wind events (8 regions), our runs, CMIP6, ERA5, TCO95 block null.
#SBATCH --job-name=cooccur_events
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/cooccur_events_%j.log
set -uo pipefail
module load python3/2025.01-gcc-13.3.0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
R=/work/ab0246/a270092/software/dark-doldrums-project
echo "=== $(date)"
NPROC=64 NPL=400 python3 -u $R/dunkelflaute/scripts/regions/cooccur_events.py || echo FAILED
echo "=== $(date) DONE"
