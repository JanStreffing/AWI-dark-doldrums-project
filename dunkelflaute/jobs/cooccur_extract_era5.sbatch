#!/bin/bash
# ERA5 hourly 10 m wind -> 6-hourly and daily onshore wind CF for the 8 European regions and the 5 deg boxes
# (co-occurrence and aggregation-scale validation).
#SBATCH --job-name=cooccur_era5
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/cooccur_era5_%j.log
set -uo pipefail
P=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z
export PATH=$P/bin:$PATH LD_LIBRARY_PATH=$P/lib OPENBLAS_NUM_THREADS=1
R=/work/ab0246/a270092/software/dark-doldrums-project
echo "=== $(date)"
N_WORKERS=96 python3 -u $R/dunkelflaute/scripts/regions/cooccur_extract_era5.py || echo FAILED
echo "=== $(date) DONE"
