#!/bin/bash
# ERA5 hourly 10 m wind -> daily (00/06/12/18 UTC) region means and station-point series (interannual test).
#SBATCH --job-name=interannual_era5_wind
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/interannual_era5_wind_%j.log
set -uo pipefail
P=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z
export PATH=$P/bin:$PATH LD_LIBRARY_PATH=$P/lib OPENBLAS_NUM_THREADS=1
R=/work/ab0246/a270092/software/dark-doldrums-project
echo "=== $(date)"
N_WORKERS=64 python3 -u $R/dunkelflaute/scripts/regimes/interannual_extract_era5_wind.py || echo FAILED
echo "=== $(date) DONE"
