#!/bin/bash
# LUMIP land-use test: daily wind region series (scripts/regimes/lumip_extract_daily.py).
# Usage: sbatch lumip_extract_daily.sbatch
#SBATCH --job-name=lumip_daily
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/lumip_extract_daily_%j.log
set -uo pipefail
P=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z
export PATH=$P/bin:$PATH LD_LIBRARY_PATH=$P/lib OPENBLAS_NUM_THREADS=1 N_WORKERS=48
R=/work/ab0246/a270092/software/dark-doldrums-project
echo "=== $(date) on $(hostname)"
python3 -u $R/dunkelflaute/scripts/regimes/lumip_extract_daily.py
echo "=== $(date) DONE"
