#!/bin/bash
# TCO95 transfer-stability diagnostics. Usage: sbatch transfer_stability_extract.sbatch RUN YEARS
#SBATCH --job-name=transfer_stability
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/transfer_stability_extract_%j.log
set -uo pipefail
P=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z
export PATH=$P/bin:$PATH LD_LIBRARY_PATH=$P/lib OPENBLAS_NUM_THREADS=1 N_WORKERS=16
R=/work/ab0246/a270092/software/dark-doldrums-project
echo "=== $(date) $1 $2"
python3 -u $R/dunkelflaute/scripts/regimes/transfer_stability_extract.py $1 $2
echo "=== $(date) DONE"
