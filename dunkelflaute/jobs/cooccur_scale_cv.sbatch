#!/bin/bash
# Variance explanation of the aggregation-scale growth (CV, p5/mean, elasticity, mean-only part per scale).
#SBATCH --job-name=cooccur_scale
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/cooccur_scale_cv_%j.log
set -uo pipefail
module load python3/2025.01-gcc-13.3.0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
R=/work/ab0246/a270092/software/dark-doldrums-project
echo "=== $(date)"
NPROC=30 python3 -u $R/dunkelflaute/scripts/regions/cooccur_scale_cv.py || echo FAILED
echo "=== $(date) DONE"
