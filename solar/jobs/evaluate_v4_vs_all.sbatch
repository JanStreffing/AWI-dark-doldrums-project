#!/bin/bash
#SBATCH --job-name=evaluate_v4_vs_all
#SBATCH --partition=compute
#SBATCH --account=ab0995
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/evaluate_v4_vs_all_%j.log
set -uo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
# Measured: sklearn's predict() on 8M rows takes 60s with unrestricted threading
# (contention overhead at very high thread counts), 10.7s at 16 threads, >100s at
# 1 thread. 4 parallel month-workers x 16 threads each = 64 cores, no oversubscription.
export OMP_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
export MKL_NUM_THREADS=16
cd /work/ab0246/a270092/software/dark-doldrums-project/solar
python3 -u scripts/evaluation/evaluate_v4_vs_all.py --workers 4
echo "=== $(date) EVALUATION COMPLETE ==="
