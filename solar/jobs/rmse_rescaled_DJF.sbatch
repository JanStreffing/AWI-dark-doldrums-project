#!/bin/bash
#SBATCH --job-name=rmse_rescaled
#SBATCH --partition=compute
#SBATCH --account=ab0995
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --time=01:30:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/rmse_rescaled_%j.log
set -euo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
export OMP_NUM_THREADS=4
cd /work/ab0246/a270092/software/dark-doldrums-project/solar
python3 -u scripts/evaluation/rmse_maps_rescaled.py --workers 16 --months 12,1,2 --tag rescaled_DJF
