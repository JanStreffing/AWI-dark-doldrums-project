#!/bin/bash
#SBATCH --job-name=tisr_era5_templates
#SBATCH --partition=compute
#SBATCH --account=ab0995
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/tisr_era5_build_%j.log

# Builds the two ERA5-grid tisr templates (365-day, 366-day), 3-hourly -- see
# build_tisr_templates_era5.py. Single process, one day (8 steps x 1.04M cells,
# ~70 MB float64) in memory at a time, so no exclusive node is needed.
#
# Submit from a login node:  sbatch jobs/build_tisr_templates_era5.sbatch

set -euo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

python3 -u scripts/preprocessing/build_tisr_templates_era5.py
