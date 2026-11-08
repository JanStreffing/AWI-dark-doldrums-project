#!/bin/bash
#SBATCH --job-name=era5_eval
#SBATCH --partition=compute
#SBATCH --account=ab0995
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=00:40:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/era5_eval_%j.log

set -euo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar

python3 -u scripts/evaluation/evaluate_era5_regional.py
