#!/bin/bash
#SBATCH --job-name=diag_training
#SBATCH --partition=compute
#SBATCH --account=ab0995
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/diag_training_%j.log
set -euo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
export OMP_NUM_THREADS=128
cd /work/ab0246/a270092/software/dark-doldrums-project/solar
python3 -u scripts/training/diag_training.py
