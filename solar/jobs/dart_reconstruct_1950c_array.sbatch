#!/bin/bash
#SBATCH --job-name=dart_recon_1950c
#SBATCH --partition=compute
#SBATCH --account=ab0995
#SBATCH --array=0-19%20
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --time=05:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/dart_recon_%x_%a_%j.log
set -euo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
export OMP_NUM_THREADS=64
cd /work/ab0246/a270092/software/dark-doldrums-project/solar
YEAR=$((1950 + SLURM_ARRAY_TASK_ID))
echo "task $SLURM_ARRAY_TASK_ID -> YEAR=$YEAR"
python3 -u scripts/evaluation/dart_reconstruct_year.py /work/ab0995/ICCP_AWI_hackthon_2025/TCo1279-DART-1950C "$YEAR"
