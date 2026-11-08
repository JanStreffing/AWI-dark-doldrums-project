#!/bin/bash
#SBATCH --job-name=v3_pipeline
#SBATCH --partition=compute
#SBATCH --account=ab0995
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/v3_pipeline_%j.log
set -euo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar

echo "##### 1. extract ONDJF training rows (2015-2025, 55 months) #####"
OMP_NUM_THREADS=2 NJOBS=16 python3 -u scripts/preprocessing/extract_ondjf.py
echo "##### 2. train v3 (plain + area weights; 3 CV folds + final) #####"
OMP_NUM_THREADS=128 python3 -u scripts/training/train_v3_model.py --rows 32000000
echo "##### 3. evaluate the full pipeline (model + exact monthly rescale) #####"
OMP_NUM_THREADS=4 python3 -u scripts/evaluation/evaluate_v3_rescaled.py --workers 16
echo "=== ALL DONE ==="
