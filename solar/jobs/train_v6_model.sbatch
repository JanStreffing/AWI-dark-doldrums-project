#!/bin/bash
#SBATCH --job-name=train_v6_model
#SBATCH --partition=compute
#SBATCH --account=ab0995
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --time=06:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/train_v6_model_%j.log
set -euo pipefail
# sklearn 1.1.1, the environment v1-v5 were trained in (their pickles do not load in 1.6.1)
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar
NJOBS=16 python3 -u scripts/preprocessing/extract_v6_6h.py
python3 -u scripts/training/train_v6_model.py
echo "=== $(date) TRAINING COMPLETE ==="
