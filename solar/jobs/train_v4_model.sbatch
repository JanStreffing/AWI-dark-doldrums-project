#!/bin/bash
#SBATCH --job-name=train_v4_model
#SBATCH --partition=compute
#SBATCH --account=ab0995
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/train_v4_model_%j.log
set -euo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar
python3 -u scripts/training/train_v4_model.py
echo "=== $(date) TRAINING COMPLETE ==="
