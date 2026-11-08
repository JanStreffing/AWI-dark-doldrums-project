#!/bin/bash
#SBATCH --job-name=extract_all_v4
#SBATCH --partition=compute
#SBATCH --account=ab0995
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/extract_all_v4_%j.log
set -uo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
export NJOBS=64
cd /work/ab0246/a270092/software/dark-doldrums-project/solar
python3 -u scripts/preprocessing/extract_all_v4.py
echo "=== $(date) EXTRACTION COMPLETE ==="
