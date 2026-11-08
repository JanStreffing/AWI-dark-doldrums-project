#!/bin/bash
#SBATCH --job-name=extract_old_proxy
#SBATCH --partition=shared
#SBATCH --account=ab0995
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/extract_old_proxy_%j.log
set -uo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar
python3 -u scripts/preprocessing/extract_old_proxy_ondjf.py
echo "=== $(date) OLD PROXY EXTRACTION COMPLETE ==="
