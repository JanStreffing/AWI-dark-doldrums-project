#!/bin/bash
#SBATCH --job-name=aggregate_hourly_3h
#SBATCH --partition=shared
#SBATCH --account=ab0995
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/aggregate_hourly_3h_%j.log
set -uo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar
python3 -u scripts/preprocessing/aggregate_hourly_to_3h.py
echo "=== $(date) AGGREGATION COMPLETE ==="
