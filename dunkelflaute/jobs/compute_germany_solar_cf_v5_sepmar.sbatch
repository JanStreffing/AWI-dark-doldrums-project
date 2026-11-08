#!/bin/bash
#SBATCH --job-name=germany_solar
#SBATCH --partition=compute
#SBATCH --account=ab0995
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/germany_solar_%x_%j.log
set -euo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project
python3 -u dunkelflaute/scripts/era5/compute_germany_solar_cf_v5_sepmar.py "$1"
echo "=== $(date) GERMANY $1 SOLAR CF COMPLETE ==="
