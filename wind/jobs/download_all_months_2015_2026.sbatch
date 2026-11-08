#!/bin/bash
#SBATCH --job-name=era5_wind_fullyear
#SBATCH --partition=shared
#SBATCH --account=ab0995
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/wind/logs/download_all_months_2015_2026_%j.log

set -uo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/wind

# Fills in Mar-Sep for 2015-2025 (currently winter-only: Jan/Feb/Oct/Nov/Dec)
# and whatever's published beyond Feb for 2026, for winds (u100,v100,u10,v10)
# and z0 (fsr). Skips existing files, safe to resubmit if interrupted.
python3 -u scripts/download/download_all_months_2015_2026.py 2015 2026

echo "=== $(date) FULL-YEAR 2015-2026 BACKFILL DONE ==="
