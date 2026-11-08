#!/bin/bash
#SBATCH --job-name=era5_hourly_backfill
#SBATCH --partition=shared
#SBATCH --account=ab0995
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/backfill_hourly_2015_2025_%j.log

set -uo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar

# Fetches the 16 hourly timestamps/day missing from era5_ssrd_3h_*.nc /
# era5_tsr_3h_*.nc for 2015-2025, merges them with the existing 8/day samples,
# and writes the full hourly (24/day) result to data/era5_1h/ -- separate from,
# and without touching, the original era5/ 3-hourly files (tisr included).
# Idempotent: already-merged months (full 24/day already in era5_1h/) are skipped,
# so this is safe to resubmit if it dies partway through.
python3 -u scripts/download/backfill_hourly_2015_2025.py

echo "=== $(date) HOURLY BACKFILL DONE ==="
