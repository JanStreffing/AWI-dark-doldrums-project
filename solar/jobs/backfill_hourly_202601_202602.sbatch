#!/bin/bash
#SBATCH --job-name=era5_hourly_2026
#SBATCH --partition=shared
#SBATCH --account=ab0995
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/backfill_hourly_202601_202602_%j.log

set -uo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar

# Same as backfill_hourly_2015_2025.py (job 27705089/27705418): fetches the
# missing hourly timestamps for ssrd/tsr, merges with the existing 3-hourly
# era5_ssrd_3h_*.nc / era5_tsr_3h_*.nc, writes to data/era5_1h/, then deletes
# the now-redundant 3-hourly originals (tisr included, per confirmed intent).
# Scoped to just Jan/Feb 2026 -- that's what already has 3-hourly source data.
python3 -u scripts/download/backfill_hourly_2015_2025.py 202601 202602

echo "=== $(date) 2026 HOURLY BACKFILL DONE ==="
