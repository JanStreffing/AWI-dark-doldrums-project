#!/bin/bash
#SBATCH --job-name=era5_download_chain5_winter
#SBATCH --partition=shared
#SBATCH --account=ab0995
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/download_chain5_winter_%j.log

# Winter-only (Jan, Feb, Oct, Nov, Dec) backfill, 2003 -> 1992, then a retry pass.
# 2004-2025 are already complete full years (winter is a subset, nothing to do there).
set -uo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar

run_year () {
  local year=$1 months=""
  for m in 01 02 10 11 12; do months="$months ${year}${m}"; done
  echo "=== $(date) ${year} (winter) ==="
  ./scripts/download/download_era5.sh $months || echo "!!! ${year}: some files failed, will retry in final pass"
}

echo "##### PASS 1: 2003 -> 1992 (winter months only) #####"
for y in $(seq 2003 -1 1992); do run_year $y; done

echo "##### PASS 2: final retry, 2003 -> 1992 #####"
for y in $(seq 2003 -1 1992); do run_year $y; done

echo "=== $(date) CHAIN 5 WINTER COMPLETE ==="
