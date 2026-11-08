#!/bin/bash
#SBATCH --job-name=era5_download_chain6_winter
#SBATCH --partition=shared
#SBATCH --account=ab0995
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/download_chain6_winter_%j.log

# Winter-only (Jan, Feb, Oct, Nov, Dec) backfill, 1991 -> 1979, then a retry pass.
# Runs in parallel with chain5 (2003 -> 1992) to cover the full 1979-2003 winter gap faster.
set -uo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar

run_year () {
  local year=$1 months=""
  for m in 01 02 10 11 12; do months="$months ${year}${m}"; done
  echo "=== $(date) ${year} (winter) ==="
  ./scripts/download/download_era5.sh $months || echo "!!! ${year}: some files failed, will retry in final pass"
}

echo "##### PASS 1: 1991 -> 1979 (winter months only) #####"
for y in $(seq 1991 -1 1979); do run_year $y; done

echo "##### PASS 2: final retry, 1991 -> 1979 #####"
for y in $(seq 1991 -1 1979); do run_year $y; done

echo "=== $(date) CHAIN 6 WINTER COMPLETE ==="
