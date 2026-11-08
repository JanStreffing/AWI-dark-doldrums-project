#!/bin/bash
#SBATCH --job-name=era5_download_chain4
#SBATCH --partition=shared
#SBATCH --account=ab0995
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/download_chain4_%j.log

# 1996-2005 (2006 is handled by chain3), then a retry pass incl. 2006.
set -uo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar

run_year () {
  local year=$1 months=""
  for m in 01 02 03 04 05 06 07 08 09 10 11 12; do months="$months ${year}${m}"; done
  echo "=== $(date) ${year} ==="
  ./scripts/download/download_era5.sh $months || echo "!!! ${year}: some files failed, will retry in final pass"
}

echo "##### PASS 1: 2005 -> 1996 #####"
for y in $(seq 2005 -1 1996); do run_year $y; done

echo "##### PASS 2: final retry, 2006 -> 1996 #####"
for y in $(seq 2006 -1 1996); do run_year $y; done

echo "=== $(date) CHAIN 4 COMPLETE ==="
