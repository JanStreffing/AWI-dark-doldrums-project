#!/bin/bash
#SBATCH --job-name=era5_download_chain3
#SBATCH --partition=shared
#SBATCH --account=ab0995
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/download_chain3_%j.log

# Self-healing: sweep 2016-2025 first (existing files are skipped, so this only
# fills gaps left by the previous job), then 2015 -> 2006, then a final retry
# pass over everything to pick up any file that failed on its first attempt.
set -uo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar

run_year () {
  local year=$1 months=""
  for m in 01 02 03 04 05 06 07 08 09 10 11 12; do months="$months ${year}${m}"; done
  echo "=== $(date) ${year} ==="
  ./scripts/download/download_era5.sh $months || echo "!!! ${year}: some files failed, will retry in final pass"
}

echo "##### PASS 0: gap-fill 2016-2025 #####"
for y in 2025 2024 2023 2022 2021 2020 2019 2018 2017 2016; do run_year $y; done

echo "##### PASS 1: 2015 -> 2006 #####"
for y in 2015 2014 2013 2012 2011 2010 2009 2008 2007 2006; do run_year $y; done

echo "##### PASS 2: final retry over 2016-2025 and 2006-2015 #####"
for y in $(seq 2025 -1 2006); do run_year $y; done

echo "=== $(date) CHAIN 3 COMPLETE ==="
