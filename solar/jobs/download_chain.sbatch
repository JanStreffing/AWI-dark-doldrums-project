#!/bin/bash
#SBATCH --job-name=era5_download_chain2
#SBATCH --partition=shared
#SBATCH --account=ab0995
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/download_chain2_%j.log

set -euo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar

for year in 2022 2021 2020 2019 2018 2017 2016; do
  echo "=== $(date) starting ${year} (full year) ==="
  months=""
  for m in 01 02 03 04 05 06 07 08 09 10 11 12; do
    months="$months ${year}${m}"
  done
  ./scripts/download/download_era5.sh $months
done

echo "=== $(date) CHAIN COMPLETE ==="
