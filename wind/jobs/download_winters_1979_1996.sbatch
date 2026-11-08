#!/bin/bash
#SBATCH --job-name=era5_wind_1979_1996
#SBATCH --partition=shared
#SBATCH --account=ab0995
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/wind/logs/download_winters_1979_1996_%j.log

set -euo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/wind

# Winters 1979-80 .. 1995-96 (17 winters), 5 parallel CDS requests at a time.
# Deliberately stops at 1995-96: winters 1996-97 .. 2005-06 are covered by the
# already-running download_all_winters.py 1996 2005 job -- this range was
# chosen to not overlap it.
python3 -u scripts/download/download_winters_parallel.py 1979 1995

echo "=== $(date) DONE ==="
