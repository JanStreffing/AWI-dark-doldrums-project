#!/bin/bash
#SBATCH --job-name=era5_wind_gap
#SBATCH --partition=shared
#SBATCH --account=ab0995
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/wind/logs/download_gap_2004_2006_%j.log

set -euo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/wind

# Replaces the login-node sequential job (download_all_winters.py 1996 2005,
# PID 63551, killed 2026-09-23 ~12:xx): that job crawled because each ~2GB
# wind file has a real chance of dying mid-transfer (seen: FileNotFoundError
# inside multiurl's own retry-resume path, and an outright truncated-download
# AssertionError), and a single-connection sequential script pays a full fresh
# CDS queue wait (10-50 min) for every such failure. Parallel workers (as
# solar's download_era5.py already does successfully) overlap that queue wait
# instead of serializing it.
#
# Scoped to winters 2004-05 and 2005-06 only -- the exact gap; 1996-2004 and
# 2006+ already exist and are skipped automatically.
python3 -u scripts/download/download_winters_parallel.py 2004 2005

echo "=== $(date) GAP FILL DONE ==="
