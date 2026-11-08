#!/bin/bash
#SBATCH --job-name=dart_ssrd_test
#SBATCH --partition=compute
#SBATCH --account=ab0995
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=00:30:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/dart_test_%j.log

# One-month DART reconstruction test -- single process, no multiprocessing.
# First attempt on `shared` with 32G was OOM-killed: at DART's full 6.6M-cell
# resolution, 5 input variables + the tisr lookup + several derived arrays
# each run ~6.5 GB uncompressed, easily exceeding 32G for one month. `compute`
# with --mem=0 (same fix as the tisr build needed) gives the whole node's
# real memory instead of a partition default that doesn't fit this workload.

set -euo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
cd /work/ab0246/a270092/software/dark-doldrums-project/solar

python3 -u scripts/core/reconstruct_ssrd_dart.py /work/ab0995/ICCP_AWI_hackthon_2025/TCo1279-DART-1950C 195001
