#!/bin/bash
# SH winter (Apr-Sep) daily wind and sea-level pressure (scripts/regions/sh_extract.py).
# Usage: sbatch sh_extract.sbatch SOURCE:NAME [SOURCE:NAME ...]
#SBATCH --job-name=sh_extract
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/sh_extract_%j.log
set -uo pipefail
P=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z
export PATH=$P/bin:$PATH LD_LIBRARY_PATH=$P/lib OPENBLAS_NUM_THREADS=1
R=/work/ab0246/a270092/software/dark-doldrums-project
for SN in "$@"; do
  echo "=== $(date) $SN"
  python3 -u $R/dunkelflaute/scripts/regions/sh_extract.py ${SN%%:*} ${SN#*:} || echo "FAILED $SN"
done
echo "=== $(date) DONE"
