#!/bin/bash
# Daily region series for several sources in sequence. Usage: sbatch cmip6_extract.sbatch SOURCE NAME [NAME ...]
#SBATCH --job-name=cmip6_extract
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/cmip6_extract_%x_%j.log
set -uo pipefail
P=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z
export PATH=$P/bin:$PATH LD_LIBRARY_PATH=$P/lib OPENBLAS_NUM_THREADS=1
R=/work/ab0246/a270092/software/dark-doldrums-project
SRC=$1; shift
for N in "$@"; do
  echo "=== $(date) $SRC $N"
  python3 -u $R/dunkelflaute/scripts/cmip6/extract_daily.py $SRC $N || echo "FAILED $SRC $N"
done
echo "=== $(date) DONE"
