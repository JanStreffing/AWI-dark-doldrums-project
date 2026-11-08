#!/bin/bash
# Run one dunkelflaute/scripts/ script once per argument, in sequence. Usage: sbatch cmip6_loop.sbatch SUBDIR/SCRIPT.py NAME [NAME ...]
#SBATCH --job-name=cmip6_loop
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/cmip6_loop_%x_%j.log
set -uo pipefail
P=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z
export PATH=$P/bin:$PATH LD_LIBRARY_PATH=$P/lib OPENBLAS_NUM_THREADS=1
R=/work/ab0246/a270092/software/dark-doldrums-project
S=$1; shift
for N in "$@"; do
  echo "=== $(date) $S $N"
  python3 -u $R/dunkelflaute/scripts/$S $N || echo "FAILED $S $N"
done
echo "=== $(date) DONE"
