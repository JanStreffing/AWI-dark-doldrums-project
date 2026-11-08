#!/bin/bash
# Run one script from dunkelflaute/scripts/ in the legacy env. Usage: sbatch regions_py.sbatch SUBDIR/SCRIPT.py [ARGS...]
#SBATCH --job-name=reg_py
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/regions_py_%x_%j.log
set -euo pipefail
P=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z
export PATH=$P/bin:$PATH LD_LIBRARY_PATH=$P/lib OPENBLAS_NUM_THREADS=1
R=/work/ab0246/a270092/software/dark-doldrums-project
S=$1; shift
echo "=== $(date) $S $* on $(hostname)"
python3 -u $R/dunkelflaute/scripts/$S "$@"
echo "=== $(date) DONE"
