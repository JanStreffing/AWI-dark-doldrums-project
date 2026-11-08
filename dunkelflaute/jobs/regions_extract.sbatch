#!/bin/bash
# Region series for the multi-region test. Usage: sbatch regions_extract.sbatch MODEL RUN [N_WORKERS]
#   MODEL: dart (1950c|2080c), tco319 (ctl1950d|ssp585d), tco95 (ctl1950g|hi1950g), eerie (hist|ctl|ssp245)
#SBATCH --job-name=reg_extract
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/regions_extract_%x_%j.log
set -euo pipefail
# legacy env: eccodes (EERIE) and sklearn 1.1.1 (v7 solar pickle); its shapely needs its own libstdc++
P=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z
export PATH=$P/bin:$PATH LD_LIBRARY_PATH=$P/lib OPENBLAS_NUM_THREADS=1
R=/work/ab0246/a270092/software/dark-doldrums-project
MODEL=$1
RUN=$2
[ $# -ge 3 ] && export N_WORKERS=$3
echo "=== $(date) $MODEL $RUN on $(hostname)"
python3 -u $R/dunkelflaute/scripts/regions/extract_${MODEL}.py $RUN
echo "=== $(date) DONE"
