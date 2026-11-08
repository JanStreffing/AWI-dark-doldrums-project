#!/bin/bash
# Inputs for the transfer artifact checks (scripts/regimes/transfer_checks*.py): global 2.5 deg msl
# for every run, then the DART raw 10 m speed. Usage: sbatch transfer_checks.sbatch
#SBATCH --job-name=transfer_checks
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/transfer_checks_%j.log
set -uo pipefail
P=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z
export PATH=$P/bin:$PATH LD_LIBRARY_PATH=$P/lib OPENBLAS_NUM_THREADS=1 N_WORKERS=16
R=/work/ab0246/a270092/software/dark-doldrums-project
S=$R/dunkelflaute/scripts/regimes
L=$R/dunkelflaute/logs
J=${SLURM_JOB_ID:-local}
run() { echo "=== $(date) start $*"; python3 -u "$@" > $L/transfer_checks_${J}_$(basename $1 .py)_${2}_${3:-}.log 2>&1 || echo "FAILED $*"; echo "=== $(date) end $*"; }
echo "=== $(date) on $(hostname)"
# group 1: small grids, in parallel
for a in "tco95 ctl1950g" "tco95 hi1950g" "tco319 ctl1950d" "tco319 ssp585d" "eerie ctl" "eerie hist" "eerie ssp245"; do
  run $S/transfer_checks_msl.py $a &
done
wait
# group 2: DART msl, both slices in parallel
run $S/transfer_checks_msl.py dart9km 1950c &
run $S/transfer_checks_msl.py dart9km 2080c &
wait
# group 3: DART raw 10 m speed, one slice at a time (memory)
run $S/transfer_checks_dart10m.py 1950c
run $S/transfer_checks_dart10m.py 2080c
echo "=== $(date) DONE"
