#!/bin/bash
#SBATCH --job-name=df_dart_mockert
#SBATCH --partition=compute
#SBATCH --account=ab0995
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/df_dart_mockert_%j.log
set -euo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
export N_WORKERS=32
cd /work/ab0246/a270092/software/dark-doldrums-project
for run in 1950c 2080c; do
  python3 -u dunkelflaute/scripts/dart/compute_germany_dunkelflaute_dart_mockert.py $run
done
echo "=== $(date) DART MOCKERT COMPLETE ==="
