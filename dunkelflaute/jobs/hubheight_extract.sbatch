#!/bin/bash
# Hub-height test (hubheight_extract.py). Array task i runs models i, i+6, i+12, i+18 of the
# inventory (pressure-level footing), and the CFday model-level footing for the 5 models that have it.
# Usage: sbatch --array=0-5 hubheight_extract.sbatch
#SBATCH --job-name=hubheight
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/hubheight_%A_%a.log
set -uo pipefail
P=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z
export PATH=$P/bin:$PATH LD_LIBRARY_PATH=$P/lib OPENBLAS_NUM_THREADS=1 N_WORKERS=48
R=/work/ab0246/a270092/software/dark-doldrums-project
S=$R/dunkelflaute/scripts/regimes/hubheight_extract.py
MODELS=($(python3 -c "import json; print(' '.join(json.load(open('$R/dunkelflaute/data/cmip6/inventory.json'))))"))
ML=" AWI-ESM-1-REcoM MPI-ESM1-2-HR MPI-ESM1-2-LR IPSL-CM6A-LR NorESM2-LM "
i=$SLURM_ARRAY_TASK_ID
for j in $i $((i+6)) $((i+12)) $((i+18)); do
  [ $j -ge ${#MODELS[@]} ] && continue
  m=${MODELS[$j]}
  mode=pl; [[ "$ML" == *" $m "* ]] && mode=both
  echo "=== $(date) $m $mode"
  python3 -u $S $m $mode || echo "!!! FAILED $m"
done
echo "=== $(date) DONE"
