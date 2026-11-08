#!/bin/bash
#SBATCH --job-name=dart_100m_2080c
#SBATCH --partition=shared
#SBATCH --account=ab0995
#SBATCH --array=0-90%20
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/wind/logs/dart_reconstruct_2080c_%A_%a.log

# 91 tasks = Sep/Oct/Nov/Dec/Jan/Feb/Mar calendar-month files across DART years
# 2080-2092 (7 months x 13 years, Sep-Mar extended winter). Same method/tuning as
# the 1950C array job: chunk=7, %20 throttle matching the account's MaxJobs limit.
set -uo pipefail
cd /work/ab0246/a270092/software/dark-doldrums-project/wind
module load analysis-toolbox/python-04.2026

MONTHS=(9 10 11 12 1 2 3)
YEAR_IDX=$(( SLURM_ARRAY_TASK_ID / 7 ))
MONTH_IDX=$(( SLURM_ARRAY_TASK_ID % 7 ))
YEAR=$(( 2080 + YEAR_IDX ))
MONTH=${MONTHS[$MONTH_IDX]}

OUTFILE="data/dart_2080c_100m/u100v100_${YEAR}$(printf '%02d' $MONTH).nc"
if [ -s "$OUTFILE" ]; then
  echo "task $SLURM_ARRAY_TASK_ID: $OUTFILE already exists, skipping"
  exit 0
fi

echo "task $SLURM_ARRAY_TASK_ID -> YEAR=$YEAR MONTH=$MONTH"
python3 scripts/methods/dart_loglaw/reconstruct_month.py $YEAR $MONTH --run 2080c --chunk 7
