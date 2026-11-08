#!/bin/bash
#SBATCH --job-name=dart_wind_calib
#SBATCH --partition=shared
#SBATCH --account=ab0995
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=03:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/logs/apply_wind_calibration_dart_%j.log

# Apply the log-law -> ERA5 100 m quantile mapping to DART 1950C and 2080C (Germany),
# then compare raw vs calibrated Mockert events. python3/2025.01 because
# analysis-toolbox/python-04.2026 is not available to a270092.
set -uo pipefail
cd /work/ab0246/a270092/software/dark-doldrums-project
module load python3/2025.01-gcc-13.3.0
export N_WORKERS=32

python3 dunkelflaute/scripts/dart/apply_wind_calibration_dart.py 1950c
python3 dunkelflaute/scripts/dart/apply_wind_calibration_dart.py 2080c
python3 dunkelflaute/scripts/dart/apply_wind_calibration_dart.py summary
