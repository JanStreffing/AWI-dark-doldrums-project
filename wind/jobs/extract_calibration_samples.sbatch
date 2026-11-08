#!/bin/bash
#SBATCH --job-name=wind_calib_samples
#SBATCH --partition=shared
#SBATCH --account=ab0995
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/wind/logs/extract_calibration_samples_%j.log

# Germany per-cell ERA5 samples (v10, log-law v100, real v100, fsr) for the
# log-law -> ERA5 100 m calibration: Sep-Mar 2015-2026 (calibration) and the
# Oct-Feb 1979-2014 winters (independent test).
# analysis-toolbox/python-04.2026 (used by the older jobs) is not available to
# a270092; python3/2025.01 has xarray, shapely 2, scipy, sklearn 1.6.
set -uo pipefail
cd /work/ab0246/a270092/software/dark-doldrums-project/wind
module load python3/2025.01-gcc-13.3.0
export N_WORKERS=16

python3 scripts/methods/calibration/extract_germany_era5_samples.py recent
python3 scripts/methods/calibration/extract_germany_era5_samples.py historic
