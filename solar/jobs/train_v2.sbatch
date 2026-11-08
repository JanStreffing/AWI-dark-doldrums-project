#!/bin/bash
#SBATCH --job-name=train_v2
#SBATCH --partition=compute
#SBATCH --account=ab0995
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --output=/work/ab0246/a270092/software/dark-doldrums-project/solar/logs/train_v2_%j.log

# v2 = same method as v1, trained on ERA5 2015-2025.
# 3 year-block CV folds (for an honest evaluation) + the final model on everything.
set -euo pipefail
export PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin:$PATH
export OMP_NUM_THREADS=64
cd /work/ab0246/a270092/software/dark-doldrums-project/solar
F=models/model_v2/cv_folds

echo "##### extract features for all 132 months #####"
python3 -u scripts/training/train_v2_model.py --years 2015-2025 --stride 8 --workers 16 --cache-only

echo "##### fold 1: test 2015-2018 #####"
python3 -u scripts/training/train_v2_model.py --years 2015-2025 --exclude 2015-2018 --stride 8 --out $F/fold1_test_2015-2018
echo "##### fold 2: test 2019-2022 #####"
python3 -u scripts/training/train_v2_model.py --years 2015-2025 --exclude 2019-2022 --stride 8 --out $F/fold2_test_2019-2022
echo "##### fold 3: test 2023-2025 #####"
python3 -u scripts/training/train_v2_model.py --years 2015-2025 --exclude 2023-2025 --stride 8 --out $F/fold3_test_2023-2025
echo "##### FINAL model: all of 2015-2025 #####"
python3 -u scripts/training/train_v2_model.py --years 2015-2025 --stride 8 --out models/model_v2
echo "=== ALL DONE ==="
