"""Redo the Germany v4/v5 reconstruction using, for each year, the CV fold model
that EXCLUDED that year from training -- fully held-out, never-seen-in-training,
for every single year -- instead of the 'final' model (trained on everything).
Answers: does using a properly held-out model change the result at all? (It
doesn't -- see solar/README.md Key Finding 3.) Promoted out of scratch/ since
its result is cited there; was previously an ad hoc check only.

Usage: python3 heldout_germany_v4v5.py v4|v5
Output: models/model_<variant>/evaluation/heldout_germany_<variant>.npz
"""
import os as _os, sys as _sys
_root = _os.path.dirname(_os.path.abspath(__file__))
while not _os.path.isdir(_os.path.join(_root, "data")):
    _root = _os.path.dirname(_root)
_os.chdir(_root)

import os, sys, time
import joblib
import numpy as np

sys.path.insert(0, "/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/scripts/era5")
from compute_germany_solar_cf_v5_sepmar import (load_month, reconstruct_solar, germany_land_mask,
                                                  month_available, MONTHS, YEARS)
from core.capacity_factor import solar_cf

TRUE3H_SECONDS = 3 * 3600.0
ROOT = "models"

FOLD_FOR_YEAR = {}
for y in range(2015, 2019): FOLD_FOR_YEAR[y] = "fold1_test_2015-2018"
for y in range(2019, 2023): FOLD_FOR_YEAR[y] = "fold2_test_2019-2022"
for y in range(2023, 2027): FOLD_FOR_YEAR[y] = "fold3_test_2023-2026"


def run(variant):
    mask = germany_land_mask()
    models = {}
    times, cf_recon, cf_real = [], [], []
    t0 = time.time()
    for y in YEARS:
        for m in MONTHS:
            if not month_available(y, m):
                continue
            fold = FOLD_FOR_YEAR[y]
            key = (variant, fold)
            if key not in models:
                models[key] = joblib.load(f"{ROOT}/model_{variant}/cv_folds/{fold}/model_kt_{variant}.joblib")["model"]
            model = models[key]
            mm = load_month(y, m)
            ssrd_recon, diag = reconstruct_solar(mm, model)
            w = np.cos(np.deg2rad(mm["lat"]))[:, None] * np.ones((1, len(mm["lon"]))) * mask
            wsum = w.sum()
            cr = (solar_cf(ssrd_recon, TRUE3H_SECONDS) * w[None, :, :]).sum(axis=(1, 2)) / wsum
            cl = (solar_cf(mm["ssrd_real"], TRUE3H_SECONDS) * w[None, :, :]).sum(axis=(1, 2)) / wsum
            times.append(mm["times"]); cf_recon.append(cr); cf_real.append(cl)
            print(f"  {variant} {y}-{m:02d} [{fold}] rescale_err={diag['max_rel_sum_err']:.1e} "
                  f"[elapsed {time.time()-t0:.0f}s]", flush=True)
    return np.concatenate(times), np.concatenate(cf_recon), np.concatenate(cf_real)


if __name__ == "__main__":
    variant = sys.argv[1]
    t, recon, real = run(variant)
    out = f"models/model_{variant}/evaluation/heldout_germany_{variant}.npz"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out, time=t, cf_recon=recon, cf_real=real)
    print(f"saved {out}, {len(t)} timesteps")
