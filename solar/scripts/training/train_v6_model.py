"""v6: v5's method (HistGradientBoostingRegressor on kt = ssrd/tisr, 250 iterations, cos(lat)
weights, no cells south of 60S or above 3000 m, Sep-Mar 2015-2026, same three year-block CV
folds) with the reduced 6-hourly feature set the TCo319/TCO95 runs can supply:
[tcc, lcc, tcc_prev, lcc_prev, mu] over 6 h windows (see core/reconstruct_ssrd_v6.py).

v7 (--variant v7) is the same with all four cloud layers at both window ends (FEATURES_V7).

Usage: python3 train_v6_model.py [--variant v6|v7] [--rows 32000000] [--only fold1_test_2015-2018]
Output: models/model_<variant>/cv_folds/<fold>/model_kt_<variant>.joblib and models/model_<variant>/model_kt_<variant>.joblib"""

import os as _os, sys as _sys
_root = _os.path.dirname(_os.path.abspath(__file__))
while not _os.path.isdir(_os.path.join(_root, "data")):
    _root = _os.path.dirname(_root)
_os.chdir(_root)
_sys.path.insert(0, _os.path.join(_root, "scripts", "core"))
_sys.path.insert(0, _os.path.join(_root, "scripts", "preprocessing"))

import argparse, glob, json, os, time
import joblib
import numpy as np
import sklearn
from sklearn.ensemble import HistGradientBoostingRegressor
from extract_v6_6h import CACHE, FEATURES, VARIANT    # both follow --variant on the command line
from reconstruct_ssrd_v6 import DT6

MONTHS = (9, 10, 11, 12, 1, 2, 3)
YEARS = range(2015, 2027)
FOLDS = {"fold1_test_2015-2018": range(2015, 2019), "fold2_test_2019-2022": range(2019, 2023),
         "fold3_test_2023-2026": range(2023, 2027), "final": range(0)}
LAT_MIN, ELEV_MAX = -60.0, 3000
ROOT = f"models/model_{VARIANT}"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=int, default=32_000_000)
    ap.add_argument("--only", default="")
    ap.add_argument("--variant", default="v6")
    args = ap.parse_args()
    all_cached = {os.path.basename(f)[:-4] for f in glob.glob(f"{CACHE}/*.npz")}
    stamps_all = sorted(s for s in all_cached if int(s[:4]) in YEARS and int(s[4:6]) in MONTHS)
    log(f"{len(stamps_all)} Sep-Mar 2015-2026 months in cache")

    for name, test_years in FOLDS.items():
        if args.only and name != args.only:
            continue
        stamps = [s for s in stamps_all if int(s[:4]) not in test_years]
        parts = [np.load(f"{CACHE}/{s}.npz") for s in stamps]
        X = np.concatenate([p["X"] for p in parts]); y = np.concatenate([p["y"] for p in parts])
        lat = np.concatenate([p["lat"] for p in parts]); elev = np.concatenate([p["elev"] for p in parts])
        n0 = len(y)
        keep = (lat >= LAT_MIN) & (elev <= ELEV_MAX)
        X, y, lat = X[keep], y[keep], lat[keep]
        log(f"[{name}] {len(stamps)} months, {n0:,} rows -> {len(y):,} after the lat/elevation mask")
        if len(y) > args.rows:
            sel = np.random.default_rng(0).choice(len(y), args.rows, replace=False)
            X, y, lat = X[sel], y[sel], lat[sel]
        log(f"[{name}] training rows: {len(y):,}")

        w = np.cos(np.deg2rad(lat)).astype("float64")
        t = time.time()
        model = HistGradientBoostingRegressor(max_iter=250, random_state=0).fit(X, y, sample_weight=w)
        out = ROOT + ("" if name == "final" else f"/cv_folds/{name}")
        os.makedirs(out, exist_ok=True)
        joblib.dump({"model": model, "features": FEATURES, "dt_seconds": DT6, "train_months": stamps, "variant": "area"},
                    f"{out}/model_kt_{VARIANT}.joblib")
        json.dump({"variant": "area", "months": "Sep-Mar 2015-2026", "train_years_excluded": sorted(test_years),
                   "n_months": len(stamps), "rows_used": int(len(y)), "mask": "lat >= -60 and elevation <= 3000 m",
                   "weights": "cos(lat)", "max_iter": 250, "n_iter_": int(model.n_iter_), "sklearn": sklearn.__version__,
                   "features": FEATURES, "dt_seconds": DT6,
                   "note": "6 h windows; clouds at window start and end only, as in 6-hourly TCo319/TCO95 output"},
                  open(f"{out}/training_config.json", "w"), indent=1)
        log(f"[{name}] fit {((time.time()-t)/60):.1f} min, n_iter_={model.n_iter_} -> {out}")


if __name__ == "__main__":
    main()
