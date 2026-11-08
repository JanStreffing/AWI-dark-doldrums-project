"""
How much do the reduced 6-hourly v6 and v7 solar models lose against v5, over Germany?

Germany land (real border polygon), ERA5 full 0.25 deg resolution, every Sep-Mar month
2015-2026 on disk. Every month is scored by the CV-fold model that never saw its year
(fold1 2015-18, fold2 2019-22, fold3 2023-26), for BOTH models, so neither number is
in-sample:
  v5  3-hourly, full features [T, tcc, hcc, mcc, lcc, mu], exact monthly rescale,
      then summed onto the same 6 h windows (ssrd(t) + ssrd(t - 3 h))
  v6  6-hourly, [tcc, lcc, tcc_prev, lcc_prev, mu], exact monthly rescale (and the raw
      pre-rescale prediction, since TCO95 has no monthly ssrd to rescale against)
  v7  6-hourly, [tcc, lcc, mcc, hcc at both window ends, mu], exact monthly rescale
  real  ERA5 ssrd on the same 6 h windows
Solar CF = clip(ssrd / 21600 / 1000, 0, 1), cos(lat)-weighted mean over German land.

Then the event test that actually matters: Mockert combined CF at 6 h (48 h = 8 steps,
Mockert 2018 weights) with real ERA5 100 m wind CF sampled at the 6 h stamps, solar
swapped between real / v5 / v6 / v6-pre-rescale, events per winter across thresholds.

Usage: python3 evaluate_v6_germany.py
Output: solar/models/model_v6/evaluation/germany_v6_vs_v5_6h.{npz,json}
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "4")

import joblib
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/solar/scripts/core")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/era5")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/analysis")
from reconstruct_ssrd_v6 import DT6, FEATURES_V6, FEATURES_V7, to_6h_fields, tisr_6h, build_features_6h  # noqa: E402
from rescale import constrained_rescale                                        # noqa: E402
from core.capacity_factor import solar_cf                                      # noqa: E402
import compute_germany_solar_cf_v5_sepmar as v5g                               # noqa: E402
from robustness_mockert_events import trailing_mean, mockert_events, W_MOCKERT  # noqa: E402

os.chdir(_ROOT)
FOLDS = {"fold1_test_2015-2018": range(2015, 2019), "fold2_test_2019-2022": range(2019, 2023),
         "fold3_test_2023-2026": range(2023, 2027)}
OUT_DIR = "solar/models/model_v6/evaluation"
COMBINED_ERA5 = "dunkelflaute/results/validation_era5_2015_2026/combined_cf_2015_2026_sepmar_mockertweights.npz"
THRESHOLDS = (0.04, 0.05, 0.06, 0.07, 0.08, 0.10)
WIN6 = 8                       # 48 h at 6 h

_models, _mask = {}, None


def fold_of(year):
    return next(n for n, yrs in FOLDS.items() if year in yrs)


def model(variant, fold):
    key = (variant, fold)
    if key not in _models:
        _models[key] = joblib.load(f"solar/models/model_{variant}/cv_folds/{fold}/model_kt_{variant}.joblib")["model"]
    return _models[key]


def _init():
    global _mask
    _mask = v5g.germany_land_mask()


def twilight_weights(times, lon, nlat, hours_back):
    hours = (times - times.astype("datetime64[D]")) / np.timedelta64(1, "h") - hours_back
    local = (hours[:, None] + np.tile(lon, nlat)[None, :] / 15.0) % 24.0
    return np.maximum(np.cos(2 * np.pi * (local - 12.0) / 24.0), 0.0) + 1e-9


def _month(ym):
    year, month = ym
    t0 = time.time()
    fold = fold_of(year)
    m = v5g.load_month(year, month)
    nt3, nlat, nlon = m["tsr"].shape

    # v5, 3-hourly, held-out fold, then onto 6 h windows
    q5, _ = v5g.reconstruct_solar(m, model("v5", fold))
    clouds = {c: m[c] for c in ("tcc", "lcc", "mcc", "hcc")}
    w6 = to_6h_fields(m["times"], clouds, m["ssrd_real"])
    v5_on6 = to_6h_fields(m["times"], clouds, q5)["ssrd"]

    # v6 and v7, 6-hourly, held-out folds
    tisr6 = tisr_6h(m["lat"], m["lon"], w6["times"])
    nt6, ncell = len(w6["times"]), nlat * nlon
    real6 = w6["ssrd"].reshape(nt6, ncell)
    tw6 = twilight_weights(w6["times"], m["lon"], nlat, 3.0)
    pq = {}
    for variant, feats in (("v6", FEATURES_V6), ("v7", FEATURES_V7)):
        X, day, _ = build_features_6h(w6, tisr6, feats)
        p = np.zeros(nt6 * ncell)
        p[day] = np.clip(model(variant, fold).predict(X), 0, 1.5) * tisr6.ravel()[day]
        p = p.reshape(nt6, ncell)
        q, diag = constrained_rescale(p, tisr6.reshape(nt6, ncell).astype("float64"), real6.sum(axis=0),
                                      kmax=1.1, twilight_w=tw6)
        pq[variant] = (p, q)

    w = (np.cos(np.deg2rad(m["lat"]))[:, None] * np.ones((1, nlon)) * _mask).ravel()

    def cf(ssrd):
        return (solar_cf(ssrd.reshape(nt6, ncell), DT6) * w).sum(axis=1) / w.sum()

    res = dict(time=w6["times"], real=cf(real6), v5=cf(v5_on6), v6=cf(pq["v6"][1]), v6_pre=cf(pq["v6"][0]),
               v7=cf(pq["v7"][1]), v7_pre=cf(pq["v7"][0]),
               rescale_err=diag["max_rel_sum_err"])
    return year, month, res, time.time() - t0


def scores(pred, real, day):
    e = pred[day] - real[day]
    return dict(r=float(np.corrcoef(pred[day], real[day])[0, 1]), rmse_cf=float(np.sqrt(np.mean(e ** 2))),
                bias_cf=float(e.mean()), rmse_wm2=float(np.sqrt(np.mean(e ** 2)) * 1000),
                mean_pred=float(pred.mean()), mean_real=float(real.mean()))


def event_test(series):
    """Mockert events at 6 h with real ERA5 100 m wind and each solar variant."""
    d = np.load(COMBINED_ERA5, allow_pickle=True)
    t3 = d["time"].astype("datetime64[s]")
    idx = {v: i for i, v in enumerate(t3)}
    t6 = series["time"]
    keep = np.array([t in idx for t in t6])
    j = np.array([idx[t] for t in t6[keep]])
    on, off, winter = d["cf_on_real"][j], d["cf_off_real"][j], d["winter"][j]
    out = {}
    for name in ("real", "v5", "v6", "v6_pre", "v7", "v7_pre"):
        sol = series[name][keep]
        comb = W_MOCKERT["solar"] * sol + W_MOCKERT["on"] * on + W_MOCKERT["off"] * off
        out[name] = {}
        for thr in THRESHOLDS:
            per = []
            for y in np.unique(winter):
                msk = winter == y
                _, ev = mockert_events(trailing_mean(comb[msk], WIN6), WIN6, thr)
                per.append(len(ev))
            out[name][f"{thr:.2f}"] = dict(events_per_winter=float(np.mean(per)), per_winter=per)
    out["n_winters"] = int(len(np.unique(winter)))
    return out


def main():
    yms = [(y, mo) for y in range(2015, 2027) for mo in (9, 10, 11, 12, 1, 2, 3) if v5g.month_available(y, mo)]
    print(f"{len(yms)} Sep-Mar months", flush=True)
    t0 = time.time()
    res = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 16)), initializer=_init) as ex:
        for y, mo, r, dt in ex.map(_month, yms):
            res[(y, mo)] = r
            print(f"  {y}-{mo:02d} {dt:.0f}s rescale_err={r['rescale_err']:.1e} [{time.time()-t0:.0f}s]", flush=True)

    series = {k: np.concatenate([res[ym][k] for ym in sorted(res)]) for k in ("time", "real", "v5", "v6", "v6_pre", "v7", "v7_pre")}
    day = series["real"] > 0.005
    summary = dict(n_6h_steps=int(len(series["time"])), n_daylight_steps=int(day.sum()),
                   scores_vs_real_daylight={k: scores(series[k], series["real"], day) for k in ("v5", "v6", "v6_pre", "v7", "v7_pre")},
                   events_6h_real_wind=event_test(series))
    os.makedirs(OUT_DIR, exist_ok=True)
    np.savez(f"{OUT_DIR}/germany_v6_vs_v5_6h.npz", **series)
    json.dump(summary, open(f"{OUT_DIR}/germany_v6_vs_v5_6h.json", "w"), indent=1)

    for k, s in summary["scores_vs_real_daylight"].items():
        print(f"{k:7s} r={s['r']:.4f} RMSE={s['rmse_wm2']:.1f} W/m2 bias={1000*s['bias_cf']:+.2f} W/m2")
    ev = summary["events_6h_real_wind"]
    for thr in THRESHOLDS:
        print(f"thr {thr:.2f}: " + "  ".join(f"{k}={ev[k][f'{thr:.2f}']['events_per_winter']:.2f}"
                                             for k in ("real", "v5", "v6", "v6_pre", "v7", "v7_pre")))
    print(f"Saved {OUT_DIR}/germany_v6_vs_v5_6h.{{npz,json}}")


if __name__ == "__main__":
    main()
