"""
Calibrate the log-law 100 m wind onto ERA5's native 100 m wind by conditional
quantile mapping, and test whether it fixes the Dunkelflaute event undercount.

Hub height stays fixed at 100 m onshore and offshore with the generic power curves
(core/capacity_factor.py); the target is ERA5's own u100/v100, not SMARD.

Mapping: for each bin, the empirical quantiles of log-law speed are mapped onto the
empirical quantiles of real speed (linear interpolation; beyond the outermost
quantiles the end offset is held constant). Bins use only quantities every model
grid has, so the same table can be applied to DART and the TCo319/TCo95 runs:
  month (Sep..Mar) x 3-hourly UTC slot (8) x surface class (4):
  offshore, onshore fsr < 0.1 m, 0.1-0.5 m, > 0.5 m (fixed edges, not terciles).

Tests:
  recent winters (Sep-Mar 2016-2026): 2-fold cross-validation by winter. Germany-mean
    onshore/offshore CF vs real, and Mockert combined-CF events with the existing v5
    ERA5 solar reconstruction, against the same events from real u100/v100.
  historic winters (Oct-Feb 1980-2014, never used for fitting): table fitted on all
    recent winters. No solar data exists there, so the event test is wind-only:
    48 h mean of capacity-weighted wind CF below thresholds.

Outputs:
  wind/models/quantile_mapping_v1/qm_table.npz        (fitted on all recent winters)
  dunkelflaute/results/wind_calibration/qm_evaluation.json

Usage: python calibrate_quantile_mapping.py
"""
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
os.chdir(_ROOT)
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts")
from core.capacity_factor import OFFSHORE, ONSHORE, wind_cf  # noqa: E402
from analysis.robustness_mockert_events import mockert_events, trailing_mean  # noqa: E402

SAMPLES = "wind/data/results/calibration/germany_era5_wind_samples_{}.npz"
VALIDATION = "dunkelflaute/results/validation_era5_2015_2026/combined_cf_2015_2026_sepmar_mockertweights.npz"
TABLE_OUT = "wind/models/quantile_mapping_v1/qm_table.npz"
EVAL_OUT = "dunkelflaute/results/wind_calibration/qm_evaluation.json"

MONTHS = (9, 10, 11, 12, 1, 2, 3)
SLOTS = 8                                   # 3-hourly UTC slots, 00..21
FSR_EDGES = (0.1, 0.5)                      # m; onshore roughness classes
N_CLASS = 4                                 # 0 offshore, 1..3 onshore by fsr
QS = np.concatenate([[0.001, 0.0025, 0.005], np.arange(0.01, 0.991, 0.01), [0.995, 0.9975, 0.999]])
W = dict(solar=0.44, on=0.50, off=0.06)     # Mockert's weights, as everywhere else
THRESHOLDS = (0.04, 0.05, 0.06, 0.07, 0.08, 0.10)
WINDOWS_H = (24, 48, 72)
WIND_THRESHOLDS = (0.05, 0.10, 0.15, 0.20)  # wind-only droughts, historic test


def load(which):
    d = np.load(SAMPLES.format(which))
    s = {k: d[k] for k in d.files}
    s["time"] = s["time"].astype("datetime64[s]")
    return s


def bins(s):
    """(T, N) integer bin index = (month_idx * SLOTS + slot) * N_CLASS + class."""
    t = s["time"]
    month = t.astype("datetime64[M]").astype(int) % 12 + 1
    lut = np.full(13, -1)
    lut[list(MONTHS)] = np.arange(len(MONTHS))
    m_idx = lut[month]
    assert (m_idx >= 0).all(), "sample outside Sep-Mar"
    slot = (t - t.astype("datetime64[D]")).astype(int) // 10800
    cls = np.where(s["land"][None, :], 1 + np.digitize(s["fsr"], FSR_EDGES), 0)
    return ((m_idx * SLOTS + slot)[:, None] * N_CLASS + cls).astype(np.int32)


def fit(s, tmask):
    b = bins(s)[tmask].ravel()
    x = s["v100_ll"][tmask].ravel()
    y = s["v100_real"][tmask].ravel()
    nb = len(MONTHS) * SLOTS * N_CLASS
    qx = np.full((nb, len(QS)), np.nan, np.float32)
    qy = np.full((nb, len(QS)), np.nan, np.float32)
    order = np.argsort(b, kind="stable")
    b, x, y = b[order], x[order], y[order]
    edges = np.searchsorted(b, np.arange(nb + 1))
    for k in range(nb):
        lo, hi = edges[k], edges[k + 1]
        if hi - lo >= 1000:
            qx[k] = np.quantile(x[lo:hi], QS)
            qy[k] = np.quantile(y[lo:hi], QS)
    return qx, qy


def apply(s, qx, qy):
    b = bins(s).ravel()
    x = s["v100_ll"].ravel()
    out = x.copy()
    order = np.argsort(b, kind="stable")
    edges = np.searchsorted(b[order], np.arange(qx.shape[0] + 1))
    for k in range(qx.shape[0]):
        idx = order[edges[k]:edges[k + 1]]
        if len(idx) == 0 or np.isnan(qx[k, 0]):
            continue
        xv = x[idx]
        y = np.interp(xv, qx[k], qy[k])
        y = np.where(xv < qx[k, 0], xv + (qy[k, 0] - qx[k, 0]), y)
        y = np.where(xv > qx[k, -1], xv + (qy[k, -1] - qx[k, -1]), y)
        out[idx] = np.maximum(y, 0.0)
    return out.reshape(s["v100_ll"].shape)


def germany_cf(s, v):
    w = np.cos(np.deg2rad(s["lat"]))
    on, off = s["onshore"], s["offshore"]
    cf_on = (wind_cf(v[:, on], **ONSHORE) * w[on]).sum(1) / w[on].sum()
    cf_off = (wind_cf(v[:, off], **OFFSHORE) * w[off]).sum(1) / w[off].sum()
    return cf_on, cf_off


def winter_of(t):
    y = t.astype("datetime64[Y]").astype(int) + 1970
    m = t.astype("datetime64[M]").astype(int) % 12 + 1
    return np.where(m >= 9, y + 1, y)


def events_per_winter(series, winter, win, thr):
    n = []
    for w_ in np.unique(winter):
        m = winter == w_
        _, ev = mockert_events(trailing_mean(series[m], win), win, thr)
        n.append(len(ev))
    return np.array(n)


def cf_stats(cf, ref):
    lo = ref < 0.05
    return dict(mean=float(cf.mean()), bias=float(cf.mean() - ref.mean()), r=float(np.corrcoef(cf, ref)[0, 1]),
                ratio_when_real_below_0p05=float(cf[lo].mean() / ref[lo].mean()))


def main():
    rec = load("recent")
    win_rec = winter_of(rec["time"])
    val = np.load(VALIDATION, allow_pickle=True)
    vt = val["time"].astype("datetime64[s]")
    keep = np.isin(rec["time"], vt)
    rec = {k: (v[keep] if k in ("time", "v10", "v100_ll", "v100_real", "fsr") else v) for k, v in rec.items()}
    win_rec = win_rec[keep]
    assert np.array_equal(rec["time"], vt), "recent samples do not cover the validation time axis"
    solar = val["cf_solar_recon"]

    # 2-fold CV by winter: each half is mapped with the table fitted on the other half
    winters = np.unique(win_rec)
    folds = [winters[::2], winters[1::2]]
    v_cal = np.empty_like(rec["v100_ll"])
    for test in folds:
        tm = np.isin(win_rec, test)
        qx, qy = fit(rec, ~tm)
        part = {k: (v[tm] if k in ("time", "v10", "v100_ll", "v100_real", "fsr") else v) for k, v in rec.items()}
        v_cal[tm] = apply(part, qx, qy)

    tracks = {"loglaw": rec["v100_ll"], "calibrated_cv": v_cal, "real": rec["v100_real"]}
    cf = {k: germany_cf(rec, v) for k, v in tracks.items()}
    comb = {k: W["solar"] * solar + W["on"] * c[0] + W["off"] * c[1] for k, c in cf.items()}

    ev = {}
    for win_h in WINDOWS_H:
        for thr in THRESHOLDS:
            ev[f"{win_h}h_{thr:.2f}"] = {k: float(events_per_winter(v, win_rec, win_h // 3, thr).mean())
                                         for k, v in comb.items()}
    pw = {k: events_per_winter(v, win_rec, 16, 0.06) for k, v in comb.items()}

    out = dict(
        method=dict(bins="month x 3h UTC slot x {offshore, onshore fsr<0.1, 0.1-0.5, >0.5 m}",
                    quantiles=len(QS), fsr_edges=FSR_EDGES, weights=W, target="ERA5 native 100 m wind"),
        recent=dict(
            winters=[int(w_) for w_ in winters], cv_folds=[[int(w_) for w_ in f] for f in folds],
            onshore_cf={k: cf_stats(c[0], cf["real"][0]) for k, c in cf.items()},
            offshore_cf={k: cf_stats(c[1], cf["real"][1]) for k, c in cf.items()},
            combined_mean={k: float(v.mean()) for k, v in comb.items()},
            events_per_winter=ev,
            per_winter_48h_0p06={k: v.tolist() for k, v in pw.items()},
            per_winter_r_vs_real={k: float(np.corrcoef(v, pw["real"])[0, 1]) for k, v in pw.items()},
        ),
    )

    # historic independent test, table fitted on all recent winters
    qx, qy = fit(rec, np.ones(len(rec["time"]), bool))
    os.makedirs(os.path.dirname(TABLE_OUT), exist_ok=True)
    np.savez(TABLE_OUT, qx=qx, qy=qy, quantiles=QS, months=np.array(MONTHS), slots=SLOTS,
             fsr_edges=np.array(FSR_EDGES), n_class=N_CLASS,
             fitted_on=np.array([str(rec["time"][0]), str(rec["time"][-1])]))
    if os.path.exists(SAMPLES.format("historic")):
        his = load("historic")
        month = his["time"].astype("datetime64[M]").astype(int) % 12 + 1
        keep = np.isin(month, MONTHS)
        his = {k: (v[keep] if k in ("time", "v10", "v100_ll", "v100_real", "fsr") else v) for k, v in his.items()}
        win_h_ = winter_of(his["time"])
        # only winters with all of Oct-Feb on disk
        full = [w_ for w_ in np.unique(win_h_) if np.sum(win_h_ == w_) >= 8 * 150]
        sel = np.isin(win_h_, full)
        his = {k: (v[sel] if k in ("time", "v10", "v100_ll", "v100_real", "fsr") else v) for k, v in his.items()}
        win_h_ = win_h_[sel]
        htracks = {"loglaw": his["v100_ll"], "calibrated": apply(his, qx, qy), "real": his["v100_real"]}
        hcf = {k: germany_cf(his, v) for k, v in htracks.items()}
        wind = {k: (W["on"] * c[0] + W["off"] * c[1]) / (W["on"] + W["off"]) for k, c in hcf.items()}
        hev = {f"48h_{thr:.2f}": {k: float(events_per_winter(v, win_h_, 16, thr).mean()) for k, v in wind.items()}
               for thr in WIND_THRESHOLDS}
        hpw = {k: events_per_winter(v, win_h_, 16, 0.10) for k, v in wind.items()}
        out["historic"] = dict(
            winters=[int(w_) for w_ in full], months="Oct-Feb only (what is on disk)",
            onshore_cf={k: cf_stats(c[0], hcf["real"][0]) for k, c in hcf.items()},
            offshore_cf={k: cf_stats(c[1], hcf["real"][1]) for k, c in hcf.items()},
            wind_only_events_per_winter=hev,
            per_winter_r_vs_real_48h_0p10={k: float(np.corrcoef(v, hpw["real"])[0, 1]) for k, v in hpw.items()},
        )

    os.makedirs(os.path.dirname(EVAL_OUT), exist_ok=True)
    with open(EVAL_OUT, "w") as f:
        json.dump(out, f, indent=1)

    r = out["recent"]
    print("RECENT (2-fold CV by winter)")
    for k in tracks:
        print(f"  {k:14s} onshore CF {r['onshore_cf'][k]['mean']:.3f} (calm ratio {r['onshore_cf'][k]['ratio_when_real_below_0p05']:.2f})"
              f"  offshore {r['offshore_cf'][k]['mean']:.3f}  events 48h<0.06 {r['events_per_winter']['48h_0.06'][k]:.2f}"
              f"  per-winter r {r['per_winter_r_vs_real'][k]:.2f}")
    print("  events/winter, 48 h:", {t: {k: round(v, 2) for k, v in r["events_per_winter"][f"48h_{t:.2f}"].items()} for t in THRESHOLDS})
    if "historic" in out:
        h = out["historic"]
        print(f"HISTORIC ({len(h['winters'])} Oct-Feb winters, wind-only)")
        for k in htracks:
            print(f"  {k:14s} onshore CF {h['onshore_cf'][k]['mean']:.3f} (calm ratio {h['onshore_cf'][k]['ratio_when_real_below_0p05']:.2f})"
                  f"  offshore {h['offshore_cf'][k]['mean']:.3f}  per-winter r {h['per_winter_r_vs_real_48h_0p10'][k]:.2f}")
        print("  wind-only events/winter, 48 h:", h["wind_only_events_per_winter"])
    print(f"Saved {TABLE_OUT}, {EVAL_OUT}")


if __name__ == "__main__":
    main()
