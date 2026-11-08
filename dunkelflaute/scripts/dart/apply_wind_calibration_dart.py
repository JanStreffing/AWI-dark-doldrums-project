"""
Does the log-law -> ERA5 100 m wind calibration change the 1950C -> 2080C Dunkelflaute signal?

Re-reads the DART log-law 100 m wind for Germany (same cells and power curves as
compute_germany_dunkelflaute_dart_mockert.py), applies the quantile-mapping table from
wind/scripts/methods/calibration/calibrate_quantile_mapping.py, and recomputes
onshore/offshore CF. Solar CF is taken unchanged from the stored results npz. Bins are
built exactly as in the calibration: month x 3-hourly UTC slot x surface class, where
the onshore roughness class uses the same year-cycled ERA5 fsr climatology the DART
log law itself used.

The raw (uncalibrated) CF is recomputed alongside and checked against the stored series,
so any difference in the output is the calibration and nothing else.

Output: dunkelflaute/results/wind_calibration/dart_calibrated_{run}.npz and
        dunkelflaute/results/wind_calibration/dart_calibration_effect.json (after both runs)

Usage: python apply_wind_calibration_dart.py {1950c|2080c|summary}
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import netCDF4 as nc
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compute_germany_dunkelflaute_dart_mockert as base  # noqa: E402  (chdirs to the repo root)
from core.capacity_factor import OFFSHORE, ONSHORE, wind_cf  # noqa: E402
from core.domain import REGRID_INDEX  # noqa: E402

TABLE = "wind/models/quantile_mapping_v1/qm_table.npz"
LAND_Z0 = "wind/data/results/tco1279_land_z0_2016_2025.npz"
OUT_DIR = "dunkelflaute/results/wind_calibration"
THRESHOLDS = (0.04, 0.05, 0.06, 0.07, 0.08, 0.10)
WINDOWS_H = (24, 48, 72)

_S = None


def _init():
    global _S
    m = base.germany_masks()
    t = np.load(TABLE)
    idx = np.load(REGRID_INDEX)
    land_rank = np.cumsum(idx["land"]) - 1
    lo, hi = m["lo"], m["hi"]
    on_cells = np.arange(lo, hi)[m["on"]]
    _S = dict(m=m, qx=t["qx"], qy=t["qy"], months=list(t["months"]), slots=int(t["slots"]),
              fsr_edges=t["fsr_edges"], n_class=int(t["n_class"]),
              z0=np.load(LAND_Z0)["land_z0"][:, :, land_rank[on_cells]])


def calibrate(v, times, month, cls):
    """v: (T, N) log-law speed; cls: (N,) surface class. Same mapping as calibrate_quantile_mapping.apply.
    month comes from the file, not the stamps: DART stamps the end of each 3 h interval, so a
    month file's last step is 00:00 on the 1st of the next month."""
    slot = (times - times.astype("datetime64[D]")).astype("timedelta64[s]").astype(int) // 10800
    m_idx = _S["months"].index(month)
    b = (m_idx * _S["slots"] + slot)[:, None] * _S["n_class"] + cls[None, :]
    out = v.copy()
    for k in np.unique(b):
        qx, qy = _S["qx"][k], _S["qy"][k]
        if np.isnan(qx[0]):
            continue
        sel = b == k
        x = v[sel]
        y = np.interp(x, qx, qy)
        y = np.where(x < qx[0], x + (qy[0] - qx[0]), y)
        y = np.where(x > qx[-1], x + (qy[-1] - qx[-1]), y)
        out[sel] = np.maximum(y, 0.0)
    return out


def _month(args):
    run, year, month = args
    t0 = time.time()
    cfg = base.RUNS[run]
    m = _S["m"]
    with nc.Dataset(f"{cfg['wind_dir']}/{base.WIND_FILE.format(year=year, month=month)}") as f:
        tv = f["time_counter"]
        times = nc.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)
        times = np.array([np.datetime64(f"{t.year:04d}-{t.month:02d}-{t.day:02d}T{t.hour:02d}:00:00") for t in times])
        u = f["u100"][:, m["lo"]:m["hi"]]
        v = f["v100"][:, m["lo"]:m["hi"]]
    spd = base.wind_speed(u, v)
    on, off = spd[:, m["on"]], spd[:, m["off"]]
    era5_year = 2016 + (year - base.RUNS[run]["years"][0]) % 10
    z0 = _S["z0"][era5_year - 2016, month - 1]
    cls_on = 1 + np.digitize(z0, _S["fsr_edges"])
    cls_off = np.zeros(off.shape[1], int)
    return year, month, dict(
        time=times,
        on_raw=wind_cf(on, **ONSHORE).mean(axis=1), off_raw=wind_cf(off, **OFFSHORE).mean(axis=1),
        on_cal=wind_cf(calibrate(on, times, month, cls_on), **ONSHORE).mean(axis=1),
        off_cal=wind_cf(calibrate(off, times, month, cls_off), **OFFSHORE).mean(axis=1)), time.time() - t0


def run_one(run):
    cfg = base.RUNS[run]
    winters = base.usable_winters(cfg["years"], cfg["spinup"])
    ym = sorted({mm for _, months in winters for mm in months})
    data = {}
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 32)), initializer=_init) as ex:
        for k, (y, mo, d, dt_s) in enumerate(ex.map(_month, [(run, y, mo) for y, mo in ym]), 1):
            data[(y, mo)] = d
            print(f"  [{k}/{len(ym)}] {y}-{mo:02d} {dt_s:.0f}s (elapsed {time.time()-t0:.0f}s)", flush=True)
    ref = np.load(f"{base.RESULTS}/{run}/combined_cf_{run}_sepmar_mockertweights.npz")
    cat = {k: np.concatenate([data[mm][k] for _, months in winters for mm in months])
           for k in ("time", "on_raw", "off_raw", "on_cal", "off_cal")}
    assert np.array_equal(cat["time"].astype("datetime64[s]"), ref["time"]), "time axis differs from stored run"
    dev = max(np.abs(cat["on_raw"] - ref["on"]).max(), np.abs(cat["off_raw"] - ref["off"]).max())
    print(f"raw CF vs stored: max abs diff {dev:.2e}", flush=True)
    assert dev < 1e-6
    os.makedirs(OUT_DIR, exist_ok=True)
    np.savez(f"{OUT_DIR}/dart_calibrated_{run}.npz", time=ref["time"], winter=ref["winter"], solar=ref["solar"],
             on_raw=cat["on_raw"], off_raw=cat["off_raw"], on_cal=cat["on_cal"], off_cal=cat["off_cal"])


def summary():
    w = base.WEIGHTS
    out = {}
    for run in ("1950c", "2080c"):
        d = np.load(f"{OUT_DIR}/dart_calibrated_{run}.npz")
        out[run] = {}
        for track in ("raw", "cal"):
            comb = w["solar"] * d["solar"] + w["onshore"] * d[f"on_{track}"] + w["offshore"] * d[f"off_{track}"]
            ev = {}
            for win_h in WINDOWS_H:
                for thr in THRESHOLDS:
                    n = [_n_events(comb[d["winter"] == y], win_h // 3, thr) for y in np.unique(d["winter"])]
                    ev[f"{win_h}h_{thr:.2f}"] = float(np.mean(n))
            out[run][track] = dict(mean_on=float(d[f"on_{track}"].mean()), mean_off=float(d[f"off_{track}"].mean()),
                                   mean_combined=float(comb.mean()), events_per_winter=ev)
    out["ratio_2080c_over_1950c"] = {
        track: {k: out["2080c"][track]["events_per_winter"][k] / out["1950c"][track]["events_per_winter"][k]
                for k in out["1950c"][track]["events_per_winter"]} for track in ("raw", "cal")}
    with open(f"{OUT_DIR}/dart_calibration_effect.json", "w") as f:
        json.dump(out, f, indent=1)
    for run in ("1950c", "2080c"):
        for track in ("raw", "cal"):
            o = out[run][track]
            print(f"{run} {track}: onshore {o['mean_on']:.3f} offshore {o['mean_off']:.3f} combined {o['mean_combined']:.3f} "
                  f"events 48h<0.06 {o['events_per_winter']['48h_0.06']:.2f}")
    for track in ("raw", "cal"):
        print(f"ratio {track} (48 h):", {t: round(out['ratio_2080c_over_1950c'][track][f'48h_{t:.2f}'], 2) for t in THRESHOLDS})


def _n_events(series, win, thr):
    """Mockert window-expansion events for any window length (base.mockert_events fixes it at 48 h)."""
    roll = np.full(len(series), np.nan)
    if len(series) >= win:
        c = np.cumsum(np.insert(series, 0, 0.0))
        roll[win - 1:] = (c[win:] - c[:-win]) / win
    flagged = np.zeros(len(series), bool)
    for i in np.where(np.nan_to_num(roll, nan=np.inf) < thr)[0]:
        flagged[max(0, i - win + 1):i + 1] = True
    return int(np.sum(flagged[1:] & ~flagged[:-1]) + flagged[0])


if __name__ == "__main__":
    if sys.argv[1] == "summary":
        summary()
    else:
        run_one(sys.argv[1].lower())
