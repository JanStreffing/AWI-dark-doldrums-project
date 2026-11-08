"""
Per-cell ERA5 samples over Germany (real land polygon onshore, EEZ offshore) for
calibrating the log-law 100 m wind against ERA5's own native 100 m wind.

For every month on disk in wind/data/era5/training_data, keeps only the Germany
onshore/offshore cells and stores, per 3-hourly step and cell:
  v10        ERA5 10 m wind speed
  v100_ll    log-law 100 m wind speed, exactly as the production pipeline computes it
             (reconstruct_wind from dunkelflaute/scripts/era5/compute_germany_wind_cf_fullyear.py:
             land z0 = ERA5 fsr at that step, ocean z0 = Charnock)
  v100_real  ERA5's native 100 m wind speed (the calibration target)
  fsr        ERA5 forecast surface roughness (land conditioning variable)

Two output files, one per sample set:
  recent:   Sep-Mar 2015-2026 (Sep 2014 .. Mar 2026 winters that are on disk) -- calibration
  historic: everything on disk before 2015 (Oct-Feb winters 1979-2014)       -- independent test

Usage: python extract_germany_era5_samples.py {recent|historic}
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import xarray as xr

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/era5")
from compute_germany_wind_cf_fullyear import (LAT_SL, LON_SL, WIND_DIR, germany_masks,  # noqa: E402
                                              reconstruct_wind)

OUT_DIR = f"{_ROOT}/wind/data/results/calibration"
SEPMAR = {9, 10, 11, 12, 1, 2, 3}

_M = None


def _init():
    global _M
    _M = germany_masks()


def months_on_disk(which):
    out = []
    for f in sorted(os.listdir(WIND_DIR)):
        if not f.startswith("era5_winds_"):
            continue
        y, m = int(f[11:15]), int(f[15:17])
        if not os.path.exists(f"{WIND_DIR}/era5_z0_{y}{m:02d}_3hourly.nc"):
            continue
        if which == "recent" and y >= 2015 and m in SEPMAR:
            out.append((y, m))
        elif which == "historic" and y < 2015:
            out.append((y, m))
    return out


def _month(ym):
    y, m = ym
    t0 = time.time()
    sub = dict(latitude=LAT_SL, longitude=LON_SL)
    with xr.open_dataset(f"{WIND_DIR}/era5_winds_{y}{m:02d}_3hourly.nc").isel(**sub) as w, \
         xr.open_dataset(f"{WIND_DIR}/era5_z0_{y}{m:02d}_3hourly.nc").isel(**sub) as z:
        times = w.valid_time.values.astype("datetime64[s]")
        u10, v10 = w.u10.values.astype("float32"), w.v10.values.astype("float32")
        u100, v100 = w.u100.values.astype("float32"), w.v100.values.astype("float32")
        fsr = z.fsr.values.astype("float32")
    land3d = np.broadcast_to(_M["true_land"], u10.shape)
    u_ll, v_ll = reconstruct_wind(u10, v10, fsr, land3d)
    cells = _M["onshore"] | _M["offshore"]
    return y, m, dict(time=times,
                      v10=np.hypot(u10, v10)[:, cells],
                      v100_ll=np.hypot(u_ll, v_ll)[:, cells],
                      v100_real=np.hypot(u100, v100)[:, cells],
                      fsr=fsr[:, cells]), time.time() - t0


def main():
    which = sys.argv[1]
    yms = months_on_disk(which)
    print(f"{which}: {len(yms)} months, {yms[0]} .. {yms[-1]}", flush=True)
    m = germany_masks()
    cells = m["onshore"] | m["offshore"]
    lat2d, lon2d = np.meshgrid(m["lat"], m["lon"], indexing="ij")
    res = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 16)), initializer=_init) as ex:
        for k, (y, mo, d, dt_s) in enumerate(ex.map(_month, yms), 1):
            res[(y, mo)] = d
            print(f"  [{k}/{len(yms)}] {y}-{mo:02d} {dt_s:.0f}s", flush=True)
    keys = sorted(res)
    os.makedirs(OUT_DIR, exist_ok=True)
    out = f"{OUT_DIR}/germany_era5_wind_samples_{which}.npz"
    np.savez(out, **{k: np.concatenate([res[ym][k] for ym in keys]) for k in ("time", "v10", "v100_ll", "v100_real", "fsr")},
             lat=lat2d[cells], lon=lon2d[cells], onshore=m["onshore"][cells], offshore=m["offshore"][cells],
             land=m["true_land"][cells])
    print(f"Saved {out}: {sum(len(res[k]['time']) for k in keys)} steps x {cells.sum()} cells")


if __name__ == "__main__":
    main()
