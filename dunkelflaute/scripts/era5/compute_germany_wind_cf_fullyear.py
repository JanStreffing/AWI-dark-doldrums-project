"""
Germany wind-only CF (onshore + offshore), ALL 12 calendar months (not just
Oct-Feb), 2015-01 through the latest complete month on disk. Wind-only sibling
of compute_germany_combined_cf_era5.py -- same log-law reconstruction, same
real-boundary Germany/EEZ masks, same wind_cf power curve -- but skips solar
entirely (no ML model, no clouds/tsr/ssrd needed), so it can cover the full
calendar year and run fast without depending on solar data coverage.

Two CF tracks per region, both from real ERA5 data for Germany:
  recon: log-law-reconstructed 100 m wind (v100 = v10*ln(100/z0)/ln(10/z0)),
         land z0 = ERA5's own real fsr for that exact month, ocean z0 by
         Charnock from that month's own u10/v10 -- same method as the
         DART/ERA5-validation pipelines.
  real:  ERA5's own genuine u100/v100.

Usage: python compute_germany_wind_cf_fullyear.py [--first-year 2015] [--last-year 2026]
"""
import datetime as dt
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import xarray as xr

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts")
from core.capacity_factor import wind_cf, ONSHORE, OFFSHORE  # noqa: E402

WIND_DIR = f"{_ROOT}/wind/data/era5/training_data"
LSM_PATH = f"{_ROOT}/wind/data/era5/era5_lsm_20260830_0000.nc"
LAND_PATH = f"{_ROOT}/boundaries/germany/germany_land.geojson"
EEZ_PATH = f"{_ROOT}/boundaries/germany/germany_eez.geojson"
OUT_PATH = f"{_ROOT}/dunkelflaute/data/germany_era5/wind_cf_era5_fullyear_2015_2026.npz"

KAPPA, G, NU, ALPHA_CH = 0.4, 9.81, 1.5e-5, 0.018


def _arg(name, default):
    return int(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


FIRST_YEAR = _arg("--first-year", 2015)
LAST_YEAR = _arg("--last-year", 2026)

GERMANY_BBOX = dict(latlo=46.75, lathi=56.25, lonlo=2.75, lonhi=15.5)


def box_slices(latlo, lathi, lonlo, lonhi):
    lat_sl = slice(int(round((90 - lathi) / 0.25)), int(round((90 - latlo) / 0.25)) + 1)
    lon_sl = slice(int(round(lonlo / 0.25)), int(round(lonhi / 0.25)) + 1)
    return lat_sl, lon_sl


LAT_SL, LON_SL = box_slices(**GERMANY_BBOX)


def germany_masks():
    ds = xr.open_dataset(LSM_PATH).isel(latitude=LAT_SL, longitude=LON_SL)
    lsm = ds.lsm.squeeze().values
    lat, lon = ds.latitude.values, ds.longitude.values
    ds.close()
    true_land = lsm >= 0.5
    ocean = ~true_land

    import json
    import shapely.geometry as sgeom
    import shapely.prepared

    with open(LAND_PATH) as f:
        germany_geom = sgeom.shape(json.load(f)["features"][0]["geometry"])
    land_prepared = shapely.prepared.prep(germany_geom)
    with open(EEZ_PATH) as f:
        eez_geom = sgeom.shape(json.load(f)["features"][0]["geometry"])
    eez_prepared = shapely.prepared.prep(eez_geom)

    lon_signed = np.where(lon > 180, lon - 360, lon)
    onshore = np.zeros(lsm.shape, dtype=bool)
    offshore = np.zeros(lsm.shape, dtype=bool)
    for i in range(len(lat)):
        for j in range(len(lon)):
            pt = sgeom.Point(lon_signed[j], lat[i])
            if true_land[i, j]:
                if land_prepared.contains(pt):
                    onshore[i, j] = True
            elif eez_prepared.contains(pt):
                offshore[i, j] = True
    return dict(true_land=true_land, onshore=onshore, offshore=offshore, lat=lat, lon=lon)


def complete(yr, m):
    import calendar
    wp = f"{WIND_DIR}/era5_winds_{yr}{m:02d}_3hourly.nc"
    zp = f"{WIND_DIR}/era5_z0_{yr}{m:02d}_3hourly.nc"
    if not (os.path.exists(wp) and os.path.exists(zp)):
        return False
    try:
        with xr.open_dataset(wp) as ds:
            return ds.sizes["valid_time"] == calendar.monthrange(yr, m)[1] * 8
    except Exception:
        return False


ALL_YM = [(y, m) for y in range(FIRST_YEAR, LAST_YEAR + 1) for m in range(1, 13)]
USABLE = [(y, m) for y, m in ALL_YM if complete(y, m)]


def reconstruct_wind(u10, v10, fsr, land3d):
    spd10 = np.hypot(u10, v10)
    spd_safe = np.maximum(spd10, 0.5)
    z0 = np.where(land3d, 1.0, 2.0e-4).astype(np.float32)
    for _ in range(8):
        ustar = KAPPA * spd_safe / np.log(10.0 / z0)
        z0_new = ALPHA_CH * ustar ** 2 / G + 0.11 * NU / ustar
        z0 = np.where(land3d, z0, z0_new)
    z0 = np.where(land3d, fsr, z0)
    ratio = np.log(100.0 / z0) / np.log(10.0 / z0)
    return (u10 * ratio).astype("float32"), (v10 * ratio).astype("float32")


_worker_masks = None


def _init_worker():
    global _worker_masks
    _worker_masks = germany_masks()


def _worker(ym):
    year, month = ym
    t0 = time.time()
    stamp = f"{year}{month:02d}"
    sub = dict(latitude=LAT_SL, longitude=LON_SL)
    winds = xr.open_dataset(f"{WIND_DIR}/era5_winds_{stamp}_3hourly.nc").isel(**sub)
    z0f = xr.open_dataset(f"{WIND_DIR}/era5_z0_{stamp}_3hourly.nc").isel(**sub)
    lat = winds.latitude.values
    times = winds.valid_time.values
    u10 = winds.u10.values.astype("float32")
    v10 = winds.v10.values.astype("float32")
    u100_real = winds.u100.values.astype("float32")
    v100_real = winds.v100.values.astype("float32")
    fsr = z0f.fsr.values.astype("float32")
    winds.close(); z0f.close()

    land3d = np.broadcast_to(_worker_masks["true_land"], u10.shape)
    u100_recon, v100_recon = reconstruct_wind(u10, v10, fsr, land3d)

    on, off = _worker_masks["onshore"], _worker_masks["offshore"]

    def spatial_mean(field2d_series, mask):
        w = np.cos(np.deg2rad(lat))[:, None] * mask
        wsum = w.sum()
        return (field2d_series * w[None, :, :]).sum(axis=(1, 2)) / wsum

    v_recon = np.hypot(u100_recon, v100_recon)
    v_real = np.hypot(u100_real, v100_real)
    cf_on_recon = spatial_mean(wind_cf(v_recon, **ONSHORE), on)
    cf_off_recon = spatial_mean(wind_cf(v_recon, **OFFSHORE), off)
    cf_on_real = spatial_mean(wind_cf(v_real, **ONSHORE), on)
    cf_off_real = spatial_mean(wind_cf(v_real, **OFFSHORE), off)

    times_dt = np.array([dt.datetime.utcfromtimestamp(t.astype("datetime64[s]").astype(int)) for t in times])
    return year, month, dict(time=times_dt, cf_on_recon=cf_on_recon, cf_off_recon=cf_off_recon,
                              cf_on_real=cf_on_real, cf_off_real=cf_off_real), time.time() - t0


def main():
    t0 = time.time()
    diag_masks = germany_masks()
    print(f"Germany wind-only CF, full year: {len(USABLE)}/{len(ALL_YM)} months usable "
          f"({FIRST_YEAR}-{LAST_YEAR}), onshore={diag_masks['onshore'].sum()} cells, "
          f"offshore={diag_masks['offshore'].sum()} cells", flush=True)

    results = {}
    with ProcessPoolExecutor(max_workers=8, initializer=_init_worker) as ex:
        for year, month, res, elapsed in ex.map(_worker, USABLE):
            results[(year, month)] = res
            print(f"  {year}-{month:02d} done ({elapsed:.1f}s)", flush=True)

    keys = sorted(results.keys())
    out = {k: np.concatenate([results[ym][k] for ym in keys]) for k in
           ["time", "cf_on_recon", "cf_off_recon", "cf_on_real", "cf_off_real"]}
    np.savez(OUT_PATH, **out)
    print(f"\nSaved {OUT_PATH}: {len(out['time'])} timesteps, {out['time'][0]} .. {out['time'][-1]}")
    print(f"Total elapsed: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
