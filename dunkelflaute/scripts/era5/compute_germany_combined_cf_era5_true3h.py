"""
Same as compute_germany_combined_cf_era5.py, except ssrd and tsr come from a
genuine 3-hour accumulation (dunkelflaute/smard_validation/data/era5_true3h/, built by
dunkelflaute/smard_validation/scripts/build_era5_true3h.py from native-hourly ERA5) instead
of the main pipeline's era5_ssrd_3h_*.nc / era5_tsr_3h_*.nc, which only ever
sampled the 1-hour accumulation ending at each 3-hourly mark
(SSRD_ACCUM_SECONDS_ERA5 = 1*3600 in core/capacity_factor.py -- see the SMARD
validation conversation for how that surfaced). This is the test of whether
that 1-hour-vs-3-hour gap actually matters for the numbers, before committing
to redoing the full global 2015-2025 training set.

Only winter 2024-2025 is supported -- that's the only period the true-3h
ssrd/tsr exist for. Clouds and wind are UNCHANGED (already correctly
instantaneous, matching DART -- see download_era5.py's docstring), so they
still come from the original global solar/wind data directories.
"""
import datetime as dt
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import joblib
import numpy as np
import xarray as xr

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/solar/scripts/core")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts")
from solar_geometry import toa_irradiance_accumulated, SOLAR_CONSTANT   # noqa: E402
from rescale import constrained_rescale                                 # noqa: E402
from core.capacity_factor import wind_cf, solar_cf, ONSHORE, OFFSHORE   # noqa: E402

SOLAR_DIR = f"{_ROOT}/solar/data/era5"                          # clouds only, here
WIND_DIR = f"{_ROOT}/wind/data/era5/training_data"
TRUE3H_DIR = f"{_ROOT}/dunkelflaute/smard_validation/data/era5_true3h"       # ssrd, tsr, genuine 3h
LSM_PATH = f"{_ROOT}/wind/data/era5/era5_lsm_20260830_0000.nc"
MODEL_PATH = f"{_ROOT}/solar/models/model_v3/area/model_kt_v3.joblib"
LAND_PATH = f"{_ROOT}/boundaries/germany/germany_land.geojson"
EEZ_PATH = f"{_ROOT}/boundaries/germany/germany_eez.geojson"

TRUE3H_SECONDS = 3 * 3600.0   # genuine 3-hour accumulation, matching DART's SSRD_ACCUM_SECONDS_DART
DAY_MU = 0.02
FINAL_FEATURES = ["T", "tcc", "hcc", "mcc", "lcc", "mu"]
KMAX = 1.1

KAPPA, G, NU, ALPHA_CH = 0.4, 9.81, 1.5e-5, 0.018

WEIGHTS = dict(solar=0.577, onshore=0.369, offshore=0.054)

WINTER_START = int(os.environ.get("GERMANY_WINTER_START", "2025"))
WINTER_END = int(os.environ.get("GERMANY_WINTER_END", "2025"))
WINTER_LABELS = range(WINTER_START, WINTER_END + 1)
MONTHS_OND = [10, 11, 12]
MONTHS_JF = [1, 2]
OUT_DIR = os.environ.get("GERMANY_OUT_DIR", f"{_ROOT}/dunkelflaute/smard_validation/data")


def box_slices(latlo, lathi, lonlo, lonhi):
    lat_sl = slice(int(round((90 - lathi) / 0.25)), int(round((90 - latlo) / 0.25)) + 1)
    lon_sl = slice(int(round(lonlo / 0.25)), int(round(lonhi / 0.25)) + 1)
    return lat_sl, lon_sl


GERMANY_BBOX = dict(latlo=46.75, lathi=56.25, lonlo=2.75, lonhi=15.5)
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
    return dict(true_land=true_land, onshore=onshore, offshore=offshore)


def winter_months(label):
    ond_year = label - 1
    return [(ond_year, m) for m in MONTHS_OND] + [(label, m) for m in MONTHS_JF]


def usable_winters():
    return [(label, winter_months(label)) for label in WINTER_LABELS]


def load_month(year, month):
    stamp = f"{year}{month:02d}"
    sub = dict(latitude=LAT_SL, longitude=LON_SL)
    clouds = xr.open_dataset(f"{SOLAR_DIR}/era5_clouds_3h_{stamp}.nc").isel(**sub)
    tsr = xr.open_dataset(f"{TRUE3H_DIR}/era5_tsr_true3h_{stamp}.nc")       # already Germany-only
    ssrd = xr.open_dataset(f"{TRUE3H_DIR}/era5_ssrd_true3h_{stamp}.nc")    # already Germany-only
    winds = xr.open_dataset(f"{WIND_DIR}/era5_winds_{stamp}_3hourly.nc").isel(**sub)
    z0f = xr.open_dataset(f"{WIND_DIR}/era5_z0_{stamp}_3hourly.nc").isel(**sub)

    lat, lon = clouds.latitude.values, clouds.longitude.values
    times = clouds.valid_time.values
    tisr = toa_irradiance_accumulated(lat, lon, times, dt_seconds=TRUE3H_SECONDS)

    out = dict(
        lat=lat, lon=lon, times=times, tisr=tisr,
        tcc=clouds.tcc.values.astype("float32"), hcc=clouds.hcc.values.astype("float32"),
        mcc=clouds.mcc.values.astype("float32"), lcc=clouds.lcc.values.astype("float32"),
        tsr=tsr.tsr.values.astype("float32"), ssrd_real=ssrd.ssrd.values.astype("float64"),
        u10=winds.u10.values.astype("float32"), v10=winds.v10.values.astype("float32"),
        u100_real=winds.u100.values.astype("float32"), v100_real=winds.v100.values.astype("float32"),
        fsr=z0f.fsr.values.astype("float32"),
    )
    for ds in (clouds, tsr, ssrd, winds, z0f):
        ds.close()
    return out


def reconstruct_solar(m, model):
    ntime, nlat, nlon = m["tsr"].shape
    ncell = nlat * nlon
    tisr = m["tisr"].reshape(ntime, ncell).astype("float64")
    tisr_flat = m["tisr"].ravel()
    mu = tisr_flat / (TRUE3H_SECONDS * SOLAR_CONSTANT)
    day = mu > DAY_MU
    feat = {
        "tcc": m["tcc"].ravel()[day], "hcc": m["hcc"].ravel()[day],
        "mcc": m["mcc"].ravel()[day], "lcc": m["lcc"].ravel()[day],
        "mu": mu[day], "T": np.clip(m["tsr"].ravel()[day] / tisr_flat[day], 0, 1.5),
    }
    X = np.column_stack([feat[c] for c in FINAL_FEATURES])
    kt_pred = np.clip(model.predict(X), 0, 1.5)

    p = np.zeros(ntime * ncell, dtype="float64")
    p[day] = kt_pred * tisr_flat[day]
    p = p.reshape(ntime, ncell)

    real = m["ssrd_real"].reshape(ntime, ncell)
    target = real.sum(axis=0)

    hours = (m["times"] - m["times"].astype("datetime64[D]")) / np.timedelta64(1, "h") - 0.5
    local = (hours[:, None] + np.tile(m["lon"], nlat)[None, :] / 15.0) % 24.0
    tw = np.maximum(np.cos(2 * np.pi * (local - 12.0) / 24.0), 0.0) + 1e-9

    q, diag = constrained_rescale(p, tisr, target, kmax=KMAX, twilight_w=tw)
    return q.reshape(ntime, nlat, nlon), diag


def reconstruct_wind(m, land):
    u10, v10, fsr = m["u10"], m["v10"], m["fsr"]
    land3d = np.broadcast_to(land, u10.shape)
    spd10 = np.hypot(u10, v10)
    spd_safe = np.maximum(spd10, 0.5)

    z0 = np.where(land3d, 1.0, 2.0e-4).astype(np.float32)
    for _ in range(8):
        ustar = KAPPA * spd_safe / np.log(10.0 / z0)
        z0_new = ALPHA_CH * ustar ** 2 / G + 0.11 * NU / ustar
        z0 = np.where(land3d, z0, z0_new)
    z0 = np.where(land3d, fsr, z0)

    ratio = np.log(100.0 / z0) / np.log(10.0 / z0)
    u100 = (u10 * ratio).astype("float32")
    v100 = (v10 * ratio).astype("float32")
    return u100, v100


def month_cf(year, month, model, masks):
    m = load_month(year, month)
    ssrd_recon, rescale_diag = reconstruct_solar(m, model)
    u100_recon, v100_recon = reconstruct_wind(m, masks["true_land"])

    on = masks["onshore"]
    off = masks["offshore"]

    def spatial_mean(field2d_series, mask):
        w = np.cos(np.deg2rad(m["lat"]))[:, None] * mask
        wsum = w.sum()
        return (field2d_series * w[None, :, :]).sum(axis=(1, 2)) / wsum

    cf_solar_recon = spatial_mean(solar_cf(ssrd_recon, TRUE3H_SECONDS), on)
    cf_solar_real = spatial_mean(solar_cf(m["ssrd_real"], TRUE3H_SECONDS), on)

    v_recon = np.hypot(u100_recon, v100_recon)
    v_real = np.hypot(m["u100_real"], m["v100_real"])
    cf_on_recon = spatial_mean(wind_cf(v_recon, **ONSHORE), on)
    cf_off_recon = spatial_mean(wind_cf(v_recon, **OFFSHORE), off)
    cf_on_real = spatial_mean(wind_cf(v_real, **ONSHORE), on)
    cf_off_real = spatial_mean(wind_cf(v_real, **OFFSHORE), off)

    times = np.array([dt.datetime.utcfromtimestamp(t.astype("datetime64[s]").astype(int))
                       for t in m["times"]])

    return dict(
        time=times,
        cf_solar_recon=cf_solar_recon, cf_on_recon=cf_on_recon, cf_off_recon=cf_off_recon,
        cf_solar_real=cf_solar_real, cf_on_real=cf_on_real, cf_off_real=cf_off_real,
        rescale_max_rel_err=rescale_diag["max_rel_sum_err"],
    )


_worker_model = None
_worker_masks = None


def _init_worker():
    global _worker_model, _worker_masks
    _worker_model = joblib.load(MODEL_PATH)["model"]
    _worker_masks = germany_masks()


def _worker(ym):
    year, month = ym
    t0 = time.time()
    res = month_cf(year, month, _worker_model, _worker_masks)
    return year, month, res, time.time() - t0


def main():
    t0 = time.time()
    winters = usable_winters()
    all_ym = sorted({ym for _, months in winters for ym in months})
    diag_masks = germany_masks()
    print(f"Germany ERA5 TRUE-3H combined-CF pipeline: {len(winters)} winters "
          f"({winters[0][0]}-{winters[-1][0]}), {len(all_ym)} month-files, "
          f"model={MODEL_PATH} "
          f"(onshore={diag_masks['onshore'].sum()} cells, offshore={diag_masks['offshore'].sum()} cells)",
          flush=True)

    results_by_ym = {}
    with ProcessPoolExecutor(max_workers=4, initializer=_init_worker) as ex:
        for year, month, res, dt_s in ex.map(_worker, all_ym):
            results_by_ym[(year, month)] = res
            print(f"  {year}-{month:02d} done in {dt_s:.1f}s "
                  f"(rescale max_rel_err={res['rescale_max_rel_err']:.2e}) "
                  f"[elapsed {time.time()-t0:.0f}s]", flush=True)

    series = {k: [] for k in
              ["time", "winter_id",
               "cf_solar_recon", "cf_on_recon", "cf_off_recon",
               "cf_solar_real", "cf_on_real", "cf_off_real"]}

    for label, months in winters:
        for ym in months:
            r = results_by_ym[ym]
            n = len(r["time"])
            series["time"].append(r["time"])
            series["winter_id"].append(np.full(n, label))
            for k in ["cf_solar_recon", "cf_on_recon", "cf_off_recon",
                      "cf_solar_real", "cf_on_real", "cf_off_real"]:
                series[k].append(r[k])

    for k in series:
        series[k] = np.concatenate(series[k])

    def combine(pfx):
        return (WEIGHTS["solar"] * series[f"cf_solar_{pfx}"] +
                WEIGHTS["onshore"] * series[f"cf_on_{pfx}"] +
                WEIGHTS["offshore"] * series[f"cf_off_{pfx}"])

    cf_combined_recon = combine("recon")
    cf_combined_real = combine("real")

    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = f"{OUT_DIR}/combined_cf_era5_{WINTER_START}_{WINTER_END}_true3h.npz"
    np.savez(out_path, **series, cf_combined_recon=cf_combined_recon, cf_combined_real=cf_combined_real)

    print("\n" + "=" * 70)
    print(f"{'':20s} {'RECON (model+log-law)':>24s} {'REAL (ERA5 direct)':>22s}")
    print(f"{'Mean CF solar':20s} {series['cf_solar_recon'].mean():>24.3f} {series['cf_solar_real'].mean():>22.3f}")
    print(f"{'Mean CF onshore':20s} {series['cf_on_recon'].mean():>24.3f} {series['cf_on_real'].mean():>22.3f}")
    print(f"{'Mean CF offshore':20s} {series['cf_off_recon'].mean():>24.3f} {series['cf_off_real'].mean():>22.3f}")
    print(f"{'Mean CF combined':20s} {cf_combined_recon.mean():>24.3f} {cf_combined_real.mean():>22.3f}")
    print("=" * 70)
    print(f"Saved: {out_path}")
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
