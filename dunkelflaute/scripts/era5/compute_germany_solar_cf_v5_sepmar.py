"""Germany solar-only CF reconstruction, September-March, years 2015-2026, using
EITHER v4 (all-12-months-trained) OR v5 (Sep-Mar-only-trained) -- both applied to
the GLOBAL true-3-hour ssrd/tsr (solar/data/era5_3h_mean/, built by
aggregate_hourly_to_3h.py), so the two runs isolate exactly one variable: which
months the model was trained on. Successor to
compute_germany_combined_cf_era5_true3h.py, which only covered a single pilot
winter (2024-2025) using pre-sliced Germany-only true3h data, before the full
global 2015-2026 backfill existed. Solar only (no wind) -- this is the direct
model-reconstruction side of the SMARD comparison (see
dunkelflaute/smard_validation/scripts/compute_smard_solar_cf_fullyear.py for the real side).

Calendar-year month grouping, matching the rest of this project (not a season
spanning two calendar years): for year Y, months are Sep-Dec of Y AND Jan-Mar
of Y -- so "2015" means Sep-Dec 2015 + Jan-Mar 2015, not a single continuous
winter. Months with no source data yet (e.g. Sep-Dec 2026) are skipped, logged,
not treated as an error.

Usage: python3 compute_germany_solar_cf_v5_sepmar.py v4|v5
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import joblib
import numpy as np
import xarray as xr

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/solar/scripts/core")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts")
from solar_geometry import toa_irradiance_accumulated, SOLAR_CONSTANT   # noqa: E402
from rescale import constrained_rescale                                 # noqa: E402
from core.capacity_factor import solar_cf                               # noqa: E402

CLOUDS_DIR = f"{_ROOT}/solar/data/era5"
MEAN3H_DIR = f"{_ROOT}/solar/data/era5_3h_mean"

TRUE3H_SECONDS = 3 * 3600.0
DAY_MU = 0.02
FINAL_FEATURES = ["T", "tcc", "hcc", "mcc", "lcc", "mu"]
KMAX = 1.1
MONTHS = (9, 10, 11, 12, 1, 2, 3)
YEARS = range(2015, 2027)


def box_slices(latlo, lathi, lonlo, lonhi):
    lat_sl = slice(int(round((90 - lathi) / 0.25)), int(round((90 - latlo) / 0.25)) + 1)
    lon_sl = slice(int(round(lonlo / 0.25)), int(round(lonhi / 0.25)) + 1)
    return lat_sl, lon_sl


GERMANY_BBOX = dict(latlo=47.3, lathi=55.1, lonlo=5.9, lonhi=15.0)
LAT_SL, LON_SL = box_slices(**GERMANY_BBOX)
LSM_PATH = f"{_ROOT}/wind/data/era5/era5_lsm_20260830_0000.nc"
LAND_PATH = f"{_ROOT}/boundaries/germany/germany_land.geojson"


def germany_land_mask():
    """True land AND inside Germany's actual border polygon -- the bounding box
    above also covers parts of France, Switzerland, Austria, Czech Republic,
    Poland, Denmark, the Netherlands and Belgium, which this excludes. Same
    method as compute_germany_combined_cf_era5_true3h.py's germany_masks()."""
    import json
    import shapely.geometry as sgeom
    import shapely.prepared

    ds = xr.open_dataset(LSM_PATH).isel(latitude=LAT_SL, longitude=LON_SL)
    lsm = ds.lsm.squeeze().values
    lat, lon = ds.latitude.values, ds.longitude.values
    ds.close()
    true_land = lsm >= 0.5

    with open(LAND_PATH) as f:
        germany_geom = sgeom.shape(json.load(f)["features"][0]["geometry"])
    land_prepared = shapely.prepared.prep(germany_geom)

    lon_signed = np.where(lon > 180, lon - 360, lon)
    onshore = np.zeros(lsm.shape, dtype=bool)
    for i in range(len(lat)):
        for j in range(len(lon)):
            if true_land[i, j] and land_prepared.contains(sgeom.Point(lon_signed[j], lat[i])):
                onshore[i, j] = True
    return onshore


def month_available(year, month):
    stamp = f"{year}{month:02d}"
    return (os.path.exists(f"{CLOUDS_DIR}/era5_clouds_3h_{stamp}.nc") and
            os.path.exists(f"{MEAN3H_DIR}/era5_tsr_3hmean_{stamp}.nc") and
            os.path.exists(f"{MEAN3H_DIR}/era5_ssrd_3hmean_{stamp}.nc"))


def load_month(year, month):
    stamp = f"{year}{month:02d}"
    sub = dict(latitude=LAT_SL, longitude=LON_SL)
    clouds = xr.open_dataset(f"{CLOUDS_DIR}/era5_clouds_3h_{stamp}.nc").isel(**sub)
    tsr_n = xr.open_dataset(f"{MEAN3H_DIR}/era5_tsr_3hmean_{stamp}.nc").isel(**sub)
    ssrd_n = xr.open_dataset(f"{MEAN3H_DIR}/era5_ssrd_3hmean_{stamp}.nc").isel(**sub)

    ct = clouds.valid_time.values.astype("datetime64[s]")
    tn = tsr_n.valid_time.values.astype("datetime64[s]")
    sn = ssrd_n.valid_time.values.astype("datetime64[s]")
    common = np.intersect1d(np.intersect1d(ct, tn), sn)
    clouds = clouds.isel(valid_time=np.isin(ct, common))
    tsr_n = tsr_n.isel(time=np.isin(tn, common))
    ssrd_n = ssrd_n.isel(time=np.isin(sn, common))

    lat, lon = clouds.latitude.values, clouds.longitude.values
    tisr = toa_irradiance_accumulated(lat, lon, common, dt_seconds=TRUE3H_SECONDS)

    out = dict(
        lat=lat, lon=lon, times=common, tisr=tisr,
        tcc=clouds.tcc.values.astype("float32"), hcc=clouds.hcc.values.astype("float32"),
        mcc=clouds.mcc.values.astype("float32"), lcc=clouds.lcc.values.astype("float32"),
        tsr=tsr_n.tsr.values.astype("float32"), ssrd_real=ssrd_n.ssrd.values.astype("float64"),
    )
    for ds in (clouds, tsr_n, ssrd_n):
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


def month_cf(year, month, model, land_mask):
    m = load_month(year, month)
    ssrd_recon, diag = reconstruct_solar(m, model)
    w = np.cos(np.deg2rad(m["lat"]))[:, None] * np.ones((1, len(m["lon"]))) * land_mask
    wsum = w.sum()

    def spatial_mean(field2d_series):
        return (field2d_series * w[None, :, :]).sum(axis=(1, 2)) / wsum

    cf_recon = spatial_mean(solar_cf(ssrd_recon, TRUE3H_SECONDS))
    cf_real = spatial_mean(solar_cf(m["ssrd_real"], TRUE3H_SECONDS))
    return dict(time=m["times"], cf_recon=cf_recon, cf_real=cf_real,
                rescale_max_rel_err=diag["max_rel_sum_err"])


_worker_model = None
_worker_model_path = None
_worker_land_mask = None


def _init_worker():
    global _worker_model, _worker_land_mask
    _worker_model = joblib.load(_worker_model_path)["model"]
    _worker_land_mask = germany_land_mask()


def _worker(ym):
    year, month = ym
    t0 = time.time()
    res = month_cf(year, month, _worker_model, _worker_land_mask)
    return year, month, res, time.time() - t0


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in ("v4", "v5"):
        sys.exit("Usage: python3 compute_germany_solar_cf_v5_sepmar.py v4|v5")
    variant = sys.argv[1]
    model_path = f"{_ROOT}/solar/models/model_{variant}/model_kt_{variant}.joblib"
    out_path = f"{_ROOT}/dunkelflaute/smard_validation/data/germany_solar_cf_{variant}_sepmar_2015_2026.npz"

    global _worker_model_path
    _worker_model_path = model_path

    t0 = time.time()
    all_ym = [(y, m) for y in YEARS for m in MONTHS]
    avail_ym = [ym for ym in all_ym if month_available(*ym)]
    skipped = sorted(set(all_ym) - set(avail_ym))
    mask = germany_land_mask()
    print(f"Germany {variant} Sep-Mar solar-CF pipeline: {len(avail_ym)}/{len(all_ym)} month-files available, "
          f"model={model_path}", flush=True)
    print(f"  bounding box: {mask.size} cells; true German land (actual border polygon): "
          f"{mask.sum()} cells ({100*mask.mean():.1f}%)", flush=True)
    if skipped:
        print(f"  skipping (no source data yet): {skipped}", flush=True)

    results = {}
    with ProcessPoolExecutor(max_workers=8, initializer=_init_worker) as ex:
        for year, month, res, dt_s in ex.map(_worker, avail_ym):
            results[(year, month)] = res
            print(f"  {year}-{month:02d} done in {dt_s:.1f}s (rescale max_rel_err={res['rescale_max_rel_err']:.2e}) "
                  f"[elapsed {time.time()-t0:.0f}s]", flush=True)

    series = {k: [] for k in ["time", "cf_recon", "cf_real"]}
    for ym in sorted(avail_ym):
        r = results[ym]
        for k in series:
            series[k].append(r[k])
    for k in series:
        series[k] = np.concatenate(series[k])

    np.savez(out_path, **series)
    print("\n" + "=" * 50)
    print(f"{'':20s} {f'RECON ({variant})':>14s} {'REAL (ERA5 direct)':>20s}")
    print(f"{'Mean CF solar':20s} {series['cf_recon'].mean():>14.3f} {series['cf_real'].mean():>20.3f}")
    print("=" * 50)
    print(f"Saved: {out_path}")
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
