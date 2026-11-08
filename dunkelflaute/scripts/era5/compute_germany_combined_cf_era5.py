"""
ERA5 validation run: apply the SAME reconstruction methods built for
TCo1279-DART (trained solar ML model + log-law wind extrapolation) to real
ERA5 data for Germany, then run the identical combined-CF / Dunkelflaute
event pipeline used on DART (dunkelflaute/scripts/compute_germany_combined_cf.py).

Because ERA5 also carries the REAL ssrd and REAL u100/v100 (unlike DART,
which has neither at 3-hourly resolution), this produces TWO parallel result
sets per winter:
  - "recon"  : solar predicted by models/model_v3/area/model_kt_v3.joblib
               (features T,tcc,hcc,mcc,lcc,mu) + exact constrained monthly
               rescale (target = sum of ERA5's own real 3-hourly ssrd);
               wind by the log-law (v100 = v10*ln(100/z0)/ln(10/z0)), land z0
               = ERA5's own real forecast_surface_roughness for that exact
               month (no decade-cycling needed here, unlike DART), ocean z0
               by the Charnock relation from that month's own u10/v10.
  - "real"   : ERA5's own genuine ssrd and u100/v100, straight through the
               same CF formulas and event thresholds.
"recon" vs "real" isolates how much of any gap to Mockert/Li's published
numbers comes from the reconstruction methods themselves vs. the combined-CF/
event-detection logic (which "real" also exercises, so a real-vs-published
gap there is NOT a reconstruction artifact).

Winters: label 2007..2025 (OND(Y-1)+JF(Y)), the range confirmed complete on
disk for both solar (ssrd/clouds since 2004) and wind (winds/z0, gap-free
from Oct 2006 on) as of 2026-09-23.
"""
import datetime as dt
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

# Must be set before numpy/sklearn are imported: each worker process otherwise
# spawns its own OpenMP thread pool (via numpy's BLAS backend and sklearn's
# HistGradientBoostingRegressor), and N_WORKERS processes x many threads each
# oversubscribes the node's thread limit -- this crashed the pool outright
# (OMP Error #34 / BrokenProcessPool) at 8 workers with default threading.
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
from core.capacity_factor import wind_cf, solar_cf, ONSHORE, OFFSHORE, SSRD_ACCUM_SECONDS_ERA5   # noqa: E402

SOLAR_DIR = f"{_ROOT}/solar/data/era5"
WIND_DIR = f"{_ROOT}/wind/data/era5/training_data"
LSM_PATH = f"{_ROOT}/wind/data/era5/era5_lsm_20260830_0000.nc"
MODEL_PATH = f"{_ROOT}/solar/models/model_v3/area/model_kt_v3.joblib"
LAND_PATH = f"{_ROOT}/boundaries/germany/germany_land.geojson"
EEZ_PATH = f"{_ROOT}/boundaries/germany/germany_eez.geojson"

DT_ERA5 = 3600.0          # ERA5 tsr/ssrd/tisr: 1h accumulation
DAY_MU = 0.02
FINAL_FEATURES = ["T", "tcc", "hcc", "mcc", "lcc", "mu"]
KMAX = 1.1

KAPPA, G, NU, ALPHA_CH = 0.4, 9.81, 1.5e-5, 0.018   # same constants as dart_loglaw

WEIGHTS = dict(solar=0.577, onshore=0.369, offshore=0.054)

WINTER_START = int(os.environ.get("GERMANY_WINTER_START", "2007"))
WINTER_END = int(os.environ.get("GERMANY_WINTER_END", "2025"))
WINTER_LABELS = range(WINTER_START, WINTER_END + 1)   # OND(Y-1)+JF(Y) per label
MONTHS_OND = [10, 11, 12]
MONTHS_JF = [1, 2]
OUT_DIR = os.environ.get("GERMANY_OUT_DIR", f"{_ROOT}/dunkelflaute/data/germany_era5")


def box_slices(latlo, lathi, lonlo, lonhi):
    """ERA5 0.25deg grid, lat 90->-90, lon 0->359.75. Germany doesn't cross
    0deg so this always returns a single lon slice here."""
    lat_sl = slice(int(round((90 - lathi) / 0.25)), int(round((90 - latlo) / 0.25)) + 1)
    lon_sl = slice(int(round(lonlo / 0.25)), int(round(lonhi / 0.25)) + 1)
    return lat_sl, lon_sl


GERMANY_BBOX = dict(latlo=46.75, lathi=56.25, lonlo=2.75, lonhi=15.5)
LAT_SL, LON_SL = box_slices(**GERMANY_BBOX)


def germany_masks(real_boundary=True):
    """real_boundary=True (default): onshore = ERA5 land cells whose centre falls inside
    Germany's actual Natural Earth admin_0_countries polygon (1:10m; extracted via cartopy
    shapereader and cached Germany-only at LAND_PATH), so neighbouring-country
    (NL/BE/FR/CH/AT/CZ/PL) land slivers in the download box are excluded. Offshore = ERA5
    ocean cells whose centre falls inside Germany's Exclusive Economic Zone (Marine
    Regions "eez" layer, which already includes the territorial sea; fetched via the
    Marine Regions WFS and cached Germany-only at EEZ_PATH) -- so only Germany's actual
    North Sea/Baltic waters count, not the whole download box's ocean.
    real_boundary=False: coarse fallback, onshore/offshore = ERA5 land mask intersected
    with the download BOUNDING BOX only (no country/EEZ polygon test) -- kept only for
    comparison against the old approximation."""
    ds = xr.open_dataset(LSM_PATH).isel(latitude=LAT_SL, longitude=LON_SL)
    lsm = ds.lsm.squeeze().values  # (nlat, nlon), fraction land 0-1
    lat, lon = ds.latitude.values, ds.longitude.values
    ds.close()
    true_land = lsm >= 0.5   # physical land/sea, used for z0 (roughness) regardless of country
    ocean = ~true_land

    if not real_boundary:
        return dict(true_land=true_land, onshore=true_land, offshore=ocean)

    import json

    import shapely.geometry as sgeom
    import shapely.prepared

    with open(LAND_PATH) as f:
        germany_geom = sgeom.shape(json.load(f)["features"][0]["geometry"])
    land_prepared = shapely.prepared.prep(germany_geom)

    with open(EEZ_PATH) as f:
        eez_geom = sgeom.shape(json.load(f)["features"][0]["geometry"])
    eez_prepared = shapely.prepared.prep(eez_geom)

    lon_signed = np.where(lon > 180, lon - 360, lon)   # Natural Earth / Marine Regions use -180..180
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
    tsr = xr.open_dataset(f"{SOLAR_DIR}/era5_tsr_3h_{stamp}.nc").isel(**sub)
    ssrd = xr.open_dataset(f"{SOLAR_DIR}/era5_ssrd_3h_{stamp}.nc").isel(**sub)
    winds = xr.open_dataset(f"{WIND_DIR}/era5_winds_{stamp}_3hourly.nc").isel(**sub)
    z0f = xr.open_dataset(f"{WIND_DIR}/era5_z0_{stamp}_3hourly.nc").isel(**sub)

    lat, lon = clouds.latitude.values, clouds.longitude.values
    times = clouds.valid_time.values
    tisr = toa_irradiance_accumulated(lat, lon, times, dt_seconds=DT_ERA5)

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
    mu = tisr_flat / (DT_ERA5 * SOLAR_CONSTANT)
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
    ntime = u10.shape[0]
    land3d = np.broadcast_to(land, u10.shape)
    spd10 = np.hypot(u10, v10)
    spd_safe = np.maximum(spd10, 0.5)

    z0 = np.where(land3d, 1.0, 2.0e-4).astype(np.float32)
    for _ in range(8):
        ustar = KAPPA * spd_safe / np.log(10.0 / z0)
        z0_new = ALPHA_CH * ustar ** 2 / G + 0.11 * NU / ustar
        z0 = np.where(land3d, z0, z0_new)
    z0 = np.where(land3d, fsr, z0)   # ERA5's own real land roughness for this exact month

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
        # field: (ntime, nlat, nlon); area-weight by cos(lat)
        w = np.cos(np.deg2rad(m["lat"]))[:, None] * mask
        wsum = w.sum()
        return (field2d_series * w[None, :, :]).sum(axis=(1, 2)) / wsum

    cf_solar_recon = spatial_mean(solar_cf(ssrd_recon, SSRD_ACCUM_SECONDS_ERA5), on)
    cf_solar_real = spatial_mean(solar_cf(m["ssrd_real"], SSRD_ACCUM_SECONDS_ERA5), on)

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
REAL_BOUNDARY = os.environ.get("GERMANY_REAL_BOUNDARY", "1") == "1"


def _init_worker():
    global _worker_model, _worker_masks
    _worker_model = joblib.load(MODEL_PATH)["model"]
    _worker_masks = germany_masks(real_boundary=REAL_BOUNDARY)


def _worker(ym):
    year, month = ym
    t0 = time.time()
    res = month_cf(year, month, _worker_model, _worker_masks)
    return year, month, res, time.time() - t0


def rolling_mean_16(x):
    if len(x) < 16:
        return np.full(len(x), np.nan)
    kernel = np.ones(16) / 16.0
    m = np.convolve(x, kernel, mode="valid")
    pad = len(x) - len(m)
    return np.concatenate([np.full(pad, np.nan), m])


def main():
    t0 = time.time()
    winters = usable_winters()
    all_ym = sorted({ym for _, months in winters for ym in months})
    diag_masks = germany_masks(real_boundary=REAL_BOUNDARY)
    print(f"Germany ERA5 combined-CF pipeline: {len(winters)} winters "
          f"({winters[0][0]}-{winters[-1][0]}), {len(all_ym)} month-files, "
          f"model={MODEL_PATH}, real_boundary={REAL_BOUNDARY} "
          f"(onshore={diag_masks['onshore'].sum()} cells, offshore={diag_masks['offshore'].sum()} cells, "
          f"true_land={diag_masks['true_land'].sum()} cells)", flush=True)

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

    winter_id = series["winter_id"]

    def per_winter_rolling(combined):
        out = []
        for label, _ in winters:
            mask = winter_id == label
            out.append(rolling_mean_16(combined[mask]))
        return np.concatenate(out)

    roll_recon = per_winter_rolling(cf_combined_recon)
    roll_real = per_winter_rolling(cf_combined_real)

    def event_stats(roll48, wind_blend, solar_cf_arr, winter_id, n_winters):
        mockert = roll48 < 0.06
        wind_blend_ok = (WEIGHTS["onshore"] * wind_blend[0] + WEIGHTS["offshore"] * wind_blend[1]) / \
                        (WEIGHTS["onshore"] + WEIGHTS["offshore"])
        li_instant = (wind_blend_ok < 0.20) & (solar_cf_arr < 0.20)

        def durations(flag, min_dur=None):
            ev, cur = [], 0
            for i in range(len(flag)):
                f = bool(flag[i]) if not (isinstance(flag[i], float) and np.isnan(flag[i])) else False
                same = (i == 0) or (winter_id[i] == winter_id[i - 1])
                if f and same:
                    cur += 1
                else:
                    if cur > 0:
                        ev.append(cur * 3 / 24.0)
                    cur = 1 if (f and not same) else 0
            if cur > 0:
                ev.append(cur * 3 / 24.0)
            ev = np.array(ev)
            return ev[ev >= min_dur] if min_dur else ev

        mock_d = durations(mockert)
        li_d = durations(li_instant, min_dur=1.0)
        return dict(
            mockert_events=len(mock_d), mockert_per_winter=len(mock_d) / n_winters,
            mockert_mean_dur=float(mock_d.mean()) if len(mock_d) else 0.0,
            li_events=len(li_d), li_per_winter=len(li_d) / n_winters,
            li_mean_dur=float(li_d.mean()) if len(li_d) else 0.0,
        )

    n_winters = len(winters)
    stats_recon = event_stats(roll_recon, (series["cf_on_recon"], series["cf_off_recon"]),
                               series["cf_solar_recon"], winter_id, n_winters)
    stats_real = event_stats(roll_real, (series["cf_on_real"], series["cf_off_real"]),
                              series["cf_solar_real"], winter_id, n_winters)

    suffix = "" if REAL_BOUNDARY else "_bboxonly"
    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = f"{OUT_DIR}/combined_cf_era5_{WINTER_START}_{WINTER_END}{suffix}.npz"
    np.savez(out_path, **series, cf_combined_recon=cf_combined_recon, cf_combined_real=cf_combined_real,
             roll_recon=roll_recon, roll_real=roll_real)

    print("\n" + "=" * 70)
    print(f"{'':20s} {'RECON (model+log-law)':>24s} {'REAL (ERA5 direct)':>22s}")
    print(f"{'Mean CF solar':20s} {series['cf_solar_recon'].mean():>24.3f} {series['cf_solar_real'].mean():>22.3f}")
    print(f"{'Mean CF onshore':20s} {series['cf_on_recon'].mean():>24.3f} {series['cf_on_real'].mean():>22.3f}")
    print(f"{'Mean CF offshore':20s} {series['cf_off_recon'].mean():>24.3f} {series['cf_off_real'].mean():>22.3f}")
    print(f"{'Mean CF combined':20s} {cf_combined_recon.mean():>24.3f} {cf_combined_real.mean():>22.3f}")
    print(f"{'Mockert events/winter':20s} {stats_recon['mockert_per_winter']:>24.2f} {stats_real['mockert_per_winter']:>22.2f}   (published ~4/yr)")
    print(f"{'Mockert mean dur (d)':20s} {stats_recon['mockert_mean_dur']:>24.2f} {stats_real['mockert_mean_dur']:>22.2f}")
    print(f"{'Li events/winter':20s} {stats_recon['li_per_winter']:>24.2f} {stats_real['li_per_winter']:>22.2f}   (published ~5-10/yr)")
    print(f"{'Li mean dur (d)':20s} {stats_recon['li_mean_dur']:>24.2f} {stats_real['li_mean_dur']:>22.2f}")
    print("=" * 70)
    print(f"Saved: {out_path}")
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
