"""
Perfect-model test, step 1: Germany capacity-factor series from TCO95 (AWI-CM3, ~100 km),
native vs reconstructed, for every year of ctl1950g (constant 1950) and hi1950g
(historical + SSP5-8.5), 1952-2099.

TCO95 writes what the 9 km DART slices lack (native 100 m wind, 3-hourly ssrd), so the same
reconstructions the DART pipeline uses can be scored against the model's own truth, in a
constant-1950 climate and in a 2080s SSP5-8.5 climate. Per 3-hourly step, Germany means of:
  native   wind CF from 100u/100v, solar CF from ssrd
  ll       log law from 10u/10v: land z0 = ERA5 fsr (2016-2025 monthly climatology, geometric
           mean of the ERA5 land cells inside each TCO95 cell), ocean z0 = Charnock from the
           model's own 10 m wind -- the DART recipe (wind/scripts/methods/dart_loglaw/)
  qm       ll, then the ERA5 quantile-mapping calibration (wind/models/quantile_mapping_v1/)
  llfsr    log law with the model's OWN monthly fsr over land, to separate the borrowed-z0 error
  v5       solar model v5 (final model, as on DART) from tsr/tisr, tcc/hcc/mcc/lcc, mu, then
           the exact monthly rescale to the model's own monthly ssrd total per cell
Grid: regular 400 x 192 (atm_remapped). Germany = land cells whose centre is inside the land
polygon, sea cells inside the EEZ polygon; cos(lat)-weighted means (regular grid).
Conventions checked against the output (2026-11-05): instantaneous fields are stamped at the
instant, accumulated fields over the 3 h ENDING at the stamp (model tisr vs analytic tisr for
that window: r = 1.0000); a stamp belongs to the month of (stamp - 1 s), as on DART.

Usage: python3 perfect_model_series.py RUN [Y0 Y1]      (default 1952 2099)
Output: dunkelflaute/data/tco95/<RUN>_germany_3h.npz
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "2")

import joblib
import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/solar/scripts/core")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts")
from solar_geometry import toa_irradiance_accumulated, SOLAR_CONSTANT   # noqa: E402
from rescale import constrained_rescale                                 # noqa: E402
from core.capacity_factor import ONSHORE, OFFSHORE, solar_cf, wind_cf   # noqa: E402

os.chdir(_ROOT)
RT = "/work/ab0995/a270092/runtime/awicm3-frontiers-xios"
OUT_DIR = "dunkelflaute/data/tco95"
LAND_PATH = "boundaries/germany/germany_land.geojson"
EEZ_PATH = "boundaries/germany/germany_eez.geojson"
ERA5_FSR = "wind/data/era5/era5_fsr_monthly_2015_2025.nc"
ERA5_LSM = "wind/data/era5/era5_lsm_20260830_0000.nc"
QM_TABLE = "wind/models/quantile_mapping_v1/qm_table.npz"
V5_MODEL = "solar/models/model_v5/model_kt_v5.joblib"
BOX = dict(lat_min=46.0, lat_max=56.5, lon_min=4.0, lon_max=16.5)
DT = 10800.0
DAY_MU = 0.02
KMAX = 1.1
V5_FEATURES = ["T", "tcc", "hcc", "mcc", "lcc", "mu"]
KAPPA, G, NU, ALPHA_CH = 0.4, 9.81, 1.5e-5, 0.018
MONTHS = (9, 10, 11, 12, 1, 2, 3)
FIELDS = ("10u", "10v", "100u", "100v", "ssrd", "tsr", "tcc", "lcc", "mcc", "hcc")

_S = None


def path(run, stream, var, year):
    return f"{RT}/{run}/outdata/oifs/atm_remapped_{stream}_{var}_{stream}_{year}-{year}.nc"


def setup(run):
    """Box slices, Germany masks, cos(lat) weights, ERA5 z0 on the onshore cells, QM table."""
    import shapely.geometry as sgeom
    import shapely.prepared
    with nc.Dataset(path(run, "1m", "lsm", 1952)) as d:
        lat, lon = d["lat"][:].astype(float), d["lon"][:].astype(float)
        lsm = d["lsm"][0].astype(float)
    lon_s = np.where(lon > 180, lon - 360, lon)
    ilat = np.where((lat >= BOX["lat_min"]) & (lat <= BOX["lat_max"]))[0]
    ilon = np.where((lon_s >= BOX["lon_min"]) & (lon_s <= BOX["lon_max"]))[0]
    assert np.all(np.diff(ilat) == 1) and np.all(np.diff(ilon) == 1), "box must be contiguous"
    lat_sl, lon_sl = slice(ilat[0], ilat[-1] + 1), slice(ilon[0], ilon[-1] + 1)
    blat, blon, blsm = lat[lat_sl], lon_s[lon_sl], lsm[lat_sl, lon_sl]
    land = shapely.prepared.prep(sgeom.shape(json.load(open(LAND_PATH))["features"][0]["geometry"]))
    eez = shapely.prepared.prep(sgeom.shape(json.load(open(EEZ_PATH))["features"][0]["geometry"]))
    on = np.zeros(blsm.shape, bool)
    off = np.zeros(blsm.shape, bool)
    for i, la in enumerate(blat):
        for j, lo in enumerate(blon):
            p = sgeom.Point(lo, la)
            on[i, j] = blsm[i, j] >= 0.5 and land.contains(p)
            off[i, j] = blsm[i, j] < 0.5 and eez.contains(p)
    w2d = np.cos(np.deg2rad(blat))[:, None] * np.ones((1, len(blon)))

    # ERA5 fsr monthly climatology 2016-2025, geometric mean over ERA5 land cells in each TCO95 cell
    dlat, dlon = abs(lat[1] - lat[0]), abs(lon[1] - lon[0])
    with nc.Dataset(ERA5_FSR) as d:
        elat, elon = d["latitude"][:], d["longitude"][:]
        ei = np.where((elat >= BOX["lat_min"] - 1) & (elat <= BOX["lat_max"] + 1))[0]
        ej = np.where((elon >= BOX["lon_min"] - 1) & (elon <= BOX["lon_max"] + 1))[0]   # box is east of 0 E
        esl, ejl = slice(ei[0], ei[-1] + 1), slice(ej[0], ej[-1] + 1)
        elat, elon = elat[esl], elon[ejl]
        tv = d["valid_time"]
        yrs = np.array([t.year for t in nc.num2date(tv[:], tv.units)])
        fsr = d["fsr"][:, esl, ejl][yrs >= 2016]
    with nc.Dataset(ERA5_LSM) as d:
        elsm = d["lsm"][:].squeeze()[esl, ejl]
    clim = np.log(np.maximum(fsr.reshape(-1, 12, *fsr.shape[1:]).mean(axis=0), 1e-5))   # (12, nlat_box, nlon_box)
    elon_s = elon
    z0_on = []
    for i, j in zip(*np.where(on)):
        mi = np.abs(elat - blat[i]) <= dlat / 2
        mj = np.abs(elon_s - blon[j]) <= dlon / 2
        sel = (elsm[np.ix_(mi, mj)] >= 0.5)
        z0_on.append(np.exp(clim[:, mi][:, :, mj][:, sel].mean(axis=1)))
    z0_on = np.array(z0_on).T                                                           # (12, n_on)
    t = np.load(QM_TABLE)
    return dict(lat_sl=lat_sl, lon_sl=lon_sl, blat=blat, blon=blon, on=on, off=off,
                w_on=w2d[on], w_off=w2d[off], z0_era5_on=z0_on, qx=t["qx"], qy=t["qy"],
                qm_months=list(t["months"]), slots=int(t["slots"]), fsr_edges=t["fsr_edges"],
                n_class=int(t["n_class"]), lon_full=lon)


def _init(run):
    global _S
    _S = setup(run)
    _S["v5"] = joblib.load(V5_MODEL)["model"]


def wmean(x, w):
    return (x * w).sum(axis=1) / w.sum()


def loglaw_speed(spd10, land_z0, n_land):
    """Speeds at 100 m for (T, n_land + n_sea) cells: land columns first with z0 = land_z0 (T or 1, n_land),
    sea columns by Charnock iterated from the 10 m speed, exactly as reconstruct_month.py."""
    spd_safe = np.maximum(spd10, 0.5)
    z0 = np.full(spd10.shape, 2.0e-4, np.float64)
    sea = slice(n_land, None)
    for _ in range(8):
        ustar = KAPPA * spd_safe[:, sea] / np.log(10.0 / z0[:, sea])
        z0[:, sea] = ALPHA_CH * ustar ** 2 / G + 0.11 * NU / ustar
    z0[:, :n_land] = land_z0
    return spd10 * np.log(100.0 / z0) / np.log(10.0 / z0)


def calibrate(v, slot, month, cls):
    """cls: (N,) one surface class per cell, or (T, N) when the class follows the month's z0."""
    m_idx = np.array([_S["qm_months"].index(m) for m in month])
    b = (m_idx * _S["slots"] + slot)[:, None] * _S["n_class"] + (cls[None, :] if cls.ndim == 1 else cls)
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


def _year(args):
    run, year = args
    t0 = time.time()
    S = _S
    f = {}
    for v in FIELDS:
        with nc.Dataset(path(run, "3h", v, year)) as d:
            if v == "10u":
                tv = d["time_instant"]
                times = np.array([np.datetime64(str(x)[:19]) for x in
                                  nc.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)])
            f[v] = d[v][:, S["lat_sl"], S["lon_sl"]].astype(np.float64)
    with nc.Dataset(path(run, "1m", "fsr", year)) as d:
        fsr_model = d["fsr"][:, S["lat_sl"], S["lon_sl"]].astype(np.float64)[:, S["on"]]   # (12, n_on)
    month = ((times - np.timedelta64(1, "s")).astype("datetime64[M]").astype(int) % 12) + 1
    keep = np.isin(month, MONTHS)
    times, month = times[keep], month[keep]
    f = {k: v[keep] for k, v in f.items()}
    slot = ((times - times.astype("datetime64[D]")) // np.timedelta64(3, "h")).astype(int) % 8

    on, off = S["on"], S["off"]
    n_on = int(on.sum())
    cells = lambda x: np.concatenate([x[:, on], x[:, off]], axis=1)   # noqa: E731
    spd10 = cells(np.hypot(f["10u"], f["10v"]))
    spd100 = cells(np.hypot(f["100u"], f["100v"]))
    res = dict(time=times)
    res["on_native"] = wmean(wind_cf(spd100[:, :n_on], **ONSHORE), S["w_on"])
    res["off_native"] = wmean(wind_cf(spd100[:, n_on:], **OFFSHORE), S["w_off"])
    ll = loglaw_speed(spd10, S["z0_era5_on"][month - 1], n_on)
    res["on_ll"] = wmean(wind_cf(ll[:, :n_on], **ONSHORE), S["w_on"])
    res["off_ll"] = wmean(wind_cf(ll[:, n_on:], **OFFSHORE), S["w_off"])
    cls = np.concatenate([1 + np.digitize(S["z0_era5_on"].mean(axis=0), S["fsr_edges"]),
                          np.zeros(int(off.sum()), int)])
    qm = calibrate(ll, slot, month, cls)
    res["on_qm"] = wmean(wind_cf(qm[:, :n_on], **ONSHORE), S["w_on"])
    res["off_qm"] = wmean(wind_cf(qm[:, n_on:], **OFFSHORE), S["w_off"])
    llf = loglaw_speed(spd10, np.maximum(fsr_model[month - 1], 1e-5), n_on)
    res["on_llfsr"] = wmean(wind_cf(llf[:, :n_on], **ONSHORE), S["w_on"])
    res["spd100_on_native"] = wmean(spd100[:, :n_on], S["w_on"])
    res["spd100_on_ll"] = wmean(ll[:, :n_on], S["w_on"])
    res["spd100_off_native"] = wmean(spd100[:, n_on:], S["w_off"])
    res["spd100_off_ll"] = wmean(ll[:, n_on:], S["w_off"])

    # solar: native, and v5 + exact monthly rescale to the model's own monthly total per cell
    ssrd = f["ssrd"][:, on]
    res["solar_native"] = wmean(solar_cf(ssrd, DT), S["w_on"])
    tisr = toa_irradiance_accumulated(S["blat"], S["blon"], times, dt_seconds=DT)[:, on].astype(np.float64)
    mu = tisr / (DT * SOLAR_CONSTANT)
    day = mu > DAY_MU
    feat = dict(T=np.clip(f["tsr"][:, on][day] / tisr[day], 0, 1.5), tcc=f["tcc"][:, on][day],
                hcc=f["hcc"][:, on][day], mcc=f["mcc"][:, on][day], lcc=f["lcc"][:, on][day], mu=mu[day])
    p = np.zeros_like(tisr)
    p[day] = np.clip(S["v5"].predict(np.column_stack([feat[c] for c in V5_FEATURES])), 0, 1.5) * tisr[day]
    q = np.zeros_like(p)
    hours = (times - times.astype("datetime64[D]")) / np.timedelta64(1, "h") - 1.5
    lon_on = np.broadcast_to(S["blon"][None, :], on.shape)[on]
    tw = np.maximum(np.cos(2 * np.pi * (((hours[:, None] + lon_on[None, :] / 15.0) % 24.0) - 12.0) / 24.0), 0.0) + 1e-9
    err = 0.0
    for m in MONTHS:
        sel = month == m
        if sel.any():
            q[sel], diag = constrained_rescale(p[sel], tisr[sel], ssrd[sel].sum(axis=0), kmax=KMAX, twilight_w=tw[sel])
            err = max(err, diag["max_rel_sum_err"])
    res["solar_v5"] = wmean(solar_cf(q, DT), S["w_on"])
    res["solar_v5_pre"] = wmean(solar_cf(p, DT), S["w_on"])
    return year, res, err, time.time() - t0


def main():
    run = sys.argv[1]
    y0, y1 = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (1952, 2099)
    S = setup(run)
    print(f"{run} {y0}-{y1}: box {len(S['blat'])}x{len(S['blon'])}, onshore {S['on'].sum()} cells, "
          f"offshore {S['off'].sum()} cells; ERA5 z0 onshore (Jan) median {np.median(S['z0_era5_on'][0]):.3f} m", flush=True)
    t0 = time.time()
    out = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 16)), initializer=_init, initargs=(run,)) as ex:
        for y, r, err, dt in ex.map(_year, [(run, y) for y in range(y0, y1 + 1)]):
            out[y] = r
            print(f"  {y} {dt:.0f}s rescale_err={err:.1e} [{time.time()-t0:.0f}s]", flush=True)
    keys = list(out[y0])
    series = {k: np.concatenate([out[y][k] for y in sorted(out)]) for k in keys}
    os.makedirs(OUT_DIR, exist_ok=True)
    fn = f"{OUT_DIR}/{run}_germany_3h.npz"
    np.savez(fn, **series, n_on=int(S["on"].sum()), n_off=int(S["off"].sum()),
             on_mask=S["on"], off_mask=S["off"], box_lat=S["blat"], box_lon=S["blon"])
    print(f"Saved {fn} ({len(series['time'])} steps, {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
