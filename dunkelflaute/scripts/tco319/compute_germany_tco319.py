"""
Germany capacity-factor series from the TCo319-DART (~31 km) d runs: ctl1950d (constant 1950,
6-hourly 1950-1969) and ssp585d (SSP5-8.5 transient, 6-hourly 2080-2099).

Same recipes as the 9 km pipeline, adapted to what these runs write (docs/methods_notes.md 1.2):
  wind   log law from 10u/10v: land z0 = ERA5 fsr monthly climatology 2016-2025 (geometric mean of
         the ERA5 land cells inside each cell), ocean z0 = Charnock from the run's own 10 m wind;
         then the ERA5 quantile-mapping calibration (wind/models/quantile_mapping_v1/), binned by
         month x 3-hourly UTC slot x surface class, the class from that month's z0 as on DART
  solar  v7 (6 h windows, tcc/lcc/mcc/hcc at both window ends, mu), then the exact monthly
         rescale to the run's own monthly ssrd. That file holds the mean of 6 h accumulations,
         so the monthly total is mean x (6 h windows in the month); the one window not built
         (00-06 UTC on 1 January, which needs the previous year's file) is night over Germany.
Land/sea from the DART-matched init file ICMGGa0iwINIT_DART (TCO319L137), nearest neighbour onto
the 1312 x 640 output grid; its sr is a flat placeholder, hence the borrowed ERA5 z0.

Usage: python3 compute_germany_tco319.py {ctl1950d|ssp585d}
Output: dunkelflaute/data/tco319/<run>_germany_6h.npz
"""
import calendar
import json
import os
import subprocess
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
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/tco95")
from reconstruct_ssrd_v6 import DT6, FEATURES_V7, CLOUDS_V7, tisr_6h, build_features_6h   # noqa: E402
from rescale import constrained_rescale                                                  # noqa: E402
from core.capacity_factor import ONSHORE, OFFSHORE, solar_cf, wind_cf                    # noqa: E402
import perfect_model_series as pm                                                        # noqa: E402

os.chdir(_ROOT)
H = "/work/ab0995/ICCP_AWI_hackthon_2025"
RUNS = dict(ctl1950d=dict(dir="TCo319-DART-ctl1950d-gibbs-charn", years=range(1950, 1970)),
            ssp585d=dict(dir="TCo319-DART-ssp585d-gibbs-charn", years=range(2080, 2100)))
INIT = "/work/ab0246/a270092/input/oifs-43r3/TCO319L137/ICMGGa0iwINIT_DART"
GRIB_GET_DATA = "/work/ab0246/HPC_libraries/intel-oneapi-compilers/2022.0.1-gcc-11.2.0/openmpi/4.1.2-gcc-11.2.0/bin/grib_get_data"
V7_MODEL = "solar/models/model_v7/model_kt_v7.joblib"
OUT_DIR = "dunkelflaute/data/tco319"
LSM_CACHE = f"{OUT_DIR}/lsm_init_dart_on_output_grid_box.npz"
MONTHS = pm.MONTHS
WIND_FIELDS = ("10u", "10v")

_S = None


def path(run, stream, var, year):
    return f"{H}/{RUNS[run]['dir']}/outdata/oifs/atm_remapped_{stream}_{var}_{stream}_{year}-{year}.nc"


def init_lsm_on(blat, blon):
    """Land-sea mask of the DART init file, nearest reduced-Gaussian point to each box cell centre."""
    if os.path.exists(LSM_CACHE):
        c = np.load(LSM_CACHE)
        if np.array_equal(c["lat"], blat) and np.array_equal(c["lon"], blon):
            return c["lsm"]
    txt = subprocess.run([GRIB_GET_DATA, "-w", "shortName=lsm", INIT], capture_output=True, text=True, check=True).stdout
    a = np.loadtxt(txt.splitlines()[1:])
    glat, glon, gval = a[:, 0], np.where(a[:, 1] > 180, a[:, 1] - 360, a[:, 1]), a[:, 2]
    near = (glat > blat.min() - 1) & (glat < blat.max() + 1) & (glon > blon.min() - 1) & (glon < blon.max() + 1)
    glat, glon, gval = glat[near], glon[near], gval[near]
    out = np.zeros((len(blat), len(blon)))
    for i, la in enumerate(blat):
        dy = np.deg2rad(glat - la)
        for j, lo in enumerate(blon):
            dx = np.deg2rad(glon - lo) * np.cos(np.deg2rad(la))
            out[i, j] = gval[np.argmin(dx * dx + dy * dy)]
    os.makedirs(OUT_DIR, exist_ok=True)
    np.savez(LSM_CACHE, lat=blat, lon=blon, lsm=out)
    return out


def setup(run):
    import shapely.geometry as sgeom
    import shapely.prepared
    y0 = RUNS[run]["years"][0]
    with nc.Dataset(path(run, "6h", "10u", y0)) as d:
        lat, lon = d["lat"][:].astype(float), d["lon"][:].astype(float)
    lon_s = np.where(lon > 180, lon - 360, lon)
    B = pm.BOX
    ilat = np.where((lat >= B["lat_min"]) & (lat <= B["lat_max"]))[0]
    ilon = np.where((lon_s >= B["lon_min"]) & (lon_s <= B["lon_max"]))[0]
    assert np.all(np.diff(ilat) == 1) and np.all(np.diff(ilon) == 1), "box must be contiguous"
    lat_sl, lon_sl = slice(ilat[0], ilat[-1] + 1), slice(ilon[0], ilon[-1] + 1)
    blat, blon = lat[lat_sl], lon_s[lon_sl]
    lsm = init_lsm_on(blat, blon)
    land = shapely.prepared.prep(sgeom.shape(json.load(open(pm.LAND_PATH))["features"][0]["geometry"]))
    eez = shapely.prepared.prep(sgeom.shape(json.load(open(pm.EEZ_PATH))["features"][0]["geometry"]))
    on = np.zeros(lsm.shape, bool)
    off = np.zeros(lsm.shape, bool)
    for i, la in enumerate(blat):
        for j, lo in enumerate(blon):
            p = sgeom.Point(lo, la)
            on[i, j] = lsm[i, j] >= 0.5 and land.contains(p)
            off[i, j] = lsm[i, j] < 0.5 and eez.contains(p)
    w2d = np.cos(np.deg2rad(blat))[:, None] * np.ones((1, len(blon)))

    # ERA5 fsr climatology: geometric mean over the ERA5 land cells inside each cell, else the nearest ERA5 land cell
    with nc.Dataset(pm.ERA5_FSR) as d:
        elat, elon = d["latitude"][:], d["longitude"][:]
        ei = np.where((elat >= B["lat_min"] - 1) & (elat <= B["lat_max"] + 1))[0]
        ej = np.where((elon >= B["lon_min"] - 1) & (elon <= B["lon_max"] + 1))[0]
        esl, ejl = slice(ei[0], ei[-1] + 1), slice(ej[0], ej[-1] + 1)
        elat, elon = elat[esl], elon[ejl]
        tv = d["valid_time"]
        yrs = np.array([t.year for t in nc.num2date(tv[:], tv.units)])
        fsr = d["fsr"][:, esl, ejl][yrs >= 2016]
    with nc.Dataset(pm.ERA5_LSM) as d:
        elsm = d["lsm"][:].squeeze()[esl, ejl]
    clim = np.log(np.maximum(fsr.reshape(-1, 12, *fsr.shape[1:]).mean(axis=0), 1e-5))
    dlat, dlon = abs(lat[1] - lat[0]), abs(lon[1] - lon[0])
    eland = np.argwhere(elsm >= 0.5)
    z0_on = []
    for i, j in zip(*np.where(on)):
        mi = np.abs(elat - blat[i]) <= dlat / 2
        mj = np.abs(elon - blon[j]) <= dlon / 2
        sel = elsm[np.ix_(mi, mj)] >= 0.5
        if sel.any():
            z0_on.append(np.exp(clim[:, mi][:, :, mj][:, sel].mean(axis=1)))
        else:
            k = np.argmin((elat[eland[:, 0]] - blat[i]) ** 2 + (elon[eland[:, 1]] - blon[j]) ** 2)
            z0_on.append(np.exp(clim[:, eland[k, 0], eland[k, 1]]))
    z0_on = np.array(z0_on).T                                                                    # (12, n_on)
    t = np.load(pm.QM_TABLE)
    return dict(lat_sl=lat_sl, lon_sl=lon_sl, blat=blat, blon=blon, on=on, off=off,
                w_on=w2d[on], w_off=w2d[off], z0_era5_on=z0_on, qx=t["qx"], qy=t["qy"],
                qm_months=list(t["months"]), slots=int(t["slots"]), fsr_edges=t["fsr_edges"],
                n_class=int(t["n_class"]))


def _init(run):
    global _S
    _S = setup(run)
    _S["v7"] = joblib.load(V7_MODEL)["model"]
    pm._S = _S                       # pm.calibrate reads its table from pm._S


def _year(args):
    run, year = args
    t0 = time.time()
    S = _S
    f = {}
    for v in WIND_FIELDS + CLOUDS_V7:
        with nc.Dataset(path(run, "6h", v, year)) as d:
            if v == "10u":
                tv = d["time_instant"]
                times = np.array([np.datetime64(str(x)[:19]) for x in
                                  nc.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)])
            f[v] = d[v][:, S["lat_sl"], S["lon_sl"]].astype(np.float64)
    with nc.Dataset(path(run, "1m", "ssrd", year)) as d:
        ssrd_month = d["ssrd"][:, S["lat_sl"], S["lon_sl"]].astype(np.float64)[:, S["on"]]     # (12, n_on), mean per 6 h

    month = ((times - np.timedelta64(1, "s")).astype("datetime64[M]").astype(int) % 12) + 1
    on, off = S["on"], S["off"]
    n_on = int(on.sum())

    # 6 h windows need the clouds 6 h earlier; build on the full year, then keep Sep-Mar
    pos = {v: i for i, v in enumerate(times)}
    prev = np.array([pos.get(t - np.timedelta64(6, "h"), -1) for t in times])
    ok = prev >= 0
    keep = ok & np.isin(month, MONTHS)
    w6 = {c: f[c][keep][:, on] for c in CLOUDS_V7}
    w6.update({f"{c}_prev": f[c][prev[keep]][:, on] for c in CLOUDS_V7})
    times_k, month_k = times[keep], month[keep]
    slot = ((times_k - times_k.astype("datetime64[D]")) // np.timedelta64(3, "h")).astype(int) % 8

    res = dict(time=times_k)
    spd = np.hypot(f["10u"], f["10v"])[keep]
    spd10 = np.concatenate([spd[:, on], spd[:, off]], axis=1)
    z0_month = S["z0_era5_on"][month_k - 1]                                                       # (T, n_on)
    ll = pm.loglaw_speed(spd10, z0_month, n_on)
    res["on_ll"] = pm.wmean(wind_cf(ll[:, :n_on], **ONSHORE), S["w_on"])
    res["off_ll"] = pm.wmean(wind_cf(ll[:, n_on:], **OFFSHORE), S["w_off"])
    cls = np.concatenate([1 + np.digitize(z0_month, S["fsr_edges"]), np.zeros((len(times_k), int(off.sum())), int)], axis=1)
    qm = pm.calibrate(ll, slot, month_k, cls)
    res["on_qm"] = pm.wmean(wind_cf(qm[:, :n_on], **ONSHORE), S["w_on"])
    res["off_qm"] = pm.wmean(wind_cf(qm[:, n_on:], **OFFSHORE), S["w_off"])
    res["spd10_on"] = pm.wmean(spd10[:, :n_on], S["w_on"])
    res["spd100_on_qm"] = pm.wmean(qm[:, :n_on], S["w_on"])
    res["spd100_off_qm"] = pm.wmean(qm[:, n_on:], S["w_off"])

    # solar: v7 on the onshore cells, exact monthly rescale to the run's own monthly total
    lat_on = np.broadcast_to(S["blat"][:, None], on.shape)[on]
    lon_on = np.broadcast_to(S["blon"][None, :], on.shape)[on]
    tisr = tisr_6h(S["blat"], S["blon"], times_k)[:, on].astype(np.float64)
    X, day, _ = build_features_6h(w6, tisr, FEATURES_V7)
    p = np.zeros(tisr.size)
    p[day] = np.clip(S["v7"].predict(X), 0, 1.5) * tisr.ravel()[day]
    p = p.reshape(tisr.shape)
    hours = (times_k - times_k.astype("datetime64[D]")) / np.timedelta64(1, "h") - 3.0
    tw = np.maximum(np.cos(2 * np.pi * (((hours[:, None] + lon_on[None, :] / 15.0) % 24.0) - 12.0) / 24.0), 0.0) + 1e-9
    q = np.zeros_like(p)
    err = 0.0
    for m in MONTHS:
        sel = month_k == m
        if not sel.any():
            continue
        n_full = 4 * calendar.monthrange(year, m)[1]
        target = ssrd_month[m - 1] * n_full
        q[sel], diag = constrained_rescale(p[sel], tisr[sel], target, kmax=pm.KMAX, twilight_w=tw[sel])
        err = max(err, diag["max_rel_sum_err"])
    res["solar_v7"] = pm.wmean(solar_cf(q, DT6), S["w_on"])
    res["solar_v7_pre"] = pm.wmean(solar_cf(p, DT6), S["w_on"])
    return year, res, err, time.time() - t0


def main():
    run = sys.argv[1]
    S = setup(run)
    print(f"{run}: box {len(S['blat'])}x{len(S['blon'])}, onshore {S['on'].sum()} cells, offshore {S['off'].sum()} cells; "
          f"ERA5 z0 onshore (Jan) median {np.median(S['z0_era5_on'][0]):.3f} m", flush=True)
    t0 = time.time()
    out = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 20)), initializer=_init, initargs=(run,)) as ex:
        for y, r, err, dt in ex.map(_year, [(run, y) for y in RUNS[run]["years"]]):
            out[y] = r
            print(f"  {y} {dt:.0f}s rescale_err={err:.1e} [{time.time()-t0:.0f}s]", flush=True)
    series = {k: np.concatenate([out[y][k] for y in sorted(out)]) for k in out[min(out)]}
    os.makedirs(OUT_DIR, exist_ok=True)
    fn = f"{OUT_DIR}/{run}_germany_6h.npz"
    np.savez(fn, **series, on_mask=S["on"], off_mask=S["off"], box_lat=S["blat"], box_lon=S["blon"])
    print(f"Saved {fn} ({len(series['time'])} steps, {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
