"""
Region series from the TCo319-DART (~31 km) d runs, Sep-Mar, 6-hourly: ctl1950d (1950-1969) and
ssp585d (2080-2099).

Same recipe as dunkelflaute/scripts/tco319/compute_germany_tco319.py, on every region's land cells:
  wind   log law from 10u/10v with the ERA5 fsr land z0 climatology (ERA5 land cells within half a
         model cell), uncalibrated ("wind") and with the Germany-fitted QM ("wind_qm"); German EEZ
         by Charnock (wind_off_DE, wind_off_DE_qm)
  solar  v7 (6 h windows, four cloud covers at both window ends, mu), exact monthly rescale to the
         run's own monthly ssrd (mean of 6 h accumulations x windows in the month)
Land/sea: lsm of the DART-matched init file ICMGGa0iwINIT_DART, nearest reduced-Gaussian point to
each output cell (3-D chord distance). Weights cos(lat).

Usage: python3 extract_tco319.py {ctl1950d|ssp585d}
Output: dunkelflaute/data/regions/tco319_<run>.npz
"""
import calendar
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
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
from reconstruct_ssrd_v6 import DT6, FEATURES_V7, CLOUDS_V7, build_features_6h   # noqa: E402
from solar_geometry import toa_irradiance_accumulated_cells                      # noqa: E402
from rescale import constrained_rescale                                          # noqa: E402
from core.capacity_factor import ONSHORE, OFFSHORE, solar_cf, wind_cf            # noqa: E402
import regions as rg                                                             # noqa: E402

H = "/work/ab0995/ICCP_AWI_hackthon_2025"
RUNS = dict(ctl1950d=dict(dir="TCo319-DART-ctl1950d-gibbs-charn", years=range(1950, 1970)),
            ssp585d=dict(dir="TCo319-DART-ssp585d-gibbs-charn", years=range(2080, 2100)))
INIT = "/work/ab0246/a270092/input/oifs-43r3/TCO319L137/ICMGGa0iwINIT_DART"
GRIB_GET_DATA = "/work/ab0246/HPC_libraries/intel-oneapi-compilers/2022.0.1-gcc-11.2.0/openmpi/4.1.2-gcc-11.2.0/bin/grib_get_data"
V7_MODEL = f"{_ROOT}/solar/models/model_v7/model_kt_v7.joblib"
CACHE = f"{rg.OUT_DIR}/static_tco319_band.npz"
KMAX = 1.1

_S = None


def path(run, stream, var, year):
    return f"{H}/{RUNS[run]['dir']}/outdata/oifs/atm_remapped_{stream}_{var}_{stream}_{year}-{year}.nc"


def _xyz(lat, lon):
    la, lo = np.deg2rad(lat), np.deg2rad(lon)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def setup(run):
    from scipy.spatial import cKDTree
    with nc.Dataset(path(run, "6h", "10u", RUNS[run]["years"][0])) as d:
        lat, lon = d["lat"][:].astype(float), d["lon"][:].astype(float)
    rows = np.where((lat >= rg.LAT_BAND[0]) & (lat <= rg.LAT_BAND[1]))[0]
    r0, r1 = int(rows[0]), int(rows[-1]) + 1
    blat = lat[r0:r1]
    flat = np.repeat(blat, len(lon))
    flon = np.tile(lon, len(blat))
    if os.path.exists(CACHE):
        c = np.load(CACHE)
        lsm, z0 = c["lsm"], None
    else:
        txt = subprocess.run([GRIB_GET_DATA, "-w", "shortName=lsm", INIT], capture_output=True, text=True, check=True).stdout
        a = np.loadtxt(txt.splitlines()[1:])
        _, k = cKDTree(_xyz(a[:, 0], a[:, 1])).query(_xyz(flat, flon))
        lsm = a[k, 2]
    land = lsm >= 0.5
    cells = rg.region_cells(flat, flon, land)
    union = np.unique(np.concatenate(list(cells.values())))
    pos = {r: np.searchsorted(union, c) for r, c in cells.items()}
    land_regions = [r for r in cells if r != "DE_off"]
    land_u = np.unique(np.concatenate([cells[r] for r in land_regions]))
    if os.path.exists(CACHE) and "z0" in np.load(CACHE).files and np.array_equal(np.load(CACHE)["land_u"], land_u):
        z0 = np.load(CACHE)["z0"]
    else:
        z0 = rg.era5_z0_on(flat[land_u], flon[land_u], abs(lat[1] - lat[0]), abs(lon[1] - lon[0]))
        os.makedirs(rg.OUT_DIR, exist_ok=True)
        np.savez(CACHE, lsm=lsm, land_u=land_u, z0=z0)
    lpos = {r: np.searchsorted(land_u, cells[r]) for r in land_regions}
    return dict(r0=r0, r1=r1, nlon=len(lon), flat=flat, flon=flon, cells=cells, union=union, pos=pos,
                land_regions=land_regions, land_u=land_u, lpos=lpos, z0=z0,
                w=np.cos(np.deg2rad(flat)))


def _init(run):
    global _S
    _S = setup(run)
    _S["v7"] = joblib.load(V7_MODEL)["model"]
    _S["qm"] = rg.QM()


def _rmean(x, idx, w):
    return np.stack([(x[:, i] * w[i]).sum(axis=1) / w[i].sum() for i in idx], axis=1)


def _year(args):
    run, year = args
    t0 = time.time()
    S = _S
    lu = S["land_u"]
    f = {}
    for v in ("10u", "10v") + CLOUDS_V7:
        with nc.Dataset(path(run, "6h", v, year)) as d:
            if v == "10u":
                tv = d["time_instant"]
                times = np.array([np.datetime64(str(x)[:19]) for x in
                                  nc.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)])
            x = d[v][:, S["r0"]:S["r1"], :].astype(np.float32).reshape(len(times), -1)
            f[v] = x[:, S["union"]].astype(np.float64)
    with nc.Dataset(path(run, "1m", "ssrd", year)) as d:
        ssrd_m = d["ssrd"][:, S["r0"]:S["r1"], :].reshape(12, -1)[:, lu].astype(np.float64)   # mean per 6 h
    month = rg.month_of(times)
    pos_t = {t: i for i, t in enumerate(times)}
    prev = np.array([pos_t.get(t - np.timedelta64(6, "h"), -1) for t in times])
    keep = (prev >= 0) & np.isin(month, rg.MONTHS)
    tk, mk = times[keep], month[keep]
    slot = ((tk - tk.astype("datetime64[D]")) // np.timedelta64(3, "h")).astype(int) % 8
    ul = np.searchsorted(S["union"], lu)              # land cells inside the union list
    uo = S["pos"]["DE_off"]
    res = dict(time=tk)

    spd10 = np.hypot(f["10u"], f["10v"])[keep]
    sl = spd10[:, ul]
    z0m = S["z0"][mk - 1]
    ll = rg.loglaw_land(sl, z0m)
    qm = S["qm"](ll, slot, mk, S["qm"].land_class(z0m))
    w = S["w"][lu]
    lp = [S["lpos"][r] for r in S["land_regions"]]
    res["wind"] = _rmean(wind_cf(ll, **ONSHORE), lp, w)
    res["wind_qm"] = _rmean(wind_cf(qm, **ONSHORE), lp, w)
    res["spd100"] = _rmean(ll, lp, w)
    res["spd10"] = _rmean(sl, lp, w)
    so = rg.loglaw_sea(spd10[:, uo])
    wo = S["w"][S["cells"]["DE_off"]]
    so_qm = S["qm"](so, slot, mk, np.zeros(so.shape[1], int))
    res["wind_off_DE"] = (wind_cf(so, **OFFSHORE) * wo).sum(axis=1) / wo.sum()
    res["wind_off_DE_qm"] = (wind_cf(so_qm, **OFFSHORE) * wo).sum(axis=1) / wo.sum()

    w6 = {c: f[c][keep][:, ul] for c in CLOUDS_V7}
    w6.update({f"{c}_prev": f[c][prev[keep]][:, ul] for c in CLOUDS_V7})
    lat_l, lon_l = S["flat"][lu], S["flon"][lu]
    tisr = toa_irradiance_accumulated_cells(lat_l, lon_l, tk, dt_seconds=DT6, n_substeps=24).astype(np.float64)
    X, day, _ = build_features_6h(w6, tisr, FEATURES_V7)
    p = np.zeros(tisr.size)
    p[day] = np.clip(S["v7"].predict(X), 0, 1.5) * tisr.ravel()[day]
    p = p.reshape(tisr.shape)
    hours = (tk - tk.astype("datetime64[D]")) / np.timedelta64(1, "h") - 3.0
    lon_s = rg.lon180(lon_l)
    tw = np.maximum(np.cos(2 * np.pi * (((hours[:, None] + lon_s[None, :] / 15.0) % 24.0) - 12.0) / 24.0), 0.0) + 1e-9
    q = np.zeros_like(p)
    err = 0.0
    for m in rg.MONTHS:
        sel = mk == m
        if not sel.any():
            continue
        target = ssrd_m[m - 1] * 4 * calendar.monthrange(year, m)[1]
        q[sel], diag = constrained_rescale(p[sel], tisr[sel], target, kmax=KMAX, twilight_w=tw[sel])
        err = max(err, diag["max_rel_sum_err"])
    res["solar"] = _rmean(solar_cf(q, DT6), lp, w)
    return year, res, err, time.time() - t0


def main():
    run = sys.argv[1]
    S = setup(run)
    print(f"tco319 {run}: band rows {S['r0']}-{S['r1']}, " + ", ".join(f"{r} {len(c)}" for r, c in S["cells"].items()), flush=True)
    t0 = time.time()
    out = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 10)), initializer=_init, initargs=(run,)) as ex:
        for y, r, err, dt in ex.map(_year, [(run, y) for y in RUNS[run]["years"]]):
            out[y] = r
            print(f"  {y} {dt:.0f}s rescale_err={err:.1e} [{time.time()-t0:.0f}s]", flush=True)
    series = {k: np.concatenate([out[y][k] for y in sorted(out)]) for k in out[min(out)]}
    fn = rg.save(f"tco319_{run}", series, {r: S["cells"][r] for r in S["land_regions"]},
                 dict(model="TCo319-DART", run=run, dt_h=6, wind="log law (ERA5 z0); wind_qm with Germany QM",
                      solar="v7 + exact monthly rescale", n_cells_DE_off=int(len(S["cells"]["DE_off"]))))
    print(f"Saved {fn} ({len(series['time'])} steps, {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
