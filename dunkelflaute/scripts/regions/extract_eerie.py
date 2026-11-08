"""
Region series from EERIE IFS-FESOM2-SR (IFS TCo1279 ~9 km + FESOM2 NG5), 0.25 deg regridded output,
6-hourly, Sep-Mar: hist-1950 (1970-2014), eerie-control-1950 (1971-2050), highres-future-ssp245
(2015-2070), member r1i1p1f1.

The 0.25 deg grid is ERA5's own (1440 x 721, 90 N to 90 S, 0 to 359.75 E), so:
  wind   log law from 10u/10v with the ERA5 fsr monthly climatology of the same cell; "wind" is
         uncalibrated, "wind_qm" applies the Germany-fitted ERA5 QM (fitted on exactly this grid);
         German EEZ by Charnock (wind_off_DE, wind_off_DE_qm)
  solar  native ssrd, accumulated over the 6 h ending at the validity time (/21600 s)
  t2m    2 m temperature, region mean (for the "darker but warmer" question)
Land/sea = ERA5 lsm. Weights cos(lat).

Data: GRIB messages located with the byte index from eerie_index.py and decoded with eccodes. Row k
of a store is normally time step k, but the control's instantaneous store has a few missing steps,
so the time of every row is resolved by bisection on decoded validity times (cached in the index
file) and every decoded message is checked against it.

Usage: python3 extract_eerie.py {hist|ctl|ssp245}       (legacy python: eccodes)
Output: dunkelflaute/data/regions/eerie_<run>.npz
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "1")
import eccodes
import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
from core.capacity_factor import ONSHORE, OFFSHORE, solar_cf, wind_cf   # noqa: E402
import regions as rg                                                    # noqa: E402

IDX = f"{rg.BASE_DIR}/eerie_index"
VARS = ("10u", "10v", "ssrd", "2t")
DT = 21600.0
H6 = np.timedelta64(6, "h")
NLAT, NLON = 721, 1440

_S = None


def _read(ix, k):
    with open(ix["path"][k], "rb") as f:
        f.seek(int(ix["offset"][k]))
        b = f.read(int(ix["size"][k]))
    return eccodes.codes_new_from_message(b)


def _vt(h):
    d, t = eccodes.codes_get(h, "validityDate"), eccodes.codes_get(h, "validityTime")
    s = str(d)
    return np.datetime64(f"{s[:4]}-{s[4:6]}-{s[6:]}T{t // 100:02d}:{t % 100:02d}")


def row_time(ix, k):
    h = _read(ix, k)
    try:
        return _vt(h)
    finally:
        eccodes.codes_release(h)


def resolve_times(run, var):
    """Validity time of every row. Missing steps make offset(k) = (t_k - t_0)/6 h - k a nondecreasing
    step function; its change points are found by bisection, then spot-checked."""
    fn = f"{IDX}/{run}_{var}.npz"
    ix = dict(np.load(fn))
    if "time" in ix:
        return ix
    n = len(ix["path"])
    t0 = row_time(ix, 0)
    cache = {}

    def off(k):
        if k not in cache:
            cache[k] = int((row_time(ix, k) - t0) // H6) - k
        return cache[k]

    jumps = []

    def split(a, b):
        if off(a) == off(b):
            return
        if b - a == 1:
            jumps.append(b)
            return
        m = (a + b) // 2
        split(a, m)
        split(m, b)

    split(0, n - 1)
    offs = np.zeros(n, np.int64)
    for j in sorted(jumps):
        offs[j:] = off(j)
    times = t0 + (np.arange(n) + offs) * H6
    for k in list(range(0, n, 997)) + [n - 1]:
        assert row_time(ix, k) == times[k], f"{run} {var}: row {k} time model fails"
    ix["time"] = times.astype("datetime64[s]")
    np.savez(fn, **ix)
    print(f"  {run} {var}: {n} rows, {t0} .. {times[-1]}, missing steps at rows {sorted(jumps)}", flush=True)
    return ix


def setup():
    with nc.Dataset(rg.ERA5_LSM) as d:
        lsm = d["lsm"][:].squeeze().astype(float)
        lat = d["latitude"][:].astype(float)
        lon = d["longitude"][:].astype(float)
    flat = np.repeat(lat, NLON)
    flon = np.tile(lon, NLAT)
    band = np.where((flat >= rg.LAT_BAND[0] - 1) & (flat <= rg.LAT_BAND[1] + 1))[0]
    cells_b = rg.region_cells(flat[band], flon[band], lsm.ravel()[band] >= 0.5)
    cells = {r: band[c] for r, c in cells_b.items()}                       # flat global indices
    union = np.unique(np.concatenate(list(cells.values())))
    land_regions = [r for r in cells if r != "DE_off"]
    land_u = np.unique(np.concatenate([cells[r] for r in land_regions]))
    cache = f"{rg.OUT_DIR}/static_eerie_era5grid.npz"
    if os.path.exists(cache) and np.array_equal(np.load(cache)["land_u"], land_u):
        z0 = np.load(cache)["z0"]
    else:
        z0 = rg.era5_z0_on(flat[land_u], flon[land_u], 0.0, 0.0)
        os.makedirs(rg.OUT_DIR, exist_ok=True)
        np.savez(cache, land_u=land_u, z0=z0)
    return dict(flat=flat, cells=cells, union=union, land_regions=land_regions, land_u=land_u, z0=z0,
                w=np.cos(np.deg2rad(flat)), ul=np.searchsorted(union, land_u),
                uo=np.searchsorted(union, cells["DE_off"]),
                lpos=[np.searchsorted(land_u, cells[r]) for r in land_regions])


def _init():
    global _S
    _S = setup()
    _S["qm"] = rg.QM()


def _rmean(x, idx, w):
    return np.stack([(x[:, i] * w[i]).sum(axis=1) / w[i].sum() for i in idx], axis=1)


def _year(args):
    run, year, rows = args          # rows: {var: row indices}, all vars on the same time axis
    t0 = time.time()
    S = _S
    f = {}
    times = None
    for v in VARS:
        ix = np.load(f"{IDX}/{run}_{v}.npz")
        out = np.empty((len(rows[v]), len(S["union"])))
        tt = []
        for n, k in enumerate(rows[v]):
            h = _read(ix, k)
            try:
                assert eccodes.codes_get(h, "Ni") == NLON and eccodes.codes_get(h, "Nj") == NLAT
                tt.append(_vt(h))
                out[n] = eccodes.codes_get_values(h)[S["union"]]
            finally:
                eccodes.codes_release(h)
        tt = np.array(tt, dtype="datetime64[s]")
        assert np.array_equal(tt, ix["time"][rows[v]]), f"{run} {year} {v}: decoded times differ from the index"
        if times is None:
            times = tt
        assert np.array_equal(tt, times), f"{run} {year}: {v} not on the common time axis"
        f[v] = out
    mk = rg.month_of(times)
    slot = ((times - times.astype("datetime64[D]")) // np.timedelta64(3, "h")).astype(int) % 8
    ul, uo, lp = S["ul"], S["uo"], S["lpos"]
    w = S["w"][S["land_u"]]
    wo = S["w"][S["cells"]["DE_off"]]
    spd10 = np.hypot(f["10u"], f["10v"])
    z0m = S["z0"][mk - 1]
    ll = rg.loglaw_land(spd10[:, ul], z0m)
    qm = S["qm"](ll, slot, mk, S["qm"].land_class(z0m))
    res = dict(time=times)
    res["wind"] = _rmean(wind_cf(ll, **ONSHORE), lp, w)
    res["wind_qm"] = _rmean(wind_cf(qm, **ONSHORE), lp, w)
    res["spd100"] = _rmean(ll, lp, w)
    res["spd10"] = _rmean(spd10[:, ul], lp, w)
    res["solar"] = _rmean(solar_cf(f["ssrd"][:, ul], DT), lp, w)
    res["t2m"] = _rmean(f["2t"][:, ul], lp, w)
    so = rg.loglaw_sea(spd10[:, uo])
    so_qm = S["qm"](so, slot, mk, np.zeros(so.shape[1], int))
    res["wind_off_DE"] = (wind_cf(so, **OFFSHORE) * wo).sum(axis=1) / wo.sum()
    res["wind_off_DE_qm"] = (wind_cf(so_qm, **OFFSHORE) * wo).sum(axis=1) / wo.sum()
    return year, res, time.time() - t0


def main():
    run = sys.argv[1]
    S = setup()
    print(f"eerie {run}: " + ", ".join(f"{r} {len(c)}" for r, c in S["cells"].items()), flush=True)
    ix = {v: resolve_times(run, v) for v in VARS}
    common = ix[VARS[0]]["time"]
    for v in VARS[1:]:
        common = np.intersect1d(common, ix[v]["time"])
    mk = rg.month_of(common)
    common = common[np.isin(mk, rg.MONTHS)]
    yr = (common - np.timedelta64(1, "s")).astype("datetime64[Y]").astype(int) + 1970
    jobs = []
    for y in np.unique(yr):
        tsel = common[yr == y]
        jobs.append((run, int(y), {v: np.searchsorted(ix[v]["time"], tsel) for v in VARS}))
    n_all = sum(len(j[2][VARS[0]]) for j in jobs)
    print(f"  {len(jobs)} years, {n_all} Sep-Mar steps on the common axis", flush=True)
    t0 = time.time()
    out = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 64)), initializer=_init) as ex:
        for y, r, dt in ex.map(_year, jobs):
            out[y] = r
            print(f"  {y} {dt:.0f}s [{time.time()-t0:.0f}s]", flush=True)
    series = {k: np.concatenate([out[y][k] for y in sorted(out)]) for k in out[min(out)]}
    fn = rg.save(f"eerie_{run}", series, {r: S["cells"][r] for r in S["land_regions"]},
                 dict(model="IFS-FESOM2-SR (EERIE)", run=run, dt_h=6, grid="0.25 deg (ERA5 grid)",
                      wind="log law (ERA5 z0 same cell); wind_qm Germany QM", solar="native ssrd 6 h",
                      n_cells_DE_off=int(len(S["cells"]["DE_off"]))))
    print(f"Saved {fn} ({len(series['time'])} steps, {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
