"""
Region series from the TCO95 (~100 km) reruns ctl1950g (constant 1950) and hi1950g (historical +
SSP5-8.5), Sep-Mar, 3-hourly, 1952-2099.

  wind         log law from 10u/10v with the ERA5 fsr land z0 climatology, uncalibrated, the recipe
               the 9 and 31 km extractors use; wind_qm with the Germany-fitted QM
  wind_native  native 100u/100v (the perfect-model truth)
  solar        native 3 h ssrd (accumulated over the 3 h ending at the stamp, /10800 s)
Land/sea from the model's own monthly lsm. Weights cos(lat).

Usage: python3 extract_tco95.py {ctl1950g|hi1950g}
Output: dunkelflaute/data/regions/tco95_<run>.npz
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "1")
import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
from core.capacity_factor import ONSHORE, OFFSHORE, solar_cf, wind_cf   # noqa: E402
import regions as rg                                                    # noqa: E402

RT = "/work/ab0995/a270092/runtime/awicm3-frontiers-xios"
YEARS = range(1952, 2100)
DT = 10800.0
FIELDS = ("10u", "10v", "100u", "100v", "ssrd")

_S = None


def path(run, stream, var, year):
    return f"{RT}/{run}/outdata/oifs/atm_remapped_{stream}_{var}_{stream}_{year}-{year}.nc"


def setup(run):
    with nc.Dataset(path(run, "1m", "lsm", 1952)) as d:
        lat, lon = d["lat"][:].astype(float), d["lon"][:].astype(float)
        lsm = d["lsm"][0].astype(float)
    rows = np.where((lat >= rg.LAT_BAND[0]) & (lat <= rg.LAT_BAND[1]))[0]
    r0, r1 = int(rows[0]), int(rows[-1]) + 1
    flat = np.repeat(lat[r0:r1], len(lon))
    flon = np.tile(lon, r1 - r0)
    land = lsm[r0:r1].ravel() >= 0.5
    cells = rg.region_cells(flat, flon, land)
    union = np.unique(np.concatenate(list(cells.values())))
    land_regions = [r for r in cells if r != "DE_off"]
    land_u = np.unique(np.concatenate([cells[r] for r in land_regions]))
    cache = f"{rg.OUT_DIR}/static_tco95_{run}_band.npz"
    if os.path.exists(cache) and np.array_equal(np.load(cache)["land_u"], land_u):
        z0 = np.load(cache)["z0"]
    else:
        z0 = rg.era5_z0_on(flat[land_u], flon[land_u], abs(lat[1] - lat[0]), abs(lon[1] - lon[0]))
        os.makedirs(rg.OUT_DIR, exist_ok=True)
        np.savez(cache, land_u=land_u, z0=z0)
    return dict(r0=r0, r1=r1, flat=flat, flon=flon, cells=cells, union=union, land_regions=land_regions,
                land_u=land_u, z0=z0, w=np.cos(np.deg2rad(flat)),
                ul=np.searchsorted(union, land_u), uo=np.searchsorted(union, cells["DE_off"]),
                lpos=[np.searchsorted(land_u, cells[r]) for r in land_regions])


def _init(run):
    global _S
    _S = setup(run)
    _S["qm"] = rg.QM()


def _rmean(x, idx, w):
    return np.stack([(x[:, i] * w[i]).sum(axis=1) / w[i].sum() for i in idx], axis=1)


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
            f[v] = d[v][:, S["r0"]:S["r1"], :].reshape(len(times), -1)[:, S["union"]].astype(np.float64)
    month = rg.month_of(times)
    keep = np.isin(month, rg.MONTHS)
    tk, mk = times[keep], month[keep]
    slot = ((tk - tk.astype("datetime64[D]")) // np.timedelta64(3, "h")).astype(int) % 8
    f = {k: v[keep] for k, v in f.items()}
    ul, uo, lp = S["ul"], S["uo"], S["lpos"]
    w = S["w"][S["land_u"]]
    wo = S["w"][S["cells"]["DE_off"]]
    spd10 = np.hypot(f["10u"], f["10v"])
    spd100 = np.hypot(f["100u"], f["100v"])
    z0m = S["z0"][mk - 1]
    ll = rg.loglaw_land(spd10[:, ul], z0m)
    qm = S["qm"](ll, slot, mk, S["qm"].land_class(z0m))
    res = dict(time=tk)
    res["wind"] = _rmean(wind_cf(ll, **ONSHORE), lp, w)
    res["wind_qm"] = _rmean(wind_cf(qm, **ONSHORE), lp, w)
    res["wind_native"] = _rmean(wind_cf(spd100[:, ul], **ONSHORE), lp, w)
    res["spd100"] = _rmean(ll, lp, w)
    res["spd100_native"] = _rmean(spd100[:, ul], lp, w)
    res["spd10"] = _rmean(spd10[:, ul], lp, w)
    res["solar"] = _rmean(solar_cf(f["ssrd"][:, ul], DT), lp, w)
    so = rg.loglaw_sea(spd10[:, uo])
    so_qm = S["qm"](so, slot, mk, np.zeros(so.shape[1], int))
    res["wind_off_DE"] = (wind_cf(so, **OFFSHORE) * wo).sum(axis=1) / wo.sum()
    res["wind_off_DE_qm"] = (wind_cf(so_qm, **OFFSHORE) * wo).sum(axis=1) / wo.sum()
    res["wind_off_DE_native"] = (wind_cf(spd100[:, uo], **OFFSHORE) * wo).sum(axis=1) / wo.sum()
    return year, res, time.time() - t0


def main():
    run = sys.argv[1]
    S = setup(run)
    print(f"tco95 {run}: band rows {S['r0']}-{S['r1']}, " + ", ".join(f"{r} {len(c)}" for r, c in S["cells"].items()), flush=True)
    t0 = time.time()
    out = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 24)), initializer=_init, initargs=(run,)) as ex:
        for y, r, dt in ex.map(_year, [(run, y) for y in YEARS]):
            out[y] = r
            print(f"  {y} {dt:.0f}s [{time.time()-t0:.0f}s]", flush=True)
    series = {k: np.concatenate([out[y][k] for y in sorted(out)]) for k in out[min(out)]}
    fn = rg.save(f"tco95_{run}", series, {r: S["cells"][r] for r in S["land_regions"]},
                 dict(model="TCO95", run=run, dt_h=3, wind="log law (ERA5 z0); wind_qm Germany QM; wind_native 100u/100v",
                      solar="native ssrd", n_cells_DE_off=int(len(S["cells"]["DE_off"]))))
    print(f"Saved {fn} ({len(series['time'])} steps, {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
