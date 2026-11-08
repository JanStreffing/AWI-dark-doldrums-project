"""
Region series from TCo1279-DART (~9 km), 1950C and 2080C, Sep-Mar, 3-hourly.

Reads the global reconstructions the Germany pipeline already uses: log-law 100 m wind
(wind/data/dart_<run>_100m, ERA5 fsr land z0 year-cycled, Charnock over sea) and the v5 solar
reconstruction with exact monthly rescale (solar/data/dart_reconstructed_<run>). No calibration:
the QM table was fitted over Germany only, and at 9 km the Germany ratio is the same with and
without it (1.73 vs 1.76). Land/sea = ERA5 lsm regridded by nearest neighbour (as for Germany).
Means over cells are unweighted (DART cells are near-equal-area).

Usage: python3 extract_dart.py {1950c|2080c}
Output: dunkelflaute/data/regions/dart9km_<run>.npz with time (T,), solar/wind/spd100 (T, R) for
regions.REGIONS, and wind_off_DE/spd100_off_DE (T,) for the German EEZ.
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

RUNS = {"1950c": dict(years=range(1950, 1970), wind=f"{_ROOT}/wind/data/dart_1950c_100m",
                      solar=f"{_ROOT}/solar/data/dart_reconstructed_1950c"),
        "2080c": dict(years=range(2080, 2093), wind=f"{_ROOT}/wind/data/dart_2080c_100m",
                      solar=f"{_ROOT}/solar/data/dart_reconstructed_2080c")}
INDEX = f"{_ROOT}/wind/data/results/tco1279_era5_regrid_index.npz"
DT = 10800.0
CHUNK = 62

_C = None


def setup():
    d = np.load(INDEX)
    cells = rg.region_cells(d["tco_lat"], d["tco_lon"], d["land"])
    allc = np.unique(np.concatenate(list(cells.values())))
    lo, hi = int(allc.min()), int(allc.max()) + 1
    # local positions inside the [lo, hi) read window, and inside the union list
    union = allc - lo
    pos = {r: np.searchsorted(allc, c) for r, c in cells.items()}
    return dict(lo=lo, hi=hi, union=union, pos=pos, cells=cells)


def _init():
    global _C
    _C = setup()


def _month(args):
    run, year, month = args
    t0 = time.time()
    cfg, C = RUNS[run], _C
    fs = f"{cfg['solar']}/{year}/ssrd_reduced_3h_{year}{month:02d}-{year}{month:02d}.nc"
    fw = f"{cfg['wind']}/u100v100_{year}{month:02d}.nc"
    if not (os.path.exists(fs) and os.path.exists(fw)):
        return run, year, month, None, 0.0
    with nc.Dataset(fs) as s, nc.Dataset(fw) as w:
        tv = s["time_counter"]
        times = np.array([np.datetime64(str(x)[:19]) for x in
                          nc.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)])
        tw = w["time_counter"]
        times_w = np.array([np.datetime64(str(x)[:19]) for x in
                            nc.num2date(tw[:], tw.units, only_use_cftime_datetimes=False)])
        assert np.array_equal(times, times_w), f"{run} {year}-{month:02d}: wind and solar time axes differ"
        n = len(times)
        out = {k: [] for k in ("solar", "wind", "spd100", "wind_off_DE", "spd100_off_DE")}
        land_regions = [r for r in C["pos"] if r != "DE_off"]
        for a in range(0, n, CHUNK):
            b = min(n, a + CHUNK)
            ss = s["ssrd"][a:b, C["lo"]:C["hi"]][:, C["union"]].astype(np.float64)
            u = w["u100"][a:b, C["lo"]:C["hi"]][:, C["union"]].astype(np.float64)
            v = w["v100"][a:b, C["lo"]:C["hi"]][:, C["union"]].astype(np.float64)
            spd = np.hypot(u, v)
            cs, cw = solar_cf(ss, DT), wind_cf(spd, **ONSHORE)
            out["solar"].append(np.stack([cs[:, C["pos"][r]].mean(axis=1) for r in land_regions], axis=1))
            out["wind"].append(np.stack([cw[:, C["pos"][r]].mean(axis=1) for r in land_regions], axis=1))
            out["spd100"].append(np.stack([spd[:, C["pos"][r]].mean(axis=1) for r in land_regions], axis=1))
            so = spd[:, C["pos"]["DE_off"]]
            out["wind_off_DE"].append(wind_cf(so, **OFFSHORE).mean(axis=1))
            out["spd100_off_DE"].append(so.mean(axis=1))
    res = {k: np.concatenate(v) for k, v in out.items()}
    res["time"] = times
    return run, year, month, res, time.time() - t0


def main():
    run = sys.argv[1]
    C = setup()
    print(f"dart9km {run}: read window {C['lo']}-{C['hi']}, union {len(C['union'])} cells; "
          + ", ".join(f"{r} {len(c)}" for r, c in C["cells"].items()), flush=True)
    jobs = [(run, y, m) for y in RUNS[run]["years"] for m in rg.MONTHS]
    t0 = time.time()
    got = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 16)), initializer=_init) as ex:
        for r, y, m, res, dt in ex.map(_month, jobs):
            if res is None:
                print(f"  {y}-{m:02d} missing", flush=True)
                continue
            got[(y, m)] = res
            print(f"  {y}-{m:02d} {dt:.0f}s [{time.time()-t0:.0f}s]", flush=True)
    keys = sorted(got, key=lambda k: got[k]["time"][0])
    series = {k: np.concatenate([got[j][k] for j in keys]) for k in got[keys[0]]}
    fn = rg.save(f"dart9km_{run}", series, {r: c for r, c in C["cells"].items() if r != "DE_off"},
                 dict(model="TCo1279-DART", run=run, dt_h=3, wind="log law (ERA5 z0), uncalibrated",
                      solar="v5 + exact monthly rescale", n_cells_DE_off=int(len(C["cells"]["DE_off"]))))
    print(f"Saved {fn} ({len(series['time'])} steps, {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
