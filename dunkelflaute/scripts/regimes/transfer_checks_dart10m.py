"""
Raw 10 m wind speed region means for TCo1279-DART (~9 km), for the transfer artifact check A4.

The region npz of extract_dart.py holds only the log-law 100 m speed (spd100); this writes the
matching 10 m speed from the native 3-hourly 10u/10v (atm_reduced_3h_10{u,v}_3h_YYYYMM-YYYYMM.nc,
reduced Gaussian, the cell order of tco1279_era5_regrid_index.npz, checked), on exactly the region
cells and with the unweighted means of extract_dart.py. Sep-Mar months of the winters the
comparisons use.

Usage: python3 transfer_checks_dart10m.py {1950c|2080c}      (legacy python)
Output: dunkelflaute/data/regions/transfer_checks_dart9km_<run>_spd10.npz  (time (T,), spd10 (T, R))
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "1")
import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("regions", "regimes"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/{p}")
import regions as rg                     # noqa: E402
import extract_dart as ed                # noqa: E402
from transfer_checks_msl import winter_range   # noqa: E402

H = "/work/ab0995/ICCP_AWI_hackthon_2025"
DIRS = {"1950c": "TCo1279-DART-1950C", "2080c": "TCo1279-DART-2080C"}
_C = None


def _init():
    global _C
    _C = ed.setup()


def _month(args):
    run, year, month = args
    t0 = time.time()
    C = _C
    f = {v: f"{H}/{DIRS[run]}/outdata/oifs/atm_reduced_3h_{v}_3h_{year}{month:02d}-{year}{month:02d}.nc" for v in ("10u", "10v")}
    if not all(os.path.exists(x) for x in f.values()):
        return year, month, None, 0.0
    land_regions = [r for r in C["pos"] if r != "DE_off"]
    with nc.Dataset(f["10u"]) as du, nc.Dataset(f["10v"]) as dv:
        tv = du["time_counter"]
        times = np.array([np.datetime64(str(x)[:19]) for x in
                          nc.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)])
        if month == 1:
            la = du["lat"][C["lo"]:C["hi"]][C["union"]]
            assert np.allclose(la, np.load(ed.INDEX)["tco_lat"][C["lo"]:C["hi"]][C["union"]])
        out = []
        for a in range(0, len(times), ed.CHUNK):
            b = min(len(times), a + ed.CHUNK)
            u = du["10u"][a:b, C["lo"]:C["hi"]][:, C["union"]].astype(np.float64)
            v = dv["10v"][a:b, C["lo"]:C["hi"]][:, C["union"]].astype(np.float64)
            spd = np.hypot(u, v)
            out.append(np.stack([spd[:, C["pos"][r]].mean(axis=1) for r in land_regions], axis=1))
    return year, month, dict(time=times, spd10=np.concatenate(out)), time.time() - t0


def main():
    run = sys.argv[1]
    w0, w1 = winter_range("dart9km", run)
    jobs = [(run, y, m) for y in range(w0 - 1, w1 + 1) for m in rg.MONTHS
            if (m >= 9 and w0 <= y + 1 <= w1) or (m <= 3 and w0 <= y <= w1)]
    t0 = time.time()
    got = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 16)), initializer=_init) as ex:
        for y, m, res, dt in ex.map(_month, jobs):
            if res is None:
                print(f"  {y}-{m:02d} missing", flush=True)
                continue
            got[(y, m)] = res
            print(f"  {y}-{m:02d} {dt:.0f}s [{time.time() - t0:.0f}s]", flush=True)
    keys = sorted(got)
    fn = f"{rg.OUT_DIR}/transfer_checks_dart9km_{run}_spd10.npz"
    np.savez(fn, time=np.concatenate([got[k]["time"] for k in keys]),
             spd10=np.concatenate([got[k]["spd10"] for k in keys]), regions=np.array(rg.REGIONS))
    print(f"Saved {fn} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
