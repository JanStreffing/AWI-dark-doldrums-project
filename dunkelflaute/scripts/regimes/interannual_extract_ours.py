"""
Daily region land-mean 2 m temperature for the 9 km and 31 km runs (the region series of
data/regions carry no 2t), for the interannual transfer test. 8 European regions
(regime_analysis.EU), Sep-Mar.

  dart9km  TCo1279-DART 1950C (1950-1969) and 2080C (2080-2092): 3-hourly reduced-Gaussian 2t,
           cells and unweighted means exactly as extract_dart.py (ERA5 lsm by nearest neighbour);
           every 2nd stamp (03, 09, 15, 21 UTC), which is enough for winter means
  tco319   TCo319-DART ctl1950d (1950-1969) and ssp585d (2080-2099): 6-hourly remapped 2t, land
           cells and cos-lat weights exactly as extract_tco319.py
Daily mean = mean of the used stamps in (d 00, d+1 00], as the wind series of geostrophic.ours.

Usage: python3 interannual_extract_ours.py {dart9km|tco319}     (legacy python; compute node)
Output: dunkelflaute/data/regimes/interannual_ours/<model>_<run>_t2.npz  (day (D,), t2 (D, 8) K)
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "1")
import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("", "regions"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/{p}")
import regions as rg                       # noqa: E402

EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
OUT = f"{_ROOT}/dunkelflaute/data/regimes/interannual_ours"
H = "/work/ab0995/ICCP_AWI_hackthon_2025"
NW = int(os.environ.get("N_WORKERS", 16))
_C = None


def daily(times, x):
    day = (times - np.timedelta64(1, "s")).astype("datetime64[D]")
    u, inv, cnt = np.unique(day, return_inverse=True, return_counts=True)
    s = np.zeros((len(u), x.shape[1]))
    np.add.at(s, inv, x)
    return u, s / cnt[:, None], cnt


# ---------------- 9 km
def _dart_init():
    global _C
    import extract_dart as xdart
    _C = xdart.setup()


def _dart_month(args):
    run, y, m = args
    C = _C
    f = f"{H}/TCo1279-DART-{run}/outdata/oifs/atm_reduced_3h_2t_3h_{y}{m:02d}-{y}{m:02d}.nc"
    if not os.path.exists(f):
        return None
    with nc.Dataset(f) as d:
        tv = d["time_counter"]
        times = np.array([np.datetime64(str(x)[:19]) for x in nc.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)])
        out = []
        for k in range(0, len(times), 2):           # every 2nd 3 h stamp: 03, 09, 15, 21 UTC
            x = d["2t"][k, C["lo"]:C["hi"]][C["union"]].astype(np.float64)
            out.append([x[C["pos"][r]].mean() for r in EU])
    return times[::2], np.array(out)


def dart():
    runs = {"1950C": range(1950, 1970), "2080C": range(2080, 2093)}
    for run, years in runs.items():
        t0 = time.time()
        jobs = [(run, y, m) for y in years for m in rg.MONTHS]
        T, X = [], []
        with ProcessPoolExecutor(NW, initializer=_dart_init) as ex:
            for r in ex.map(_dart_month, jobs):
                if r is not None:
                    T.append(r[0])
                    X.append(r[1])
        t, x = np.concatenate(T), np.concatenate(X)
        o = np.argsort(t)
        day, t2, cnt = daily(t[o], x[o])
        full = cnt == 4
        fn = f"{OUT}/dart9km_{run.lower()}_t2.npz"
        np.savez(fn, day=day[full], t2=t2[full], regions=np.array(EU))
        print(f"dart9km {run}: {full.sum()} days, {time.time()-t0:.0f}s -> {fn}", flush=True)


# ---------------- 31 km
def _tco319_init(run):
    global _C
    import extract_tco319 as xt
    _C = xt.setup(run)
    _C["path"] = xt.path


def _tco319_year(args):
    run, y = args
    S = _C
    with nc.Dataset(S["path"](run, "6h", "2t", y)) as d:
        tv = d["time_instant"]
        times = np.array([np.datetime64(str(x)[:19]) for x in nc.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)])
        month = rg.month_of(times)
        keep = np.isin(month, rg.MONTHS)
        idx = np.where(keep)[0]
        cols = []
        for a in range(0, len(idx), 120):
            ii = idx[a:a + 120]
            x = d["2t"][ii, S["r0"]:S["r1"], :].astype(np.float64).reshape(len(ii), -1)
            cols.append(np.stack([(x[:, S["cells"][r]] * S["w"][S["cells"][r]]).sum(1) / S["w"][S["cells"][r]].sum() for r in EU], axis=1))
    return times[keep], np.concatenate(cols)


def tco319():
    for run, years in (("ctl1950d", range(1950, 1970)), ("ssp585d", range(2080, 2100))):
        t0 = time.time()
        T, X = [], []
        with ProcessPoolExecutor(min(NW, len(years)), initializer=_tco319_init, initargs=(run,)) as ex:
            for t, x in ex.map(_tco319_year, [(run, y) for y in years]):
                T.append(t)
                X.append(x)
        t, x = np.concatenate(T), np.concatenate(X)
        o = np.argsort(t)
        day, t2, cnt = daily(t[o], x[o])
        full = cnt == 4
        fn = f"{OUT}/tco319_{run}_t2.npz"
        np.savez(fn, day=day[full], t2=t2[full], regions=np.array(EU))
        print(f"tco319 {run}: {full.sum()} days, {time.time()-t0:.0f}s -> {fn}", flush=True)


def main():
    os.makedirs(OUT, exist_ok=True)
    {"dart9km": dart, "tco319": tco319}[sys.argv[1]]()


if __name__ == "__main__":
    main()
