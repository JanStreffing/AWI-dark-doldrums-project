"""
Daily-mean msl of our runs on the CMIP6 psl grid, for the transfer artifact checks (A2, A3).

Grid and interpolation exactly as dunkelflaute/scripts/cmip6/extract_psl.py: lat 15.0, 17.5, ..., 80.0,
lon -180.0, ..., 177.5, bilinear, periodic in longitude. This replaces the 2.5 deg box mean of
extract_msl.py (Euro-Atlantic sector only) by the construction CMIP6 gets, over the whole 15-80 N band.

Daily means as extract_msl.py:
  tco95    ctl1950g / hi1950g, daily file atm_remapped_1d_msl (400 x 192)
  tco319   ctl1950d / ssp585d, mean of the 06/12/18/24 UTC instants (1312 x 640)
  dart9km  1950c / 2080c, mean of the 06/12/18/24 UTC instants (5136 x 2560)
  eerie    ctl / hist / ssp245, daily average mmsl (0.25 deg), GRIB via the byte index
Only the winters the comparisons use (union over region_events.COMPARISONS) are written.

Usage: python3 transfer_checks_msl.py SOURCE RUN      (legacy python)
Output: dunkelflaute/data/regimes/transfer_checks_msl_<SOURCE>_<RUN>.npz
        (dates (D,), msl (D, 27, 144) float32 Pa, lat, lon)
"""
import glob
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
import regions as rg            # noqa: E402
import region_events as ev      # noqa: E402
import extract_msl as em        # noqa: E402

OUT = f"{_ROOT}/dunkelflaute/data/regimes"
LAT = np.arange(15.0, 80.01, 2.5)
LON = np.arange(-180.0, 180.0, 2.5)
MONTHS = rg.MONTHS


def winter_range(src, run):
    ys = [(a, b) for c in ev.COMPARISONS.values() for n, a, b in c if n == f"{src}_{run}"]
    assert ys, f"{src}_{run} is in no comparison"
    return min(a for a, _ in ys), max(b for _, b in ys)


def in_winters(day, w0, w1):
    y = int(str(day)[:4])
    m = int(str(day)[5:7])
    w = y + 1 if m >= 9 else y
    return m in MONTHS and w0 <= w <= w1


def bilinear_setup(lat, lon):
    """Source axes (any order) -> rows to read and the 4-point weights onto (LAT, LON)."""
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float) % 360
    oy = np.argsort(lat)
    ls = lat[oy]
    j = np.clip(np.searchsorted(ls, LAT) - 1, 0, len(ls) - 2)
    wy = (LAT - ls[j]) / (ls[j + 1] - ls[j])
    assert ((wy >= 0) & (wy <= 1)).all()
    ra, rb = oy[j], oy[j + 1]                       # source row indices
    rows = np.unique(np.concatenate([ra, rb]))
    ia, ib = np.searchsorted(rows, ra), np.searchsorted(rows, rb)
    ox = np.argsort(lon)
    lo = lon[ox]
    lop = np.concatenate([lo[-1:] - 360, lo, lo[:1] + 360])
    oxp = np.concatenate([ox[-1:], ox, ox[:1]])
    tl = LON % 360
    k = np.clip(np.searchsorted(lop, tl) - 1, 0, len(lop) - 2)
    wx = (tl - lop[k]) / (lop[k + 1] - lop[k])
    assert ((wx >= 0) & (wx <= 1)).all()
    return dict(rows=rows, ia=ia, ib=ib, wy=wy, ja=oxp[k], jb=oxp[k + 1], wx=wx)


def bilinear(x, B):
    """x (..., len(rows), nlon) -> (..., 27, 144)."""
    a, b = x[..., B["ia"], :], x[..., B["ib"], :]
    wy, wx = B["wy"][:, None], B["wx"][None, :]
    top = a[..., B["ja"]] * (1 - wx) + a[..., B["jb"]] * wx
    bot = b[..., B["ja"]] * (1 - wx) + b[..., B["jb"]] * wx
    return (top * (1 - wy) + bot * wy).astype(np.float32)


def read_rows(v, rows, tsel=slice(None)):
    """v (T, lat, lon) netCDF variable -> (T', len(rows), lon), reading contiguous row runs."""
    runs, s = [], 0
    for i in range(1, len(rows) + 1):
        if i == len(rows) or rows[i] != rows[i - 1] + 1:
            runs.append((s, i))
            s = i
    parts = [np.asarray(v[tsel, int(rows[a]):int(rows[b - 1]) + 1, :], dtype=np.float64) for a, b in runs]
    return np.concatenate(parts, axis=1)


# ---------------- regular lat/lon NetCDF
def _reg_file(args):
    f, mode, w0, w1 = args
    with nc.Dataset(f) as d:
        B = bilinear_setup(d["lat"][:], d["lon"][:])
        t = em._times(d)
        day = (t - np.timedelta64(1, "s")).astype("datetime64[D]")
        out_d, out_x = [], []
        if mode == "daily":
            k = np.array([i for i in range(len(t)) if in_winters(day[i], w0, w1)])
            if len(k):
                x = read_rows(d["msl"], B["rows"], slice(int(k[0]), int(k[-1]) + 1))[k - k[0]]
                out_d, out_x = list(day[k]), list(bilinear(x, B))
        else:
            hour = (t - t.astype("datetime64[D]")) // np.timedelta64(1, "h")
            days = [dd for dd in np.unique(day) if in_winters(dd, w0, w1)]
            if days:
                ks = {dd: np.where((day == dd) & (hour % 6 == 0))[0] for dd in days}
                days = [dd for dd in days if len(ks[dd]) == 4]
            if days:
                k0 = min(int(ks[dd][0]) for dd in days)
                k1 = max(int(ks[dd][-1]) for dd in days) + 1
                x = read_rows(d["msl"], B["rows"], slice(k0, k1))
                for dd in days:
                    out_d.append(dd)
                    out_x.append(bilinear(x[ks[dd] - k0].mean(axis=0), B))
    return f, out_d, out_x


def regular(files, mode, w0, w1, workers):
    res = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for f, dd, xx in ex.map(_reg_file, [(f, mode, w0, w1) for f in files]):
            for a, b in zip(dd, xx):
                res[a] = b
            print(f"  {os.path.basename(f)}: {len(dd)} days", flush=True)
    return res


# ---------------- EERIE GRIB via the byte index
def _eerie_chunk(args):
    import eccodes
    run, rows, w0, w1 = args
    ix = np.load(f"{em.EIDX}/{run}_mmsl.npz")
    out, B = [], None
    for k in rows:
        if ix["path"][k] == "":
            continue
        with open(ix["path"][k], "rb") as fh:
            fh.seek(int(ix["offset"][k]))
            h = eccodes.codes_new_from_message(rg.grib2_strip_local(fh.read(int(ix["size"][k]))))
        try:
            d = eccodes.codes_get(h, "dataDate")
            day = np.datetime64(f"{d // 10000:04d}-{d // 100 % 100:02d}-{d % 100:02d}")
            if not in_winters(day, w0, w1):
                continue
            ni, nj = eccodes.codes_get(h, "Ni"), eccodes.codes_get(h, "Nj")
            if B is None:
                la = eccodes.codes_get_array(h, "latitudes").reshape(nj, ni)
                lo = eccodes.codes_get_array(h, "longitudes").reshape(nj, ni)
                B = bilinear_setup(la[:, 0], lo[0, :])
            x = eccodes.codes_get_values(h).reshape(nj, ni)[B["rows"]]
            out.append((day, bilinear(x, B)))
        finally:
            eccodes.codes_release(h)
    return out


def eerie(run, w0, w1, workers):
    n = len(np.load(f"{em.EIDX}/{run}_mmsl.npz")["path"])
    chunks = [(run, list(range(a, min(n, a + 365))), w0, w1) for a in range(0, n, 365)]
    res = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for out in ex.map(_eerie_chunk, chunks):
            for a, b in out:
                res[a] = b
            print(f"  chunk: {len(out)} days", flush=True)
    return res


def main():
    src, run = sys.argv[1], sys.argv[2]
    w = int(os.environ.get("N_WORKERS", 16))
    w0, w1 = winter_range(src, run)
    t0 = time.time()
    print(f"{src} {run}: winters {w0}-{w1}", flush=True)
    yrs = range(w0 - 1, w1 + 1)
    if src == "dart9km":
        D = {"1950c": "TCo1279-DART-1950C", "2080c": "TCo1279-DART-2080C"}[run]
        files = [f for f in sorted(glob.glob(f"{em.H}/{D}/outdata/oifs/atm_remapped_3h_msl_3h_*.nc"))
                 if int(f.split("_")[-1][4:6]) in MONTHS and int(f.split("_")[-1][:4]) in yrs]
        res = regular(files, "instant", w0, w1, w)
    elif src == "tco319":
        D = {"ctl1950d": "TCo319-DART-ctl1950d-gibbs-charn", "ssp585d": "TCo319-DART-ssp585d-gibbs-charn"}[run]
        files = [f for f in sorted(glob.glob(f"{em.H}/{D}/outdata/oifs/atm_remapped_6h_msl_6h_*.nc"))
                 if int(os.path.basename(f).split("_")[-1][:4]) in yrs]
        res = regular(files, "instant", w0, w1, w)
    elif src == "tco95":
        files = [f for f in sorted(glob.glob(f"{em.RT}/{run}/outdata/oifs/atm_remapped_1d_msl_1d_*.nc"))
                 if int(os.path.basename(f).split("_")[-1][:4]) in yrs]
        res = regular(files, "daily", w0, w1, w)
    elif src == "eerie":
        res = eerie(run, w0, w1, w)
    else:
        raise ValueError(src)
    dates = np.array(sorted(res))
    os.makedirs(OUT, exist_ok=True)
    fn = f"{OUT}/transfer_checks_msl_{src}_{run}.npz"
    np.savez(fn, dates=dates, msl=np.stack([res[d] for d in dates]), lat=LAT, lon=LON)
    print(f"Saved {fn}: {len(dates)} days, {dates[0]} .. {dates[-1]} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
