"""
Daily-mean mean-sea-level pressure on a common 2.5 deg Euro-Atlantic grid (20-80 N, 80 W-40 E,
25 x 49 box centres), Sep-Mar, for the weather-regime analysis.

Regridding is a box mean: every source grid point whose centre falls inside a 2.5 deg box
(centre +-1.25 deg) is averaged with equal weight (all source grids are fine compared to 2.5 deg,
and their cells are near-equal-area within a box). This works the same for regular and reduced
Gaussian grids.

Sources and daily means:
  era5        /pool/data/ERA5/E5/sf/an/1D/151 (daily average, reduced Gaussian N320), 1979-2025
  dart9km     TCo1279-DART 1950C / 2080C, 3-hourly msl on the remapped 5136 x 2560 grid: mean of the
              06/12/18/24 UTC instants of each day
  tco319      ctl1950d / ssp585d, 6-hourly: mean of the 06/12/18/24 UTC instants of each day
  tco95       ctl1950g / hi1950g, daily file atm_remapped_1d_msl
  eerie       hist / ctl / ssp245, daily average mmsl (0.25 deg), GRIB via the byte index

Usage: python3 extract_msl.py SOURCE RUN     (legacy python; RUN '-' for era5)
Output: dunkelflaute/data/regimes/msl_<SOURCE>_<RUN>.npz  (dates (D,), msl (D, 25, 49) float32, Pa)
"""
import glob
import os
import sys
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "1")
import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
import regions as rg   # noqa: E402
OUT = f"{_ROOT}/dunkelflaute/data/regimes"
LATC = np.arange(20.0, 80.01, 2.5)
LONC = np.arange(-80.0, 40.01, 2.5)
NB = len(LATC) * len(LONC)
MONTHS = (9, 10, 11, 12, 1, 2, 3)
H = "/work/ab0995/ICCP_AWI_hackthon_2025"
RT = "/work/ab0995/a270092/runtime/awicm3-frontiers-xios"
ERA5 = "/pool/data/ERA5/E5/sf/an/1D/151"
EIDX = f"{_ROOT}/dunkelflaute/data/regions/eerie_index"


def boxes(lat, lon):
    """Box id per point (or -1 outside the sector)."""
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    lon = np.where(lon > 180, lon - 360, lon)
    i = np.rint((lat - LATC[0]) / 2.5).astype(int)
    j = np.rint((lon - LONC[0]) / 2.5).astype(int)
    ok = (i >= 0) & (i < len(LATC)) & (j >= 0) & (j < len(LONC))
    return np.where(ok, i * len(LONC) + j, -1)


def box_mean(values, bid, cnt):
    s = np.bincount(bid, weights=values, minlength=NB)
    return (s / cnt).reshape(len(LATC), len(LONC)).astype(np.float32)


def _date(t):
    return np.datetime64(str(t)[:10])


# ---------------- regular lat/lon NetCDF (dart9km remapped, tco319, tco95)
_G = {}


def _reg_setup(f):
    with nc.Dataset(f) as d:
        lat = d["lat"][:].astype(float)
        lon = d["lon"][:].astype(float)
    rows = np.where((lat >= LATC[0] - 1.25) & (lat <= LATC[-1] + 1.25))[0]
    r0, r1 = int(rows[0]), int(rows[-1]) + 1
    bid2 = boxes(np.repeat(lat[r0:r1], len(lon)), np.tile(lon, r1 - r0))
    sel = np.where(bid2 >= 0)[0]
    bid = bid2[sel]
    cnt = np.bincount(bid, minlength=NB).astype(float)
    assert (cnt > 0).all(), "every 2.5 deg box needs source points"
    return dict(r0=r0, r1=r1, sel=sel, bid=bid, cnt=cnt)


def _times(d):
    tv = d["time_instant"] if "time_instant" in d.variables else d["time_counter"]
    return np.array([np.datetime64(str(x)[:19]) for x in nc.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)])


def _reg_file(args):
    f, mode = args
    if f not in _G:
        _G[f] = _reg_setup(f)
    G = _G[f]
    out_d, out_x = [], []
    with nc.Dataset(f) as d:
        t = _times(d)
        v = d["msl"]
        if mode == "daily":
            # a daily-mean file: one record per day, stamped inside (or at the end of) the day
            day = (t - np.timedelta64(1, "s")).astype("datetime64[D]")
            for k in range(len(t)):
                m = int(str(day[k])[5:7])
                if m in MONTHS:
                    x = v[k, G["r0"]:G["r1"], :].astype(np.float64).ravel()[G["sel"]]
                    out_d.append(day[k])
                    out_x.append(box_mean(x, G["bid"], G["cnt"]))
        else:
            # a day is the window (d 00 UTC, d+1 00 UTC]: instants 06, 12, 18 and 24 UTC, so every
            # day lies inside one monthly (DART) or yearly (TCo319) file
            hour = (t - t.astype("datetime64[D]")) // np.timedelta64(1, "h")
            day = (t - np.timedelta64(1, "s")).astype("datetime64[D]")
            for dd in np.unique(day):
                if int(str(dd)[5:7]) not in MONTHS:
                    continue
                k = np.where((day == dd) & (hour % 6 == 0))[0]
                if len(k) != 4:
                    continue
                acc = 0.0
                for i in k:
                    acc = acc + v[i, G["r0"]:G["r1"], :].astype(np.float64).ravel()[G["sel"]]
                out_d.append(dd)
                out_x.append(box_mean(acc / 4, G["bid"], G["cnt"]))
    return f, out_d, out_x


def regular(files, mode, workers):
    res = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for f, dd, xx in ex.map(_reg_file, [(f, mode) for f in files]):
            for a, b in zip(dd, xx):
                res[a] = b
            print(f"  {os.path.basename(f)}: {len(dd)} days", flush=True)
    return res


# ---------------- GRIB (era5 reduced Gaussian files, eerie byte index)
_E = {}


def _grib_box(h, key):
    import eccodes
    if key not in _E:
        lat = eccodes.codes_get_array(h, "latitudes")
        lon = eccodes.codes_get_array(h, "longitudes")
        b = boxes(lat, lon)
        sel = np.where(b >= 0)[0]
        cnt = np.bincount(b[sel], minlength=NB).astype(float)
        assert (cnt > 0).all()
        _E[key] = (sel, b[sel], cnt)
    return _E[key]


def _era5_month(f):
    import eccodes
    out = []
    with open(f, "rb") as fh:
        while True:
            h = eccodes.codes_grib_new_from_file(fh)
            if h is None:
                break
            try:
                d = eccodes.codes_get(h, "dataDate")
                sel, bid, cnt = _grib_box(h, "era5")
                x = eccodes.codes_get_values(h)[sel]
                out.append((np.datetime64(f"{d // 10000:04d}-{d // 100 % 100:02d}-{d % 100:02d}"), box_mean(x, bid, cnt)))
            finally:
                eccodes.codes_release(h)
    return f, out


def era5(workers):
    files = []
    for f in sorted(glob.glob(f"{ERA5}/E5sf00_1D_*_151.grb")):
        y, m = map(int, os.path.basename(f).split("_")[2].split("-"))
        if 1979 <= y <= 2025 and m in MONTHS:
            files.append(f)
    res = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for f, out in ex.map(_era5_month, files):
            for a, b in out:
                res[a] = b
            print(f"  {os.path.basename(f)}: {len(out)} days", flush=True)
    return res


def _eerie_chunk(args):
    import eccodes
    run, rows = args
    ix = np.load(f"{EIDX}/{run}_mmsl.npz")
    out = []
    for k in rows:
        with open(ix["path"][k], "rb") as fh:
            fh.seek(int(ix["offset"][k]))
            h = eccodes.codes_new_from_message(rg.grib2_strip_local(fh.read(int(ix["size"][k]))))
        try:
            d = eccodes.codes_get(h, "dataDate")
            if (d // 100 % 100) not in MONTHS:
                continue
            sel, bid, cnt = _grib_box(h, "eerie")
            x = eccodes.codes_get_values(h)[sel]
            out.append((np.datetime64(f"{d // 10000:04d}-{d // 100 % 100:02d}-{d % 100:02d}"), box_mean(x, bid, cnt)))
        finally:
            eccodes.codes_release(h)
    return out


def eerie(run, workers):
    n = len(np.load(f"{EIDX}/{run}_mmsl.npz")["path"])
    chunks = [(run, list(range(a, min(n, a + 365)))) for a in range(0, n, 365)]
    res = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for out in ex.map(_eerie_chunk, chunks):
            for a, b in out:
                res[a] = b
            print(f"  chunk: {len(out)} Sep-Mar days", flush=True)
    return res


def main():
    src, run = sys.argv[1], sys.argv[2]
    w = int(os.environ.get("N_WORKERS", 32))
    if src == "era5":
        res = era5(w)
    elif src == "dart9km":
        D = {"1950c": "TCo1279-DART-1950C", "2080c": "TCo1279-DART-2080C"}[run]
        files = [f for f in sorted(glob.glob(f"{H}/{D}/outdata/oifs/atm_remapped_3h_msl_3h_*.nc"))
                 if int(f.split("_")[-1][4:6]) in MONTHS]
        res = regular(files, "instant", w)
    elif src == "tco319":
        D = {"ctl1950d": "TCo319-DART-ctl1950d-gibbs-charn", "ssp585d": "TCo319-DART-ssp585d-gibbs-charn"}[run]
        res = regular(sorted(glob.glob(f"{H}/{D}/outdata/oifs/atm_remapped_6h_msl_6h_*.nc")), "instant", w)
    elif src == "tco95":
        res = regular(sorted(glob.glob(f"{RT}/{run}/outdata/oifs/atm_remapped_1d_msl_1d_*.nc")), "daily", w)
    elif src == "eerie":
        res = eerie(run, w)
    else:
        raise ValueError(src)
    dates = np.array(sorted(res))
    os.makedirs(OUT, exist_ok=True)
    fn = f"{OUT}/msl_{src}_{run}.npz"
    np.savez(fn, dates=dates, msl=np.stack([res[d] for d in dates]), lat=LATC, lon=LONC)
    print(f"Saved {fn}: {len(dates)} days, {dates[0]} .. {dates[-1]}")


if __name__ == "__main__":
    main()
