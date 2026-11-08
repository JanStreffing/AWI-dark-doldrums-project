"""
Daily CMIP6 sea-level pressure (day/psl) on a global 2.5 deg grid, 15-80 N, Sep-Mar, for the
geostrophic-wind test and for assigning CMIP6 days to the ERA5 weather regimes.

Grid: lat 15.0, 17.5, ..., 80.0; lon -180.0, -177.5, ..., 177.5. Its 20-80 N, 80 W-40 E part is
the regime grid of dunkelflaute/scripts/regimes/extract_msl.py (box centres on the same 2.5 deg
points). Regridding is bilinear (periodic in longitude), not the box mean extract_msl.py uses for
the much finer model grids: most CMIP6 grids are about as coarse as 2.5 deg, so a box mean would
leave boxes empty. msl is smooth on these scales, so the two agree closely.

Winters as daily_events.py: historical 1951-1980, ssp585 2041-2100 (both future windows).

Usage: python3 extract_psl.py MODEL      (legacy python)
Output: dunkelflaute/data/cmip6/psl_<MODEL>.npz  ({exp}|ymd (D,3), {exp}|psl (D, 27, 144) float32 Pa)
"""
import glob
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "1")
import netCDF4 as nc
import numpy as np
from scipy.interpolate import RegularGridInterpolator

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
OUT = f"{_ROOT}/dunkelflaute/data/cmip6"
INV = f"{OUT}/inventory.json"
LAT = np.arange(15.0, 80.01, 2.5)
LON = np.arange(-180.0, 180.0, 2.5)
MONTHS = (9, 10, 11, 12, 1, 2, 3)
WINTERS = {"historical": (1951, 1980), "ssp585": (2041, 2100)}


def _span(fn):
    s = os.path.basename(fn).rsplit("_", 1)[-1].replace(".nc", "").split("-")
    return int(s[0][:4]), int(s[1][:4])


def _job(args):
    fn, year, w0, w1 = args
    with nc.Dataset(fn) as a:
        t = a["time"]
        ds = nc.num2date(t[:], t.units, calendar=getattr(t, "calendar", "standard"))
        ymd = np.array([(x.year, x.month, x.day) for x in ds])
        win = np.where(ymd[:, 1] >= 9, ymd[:, 0] + 1, ymd[:, 0])
        keep = (ymd[:, 0] == year) & np.isin(ymd[:, 1], MONTHS) & (win >= w0) & (win <= w1)
        if not keep.any():
            return np.zeros((0, 3), int), np.zeros((0, len(LAT), len(LON)), np.float32)
        lat = a["lat"][:].astype(float)
        lon = a["lon"][:].astype(float) % 360
        rows = np.where((lat >= LAT[0] - 5) & (lat <= LAT[-1] + 5))[0]
        r0, r1 = int(rows.min()), int(rows.max()) + 1
        idx = np.where(keep)[0]
        sl = slice(int(idx[0]), int(idx[-1]) + 1)
        x = a["psl"][sl, r0:r1, :].astype(np.float64)[keep[sl]]
    lat = lat[r0:r1]
    if lat[0] > lat[-1]:
        lat, x = lat[::-1], x[:, ::-1]
    o = np.argsort(lon)
    lon, x = lon[o], x[:, :, o]
    lonp = np.concatenate([lon[-1:] - 360, lon, lon[:1] + 360])          # periodic padding
    xp = np.concatenate([x[:, :, -1:], x, x[:, :, :1]], axis=2)
    LA, LO = np.meshgrid(LAT, LON % 360, indexing="ij")
    pts = np.column_stack([LA.ravel(), LO.ravel()])
    out = np.empty((len(x), len(LAT), len(LON)), np.float32)
    for k in range(len(x)):
        f = RegularGridInterpolator((lat, lonp), xp[k], bounds_error=True)
        out[k] = f(pts).reshape(len(LAT), len(LON))
    return ymd[keep], out


def main():
    model = sys.argv[1]
    inv = json.load(open(INV))[model]
    res = {}
    for exp, (w0, w1) in WINTERS.items():
        base = inv["paths"][exp]["day/sfcWind"].split("/day/")[0]
        d = sorted(x for x in glob.glob(f"{base}/day/psl/*/*") if glob.glob(f"{x}/*.nc"))[-1]
        jobs = []
        for f in sorted(glob.glob(f"{d}/*.nc")):
            a, b = _span(f)
            for y in range(max(a, w0 - 1), min(b, w1) + 1):
                jobs.append((f, y, w0, w1))
        got = {}
        with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 32))) as ex:
            for ymd, x in ex.map(_job, jobs):
                for k in range(len(ymd)):
                    got[tuple(ymd[k])] = x[k]
        keys = sorted(got)
        res[f"{exp}|ymd"] = np.array(keys)
        res[f"{exp}|psl"] = np.stack([got[k] for k in keys])
        print(f"  {model} {exp}: {len(keys)} days {keys[0]} .. {keys[-1]}", flush=True)
    fn = f"{OUT}/psl_{model}.npz"
    np.savez(fn, **res, lat=LAT, lon=LON)
    print(f"Saved {fn}")


if __name__ == "__main__":
    main()
