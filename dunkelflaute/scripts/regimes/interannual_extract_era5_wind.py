"""
ERA5 10 m wind speed for a method check of the observed interannual transfer test: the same daily
definition as the stations (mean of the 00, 06, 12, 18 UTC speeds), Sep-Mar of the winters
1980-2025, (a) as region land means (land as interannual_extract_era5.py, cos-lat weights) and
(b) at the nearest ERA5 grid point of every selected HadISD station. (b) run through the station
pipeline tells whether a difference between the station result and ERA5 comes from the sampling
(points, normalisation) or from the wind itself.

  source  /pool/data/ERA5/E5/sf/an/1H/165 and 166 (10u, 10v; daily files of 24 hourly GRIB messages,
          reduced Gaussian N320). Only the 4 needed messages per file are decoded.

Usage: python3 interannual_extract_era5_wind.py      (legacy python; compute node, N_WORKERS)
Output: dunkelflaute/data/regimes/interannual_era5_wind.npz
  day (D,), s_region (D, 8), s_station (D, N), station ids (N,)
"""
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regimes")
import regions as rg                       # noqa: E402
import interannual_extract_era5 as xe      # noqa: E402

EU = xe.EU
SRC = "/pool/data/ERA5/E5/sf/an/1H/{v}/E5sf00_1H_{d}_{v}.grb"
SEL = f"{_ROOT}/dunkelflaute/data/obs/hadisd_selection.json"
OUT = f"{_ROOT}/dunkelflaute/data/regimes/interannual_era5_wind.npz"
HOURS = (0, 6, 12, 18)
_S = {}


def messages(fn, want):
    """Raw GRIB messages number `want` (0-based) of a file, skipping the others by their length."""
    out = {}
    with open(fn, "rb") as fh:
        k = 0
        while True:
            p = fh.tell()
            h = fh.read(16)
            if len(h) < 16:
                break
            assert h[:4] == b"GRIB", fn
            n = int.from_bytes(h[4:7], "big") if h[7] == 1 else int.from_bytes(h[8:16], "big")
            if k in want:
                fh.seek(p)
                out[k] = fh.read(n)
            else:
                fh.seek(p + n)
            k += 1
    return out


def setup(h):
    import eccodes
    if "cells" in _S:
        return
    lat = eccodes.codes_get_array(h, "latitudes")
    lon = rg.lon180(eccodes.codes_get_array(h, "longitudes"))
    cells, w = xe.cells_of(h)
    _S["cells"], _S["w"] = cells, w
    from scipy.spatial import cKDTree
    sel = json.load(open(SEL))
    ids = sorted(sel["files"])
    sl = np.array([sel["info"][i][0] for i in ids])
    so = np.array([sel["info"][i][1] for i in ids])
    xyz = lambda a, b: np.column_stack([np.cos(np.deg2rad(a)) * np.cos(np.deg2rad(b)), np.cos(np.deg2rad(a)) * np.sin(np.deg2rad(b)), np.sin(np.deg2rad(a))])  # noqa: E731
    _S["near"] = cKDTree(xyz(lat, lon)).query(xyz(sl, so))[1]
    _S["ids"] = ids


def _day(d):
    import eccodes
    vals = {}
    for v in ("165", "166"):
        ms = messages(SRC.format(v=v, d=d), set(HOURS))
        xs = []
        for kk in HOURS:
            h = eccodes.codes_new_from_message(ms[kk])
            try:
                assert eccodes.codes_get(h, "dataTime") == kk * 100, (d, v, kk)
                setup(h)
                xs.append(eccodes.codes_get_values(h))
            finally:
                eccodes.codes_release(h)
        vals[v] = np.array(xs)
    spd = np.hypot(vals["165"], vals["166"]).mean(axis=0)
    reg = [float((spd[c] * _S["w"][c]).sum() / _S["w"][c].sum()) for c in _S["cells"]]
    return d, np.array(reg), spd[_S["near"]].astype(np.float32), _S["ids"]


def main():
    days = np.arange(np.datetime64("1979-09-01"), np.datetime64("2025-04-01"))
    days = days[np.isin(days.astype("datetime64[M]").astype(int) % 12 + 1, rg.MONTHS)]
    R, S, ids = [], [], None
    with ProcessPoolExecutor(int(os.environ.get("N_WORKERS", 64))) as ex:
        for n, (d, reg, st, ii) in enumerate(ex.map(_day, [str(x) for x in days], chunksize=8)):
            R.append(reg)
            S.append(st)
            ids = ii
            if n % 1000 == 0:
                print(f"  {d}", flush=True)
    np.savez(OUT, day=days, s_region=np.array(R), s_station=np.array(S), ids=np.array(ids), regions=np.array(EU))
    print(f"Saved {OUT}: {len(days)} days")


if __name__ == "__main__":
    main()
