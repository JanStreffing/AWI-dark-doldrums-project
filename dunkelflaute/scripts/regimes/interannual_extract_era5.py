"""
ERA5 monthly region land means of 2 m temperature and 850 hPa temperature, Sep-Mar of the winters
1979-2025, for the observed interannual transfer test (ERA5 t2m as a check on the station
temperature; t850 - t2m as a winter inversion index, ERA5 temperature profiles being
radiosonde-constrained).

  sources  /pool/data/ERA5/E5/sf/an/1M/167 (2t) and /pool/data/ERA5/E5/pl/an/1M/130 (t, level 850),
           reduced Gaussian N320
  land     ERA5 lsm (regions.ERA5_LSM, 0.25 deg) at the nearest 0.25 deg point >= 0.5;
           region cells by regions.region_cells (DE = Germany polygon); cos-lat weights

Usage: python3 interannual_extract_era5.py      (legacy python, login node, a few minutes)
Output: dunkelflaute/data/regimes/interannual_era5_t.npz  (ym (M, 2), t2m (M, 8), t850 (M, 8), K)
"""
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
import regions as rg                       # noqa: E402

EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
SF = "/pool/data/ERA5/E5/sf/an/1M/167/E5sf00_1M_{y}_167.grb"
PL = "/pool/data/ERA5/E5/pl/an/1M/130/E5pl00_1M_{y}_130.grb"
OUT = f"{_ROOT}/dunkelflaute/data/regimes/interannual_era5_t.npz"
_C = {}


def cells_of(h):
    import eccodes
    if "cells" not in _C:
        lat = eccodes.codes_get_array(h, "latitudes")
        lon = rg.lon180(eccodes.codes_get_array(h, "longitudes"))
        with nc.Dataset(rg.ERA5_LSM) as d:
            elat = d["latitude"][:].astype(float)
            elon = rg.lon180(d["longitude"][:].astype(float))
            lsm = d["lsm"][:].squeeze().astype(float)
        i = np.clip(np.rint((elat[0] - lat) / abs(elat[1] - elat[0])).astype(int), 0, len(elat) - 1)
        oj = np.argsort(elon)
        j = oj[np.clip(np.searchsorted(elon[oj], lon), 0, len(elon) - 1)]
        land = lsm[i, j] >= 0.5
        cells = rg.region_cells(lat, lon, land)
        _C["cells"] = [cells[r] for r in EU]
        _C["w"] = np.cos(np.deg2rad(lat))
    return _C["cells"], _C["w"]


def _year(y):
    import eccodes
    res = {}
    for f, var, lev in ((SF.format(y=y), "t2m", None), (PL.format(y=y), "t850", 850)):
        with open(f, "rb") as fh:
            while True:
                h = eccodes.codes_grib_new_from_file(fh)
                if h is None:
                    break
                try:
                    d = eccodes.codes_get(h, "dataDate")
                    m = d // 100 % 100
                    if m not in rg.MONTHS or (lev is not None and eccodes.codes_get(h, "level") != lev):
                        continue
                    cells, w = cells_of(h)
                    x = eccodes.codes_get_values(h)
                    res[(var, y, m)] = [float((x[c] * w[c]).sum() / w[c].sum()) for c in cells]
                finally:
                    eccodes.codes_release(h)
    return res


def main():
    got = {}
    with ProcessPoolExecutor(8) as ex:
        for r in ex.map(_year, range(1978, 2026)):
            got.update(r)
    ym = sorted({(y, m) for _, y, m in got if (y + 1 if m >= 9 else y) in range(1979, 2026)})
    ym = [k for k in ym if ("t2m", *k) in got and ("t850", *k) in got]
    np.savez(OUT, ym=np.array(ym), t2m=np.array([got[("t2m", *k)] for k in ym]),
             t850=np.array([got[("t850", *k)] for k in ym]), regions=np.array(EU))
    print(f"Saved {OUT}: {len(ym)} months {ym[0]} .. {ym[-1]}")


if __name__ == "__main__":
    main()
