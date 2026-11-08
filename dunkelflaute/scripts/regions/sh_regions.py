"""
Southern-Hemisphere land regions for the SH test of the geostrophic mechanism (sh_extract.py,
sh_geostrophic.py). Counterpart of regions.py, which holds the Northern-Hemisphere boxes.

Choice of regions (2026-11-07):
  * Seven named lat/lon boxes over the SH mid-latitude land masses that have wind-power interest
    and lie in the westerly or subtropical-high belt. Land cells (land fraction >= 0.5) only, as in
    regions.py. Every box keeps its equatorward edge at or poleward of 28 S, so the Coriolis
    parameter is not small (geostrophic.py drops NH boxes with a southern edge south of 29 N;
    sin 28 deg / sin 29 deg = 0.97).
      AU_SE  southern and eastern Australia (SA south-east, Victoria, NSW, south Queensland)
      AU_SW  south-western Australia
      ZA     southern Africa (South Africa south of 28 S, Lesotho)
      AR_C   central Chile and Argentina (Chilean central valley, Cuyo, Pampas)
      PA     Patagonia (Argentina and Chile south of 38 S)
      BR_S   southern Brazil and Uruguay
      NZ     New Zealand
  * A 5 x 5 deg land grid over 25-55 S (half-open boxes, no cell in two boxes), kept where the
    ERA5 land fraction of the box is >= 0.25 (as the NH grid5 of regions.py). Equatorward edge 25 S
    keeps f at >= 0.86 of its 29-degree value.
SH winter = April-September of one calendar year, labelled by that year.

Helpers:
  sh_boxes()                    {name: (lat_min, lat_max, lon_min, lon_max)}, named + grid
  sh_cells(lat, lon, land)      {name: index array into a flat cell list}
  era5_z0_sh(lat, lon, dl, dl)  (12, n) land z0, the ERA5 fsr climatology of regions.era5_z0_on
"""
import json
import os
import sys

import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
import regions as rg   # noqa: E402

NAMED = {
    "AU_SE": (-39.0, -28.0, 138.0, 154.0),
    "AU_SW": (-35.5, -28.0, 114.0, 124.0),
    "ZA":    (-35.0, -28.0,  17.0,  32.0),
    "AR_C":  (-38.0, -30.0, -72.0, -57.0),
    "PA":    (-52.0, -38.0, -73.0, -63.0),
    "BR_S":  (-35.0, -28.0, -58.0, -48.0),
    "NZ":    (-47.5, -34.0, 166.0, 179.0),
}
LONG = {"AU_SE": "southern and eastern Australia", "AU_SW": "south-western Australia",
        "ZA": "southern Africa", "AR_C": "central Chile and Argentina", "PA": "Patagonia",
        "BR_S": "southern Brazil and Uruguay", "NZ": "New Zealand"}
GRID_DEG, GRID_LAT, GRID_MIN_LAND = 5.0, (-55.0, -25.0), 0.25
LAT_BAND = (-56.0, -24.0)                   # every SH box lies inside
MONTHS = (4, 5, 6, 7, 8, 9)
OUT_DIR = f"{_ROOT}/dunkelflaute/data/sh"
# 2.5 deg SH sea-level-pressure grid (global ring, so periodic in longitude)
PLAT = np.arange(-65.0, -14.99, 2.5)
PLON = np.arange(-180.0, 180.0, 2.5)


def _grid_boxes():
    cache = f"{OUT_DIR}/sh_grid5_boxes.json"
    if os.path.exists(cache):
        return {k: tuple(v) for k, v in json.load(open(cache)).items()}
    with nc.Dataset(rg.ERA5_LSM) as d:
        elat = d["latitude"][:].astype(float)
        elon = rg.lon180(d["longitude"][:].astype(float))
        elsm = d["lsm"][:].squeeze().astype(float)
    out = {}
    eps = 1e-6
    for la0 in np.arange(GRID_LAT[0], GRID_LAT[1], GRID_DEG):
        mi = (elat >= la0) & (elat < la0 + GRID_DEG)
        for lo0 in np.arange(-180.0, 180.0, GRID_DEG):
            mj = (elon >= lo0) & (elon < lo0 + GRID_DEG)
            if elsm[np.ix_(mi, mj)].mean() >= GRID_MIN_LAND:
                out[f"g{la0:+05.1f}{lo0:+06.1f}"] = (float(la0), float(la0 + GRID_DEG - eps),
                                                    float(lo0), float(lo0 + GRID_DEG - eps))
    os.makedirs(OUT_DIR, exist_ok=True)
    json.dump(out, open(cache, "w"), indent=0)
    return out


def sh_boxes():
    return dict(NAMED, **_grid_boxes())


def sh_cells(lat, lon, land):
    lat = np.asarray(lat, float)
    lon = rg.lon180(lon)
    land = np.asarray(land, bool)
    return {k: np.where(land & (lat >= a) & (lat <= b) & (lon >= c) & (lon <= d))[0]
            for k, (a, b, c, d) in sh_boxes().items()}


def era5_z0_sh(lat, lon, dlat, dlon):
    """regions.era5_z0_on with the SH latitude band."""
    old = rg.LAT_BAND
    rg.LAT_BAND = LAT_BAND
    try:
        return rg.era5_z0_on(lat, lon, dlat, dlon)
    finally:
        rg.LAT_BAND = old


def era5_z0_clim_global():
    """(12, 721, 1440) ERA5 fsr climatology 2016-2025 (geometric mean per calendar month)."""
    with nc.Dataset(rg.ERA5_FSR) as d:
        tv = d["valid_time"]
        yrs = np.array([t.year for t in nc.num2date(tv[:], tv.units)])
        fsr = d["fsr"][yrs >= 2016].astype(float)
    return np.exp(np.log(np.maximum(fsr.reshape(-1, 12, *fsr.shape[1:]).mean(axis=0), 1e-5)))


def psl_boxes(lat, lon):
    """2.5 deg box id (row-major on PLAT x PLON) per source point, -1 outside; longitude wraps."""
    lat = np.asarray(lat, float)
    lon = rg.lon180(lon)
    i = np.rint((lat - PLAT[0]) / 2.5).astype(int)
    j = np.rint((lon - PLON[0]) / 2.5).astype(int) % len(PLON)
    ok = (i >= 0) & (i < len(PLAT))
    return np.where(ok, i * len(PLON) + j, -1)
