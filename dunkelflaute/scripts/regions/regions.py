"""
Region definitions and grid helpers for the multi-region test of the 9 km Dunkelflaute increase.

Northern-Hemisphere land boxes only. The 9 km DART wind and solar reconstructions exist for Sep-Mar
alone, so a Southern-Hemisphere winter (Apr-Sep) has no 9 km data. Every region is land cells only
(lsm >= 0.5) inside a lat/lon box, except DE, which is the Germany land polygon used everywhere else
in the project. Offshore Germany is kept as its own entry (DE_off: sea cells inside the German EEZ)
so the Mockert Germany mix can be rebuilt where a model has not been run through the Germany
pipeline (EERIE).

Shared pieces used by every extractor:
  region_cells(lat, lon, land)  -> {name: index array into the flat cell list}
  era5_z0_on(lat, lon, ...)     -> (12, n) land z0 from the ERA5 fsr climatology 2016-2025,
                                   geometric mean of the ERA5 land cells nearest to each model cell
  loglaw_land / loglaw_sea      -> 100 m speed from 10 m speed (fixed z0 / Charnock)
"""
import json
import os

import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
ERA5_FSR = f"{_ROOT}/wind/data/era5/era5_fsr_monthly_2015_2025.nc"
ERA5_LSM = f"{_ROOT}/wind/data/era5/era5_lsm_20260830_0000.nc"
DE_LAND = f"{_ROOT}/boundaries/germany/germany_land.geojson"
DE_EEZ = f"{_ROOT}/boundaries/germany/germany_eez.geojson"
QM_TABLE = f"{_ROOT}/wind/models/quantile_mapping_v1/qm_table.npz"
OUT_DIR = f"{_ROOT}/dunkelflaute/data/regions"

# (lat_min, lat_max, lon_min, lon_max), lon in -180..180
BOXES = {
    "FR":    (43.0, 51.0,   -4.5,   7.5),
    "IB":    (36.5, 43.5,   -9.5,   3.0),
    "UK":    (50.0, 58.5,  -10.5,   1.8),
    "PL":    (49.5, 54.5,   14.5,  24.0),
    "SC":    (56.0, 64.0,    5.0,  25.0),
    "IT":    (37.0, 46.5,    7.0,  18.5),
    "BK":    (42.0, 48.5,   16.0,  28.0),
    "UA":    (46.0, 52.0,   24.0,  40.0),
    "RU":    (52.0, 60.0,   30.0,  50.0),
    "TR":    (37.0, 41.5,   27.0,  44.0),
    "KZ":    (45.0, 53.0,   60.0,  80.0),
    "US_NP": (40.0, 49.0, -104.0, -94.0),
    "US_SP": (29.0, 37.0, -103.0, -95.0),
    "US_MW": (37.0, 44.0,  -94.0, -82.0),
    "US_NE": (39.0, 45.0,  -80.0, -70.0),
    "US_CA": (33.0, 41.0, -123.0, -115.0),
    "CA_PR": (49.0, 55.0, -114.0, -98.0),
    "CN_NE": (40.0, 48.0,  115.0, 128.0),
    "CN_N":  (34.0, 40.0,  110.0, 120.0),
    "CN_E":  (26.0, 33.0,  110.0, 121.0),
    "KR":    (34.0, 38.7,  126.0, 129.6),
    "JP":    (33.0, 41.0,  130.0, 142.0),
    "IN_N":  (22.0, 30.0,   72.0,  88.0),
}
DE_BOX = (46.0, 56.5, 4.0, 16.5)
LAT_BAND = (20.0, 67.0)                 # every region lies inside; extractors read only this band
BASE_DIR = OUT_DIR                      # region-set independent inputs (EERIE byte index, global 2t)

# Grid mode for the per-K maps: REGION_SET=grid5 replaces the boxes by 5 x 5 deg land boxes over
# 20-70 N (half-open, so no cell falls in two boxes), kept where the ERA5 land fraction of the box
# is >= GRID_MIN_LAND. DE (polygon) and DE_off stay, so every extractor runs unchanged; outputs go
# to data/regions_grid5 (and data/cmip6_grid5 for the daily extractor).
REGION_SET = os.environ.get("REGION_SET", "regions")
GRID_DEG, GRID_LAT, GRID_MIN_LAND = 5.0, (20.0, 70.0), 0.25


def _grid_boxes():
    cache = f"{_ROOT}/dunkelflaute/data/regions_grid5/grid5_boxes.json"
    if os.path.exists(cache):
        return {k: tuple(v) for k, v in json.load(open(cache)).items()}
    with nc.Dataset(ERA5_LSM) as d:
        elat = d["latitude"][:].astype(float)
        elon = d["longitude"][:].astype(float)
        elon = np.where(elon > 180, elon - 360, elon)
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
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    json.dump(out, open(cache, "w"), indent=0)
    return out


REGION_BOXES = dict(BOXES)              # the 23 named boxes, whatever REGION_SET is
if REGION_SET == "grid5":
    BOXES = _grid_boxes()
    LAT_BAND = GRID_LAT
    OUT_DIR = f"{_ROOT}/dunkelflaute/data/regions_grid5"
elif REGION_SET != "regions":
    raise ValueError(f"REGION_SET={REGION_SET}")
REGIONS = ["DE"] + list(BOXES)          # land regions, in output column order
MONTHS = (9, 10, 11, 12, 1, 2, 3)
KAPPA, G, NU, ALPHA_CH = 0.4, 9.81, 1.5e-5, 0.018


def grib2_strip_local(b):
    """Drop the local-use section 2 from a GRIB2 message. The EERIE daily and monthly statistics
    messages carry a 52-byte section 2 where eccodes (2.43 and 2.46) expects 57 bytes for that
    local definition, which breaks every later key; section 2 is optional in GRIB2."""
    import struct
    out, p = [b[:16]], 16
    while p < len(b) - 4:
        n = struct.unpack(">I", b[p:p + 4])[0]
        if b[p + 4] != 2:
            out.append(b[p:p + n])
        p += n
    out.append(b"7777")
    m = b"".join(out)
    return m[:8] + struct.pack(">Q", len(m)) + m[16:]


def lon180(lon):
    lon = np.asarray(lon, dtype=float)
    return np.where(lon > 180, lon - 360, lon)


def _polygon(path):
    import shapely.geometry as sgeom
    import shapely.prepared
    return shapely.prepared.prep(sgeom.shape(json.load(open(path))["features"][0]["geometry"]))


def region_cells(lat, lon, land):
    """Flat cell lists (lat, lon, land of equal length) -> {region: indices}; also 'DE_off'."""
    import shapely.geometry as sgeom
    lat = np.asarray(lat, float)
    lon = lon180(lon)
    land = np.asarray(land, bool)
    out = {}
    la0, la1, lo0, lo1 = DE_BOX
    cand = np.where((lat >= la0) & (lat <= la1) & (lon >= lo0) & (lon <= lo1))[0]
    pl, pe = _polygon(DE_LAND), _polygon(DE_EEZ)
    on = [i for i in cand if land[i] and pl.contains(sgeom.Point(lon[i], lat[i]))]
    off = [i for i in cand if (not land[i]) and pe.contains(sgeom.Point(lon[i], lat[i]))]
    out["DE"] = np.array(on, dtype=np.int64)
    for k, (la0, la1, lo0, lo1) in BOXES.items():
        out[k] = np.where(land & (lat >= la0) & (lat <= la1) & (lon >= lo0) & (lon <= lo1))[0]
    out["DE_off"] = np.array(off, dtype=np.int64)
    return out


def _era5_static():
    with nc.Dataset(ERA5_FSR) as d:
        elat = d["latitude"][:].astype(float)
        elon = d["longitude"][:].astype(float)
        tv = d["valid_time"]
        yrs = np.array([t.year for t in nc.num2date(tv[:], tv.units)])
        b = np.where((elat >= LAT_BAND[0] - 2) & (elat <= LAT_BAND[1] + 2))[0]
        fsr = d["fsr"][yrs >= 2016][:, b[0]:b[-1] + 1].astype(float)
    with nc.Dataset(ERA5_LSM) as d:
        elsm = d["lsm"][:].squeeze()[b[0]:b[-1] + 1].astype(float)
    clim = np.log(np.maximum(fsr.reshape(-1, 12, *fsr.shape[1:]).mean(axis=0), 1e-5))   # (12, ny, nx)
    return elat[b[0]:b[-1] + 1], elon, elsm, clim


def era5_z0_on(lat, lon, dlat, dlon):
    """(12, n) z0 for model land cells at (lat, lon): geometric mean of the ERA5 fsr climatology
    (mean over 2016-2025 per calendar month) over the ERA5 land cells within half a model cell,
    else the nearest ERA5 land cell. dlat/dlon: model cell size in degrees (0 = ERA5 grid itself)."""
    elat, elon, elsm, clim = _era5_static()
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float) % 360.0
    out = np.empty((12, len(lat)))
    eland = np.argwhere(elsm >= 0.5)
    for k, (la, lo) in enumerate(zip(lat, lon)):
        mi = np.abs(elat - la) <= max(dlat / 2, 0.125)
        dl = np.abs(((elon - lo) + 180.0) % 360.0 - 180.0)
        mj = dl <= max(dlon / 2, 0.125)
        sel = elsm[np.ix_(mi, mj)] >= 0.5
        if sel.any():
            out[:, k] = np.exp(clim[:, mi][:, :, mj][:, sel].mean(axis=1))
        else:
            dd = (elat[eland[:, 0]] - la) ** 2 + (((elon[eland[:, 1]] - lo) + 180.0) % 360.0 - 180.0) ** 2
            i, j = eland[np.argmin(dd)]
            out[:, k] = np.exp(clim[:, i, j])
    return out


def loglaw_land(spd10, z0):
    return spd10 * np.log(100.0 / z0) / np.log(10.0 / z0)


def loglaw_sea(spd10):
    """Charnock, 8 fixed-point iterations from the step's own 10 m speed (as reconstruct_month.py)."""
    s = np.maximum(spd10, 0.5)
    z0 = np.full(spd10.shape, 2.0e-4)
    for _ in range(8):
        ustar = KAPPA * s / np.log(10.0 / z0)
        z0 = ALPHA_CH * ustar ** 2 / G + 0.11 * NU / ustar
    return spd10 * np.log(100.0 / z0) / np.log(10.0 / z0)


class QM:
    """ERA5 quantile-mapping calibration v1 (fitted over Germany), month x 3-hourly slot x class."""

    def __init__(self):
        t = np.load(QM_TABLE)
        self.qx, self.qy = t["qx"], t["qy"]
        self.months = list(t["months"])
        self.slots, self.n_class = int(t["slots"]), int(t["n_class"])
        self.edges = t["fsr_edges"]

    def land_class(self, z0):
        return 1 + np.digitize(z0, self.edges)

    def __call__(self, v, slot, month, cls):
        """v (T, N); slot, month (T,); cls (T, N) or (N,)."""
        m_idx = np.array([self.months.index(m) for m in month])
        b = (m_idx * self.slots + slot)[:, None] * self.n_class + (cls[None, :] if np.ndim(cls) == 1 else cls)
        out = v.copy()
        for k in np.unique(b):
            qx, qy = self.qx[k], self.qy[k]
            if np.isnan(qx[0]):
                continue
            sel = b == k
            x = v[sel]
            y = np.interp(x, qx, qy)
            y = np.where(x < qx[0], x + (qy[0] - qx[0]), y)
            y = np.where(x > qx[-1], x + (qy[-1] - qx[-1]), y)
            out[sel] = np.maximum(y, 0.0)
        return out


def month_of(times):
    """A stamp belongs to the month of (stamp - 1 s)."""
    return ((times - np.timedelta64(1, "s")).astype("datetime64[M]").astype(int) % 12) + 1


def region_means(x, cells, w):
    """x (T, N) over the flat cell list -> (T, R) weighted means for REGIONS (+ DE_off if asked)."""
    return np.stack([(x[:, cells[r]] * w[cells[r]]).sum(axis=1) / w[cells[r]].sum() for r in cells], axis=1)


def save(name, series, cells, meta):
    os.makedirs(OUT_DIR, exist_ok=True)
    fn = f"{OUT_DIR}/{name}.npz"
    np.savez(fn, **series, regions=np.array(list(cells)), n_cells=np.array([len(cells[r]) for r in cells]),
             meta=json.dumps(meta))
    return fn
