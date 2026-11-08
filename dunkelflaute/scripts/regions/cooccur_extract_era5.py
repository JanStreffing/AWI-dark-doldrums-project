"""
ERA5 onshore wind CF for the co-occurrence and aggregation-scale validation (paper section on
Europe-wide low-wind spells): the same recipe as the model region series, on both footings.

  source  /pool/data/ERA5/E5/sf/an/1H/165 and 166 (10u, 10v; daily files of 24 hourly GRIB messages,
          reduced Gaussian N320). Only the 00, 06, 12, 18 UTC messages are decoded.
  land    ERA5 lsm (regions.ERA5_LSM) at the nearest 0.25 deg point >= 0.5 (as
          interannual_extract_era5.py); z0 = ERA5 fsr climatology 2016-2025 of the calendar month
          (regions.era5_z0_on with dlat = dlon = 0, i.e. the ERA5 land cells within 0.125 deg)
  wind    log law to 100 m, onshore curve (core.capacity_factor.ONSHORE), cos-lat weighted means
          6h    CF of the instantaneous speed at 00/06/12/18 UTC (as the 6-hourly model footing)
          day   CF of the mean speed of the four instants of the calendar day (as the CMIP6 daily
                footing, which uses the daily-mean sfcWind)
  regions the 8 European named regions (DE polygon, FR, IB, UK, PL, SC, IT, BK boxes of regions.py)
          and the 420 5 x 5 deg land boxes of REGION_SET=grid5 (same cells rule as the models)
  period  Sep-Mar, 1979-09-01 to 2026-03-31 (winters 1980-2026)

Usage: python3 cooccur_extract_era5.py      (legacy python with eccodes; compute node, N_WORKERS)
Output: dunkelflaute/data/regions/cooccur_era5_wind.npz
  time6 (T,) datetime64[s], wind6 (T, R); day (D,), wind_day (D, R); regions (R,) = EU8 + grid boxes
"""
import os
import sys
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ["REGION_SET"] = "grid5"
import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
from core.capacity_factor import ONSHORE, wind_cf   # noqa: E402
import regions as rg                                # noqa: E402

EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
SRC = "/pool/data/ERA5/E5/sf/an/1H/{v}/E5sf00_1H_{d}_{v}.grb"
OUT = f"{_ROOT}/dunkelflaute/data/regions/cooccur_era5_wind.npz"
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


def build(lat, lon):
    """Region weights on the N320 grid -> (names, used cell index (n,), z0 (12, n), W (R, n) rows sum 1)."""
    lon = rg.lon180(lon)
    with nc.Dataset(rg.ERA5_LSM) as d:
        elat = d["latitude"][:].astype(float)
        elon = rg.lon180(d["longitude"][:].astype(float))
        lsm = d["lsm"][:].squeeze().astype(float)
    i = np.clip(np.rint((elat[0] - lat) / abs(elat[1] - elat[0])).astype(int), 0, len(elat) - 1)
    oj = np.argsort(elon)
    j = oj[np.clip(np.searchsorted(elon[oj], lon), 0, len(elon) - 1)]
    land = lsm[i, j] >= 0.5
    cells = rg.region_cells(lat, lon, land)          # DE polygon + grid5 boxes (+ DE_off)
    named = {"DE": cells["DE"]}
    for r in EU[1:]:
        a, b, c, e = rg.REGION_BOXES[r]
        named[r] = np.where(land & (lat >= a) & (lat <= b) & (lon >= c) & (lon <= e))[0]
    grid = [r for r in rg.REGIONS if r != "DE"]
    names = EU + grid
    lists = [named[r] for r in EU] + [cells[r] for r in grid]
    used = np.unique(np.concatenate([x for x in lists if len(x)]))
    w = np.cos(np.deg2rad(lat))
    W = np.zeros((len(names), len(used)), np.float32)
    for k, x in enumerate(lists):
        if len(x):
            p = np.searchsorted(used, x)
            W[k, p] = w[x] / w[x].sum()
    z0 = rg.era5_z0_on(lat[used], lon[used], 0.0, 0.0)
    return names, used, z0, W


def init():
    import eccodes
    ms = messages(SRC.format(v="165", d="2000-01-01"), {0})
    h = eccodes.codes_new_from_message(ms[0])
    try:
        lat = eccodes.codes_get_array(h, "latitudes")
        lon = eccodes.codes_get_array(h, "longitudes")
    finally:
        eccodes.codes_release(h)
    names, used, z0, W = build(lat, lon)
    _S.update(names=names, used=used, z0=z0, W=W)


def _day(d):
    import eccodes
    if "W" not in _S:
        init()
    vals = {}
    for v in ("165", "166"):
        ms = messages(SRC.format(v=v, d=d), set(HOURS))
        xs = []
        for kk in HOURS:
            h = eccodes.codes_new_from_message(ms[kk])
            try:
                assert eccodes.codes_get(h, "dataTime") == kk * 100, (d, v, kk)
                xs.append(eccodes.codes_get_values(h)[_S["used"]])
            finally:
                eccodes.codes_release(h)
        vals[v] = np.array(xs)
    spd = np.hypot(vals["165"], vals["166"])                       # (4, n)
    m = int(d[5:7])
    z0 = _S["z0"][m - 1]
    cf6 = wind_cf(rg.loglaw_land(spd, z0), **ONSHORE).astype(np.float64) @ _S["W"].T          # (4, R)
    cfd = wind_cf(rg.loglaw_land(spd.mean(axis=0), z0), **ONSHORE).astype(np.float64) @ _S["W"].T
    return d, cf6.astype(np.float32), cfd.astype(np.float32)


def main():
    days = np.arange(np.datetime64("1979-09-01"), np.datetime64("2026-04-01"))
    days = days[np.isin(days.astype("datetime64[M]").astype(int) % 12 + 1, rg.MONTHS)]
    init()
    print(f"{len(_S['names'])} regions, {len(_S['used'])} land cells; cells per EU region: "
          + ", ".join(f"{r} {int((_S['W'][k] > 0).sum())}" for k, r in enumerate(_S['names'][:8])), flush=True)
    T6, X6, XD = [], [], []
    with ProcessPoolExecutor(int(os.environ.get("N_WORKERS", 64))) as ex:          # fork: workers inherit _S
        for n, (d, cf6, cfd) in enumerate(ex.map(_day, [str(x) for x in days], chunksize=8)):
            T6.extend(np.datetime64(d) + np.array(HOURS) * np.timedelta64(1, "h"))
            X6.append(cf6)
            XD.append(cfd)
            if n % 1000 == 0:
                print(f"  {d}", flush=True)
    np.savez(OUT, time6=np.array(T6).astype("datetime64[s]"), wind6=np.concatenate(X6), day=days,
             wind_day=np.array(XD), regions=np.array(_S["names"]))
    print(f"Saved {OUT}: {len(days)} days")


if __name__ == "__main__":
    main()
