"""
HadISD v3.4.3.2025f station files -> daily mean 10 m wind speed and 2 m temperature, Sep-Mar of the
winters 1979-2025, for the stations of interannual_hadisd_download.py.

Daily mean = mean of the 00, 06, 12 and 18 UTC observations, valid only if all four are present
(HadISD values are already QC'd; flagged values are set to the missing value). Using the four main
synoptic hours keeps the diurnal sampling the same whether a station reports hourly, 3-hourly or
6-hourly, and through changes of reporting frequency. A second wind series uses every 3 h stamp
(00, 03, ..., 21 UTC), valid with at least 6 of the 8.

The gzipped netCDF files are opened in memory; nothing in them is executed.

Usage: python3 -I interannual_hadisd_process.py      (login node, a few minutes)
Output: dunkelflaute/data/obs/hadisd_daily.npz
  day (D,) datetime64[D], ids (N,), lat/lon/elev (N,), ws/ws3h/ta (N, D) float32 (NaN = invalid)
"""
import gzip
import json
import os
from concurrent.futures import ProcessPoolExecutor

import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
RAW = f"{_ROOT}/dunkelflaute/data/obs/hadisd_raw"
SEL = f"{_ROOT}/dunkelflaute/data/obs/hadisd_selection.json"
OUT = f"{_ROOT}/dunkelflaute/data/obs/hadisd_daily.npz"
MONTHS = (9, 10, 11, 12, 1, 2, 3)
T0 = np.datetime64("1978-09-01")
T1 = np.datetime64("2025-04-01")
DAYS = np.arange(T0, T1, dtype="datetime64[D]")
DAYS = DAYS[np.isin(DAYS.astype("datetime64[M]").astype(int) % 12 + 1, MONTHS)]
BASE = np.datetime64("1931-01-01T00", "h")


def _station(args):
    sid, fn = args
    with gzip.open(f"{RAW}/{fn}") as fh:
        b = fh.read()
    with nc.Dataset("mem", memory=b) as d:
        t = BASE + d["time"][:].astype(np.int64).astype("timedelta64[h]")
        keep = (t >= T0.astype("datetime64[h]")) & (t < T1.astype("datetime64[h]"))
        t = t[keep]
        ws = np.ma.filled(d["windspeeds"][:][keep].astype(np.float64), np.nan)
        ta = np.ma.filled(d["temperatures"][:][keep].astype(np.float64), np.nan)
    ws[ws < -100] = np.nan
    ta[ta < -100] = np.nan
    day = t.astype("datetime64[D]")
    hr = (t - day.astype("datetime64[h]")).astype(int)
    di = np.searchsorted(DAYS, day)
    ok = (di < len(DAYS)) & (DAYS[np.minimum(di, len(DAYS) - 1)] == day)
    out = {}
    for name, x, hours, need in (("ws", ws, (0, 6, 12, 18), 4), ("ws3h", ws, tuple(range(0, 24, 3)), 6), ("ta", ta, (0, 6, 12, 18), 4)):
        m = ok & np.isin(hr, hours) & np.isfinite(x)
        S = np.bincount(di[m], weights=x[m], minlength=len(DAYS))
        N = np.bincount(di[m], minlength=len(DAYS))
        v = np.where(N >= need, S / np.maximum(N, 1), np.nan)
        out[name] = v.astype(np.float32)
    return sid, out


def main():
    sel = json.load(open(SEL))
    ids = sorted(sel["files"])
    res = {}
    with ProcessPoolExecutor(8) as ex:
        for sid, o in ex.map(_station, [(i, sel["files"][i]) for i in ids]):
            res[sid] = o
    info = sel["info"]
    np.savez(OUT, day=DAYS, ids=np.array(ids), lat=np.array([info[i][0] for i in ids]), lon=np.array([info[i][1] for i in ids]),
             elev=np.array([info[i][2] for i in ids]), **{k: np.stack([res[i][k] for i in ids]) for k in ("ws", "ws3h", "ta")})
    fr = {k: float(np.mean([np.isfinite(res[i][k]).mean() for i in ids])) for k in ("ws", "ws3h", "ta")}
    print(f"Saved {OUT}: {len(ids)} stations, {len(DAYS)} days; mean valid fraction {fr}")


if __name__ == "__main__":
    main()
