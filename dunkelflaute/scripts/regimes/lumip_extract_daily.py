"""
LUMIP land-use test, daily wind. Where day/sfcWind exists in both experiments of a pair, the region
series of extract_daily.py on the same footing: daily-mean 10 m speed, log law to 100 m with the ERA5
fsr land z0 climatology (fixed, so the land-use roughness change is NOT in the 10 -> 100 m step),
onshore power curve, cos-lat weighted region means over the ERA5-land-fraction >= 0.5 cells.
Grid statics (land mask, z0) are computed with extract_daily.setup_grid into lumip_parts/ (tag
lumip_<model>), so nothing in data/cmip6 is touched.

Pairs and winters (winter = year of its January):
  historical vs hist-noLu          1985-2014 (IPSL: hist-noLu daily ends 2009, so 1985-2009)
  ssp370 vs ssp370-ssp126Lu        2071-2100  (UKESM only)
  ssp126-ssp370Lu vs ssp126        2071-2100  (UKESM only)
Members restricted to the i/p/f variant of the LUMIP members.

Usage: python3 lumip_extract_daily.py        (legacy python; N_WORKERS, default 32)
Output: dunkelflaute/data/regimes/lumip_parts/lumip_day_<model>_<exp>.npz
        <member>|ymd (D, 3), <member>|wind (D, R) CF, <member>|spd10 (D, R) m/s, float32
"""
import glob
import json
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "1")
import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("", "/regions", "/cmip6"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts{p}")
import regions as rg            # noqa: E402
import extract_daily as xd      # noqa: E402

OUTD = f"{_ROOT}/dunkelflaute/data/regimes/lumip_parts"
xd.OUT = OUTD                   # statics go to lumip_parts/static_lumip_<model>.npz
INVL = f"{_ROOT}/dunkelflaute/results/regimes/lumip_inventory.json"
PAIRS = [("hist-noLu", "historical", (1985, 2014)), ("ssp370-ssp126Lu", "ssp370", (2071, 2100)),
         ("ssp126-ssp370Lu", "ssp126", (2071, 2100))]
_S = {}


def variant(m):
    return re.sub(r"^r\d+", "", m)


def span(fn):
    s = os.path.basename(fn).rsplit("_", 1)[-1].replace(".nc", "").split("-")
    return int(s[0][:4]), int(s[1][:4])


def plan():
    inv = json.load(open(INVL))
    out = []
    for model, ex in inv.items():
        for lu, ctl, (w0, w1) in PAIRS:
            a, b = ex.get(lu, {}).get("day/sfcWind"), ex.get(ctl, {}).get("day/sfcWind")
            if not a or not b:
                continue
            wl1 = min(w1, min(x["years"][1] for x in a.values()))       # IPSL hist-noLu ends 2009
            vset = {variant(m) for m in a}
            for exp, mems in ((lu, a), (ctl, b)):
                for m, x in sorted(mems.items()):
                    if variant(m) in vset:
                        out.append((model, exp, m, x["dir"], w0, wl1))
    return out


def setup(model, f):
    if model not in _S:
        with nc.Dataset(f) as d:
            _S[model] = xd.setup_grid(d["lat"][:], d["lon"][:], f"lumip_{model}")
    return _S[model]


def job(a):
    model, exp, mem, fn, y0, y1, w0, w1 = a
    S = setup(model, fn)
    with nc.Dataset(fn) as d:
        t = d["time"]
        dates = nc.num2date(t[:], t.units, calendar=getattr(t, "calendar", "standard"))
        ymd = np.array([(x.year, x.month, x.day) for x in dates])
        win = np.where(ymd[:, 1] >= 9, ymd[:, 0] + 1, ymd[:, 0])
        keep = np.isin(ymd[:, 1], rg.MONTHS) & (win >= max(w0, y0)) & (win <= min(w1, y1))
        if not keep.any():
            return model, exp, mem, np.zeros((0, 3), int), None
        idx = np.where(keep)[0]
        sl = slice(int(idx[0]), int(idx[-1]) + 1)
        x = d["sfcWind"][sl, S["r0"]:S["r1"], :]
        x = np.ma.filled(np.ma.masked_invalid(x).astype(np.float64), np.nan)
        x = x.reshape(x.shape[0], -1)[:, S["lu"]][keep[sl]]
    yk = ymd[sl][keep[sl]]
    r = xd.region_series(S, x, np.full_like(x, np.nan), yk[:, 1])
    return model, exp, mem, yk, dict(wind=r["wind"].astype(np.float32), spd10=r["spd10"].astype(np.float32))


def main():
    P = plan()
    jobs = []
    for model, exp, mem, d, w0, w1 in P:
        for f in sorted(glob.glob(d + "/*.nc")):
            a, b = span(f)
            for y in range(w0, w1 + 1, 5):                       # 5-winter chunks keep the reads small
                y_end = y + 4
                if b < y - 1 or a > y_end:
                    continue
                jobs.append((model, exp, mem, f, y, y_end, w0, w1))
    print(f"{len(P)} member-experiments, {len(jobs)} chunks", flush=True)
    # build the statics once per model before the pool forks
    for model in sorted({p[0] for p in P}):
        f = sorted(glob.glob([p for p in P if p[0] == model][0][3] + "/*.nc"))[0]
        setup(model, f)
        print(f"  statics {model}: " + ", ".join(f"{r} {len(_S[model]['cells'][r])}" for r in rg.REGIONS), flush=True)
    got = {}
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 32))) as ex:
        for model, exp, mem, ymd, r in ex.map(job, jobs):
            if r is None:
                continue
            g = got.setdefault((model, exp), {}).setdefault(mem, {"ymd": [], "wind": [], "spd10": []})
            g["ymd"].append(ymd)
            g["wind"].append(r["wind"])
            g["spd10"].append(r["spd10"])
    for (model, exp), mems in got.items():
        z = {}
        for m, g in mems.items():
            ymd = np.concatenate(g["ymd"])
            k = ymd @ np.array([10000, 100, 1])
            _, o = np.unique(k, return_index=True)
            z[f"{m}|ymd"] = ymd[o]
            z[f"{m}|wind"] = np.concatenate(g["wind"])[o]
            z[f"{m}|spd10"] = np.concatenate(g["spd10"])[o]
            print(f"  {model:14s} {exp:18s} {m:10s} {len(o)} days {tuple(ymd[o][0])} .. {tuple(ymd[o][-1])}", flush=True)
        np.savez_compressed(f"{OUTD}/lumip_day_{model}_{exp}.npz", **z, regions=np.array(rg.REGIONS), members=np.array(sorted(mems)))
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
