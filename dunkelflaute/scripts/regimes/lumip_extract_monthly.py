"""
LUMIP land-use test, monthly extraction. Per model, experiment and member: winter (Sep-Mar) monthly
region means over the region land cells of stability_cmip6_extract.Grid (ERA5 land fraction >= 0.5
within half a model cell, cos-lat weights; the cached land mask of data/cmip6/static_cmip6_<model>.npz
is reused read-only where it matches the grid), and the monthly geostrophic speed.

  sfcWind   Amon/sfcWind (monthly mean of the 10 m speed)
  tas       Amon/tas
  G         |grad p| / (rho f) of the MONTHLY-mean psl, bilinear to the global 2.5 deg 15-80 N grid of
            extract_psl.py, region boxes as geostrophic.region_G (so not the mean of daily G)
  cropFrac, treeFrac, grassFrac, pastureFrac, lai   Lmon (winter means); Eyr/cropFrac and
            Eyr/treeFrac as an annual fallback (value of the year of the January)
Members are restricted to the i/p/f variant of the model's LUMIP members (CanESM5 LUMIP runs are p2,
UKESM f2), so control and LUMIP experiments share physics and forcing variant.
Winters: historical-type 1951-2014, scenarios 2016-2100 (winter = year of its January).

Usage: python3 lumip_extract_monthly.py [MODEL ...]      (legacy python; N_WORKERS, default 32)
Output: dunkelflaute/data/regimes/lumip_parts/lumip_mon_<model>_<exp>.npz
        <member>|<var> (W, 7, R) float32 (land cover: (W, R)), winters (W,), regions
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
from scipy.interpolate import RegularGridInterpolator

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("", "/regions", "/cmip6", "/regimes", "/analysis"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts{p}")
import regions as rg                       # noqa: E402
import stability_cmip6_extract as sx       # noqa: E402

INVL = f"{_ROOT}/dunkelflaute/results/regimes/lumip_inventory.json"
OUTD = f"{_ROOT}/dunkelflaute/data/regimes/lumip_parts"
MON = [9, 10, 11, 12, 1, 2, 3]
WIN = {"hist": (1951, 2014), "scen": (2016, 2100)}
EXPKIND = {"historical": "hist", "hist-noLu": "hist", "ssp370": "scen", "ssp370-ssp126Lu": "scen", "ssp126": "scen",
           "ssp126-ssp370Lu": "scen", "esm-ssp585": "scen", "esm-ssp585-ssp126Lu": "scen"}
LUEXP = ("hist-noLu", "ssp370-ssp126Lu", "ssp126-ssp370Lu", "esm-ssp585-ssp126Lu")
CTL = {"hist-noLu": "historical", "ssp370-ssp126Lu": "ssp370", "ssp126-ssp370Lu": "ssp126", "esm-ssp585-ssp126Lu": "esm-ssp585"}
MONVARS = ("Amon/sfcWind", "Amon/tas")
LCVARS = ("Lmon/cropFrac", "Lmon/treeFrac", "Lmon/grassFrac", "Lmon/pastureFrac", "Lmon/lai")
EYR = ("Eyr/cropFrac", "Eyr/treeFrac")
LAT = np.arange(15.0, 80.01, 2.5)
LON = np.arange(-180.0, 180.0, 2.5)
RHO, OMEGA, RE = 1.25, 7.292e-5, 6.371e6
BOX = dict(rg.BOXES, DE=rg.DE_BOX)
_GRIDS = {}


def variant(m):
    return re.sub(r"^r\d+", "", m)


def tasks():
    """(model, exp, member, {var: dir}) for every experiment of a pair with a wind field."""
    inv = json.load(open(INVL))
    out = []
    for model, ex in inv.items():
        for lu in LUEXP:
            if lu not in ex or CTL[lu] not in ex:
                continue
            if "Amon/sfcWind" not in ex[lu] or "Amon/sfcWind" not in ex[CTL[lu]]:
                continue
            vars_ = set(ex[lu]["Amon/sfcWind"])
            vset = {variant(m) for m in vars_}
            for exp in (lu, CTL[lu]):
                mems = sorted(m for m in ex[exp]["Amon/sfcWind"] if variant(m) in vset)
                for m in mems:
                    d = {v: ex[exp][v][m]["dir"] for v in ex[exp] if m in ex[exp][v]}
                    out.append((model, exp, m, d))
    seen, uniq = set(), []
    for t in out:
        if t[:3] not in seen:
            seen.add(t[:3])
            uniq.append(t)
    return uniq


def grid(model, lat, lon):
    key = (model, len(lat), len(lon), float(lat[0]), float(lon[0]))
    if key not in _GRIDS:
        _GRIDS[key] = sx.Grid(model, lat, lon)
    return _GRIDS[key]


def read_monthly(dd, var, model, w0, w1, annual=False):
    """-> (W, 7, R) winter-month region means (annual: (W, R)), or None."""
    fs = [f for f in sorted(glob.glob(dd + "/*.nc")) if not (sx.span(f)[1] < w0 - 1 or sx.span(f)[0] > w1)]
    if not fs:
        return None
    W = w1 - w0 + 1
    arr = None
    for f in fs:
        with nc.Dataset(f) as d:
            d.set_auto_mask(True)
            lat, lon = d["lat"][:], d["lon"][:]
            G = grid(model, lat, lon)
            t = d["time"]
            dates = nc.num2date(t[:], t.units, calendar=getattr(t, "calendar", "standard"))
            ym = np.array([(x.year, x.month) for x in dates])
            if annual:
                win = ym[:, 0]
                keep = (win >= w0) & (win <= w1)
            else:
                win = np.where(ym[:, 1] >= 9, ym[:, 0] + 1, ym[:, 0])
                keep = np.isin(ym[:, 1], MON) & (win >= w0) & (win <= w1)
            if not keep.any():
                continue
            idx = np.where(keep)[0]
            sl = slice(int(idx[0]), int(idx[-1]) + 1)
            x = d[var][sl, G.r0:G.r1, :]
            x = np.ma.filled(np.ma.masked_invalid(x).astype(np.float64), np.nan)
            x = x.reshape(x.shape[0], -1)[:, G.lu]
            x[np.abs(x) > 1e15] = np.nan
            r = G.rmean(x)                                              # (T, R)
            if arr is None:
                arr = np.full((W, R_N) if annual else (W, 7, R_N), np.nan)
            kk = keep[sl]
            for ri, (yy, mm), ww in zip(r[kk], ym[sl][kk], win[sl][kk]):
                if annual:
                    arr[ww - w0] = ri
                else:
                    arr[ww - w0, MON.index(mm)] = ri
    return arr


R_N = len(rg.REGIONS)


def geostrophic(p, lat, lon):
    phi = np.deg2rad(lat)[None, :, None]
    dphi = np.deg2rad(lat[1] - lat[0])
    dlam = np.deg2rad(lon[1] - lon[0])
    dpy = np.full_like(p, np.nan)
    dpy[:, 1:-1] = (p[:, 2:] - p[:, :-2]) / (2 * dphi * RE)
    dpx = (np.roll(p, -1, axis=2) - np.roll(p, 1, axis=2)) / (2 * dlam * RE * np.cos(phi))
    return np.hypot(dpx, dpy) / (RHO * 2 * OMEGA * np.sin(phi))


def region_G(vg, regions):
    LA, LO = np.meshgrid(LAT, LON, indexing="ij")
    out = []
    for r in regions:
        la0, la1, lo0, lo1 = BOX[r]
        sel = (LA >= la0) & (LA <= la1) & (LO >= lo0) & (LO <= lo1)
        if not sel.any():
            sel = (LA >= la0 - 1.25) & (LA <= la1 + 1.25) & (LO >= lo0 - 1.25) & (LO <= lo1 + 1.25)
        out.append(np.nanmean(vg[:, sel], axis=1))
    return np.stack(out, axis=1)


def read_G(dd, w0, w1):
    """monthly psl -> (W, 7, R) region G of the monthly-mean field (geostrophic.py recipe, global 2.5 deg)."""
    fs = [f for f in sorted(glob.glob(dd + "/*.nc")) if not (sx.span(f)[1] < w0 - 1 or sx.span(f)[0] > w1)]
    if not fs:
        return None
    W = w1 - w0 + 1
    arr = np.full((W, 7, R_N), np.nan)
    LA, LO = np.meshgrid(LAT, LON % 360, indexing="ij")
    pts = np.column_stack([LA.ravel(), LO.ravel()])
    for f in fs:
        with nc.Dataset(f) as a:
            t = a["time"]
            dates = nc.num2date(t[:], t.units, calendar=getattr(t, "calendar", "standard"))
            ym = np.array([(x.year, x.month) for x in dates])
            win = np.where(ym[:, 1] >= 9, ym[:, 0] + 1, ym[:, 0])
            keep = np.isin(ym[:, 1], MON) & (win >= w0) & (win <= w1)
            if not keep.any():
                continue
            lat = a["lat"][:].astype(float)
            lon = a["lon"][:].astype(float) % 360
            rows = np.where((lat >= LAT[0] - 5) & (lat <= LAT[-1] + 5))[0]
            r0, r1 = int(rows.min()), int(rows.max()) + 1
            idx = np.where(keep)[0]
            sl = slice(int(idx[0]), int(idx[-1]) + 1)
            x = np.ma.filled(a["psl"][sl, r0:r1, :].astype(np.float64), np.nan)[keep[sl]]
            ymk, wk = ym[sl][keep[sl]], win[sl][keep[sl]]
        la = lat[r0:r1]
        if la[0] > la[-1]:
            la, x = la[::-1], x[:, ::-1]
        o = np.argsort(lon)
        lo, x = lon[o], x[:, :, o]
        lonp = np.concatenate([lo[-1:] - 360, lo, lo[:1] + 360])
        xp = np.concatenate([x[:, :, -1:], x, x[:, :, :1]], axis=2)
        out = np.empty((len(x), len(LAT), len(LON)))
        for k in range(len(x)):
            out[k] = RegularGridInterpolator((la, lonp), xp[k], bounds_error=True)(pts).reshape(len(LAT), len(LON))
        g = region_G(geostrophic(out, LAT, LON), rg.REGIONS)
        for gi, (yy, mm), ww in zip(g, ymk, wk):
            arr[ww - w0, MON.index(mm)] = gi
    return arr


def one(t):
    model, exp, mem, dirs = t
    t0 = time.time()
    w0, w1 = WIN[EXPKIND[exp]]
    out = {}
    for v in MONVARS + LCVARS:
        if v in dirs:
            a = read_monthly(dirs[v], v.split("/")[1], model, w0, w1)
            if a is not None:
                out[v.split("/")[1]] = a.astype(np.float32) if v in MONVARS else np.nanmean(a, axis=1).astype(np.float32)
    for v in EYR:
        name = v.split("/")[1]
        if v in dirs and name not in out:
            a = read_monthly(dirs[v], name, model, w0, w1, annual=True)
            if a is not None:
                out[name] = a.astype(np.float32)
                out[f"{name}_source"] = np.array("Eyr")
    if "Amon/psl" in dirs:
        a = read_G(dirs["Amon/psl"], w0, w1)
        if a is not None:
            out["G"] = a.astype(np.float32)
    return model, exp, mem, out, time.time() - t0


def main():
    want = set(sys.argv[1:])
    T = [t for t in tasks() if not want or t[0] in want]
    print(f"{len(T)} member-experiments", flush=True)
    os.makedirs(OUTD, exist_ok=True)
    res = {}
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 32))) as ex:
        for model, exp, mem, out, dt in ex.map(one, T):
            res.setdefault((model, exp), {})[mem] = out
            print(f"  {model:15s} {exp:20s} {mem:10s} {dt:5.0f}s {','.join(sorted(out))}  [{time.time() - t0:.0f}s]", flush=True)
    for (model, exp), mems in res.items():
        w0, w1 = WIN[EXPKIND[exp]]
        z = {f"{m}|{k}": v for m, o in mems.items() for k, v in o.items()}
        np.savez_compressed(f"{OUTD}/lumip_mon_{model}_{exp}.npz", **z, winters=np.arange(w0, w1 + 1),
                            regions=np.array(rg.REGIONS), members=np.array(sorted(mems)))
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
