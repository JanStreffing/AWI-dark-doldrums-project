"""
CMIP6 inputs for the interannual transfer test, historical winters 1951-2014 (Sep 1950 - Mar 2014),
plus the regional winter temperature of 1951-1980 and ssp585 2071-2100 for the forced warming.

Per model, for the 8 European regions (regime_analysis.EU):
  s        daily region land mean day/sfcWind, cells and cos-lat weights exactly as extract_daily.py
           (cmip6_setup: ERA5 land fraction >= 0.5 within half a model cell)
  G        daily geostrophic speed from day/psl bilinearly interpolated to the global 2.5 deg grid
           (extract_psl._job, same psl directory choice), geostrophic.geostrophic / region_G;
           also the region-mean geostrophic components ug (eastward) and vg (northward)
  T        monthly region land mean Amon/tas (Sep-Mar), same land definition on the tas grid
Daily data only Sep-Mar.

Usage: python3 interannual_extract_cmip6.py MODEL      (legacy python; compute node, N_WORKERS)
Output: dunkelflaute/data/regimes/interannual_cmip6/<MODEL>.npz
"""
import glob
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "1")
import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("regions", "regimes", "cmip6"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/{p}")
import regions as rg                       # noqa: E402
import extract_daily as xd                 # noqa: E402
import extract_psl as xp                   # noqa: E402
import geostrophic as gs                   # noqa: E402

EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
OUT = f"{_ROOT}/dunkelflaute/data/regimes/interannual_cmip6"
W0, W1 = 1951, 2014
NW = int(os.environ.get("N_WORKERS", 32))
RHO, OMEGA, RE = 1.25, 7.292e-5, 6.371e6
BOX = dict(rg.BOXES, DE=rg.DE_BOX)


def geo_components(p, lat, lon):
    """p (D, ny, nx) Pa on the periodic global grid -> ug, vg (D, ny, nx); NaN at the lat edges."""
    phi = np.deg2rad(lat)[None, :, None]
    dphi, dlam = np.deg2rad(lat[1] - lat[0]), np.deg2rad(lon[1] - lon[0])
    dpy = np.full_like(p, np.nan)
    dpy[:, 1:-1] = (p[:, 2:] - p[:, :-2]) / (2 * dphi * RE)
    dpx = (np.roll(p, -1, axis=2) - np.roll(p, 1, axis=2)) / (2 * dlam * RE * np.cos(phi))
    f = 2 * OMEGA * np.sin(phi)
    return -dpy / (RHO * f), dpx / (RHO * f)


def region_mean(x, lat, lon):
    """Same point selection as geostrophic.region_G."""
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    out = []
    for r in EU:
        la0, la1, lo0, lo1 = BOX[r]
        sel = (LA >= la0) & (LA <= la1) & (LO >= lo0) & (LO <= lo1)
        if not sel.any():
            sel = (LA >= la0 - 1.25) & (LA <= la1 + 1.25) & (LO >= lo0 - 1.25) & (LO <= lo1 + 1.25)
        out.append(np.nanmean(x[:, sel], axis=1))
    return np.stack(out, axis=1)


def _psl_job(args):
    ymd, p = xp._job(args)
    if len(ymd) == 0:
        return ymd, np.zeros((0, 3, len(EU)))
    p = p.astype(np.float64)
    G = gs.region_G(gs.geostrophic(p, xp.LAT, xp.LON, True), xp.LAT, xp.LON, EU)
    ug, vg = geo_components(p, xp.LAT, xp.LON)
    return ymd, np.stack([G, region_mean(ug, xp.LAT, xp.LON), region_mean(vg, xp.LAT, xp.LON)], axis=1)


def psl_series(inv):
    base = inv["paths"]["historical"]["day/sfcWind"].split("/day/")[0]
    d = sorted(x for x in glob.glob(f"{base}/day/psl/*/*") if glob.glob(f"{x}/*.nc"))[-1]
    jobs = []
    for f in sorted(glob.glob(f"{d}/*.nc")):
        a, b = xp._span(f)
        for y in range(max(a, W0 - 1), min(b, W1) + 1):
            jobs.append((f, y, W0, W1))
    got = {}
    with ProcessPoolExecutor(max_workers=NW) as ex:
        for ymd, x in ex.map(_psl_job, jobs):
            for k in range(len(ymd)):
                got[tuple(ymd[k])] = x[k]
    keys = sorted(got)
    return np.array(keys), np.stack([got[k] for k in keys]), d


def wind_series(model, inv):
    jobs = []
    for f in sorted(glob.glob(inv["paths"]["historical"]["day/sfcWind"] + "/*.nc")):
        a, b = xd._span(f)
        if not (b < W0 - 1 or a > W1):
            jobs.append((model, "sfcWind", f, W0 - 1, W1))
    S = xd.cmip6_setup(model)
    got = {}
    with ProcessPoolExecutor(max_workers=NW, initializer=xd._cmip6_init, initargs=(model,)) as ex:
        for var, ymd, x in ex.map(xd._cmip6_file, jobs):
            for k in range(len(ymd)):
                got[tuple(ymd[k])] = x[k]
    keys = [k for k in sorted(got) if W0 <= (k[0] + 1 if k[1] >= 9 else k[0]) <= W1]
    u = np.stack([got[k] for k in keys])
    w = S["w"]
    s = np.stack([(u[:, i] * w[i]).sum(1) / w[i].sum() for i in (S["lpos"][rg.REGIONS.index(r)] for r in EU)], axis=1)
    return np.array(keys), s, S


def tas_cells(lat, lon, S_wind, lat_w, lon_w):
    """Land cells of the tas grid: the wind grid's if identical, else the extract_daily recipe (no cache)."""
    if len(lat) == len(lat_w) and len(lon) == len(lon_w) and np.allclose(lat, lat_w) and np.allclose(lon, lon_w):
        return S_wind["r0"], S_wind["r1"], [S_wind["lu"][S_wind["lpos"][rg.REGIONS.index(r)]] for r in EU], "wind grid"
    rows = np.where((lat >= rg.LAT_BAND[0] - 1) & (lat <= rg.LAT_BAND[1] + 1))[0]
    r0, r1 = int(rows.min()), int(rows.max()) + 1
    flat, flon = np.repeat(lat[r0:r1], len(lon)), np.tile(lon, r1 - r0)
    dlat, dlon = float(np.median(np.abs(np.diff(lat)))), float(np.median(np.abs(np.diff(lon))))
    cand = np.zeros(len(flat), bool)
    L = rg.lon180(flon)
    for la0, la1, lo0, lo1 in [BOX[r] for r in EU]:
        cand |= (flat >= la0 - dlat) & (flat <= la1 + dlat) & (L >= lo0 - dlon) & (L <= lo1 + dlon)
    frac = np.zeros(len(flat))
    frac[cand] = xd.land_fraction(flat[cand], flon[cand], dlat, dlon)
    cells = rg.region_cells(flat, flon, frac >= 0.5)
    return r0, r1, [cells[r] for r in EU], "own land mask"


def tas_series(inv, exp, y0, y1, S_wind, lat_w, lon_w):
    """Monthly region land-mean tas for Sep-Mar months of the winters y0..y1 -> (years, months, R)."""
    out, how = {}, None
    cells = None
    for f in sorted(glob.glob(inv["paths"][exp]["Amon/tas"] + "/*.nc")):
        a, b = xd._span(f)
        if b < y0 - 1 or a > y1:
            continue
        with nc.Dataset(f) as d:
            lat, lon = d["lat"][:].astype(float), d["lon"][:].astype(float)
            if cells is None:
                r0, r1, cells, how = tas_cells(lat, lon, S_wind, lat_w, lon_w)
                wl = np.cos(np.deg2rad(np.repeat(lat[r0:r1], len(lon))))
            t = d["time"]
            ds = nc.num2date(t[:], t.units, calendar=getattr(t, "calendar", "standard"))
            ym = np.array([(x.year, x.month) for x in ds])
            win = np.where(ym[:, 1] >= 9, ym[:, 0] + 1, ym[:, 0])
            keep = np.isin(ym[:, 1], rg.MONTHS) & (win >= y0) & (win <= y1)
            if not keep.any():
                continue
            x = d["tas"][np.where(keep)[0], r0:r1, :].reshape(int(keep.sum()), -1).astype(np.float64)
            for k, (y, m) in enumerate(ym[keep]):
                out[(int(y), int(m))] = [float((x[k, c] * wl[c]).sum() / wl[c].sum()) for c in cells]
    keys = sorted(out)
    return np.array(keys), np.array([out[k] for k in keys]), how


def main():
    model = sys.argv[1]
    inv = json.load(open(xd.INV))[model]
    os.makedirs(OUT, exist_ok=True)
    kw, s, S = wind_series(model, inv)
    print(f"  {model} sfcWind: {len(kw)} days {tuple(kw[0])} .. {tuple(kw[-1])}", flush=True)
    kp, g, pdir = psl_series(inv)
    print(f"  {model} psl: {len(kp)} days {tuple(kp[0])} .. {tuple(kp[-1])} ({pdir})", flush=True)
    f0 = sorted(glob.glob(inv["paths"]["ssp585"]["day/sfcWind"] + "/*.nc"))[0]
    with nc.Dataset(f0) as d:
        lat_w, lon_w = d["lat"][:].astype(float), d["lon"][:].astype(float)
    res = {}
    for exp, (y0, y1) in (("historical", (W0, W1)), ("ssp585", (2071, 2100))):
        ym, t, how = tas_series(inv, exp, y0, y1, S, lat_w, lon_w)
        res[f"{exp}|tas_ym"], res[f"{exp}|tas"] = ym, t
        print(f"  {model} tas {exp}: {len(ym)} months {tuple(ym[0])} .. {tuple(ym[-1])} ({how})", flush=True)
    np.savez(f"{OUT}/{model}.npz", wind_ymd=kw, s=s, psl_ymd=kp, G=g[:, 0], ug=g[:, 1], vg=g[:, 2], regions=np.array(EU),
             psl_dir=pdir, **res)
    print(f"Saved {OUT}/{model}.npz", flush=True)


if __name__ == "__main__":
    main()
