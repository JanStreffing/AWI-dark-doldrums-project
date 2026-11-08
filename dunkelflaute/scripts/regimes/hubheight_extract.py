"""
Hub-height test of the CMIP6 transfer drop (methods notes 13e): is the forced decline of the 10 m
wind larger than the decline at ~100 m, as a near-surface roughness effect of prescribed land use
would make it?

Per CMIP6 model, historical winters 1951-1980 and ssp585 winters 2071-2100 (Sep-Mar, winter = year
of its January), on the region land cells of extract_daily.py (cached static_cmip6_<model>.npz,
read only; cos-lat weights), daily region means of:

  pressure-level footing (all models with day/ua or Eday/ua):
    sfc       day/sfcWind (daily mean speed)
    s10v      |daily-mean (uas, vas)|          (same footing as the |daily-mean (ua, va)| below)
    V1000     |daily-mean (ua, va)| at 1000 hPa (day/ua plev8; Eday plev19 if no day/ua)
    V925      at 925 hPa (Eday plev19 only)
    V850      at 850 hPa (day/ua plev8)
    h1000     height of 1000 hPa above ground, hypsometric: (Rd T / g) ln(ps / 1000 hPa) with
              daily ps estimated as day/psl x (Amon ps / Amon psl) of that month and cell, and
              T = Amon tas minus half a standard-atmosphere lapse over the layer
    zgh1000   day/zg(1000 hPa) - fx/orog, where both exist (cross-check of h1000)
  masks: a cell-day counts for "m1000" quantities if h1000 >= HMIN (50 m) and V1000 is valid; for
  "m925" if h925 >= HMIN and V925 valid. sfc and s10v are also given on the m1000 cell-days, so the
  ratios s10/V1000 are like-for-like samples.

  model-level footing (MODE ml; CFday ua, va, ps on hybrid levels, daily means):
    V100ml    |V| interpolated linearly in ln z to 100 m above ground; level heights
              z_k = (Rd T / g) ln(ps / p_k), p_k = ap_k + b_k ps (or a_k p0 + b_k ps; IPSL: mean of
              the bounding half levels), T as above
    Vlow, zlow  lowest model level speed and height

Each region-day value is the weighted mean over the valid cells of that day; frac_* is the valid
weight fraction, so cell-day means over a winter are sum(x frac) / sum(frac).

Usage: python3 hubheight_extract.py MODEL [pl|ml|both]      (legacy python, compute node)
Output: dunkelflaute/data/regimes/hubheight_parts/hubheight_<pl|ml>_<MODEL>.npz
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

INV = f"{_ROOT}/dunkelflaute/data/cmip6/inventory.json"
OUTD = f"{_ROOT}/dunkelflaute/data/regimes/hubheight_parts"
PERIODS = {"historical": (1951, 1980), "ssp585": (2071, 2100)}
if os.environ.get("HH_TEST"):                       # quick test: two winters per period
    PERIODS = {"historical": (1951, 1952), "ssp585": (2071, 2072)}
MON = (9, 10, 11, 12, 1, 2, 3)
RD, GRAV, LAPSE = 287.05, 9.80665, 0.0065
HMIN = float(os.environ.get("HMIN", 50.0))
ZHUB = 100.0
NW = int(os.environ.get("N_WORKERS", 32))
SUF = "_test" if os.environ.get("HH_TEST") else ""


def vdir(member, table, var):
    ds = [d for d in glob.glob(f"{member}/{table}/{var}/*/v*") if glob.glob(f"{d}/*.nc")]
    if not ds:
        return None
    return sorted(ds, key=lambda d: d.split("/")[-1])[-1]


def span(fn):
    s = re.findall(r"_(\d{4})\d*-(\d{4})\d*\.nc$", os.path.basename(fn))
    return (int(s[0][0]), int(s[0][1])) if s else (-9999, 9999)


def jobs_for(d, w0, w1):
    """(file, calendar year) pairs covering winters w0..w1 (Sep of w0-1 .. Mar of w1)."""
    out = []
    for f in sorted(glob.glob(f"{d}/*.nc")):
        a, b = span(f)
        for y in range(max(a, w0 - 1), min(b, w1) + 1):
            out.append((f, y))
    return out


def _clean(x):
    x = np.ma.filled(np.ma.masked_invalid(np.ma.asarray(x)).astype(np.float32), np.nan)
    x[np.abs(x) > 1e15] = np.nan
    return x


def _keep(t, year, w0, w1):
    dates = nc.num2date(t[:], t.units, calendar=getattr(t, "calendar", "standard"))
    ymd = np.array([(x.year, x.month, x.day) for x in dates])
    win = np.where(ymd[:, 1] >= 9, ymd[:, 0] + 1, ymd[:, 0])
    keep = (ymd[:, 0] == year) & np.isin(ymd[:, 1], MON) & (win >= w0) & (win <= w1)
    return ymd, keep


_G = None


def _to_lu(a, G, get):
    """Values on the region land cells of the target (day/sfcWind) grid. get(row_slice) -> (..., rows, lon)
    of the source grid. Same grid: direct indexing; otherwise (e.g. HadGEM3/UKESM staggered u, v points)
    bilinear interpolation, periodic in longitude."""
    slat = np.asarray(a["lat"][:], float)
    slon = np.asarray(a["lon"][:], float) % 360
    if len(slat) == G["nlat"] and len(slon) == G["nlon"] and np.allclose(slat, G["lat"]) and np.allclose(slon, G["lon"] % 360):
        x = get(slice(G["r0"], G["r1"]))
        return x.reshape(x.shape[:-2] + (-1,))[..., G["lu"]]
    tla, tlo = G["tlat"], G["tlon"] % 360
    asc = slat[1] > slat[0]
    rows = np.where((slat >= tla.min() - 5) & (slat <= tla.max() + 5))[0]
    r0, r1 = int(rows.min()), int(rows.max()) + 1
    x = get(slice(r0, r1))
    la = slat[r0:r1]
    if not asc:
        la, x = la[::-1], x[..., ::-1, :]
    i = np.clip(np.searchsorted(la, tla) - 1, 0, len(la) - 2)
    wy = np.clip((tla - la[i]) / (la[i + 1] - la[i]), 0, 1)
    o = np.argsort(slon)
    lo = slon[o]
    x = x[..., o]
    lop = np.concatenate([lo[-1:] - 360, lo, lo[:1] + 360])
    j = np.clip(np.searchsorted(lop, tlo) - 1, 0, len(lop) - 2)
    wx = (tlo - lop[j]) / (lop[j + 1] - lop[j])
    jj0 = (j - 1) % len(lo)
    jj1 = j % len(lo)
    return ((1 - wy) * ((1 - wx) * x[..., i, jj0] + wx * x[..., i, jj1]) +
            wy * ((1 - wx) * x[..., i + 1, jj0] + wx * x[..., i + 1, jj1])).astype(np.float32)


def _init(g):
    global _G
    _G = g


def _read(args):
    """One (file, year): -> ymd (n,3), x (n, [L,] n_lu) float32.  lev: None, list of Pa levels."""
    fn, year, var, w0, w1, levs = args
    G = _G
    with nc.Dataset(fn) as a:
        a.set_auto_mask(True)
        ymd, keep = _keep(a["time"], year, w0, w1)
        if not keep.any():
            return np.zeros((0, 3), int), None
        idx = np.where(keep)[0]
        sl = slice(int(idx[0]), int(idx[-1]) + 1)
        v = a[var]
        if levs is None:
            x = _to_lu(a, G, lambda rs: _clean(v[sl, rs, :]))
        else:
            pl = a["plev"][:].astype(float)
            if np.nanmax(pl) < 2000:
                pl = pl * 100
            ks = []
            for L in levs:
                k = int(np.argmin(np.abs(pl - L)))
                assert abs(pl[k] - L) < 1, f"{fn}: no {L} Pa"
                ks.append(k)
            x = np.stack([_to_lu(a, G, lambda rs, k=k: _clean(v[sl, k, rs, :])) for k in ks], axis=1)
        x = x[keep[sl]]
    return ymd[keep], x


def _read_ml(args):
    """CFday ua/va lowest K levels + ps -> ymd, |V| (n, K, n_lu), p_k (n, K, n_lu)."""
    fn_u, fn_v, fn_ps, year, w0, w1, K = args
    G = _G
    out = {}
    for var, fn in (("ua", fn_u), ("va", fn_v)):
        with nc.Dataset(fn) as a:
            a.set_auto_mask(True)
            ymd, keep = _keep(a["time"], year, w0, w1)
            if not keep.any():
                return np.zeros((0, 3), int), None, None, None
            idx = np.where(keep)[0]
            sl = slice(int(idx[0]), int(idx[-1]) + 1)
            out[var] = _to_lu(a, G, lambda rs: _clean(a[var][sl, 0:K, rs, :]))[keep[sl]]
            out["ymd_" + var] = ymd[keep]
            if var == "ua":
                hy = {k: np.ma.filled(np.ma.asarray(a[k][:]).astype(float), np.nan) for k in ("ap", "b", "a", "p0") if k in a.variables}
                in_file_ps = "ps" in a.variables
                if in_file_ps and fn_ps is None:
                    out["ps"] = _to_lu(a, G, lambda rs: _clean(a["ps"][sl, rs, :]))[keep[sl]]
                    out["ymd_ps"] = ymd[keep]
    if "ps" not in out:
        with nc.Dataset(fn_ps) as a:
            ymd, keep = _keep(a["time"], year, w0, w1)
            idx = np.where(keep)[0]
            sl = slice(int(idx[0]), int(idx[-1]) + 1)
            out["ps"] = _to_lu(a, G, lambda rs: _clean(a["ps"][sl, rs, :]))[keep[sl]]
            out["ymd_ps"] = ymd[keep]
    assert np.array_equal(out["ymd_ua"], out["ymd_va"]) and np.array_equal(out["ymd_ua"], out["ymd_ps"]), fn_u
    ps = out["ps"].astype(np.float64)
    # full-level pressure of the lowest K levels
    if "a" in hy and "p0" in hy:
        pk = hy["a"][:K, None, None] * hy["p0"] + hy["b"][:K, None, None] * ps[None]
    elif len(hy["ap"]) == len(hy["b"]) and len(hy["ap"]) > K and hy["b"][0] < 1.0 - 1e-9 and hy["b"][0] > 0.5:
        pk = hy["ap"][:K, None, None] + hy["b"][:K, None, None] * ps[None]
    else:                                   # IPSL: ap, b on half levels (klevp1), bottom first, b[0] = 1
        ph = hy["ap"][:K + 1, None, None] + hy["b"][:K + 1, None, None] * ps[None]
        pk = 0.5 * (ph[:-1] + ph[1:])
    pk = np.moveaxis(pk, 0, 1).astype(np.float32)            # (n, K, n_lu)
    spd = np.hypot(out["ua"], out["va"]).astype(np.float32)
    return out["ymd_ua"], spd, pk, out["ps"]


def amon(member, var, G, w0, w1):
    """Amon field -> dict (year, month) -> (n_lu,) float64, for the months of winters w0..w1."""
    d = vdir(member, "Amon", var)
    res = {}
    for f in sorted(glob.glob(f"{d}/*.nc")):
        a_, b_ = span(f)
        if b_ < w0 - 1 or a_ > w1:
            continue
        with nc.Dataset(f) as a:
            t = a["time"]
            dates = nc.num2date(t[:], t.units, calendar=getattr(t, "calendar", "standard"))
            ym = np.array([(x.year, x.month) for x in dates])
            win = np.where(ym[:, 1] >= 9, ym[:, 0] + 1, ym[:, 0])
            keep = np.isin(ym[:, 1], MON) & (win >= w0) & (win <= w1)
            if not keep.any():
                continue
            idx = np.where(keep)[0]
            sl = slice(int(idx[0]), int(idx[-1]) + 1)
            x = _to_lu(a, G, lambda rs: _clean(a[var][sl, rs, :]))
            for xi, (yy, mm), k in zip(x, ym[sl], keep[sl]):
                if k:
                    res[(int(yy), int(mm))] = xi.astype(np.float64)
    return res


def daily(member, table, var, G, w0, w1, levs=None):
    d = vdir(member, table, var)
    if d is None:
        return None
    jobs = [(f, y, var, w0, w1, levs) for f, y in jobs_for(d, w0, w1)]
    got = {}
    with ProcessPoolExecutor(max_workers=NW, initializer=_init, initargs=(G,)) as ex:
        for ymd, x in ex.map(_read, jobs):
            for k in range(len(ymd)):
                got[tuple(int(v) for v in ymd[k])] = x[k]
    return got if got else None


def hyps(ps, p, tas):
    """Height (m) of pressure p above a surface at ps, layer T = tas - LAPSE h / 2."""
    L = np.log(ps / p)
    return RD * tas * L / GRAV / (1.0 + 0.5 * LAPSE * RD * L / GRAV)


def setup(model):
    S = xd.cmip6_setup(model)
    nlon = len(S["flat"]) // (S["r1"] - S["r0"])
    inv = json.load(open(INV))[model]
    f = sorted(glob.glob(inv["paths"]["ssp585"]["day/sfcWind"] + "/*.nc"))[0]
    with nc.Dataset(f) as d:
        lat, lon = np.asarray(d["lat"][:], float), np.asarray(d["lon"][:], float)
    assert len(lon) == nlon
    flon = np.tile(lon, S["r1"] - S["r0"])
    G = dict(r0=S["r0"], r1=S["r1"], lu=S["lu"], nlat=len(lat), nlon=nlon, lat=lat, lon=lon,
             tlat=S["flat"][S["lu"]], tlon=flon[S["lu"]])
    R = len(rg.REGIONS)
    Wm = np.zeros((len(S["lu"]), R))
    for r, pos in enumerate(S["lpos"]):
        Wm[pos, r] = S["w"][pos]
    return S, G, Wm


def rmean(x, m, Wm):
    """x, m (D, n) -> region means over valid cells (D, R), valid weight fraction (D, R)."""
    m = m & np.isfinite(x)
    num = np.where(m, x, 0.0) @ Wm
    den = m.astype(float) @ Wm
    tot = Wm.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (num / den).astype(np.float32), (den / np.where(tot > 0, tot, np.nan)).astype(np.float32)


def stack(got, keys):
    return np.stack([got[k] for k in keys]).astype(np.float64)


def members(model):
    inv = json.load(open(INV))[model]
    return {e: inv["paths"][e]["day/sfcWind"].split("/day/")[0] for e in PERIODS}


def run_pl(model):
    t0 = time.time()
    S, G, Wm = setup(model)
    mem = members(model)
    out, info = {}, {"model": model, "HMIN": HMIN}
    # orography (fx, this member or any experiment of the model)
    orog = None
    for e in PERIODS:
        d = vdir(mem[e], "fx", "orog")
        if d:
            with nc.Dataset(sorted(glob.glob(f"{d}/*.nc"))[0]) as a:
                orog = _to_lu(a, G, lambda rs: _clean(a["orog"][rs, :])).astype(np.float64)
                info["orog"] = d
            break
    for exp, (w0, w1) in PERIODS.items():
        m = mem[exp]
        src = {}
        src["sfc"] = daily(m, "day", "sfcWind", G, w0, w1)
        src["uas"] = daily(m, "day", "uas", G, w0, w1)
        src["vas"] = daily(m, "day", "vas", G, w0, w1)
        src["psl"] = daily(m, "day", "psl", G, w0, w1)
        tab = "day" if vdir(m, "day", "ua") else "Eday"
        src["u18"] = daily(m, tab, "ua", G, w0, w1, [100000.0, 85000.0])
        src["v18"] = daily(m, tab, "va", G, w0, w1, [100000.0, 85000.0])
        info[f"{exp}|uv1000_table"] = tab
        e9 = vdir(m, "Eday", "ua") and vdir(m, "Eday", "va")
        src["u9"] = daily(m, "Eday", "ua", G, w0, w1, [92500.0]) if e9 else None
        src["v9"] = daily(m, "Eday", "va", G, w0, w1, [92500.0]) if e9 else None
        src["zg"] = daily(m, "day", "zg", G, w0, w1, [100000.0]) if (orog is not None and vdir(m, "day", "zg")) else None
        ps_m, psl_m, tas_m = (amon(m, v, G, w0, w1) for v in ("ps", "psl", "tas"))
        need = ["sfc", "psl", "u18", "v18"]
        keys = sorted(set.intersection(*[set(src[k]) for k in need]) & {k for k in src["sfc"] if (k[0], k[1]) in ps_m and (k[0], k[1]) in psl_m and (k[0], k[1]) in tas_m})
        info[f"{exp}|n_days"] = len(keys)
        info[f"{exp}|have"] = {k: (v is not None) for k, v in src.items()}
        ymd = np.array(keys)
        ym = [(k[0], k[1]) for k in keys]
        sfc = stack(src["sfc"], keys)
        psl = stack(src["psl"], keys)
        ratio = np.stack([ps_m[k] / psl_m[k] for k in ym])
        tas = np.stack([tas_m[k] for k in ym])
        ps = psl * ratio
        h1000 = hyps(ps, 1e5, tas)
        h925 = hyps(ps, 92500.0, tas)
        u18, v18 = stack(src["u18"], keys), stack(src["v18"], keys)
        V1000, V850 = np.hypot(u18[:, 0], v18[:, 0]), np.hypot(u18[:, 1], v18[:, 1])
        if src["uas"] is not None and src["vas"] is not None and all(k in src["uas"] and k in src["vas"] for k in keys):
            s10v = np.hypot(stack(src["uas"], keys), stack(src["vas"], keys))
        else:
            s10v = np.full_like(sfc, np.nan)
        if src["u9"] is not None and src["v9"] is not None and all(k in src["u9"] and k in src["v9"] for k in keys):
            V925 = np.hypot(stack(src["u9"], keys)[:, 0], stack(src["v9"], keys)[:, 0])
        else:
            V925 = np.full_like(sfc, np.nan)
            info[f"{exp}|V925"] = "missing or incomplete Eday"
        zgh = (stack(src["zg"], keys)[:, 0] - orog[None]) if (src["zg"] is not None and all(k in src["zg"] for k in keys)) else np.full_like(sfc, np.nan)
        allm = np.ones_like(sfc, bool)
        m1 = (h1000 >= HMIN) & np.isfinite(V1000) & np.isfinite(sfc)
        m9 = (h925 >= HMIN) & np.isfinite(V925) & np.isfinite(sfc)
        info[f"{exp}|V1000_missing_frac_where_h>=HMIN"] = float(np.mean(~np.isfinite(V1000[h1000 >= HMIN]))) if (h1000 >= HMIN).any() else None
        info[f"{exp}|V1000_missing_frac_where_h<0"] = float(np.mean(~np.isfinite(V1000[h1000 < 0]))) if (h1000 < 0).any() else None
        o = {}
        o["sfc_all"], _ = rmean(sfc, allm, Wm)
        o["s10v_all"], _ = rmean(s10v, allm, Wm)
        o["V850_all"], o["frac850"] = rmean(V850, hyps(ps, 85000.0, tas) >= HMIN, Wm)
        o["h1000_all"], _ = rmean(h1000, allm, Wm)
        o["zgh1000_all"], _ = rmean(zgh, allm, Wm)
        o["ps_all"], _ = rmean(ps, allm, Wm)
        for tag, mk in (("m1000", m1), ("m925", m9)):
            o[f"frac_{tag}"] = rmean(sfc, mk, Wm)[1]
            for nm, x in (("sfc", sfc), ("s10v", s10v), ("V1000", V1000), ("V925", V925), ("V850", V850),
                          ("h1000", h1000), ("h925", h925), ("zgh1000", zgh)):
                if tag == "m925" and nm in ("V1000", "h1000", "zgh1000"):
                    continue
                o[f"{nm}_{tag}"], _ = rmean(x, mk, Wm)
        # cell climatology of h1000 (control) for the height report
        for k, v in o.items():
            out[f"{exp}|{k}"] = v
        out[f"{exp}|ymd"] = ymd
        if exp == "historical":
            out["cell_h1000_ctl"] = np.nanmean(h1000, 0).astype(np.float32)
            out["cell_fracabove_ctl"] = np.mean(h1000 >= HMIN, 0).astype(np.float32)
            if orog is not None:
                out["cell_orog"] = orog.astype(np.float32)
        print(f"  {model} {exp}: {len(keys)} days, V925 {'yes' if np.isfinite(V925).any() else 'no'}, s10v {'yes' if np.isfinite(s10v).any() else 'no'}, "
              f"zg-orog {'yes' if np.isfinite(zgh).any() else 'no'} [{time.time() - t0:.0f}s]", flush=True)
    out["regions"] = np.array(rg.REGIONS)
    out["info"] = np.array(json.dumps(info))
    os.makedirs(OUTD, exist_ok=True)
    np.savez_compressed(f"{OUTD}/hubheight_pl_{model}{SUF}.npz", **out)
    print(f"Saved {OUTD}/hubheight_pl_{model}{SUF}.npz [{time.time() - t0:.0f}s]", flush=True)


def interp_z(spd, z, zt):
    """spd, z (D, K, n) -> speed at height zt, linear in ln z (lowest two levels extrapolate below)."""
    lz, lt = np.log(np.maximum(z, 0.1)), np.log(zt)
    j = np.clip((z < zt).sum(1), 1, z.shape[1] - 1)               # upper bracketing index
    jj = j[:, None]
    z0, z1 = np.take_along_axis(lz, jj - 1, 1)[:, 0], np.take_along_axis(lz, jj, 1)[:, 0]
    s0, s1 = np.take_along_axis(spd, jj - 1, 1)[:, 0], np.take_along_axis(spd, jj, 1)[:, 0]
    return s0 + (s1 - s0) * (lt - z0) / (z1 - z0)


def run_ml(model):
    t0 = time.time()
    S, G, Wm = setup(model)
    mem = members(model)
    out, info = {}, {"model": model, "ZHUB": ZHUB}
    for exp, (w0, w1) in PERIODS.items():
        m = mem[exp]
        du, dv, dp = vdir(m, "CFday", "ua"), vdir(m, "CFday", "va"), vdir(m, "CFday", "ps")
        f0 = sorted(glob.glob(f"{du}/*.nc"))[0]
        with nc.Dataset(f0) as a:
            in_ps = "ps" in a.variables
            if "b" in a.variables:
                b = np.ma.filled(np.ma.asarray(a["b"][:]).astype(float), np.nan)
            ldim = a["ua"].dimensions[1]
            # number of levels with b above 0.9 (lowest ~800 m), bottom-first ordering checked
            assert b[0] > b[1], f"{model}: levels not bottom first"
            K = int((b > 0.9).sum()) + 1
            info[f"{exp}|K"] = K
            info[f"{exp}|levdim"] = ldim
        jobs = []
        psfiles = sorted(glob.glob(f"{dp}/*.nc")) if (dp and not in_ps) else []
        for fu, y in jobs_for(du, w0, w1):
            fv = fu.replace("/ua/", "/va/").replace("ua_CFday", "va_CFday")
            fv = fv if os.path.exists(fv) else None
            if fv is None:
                cand = [f for f in glob.glob(f"{dv}/*.nc") if span(f)[0] <= y <= span(f)[1]]
                fv = cand[0]
            fp = None
            if not in_ps:
                fp = [f for f in psfiles if span(f)[0] <= y <= span(f)[1]][0]
            jobs.append((fu, fv, fp, y, w0, w1, K))
        got = {}
        with ProcessPoolExecutor(max_workers=NW, initializer=_init, initargs=(G,)) as ex:
            for r in ex.map(_read_ml, jobs):
                ymd = r[0]
                for k in range(len(ymd)):
                    got[tuple(int(v) for v in ymd[k])] = (r[1][k], r[2][k], r[3][k])
        tas_m = amon(m, "tas", G, w0, w1)
        keys = sorted(k for k in got if (k[0], k[1]) in tas_m)
        spd = np.stack([got[k][0] for k in keys]).astype(np.float64)
        pk = np.stack([got[k][1] for k in keys]).astype(np.float64)
        ps = np.stack([got[k][2] for k in keys]).astype(np.float64)
        tas = np.stack([tas_m[(k[0], k[1])] for k in keys])
        z = hyps(ps[:, None], pk, tas[:, None])
        V100 = interp_z(spd, z, ZHUB)
        o = {}
        allm = np.ones(V100.shape, bool)
        o["V100ml_all"], o["frac_ml"] = rmean(V100, allm, Wm)
        o["Vlow_all"], _ = rmean(spd[:, 0], allm, Wm)
        o["zlow_all"], _ = rmean(z[:, 0], allm, Wm)
        o["z2_all"], _ = rmean(z[:, 1], allm, Wm)
        o["ps_ml_all"], _ = rmean(ps, allm, Wm)
        for k, v in o.items():
            out[f"{exp}|{k}"] = v
        out[f"{exp}|ymd"] = np.array(keys)
        info[f"{exp}|n_days"] = len(keys)
        info[f"{exp}|zlow_median"] = float(np.nanmedian(z[:, 0]))
        info[f"{exp}|z_levels_median"] = [float(np.nanmedian(z[:, k])) for k in range(K)]
        print(f"  {model} {exp} ml: {len(keys)} days, K={K}, level heights (median) "
              + ", ".join(f"{v:.0f}" for v in info[f"{exp}|z_levels_median"]) + f" m [{time.time() - t0:.0f}s]", flush=True)
        del spd, pk, ps, z, got
    out["regions"] = np.array(rg.REGIONS)
    out["info"] = np.array(json.dumps(info))
    os.makedirs(OUTD, exist_ok=True)
    np.savez_compressed(f"{OUTD}/hubheight_ml_{model}{SUF}.npz", **out)
    print(f"Saved {OUTD}/hubheight_ml_{model}{SUF}.npz [{time.time() - t0:.0f}s]", flush=True)


def main():
    model = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) > 2 else "pl"
    if mode in ("pl", "both"):
        run_pl(model)
    if mode in ("ml", "both"):
        run_ml(model)


if __name__ == "__main__":
    main()
