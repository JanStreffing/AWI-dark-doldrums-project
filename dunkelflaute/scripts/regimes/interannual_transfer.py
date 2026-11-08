"""
Interannual test of the surface-wind transfer s / G at fixed G: in warm against cold winters, does
the transfer rise (as in the AWI-CM3 / EERIE runs under warming) or fall (as the CMIP6 forced
change at fixed G does)? And does each model's interannual sensitivity predict its forced change?

Definitions as geostrophic.py / transfer_stability.py: s daily region land-mean 10 m speed, G daily
region-mean geostrophic speed from daily-mean msl on 2.5 deg, 8 European regions (regime_analysis.EU),
Sep-Mar winters labelled by the Jan-Mar year.

Per dataset, region and winter w:
  y_q5   ln( sum_k sbar_wk / sum_k Gbar_wk ) over the 5 G quintile classes k of the region (edges from
         all days of the analysed winters): the winter transfer at climatological G-class weights
  y_res  winter mean of the residual of a daily fit ln s = a_m + b_m ln G + c_m (ln G)^2 per
         calendar month m (y_dir: the same plus cos and sin of the geostrophic direction)
  y_rat  ln( mean s / mean G ) (plain ratio)
  T      winter-mean region land 2 m temperature; lnG = ln(winter-mean G)
  ug, vg winter-mean region-mean geostrophic components (direction control)
All series are detrended (linear by default; quadratic and first differences as variants). The
sensitivity is the coefficient b of  y' = a + b T' + c lnG'  (y_q5 primary), in % per K (100 b).
EU value = mean over the 8 regions; 90 % ranges from 1000 bootstrap draws of winters (same draw in
every region, so cross-region correlation is kept).

Forced change at fixed G (for the cross-model test): within-G-quintile part of the transfer change
(transfer_stability.run_decomp, as transfer_stability_cmip6.py), per region; per K of the regional
winter land warming of the same windows; EU = mean within / mean warming.
  CMIP6      historical 1951-1980 vs ssp585 2071-2100 (psl_*.npz, daily_cmip6_*.npz); interannual on
             historical 1951-2014 (interannual_extract_cmip6.py)
  100 km     hi1950g 1951-1980 vs 2071-2099 (s10); interannual on ctl1950g 1953-2099
  31 km      ctl1950d 1953-1969 vs ssp585d 2081-2099 (spd10); interannual on ctl1950d 1953-1969
  9 km       1950c 1953-1969 vs 2080c 2083-2092 (raw 10 m speed); interannual on 1950c 1953-1969
  EERIE      ctl 1972-2050 vs ssp245 2041-2070 (spd10); interannual on ctl 1972-2050
Observations: HadISD v3.4.3.2025f stations (interannual_hadisd_process.py), winters 1980-2025,
stations with >= 30 winters of >= 80 % valid days; per-station speed normalised by its calendar-month
mean, region index = mean over stations; G from ERA5 daily msl (data/regimes/msl_era5_-.npz);
T from the stations (mean of per-station winter anomalies) and from ERA5 t2m; ERA5 t850 - t2m as
inversion index. Variant B: per-station y_q5, detrended per station, then the region mean.
Method check: ERA5 10 m wind (interannual_extract_era5_wind.py) as region land means (the model
pipeline) and at the station points on the station days (the station pipeline).

Usage: python3 interannual_transfer.py      (legacy python, login node, a few minutes)
Output: dunkelflaute/results/regimes/interannual_transfer.json, interannual_transfer_*.png
"""
import glob
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("regions", "regimes"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/{p}")
os.environ.setdefault("NB", "200")
import regions as rg                       # noqa: E402
import extract_msl as xm                   # noqa: E402
import geostrophic as gs                   # noqa: E402
import transfer_stability as ts            # noqa: E402

EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
D = f"{_ROOT}/dunkelflaute/data"
DC = f"{D}/cmip6"
DI = f"{D}/regimes/interannual_cmip6"
DO = f"{D}/regimes/interannual_ours"
OUT = f"{_ROOT}/dunkelflaute/results/regimes"
NBOOT = 1000
RHO, OMEGA, RE = 1.25, 7.292e-5, 6.371e6
BOX = dict(rg.BOXES, DE=rg.DE_BOX)
MDAYS = {9: 30, 10: 31, 11: 30, 12: 31, 1: 31, 2: 28, 3: 31}
RNG = np.random.default_rng(7)


# ------------------------------------------------------------------ geostrophy helpers
def geo_uv(p, lat, lon, periodic):
    phi = np.deg2rad(lat)[None, :, None]
    dphi, dlam = np.deg2rad(lat[1] - lat[0]), np.deg2rad(lon[1] - lon[0])
    dpy = np.full_like(p, np.nan)
    dpy[:, 1:-1] = (p[:, 2:] - p[:, :-2]) / (2 * dphi * RE)
    if periodic:
        dpx = (np.roll(p, -1, axis=2) - np.roll(p, 1, axis=2)) / (2 * dlam * RE * np.cos(phi))
    else:
        dpx = np.full_like(p, np.nan)
        dpx[:, :, 1:-1] = (p[:, :, 2:] - p[:, :, :-2]) / (2 * dlam * RE * np.cos(phi))
    f = 2 * OMEGA * np.sin(phi)
    return -dpy / (RHO * f), dpx / (RHO * f)


def G_uv(msl, lat, lon, periodic):
    p = msl.astype(np.float64)
    G = gs.region_G(gs.geostrophic(p, lat, lon, periodic), lat, lon, EU)
    u, v = geo_uv(p, lat, lon, periodic)
    return G, gs.region_G(u, lat, lon, EU), gs.region_G(v, lat, lon, EU)


def daily_mean(times, x):
    day = (times - np.timedelta64(1, "s")).astype("datetime64[D]")
    u, inv, cnt = np.unique(day, return_inverse=True, return_counts=True)
    s = np.zeros((len(u),) + x.shape[1:])
    np.add.at(s, inv, x)
    full = cnt == cnt.max()
    return u[full], s[full] / cnt[full].reshape((-1,) + (1,) * (x.ndim - 1))


def month_of(k):
    return k // 100 % 100


# ------------------------------------------------------------------ winter tables
def winter_table(k, s, G, ug, vg, winters, min_days):
    """k (D,) ymd ints; s, G, ug, vg (D, R). -> dict of (W, R) arrays."""
    w = gs.winter_of(k)
    mon = month_of(k)
    sel = np.isin(w, winters) & np.isin(mon, rg.MONTHS)
    W, R = len(winters), s.shape[1]
    out = {n: np.full((W, R), np.nan) for n in ("y_q5", "y_res", "y_dir", "y_rat", "lnG", "ug", "vg", "n")}
    for r in range(R):
        ok = sel & np.isfinite(s[:, r]) & np.isfinite(G[:, r]) & (s[:, r] > 0)
        if ok.sum() < 100:
            continue
        e = np.nanpercentile(G[ok, r], [20, 40, 60, 80])
        c = np.digitize(G[:, r], e)
        lns, lng = np.log(np.where(ok, s[:, r], np.nan)), np.log(np.where(ok, G[:, r], np.nan))
        th = np.arctan2(vg[:, r], ug[:, r]) if ug is not None else None
        res = np.full(len(k), np.nan)
        rdir = np.full(len(k), np.nan)
        for m in rg.MONTHS:
            mm = ok & (mon == m)
            X = np.column_stack([np.ones(mm.sum()), lng[mm], lng[mm] ** 2])
            b = np.linalg.lstsq(X, lns[mm], rcond=None)[0]
            res[mm] = lns[mm] - X @ b
            if th is not None:
                X2 = np.column_stack([X, np.cos(th[mm]), np.sin(th[mm])])
                b2 = np.linalg.lstsq(X2, lns[mm], rcond=None)[0]
                rdir[mm] = lns[mm] - X2 @ b2
        for j, y in enumerate(winters):
            d = ok & (w == y)
            if d.sum() < min_days:
                continue
            cs = [s[d & (c == q), r] for q in range(5)]
            cg = [G[d & (c == q), r] for q in range(5)]
            if min(len(x) for x in cs) >= 3:
                out["y_q5"][j, r] = np.log(sum(x.mean() for x in cs) / sum(x.mean() for x in cg))
            out["y_res"][j, r] = res[d].mean()
            out["y_dir"][j, r] = rdir[d].mean()
            out["y_rat"][j, r] = np.log(s[d, r].mean() / G[d, r].mean())
            out["lnG"][j, r] = np.log(G[d, r].mean())
            if ug is not None:
                out["ug"][j, r] = ug[d, r].mean()
                out["vg"][j, r] = vg[d, r].mean()
            out["n"][j, r] = d.sum()
    return out


def winter_T_daily(k, t, winters, min_days):
    w = gs.winter_of(k)
    sel = np.isin(month_of(k), rg.MONTHS)
    out = np.full((len(winters), t.shape[1]), np.nan)
    for j, y in enumerate(winters):
        d = sel & (w == y)
        if d.sum() >= min_days:
            out[j] = np.nanmean(t[d], axis=0)
    return out


def winter_T_monthly(ym, t, winters):
    w = np.where(ym[:, 1] >= 9, ym[:, 0] + 1, ym[:, 0])
    out = np.full((len(winters), t.shape[1]), np.nan)
    for j, y in enumerate(winters):
        d = (w == y) & np.isin(ym[:, 1], rg.MONTHS)
        if d.sum() == 7:
            wt = np.array([MDAYS[m] for m in ym[d, 1]], float)
            out[j] = (t[d] * wt[:, None]).sum(0) / wt.sum()
    return out


# ------------------------------------------------------------------ regression
def detrend(x, yrs, order):
    x = np.asarray(x, float)
    ok = np.isfinite(x)
    if ok.sum() < order + 3:
        return np.full_like(x, np.nan)
    t = (yrs - yrs.mean()) / 10.0
    c = np.polyfit(t[ok], x[ok], order)
    return x - np.polyval(c, t)


def ols_b(y, X):
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    if ok.sum() < X.shape[1] + 4:
        return np.full(X.shape[1], np.nan)
    A = np.column_stack([np.ones(ok.sum()), X[ok]])
    return np.linalg.lstsq(A, y[ok], rcond=None)[0][1:]


def prepare(tab, T, winters, ycol, xcols, order):
    """-> y' (W, R), X' (W, R, p) detrended (order) or first-differenced (order 'diff')."""
    yrs = np.asarray(winters, float)
    cols = [T] + [tab[c] for c in xcols if c != "T"]
    Y = tab[ycol] if isinstance(ycol, str) else ycol
    R = Y.shape[1]
    if order == "diff":
        cons = np.diff(yrs) == 1
        f = lambda a: np.where(cons[:, None], np.diff(a, axis=0), np.nan)    # noqa: E731
        return f(Y), np.stack([f(c) for c in cols], axis=2)
    y = np.stack([detrend(Y[:, r], yrs, order) for r in range(R)], axis=1)
    X = np.stack([np.stack([detrend(c[:, r], yrs, order) for r in range(R)], axis=1) for c in cols], axis=2)
    return y, X


def sensitivity(tab, T, winters, ycol="y_q5", controls=("lnG",), order=1, nboot=NBOOT, Xextra=None):
    """b_T per region and EU mean (% per K) with a winter bootstrap."""
    y, X = prepare(tab, T, winters, ycol, controls, order)
    if Xextra is not None:                # replace T by another regressor (e.g. inversion)
        X[:, :, 0] = Xextra
    R = y.shape[1]
    b = np.array([ols_b(y[:, r], X[:, r]) [0] for r in range(R)]) * 100
    W = y.shape[0]
    bs = np.empty(nboot)
    for i in range(nboot):
        ix = RNG.integers(0, W, W)
        bs[i] = np.nanmean([ols_b(y[ix, r], X[ix, r])[0] for r in range(R)]) * 100
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=2)
    sdT = np.nanmean([np.nanstd(X[ok[:, r], r, 0]) for r in range(R)])
    return dict(eu=float(np.nanmean(b)), ci90=[float(v) for v in np.nanpercentile(bs, [5, 95])], se=float(np.nanstd(bs)),
                per_region={r_: float(v) for r_, v in zip(EU, b)}, n_winters=int(ok.any(axis=1).sum()), sd_T=float(sdT))


def all_sens(tab, T, winters, extra=None, order=1):
    out = dict(primary=sensitivity(tab, T, winters, "y_q5", ("lnG",), order))
    out["q5_noG"] = sensitivity(tab, T, winters, "y_q5", (), order, nboot=200)
    out["res"] = sensitivity(tab, T, winters, "y_res", ("lnG",), order, nboot=200)
    out["ratio"] = sensitivity(tab, T, winters, "y_rat", ("lnG",), order, nboot=200)
    if np.isfinite(tab["ug"]).any():
        out["dir_daily"] = sensitivity(tab, T, winters, "y_dir", ("lnG",), order, nboot=200)
        out["q5_dir"] = sensitivity(tab, T, winters, "y_q5", ("lnG", "ug", "vg"), order)
    alt = "diff" if order != "diff" else 1
    out["detrend_quadratic"] = sensitivity(tab, T, winters, "y_q5", ("lnG",), 2, nboot=200)
    out["first_diff"] = sensitivity(tab, T, winters, "y_q5", ("lnG",), alt, nboot=200)
    # how much of the T variance the circulation carries
    y, X = prepare(tab, T, winters, "y_q5", ("lnG", "ug", "vg") if np.isfinite(tab["ug"]).any() else ("lnG",), order)
    r2 = []
    for r in range(y.shape[1]):
        ok = np.all(np.isfinite(X[:, r]), axis=1)
        if ok.sum() > 8 and X.shape[2] > 1:
            A = np.column_stack([np.ones(ok.sum()), X[ok, r, 1:]])
            c = np.linalg.lstsq(A, X[ok, r, 0], rcond=None)[0]
            r2.append(1 - np.var(X[ok, r, 0] - A @ c) / np.var(X[ok, r, 0]))
    out["r2_T_on_circulation"] = float(np.mean(r2)) if r2 else None
    return out


# ------------------------------------------------------------------ forced change at fixed G
def forced(kc, sc, Gc, kf, sf, Gf, Tc, Tf, rng):
    dc = dict(G=Gc, s=sc, winter=gs.winter_of(kc))
    df = dict(G=Gf, s=sf, winter=gs.winter_of(kf))
    selc, self_ = np.ones(len(Gc), bool), np.ones(len(Gf), bool)
    cc, cf, _ = ts.quintile_classes(Gc, Gf, selc)
    r = ts.run_decomp(dc, df, selc, self_, cc, cf, 5, "s", rng)
    within = np.array([r["within"]["per_region"][x] for x in EU])
    tot = ts.tr_change(dc, df, selc, self_, "s")
    dT = Tf - Tc
    return dict(within_G=float(within.mean()), within_ci90=r["within"]["ci90"], total=float(tot.mean()),
                between_G=r["between"]["eu_mean"], dT_region=float(dT.mean()),
                within_per_K=float(within.mean() / dT.mean()),
                per_region={x: dict(within_G=float(a), dT=float(b)) for x, a, b in zip(EU, within, dT)})


# ------------------------------------------------------------------ loaders: our runs
def load_regions_run(run, svar, path=None):
    z = np.load(path or f"{rg.OUT_DIR}/{run}.npz")
    t = z["time"].astype("datetime64[s]")
    ri = [list(z["regions"]).index(r) for r in EU]
    day, s = daily_mean(t, z[svar][:, ri].astype(np.float64))
    m = np.load(f"{xm.OUT}/msl_{run}.npz")
    md = m["dates"].astype("datetime64[D]")
    k_s, k_m = gs.ymd_of(day), gs.ymd_of(md)
    common, i_s, i_m = np.intersect1d(k_s, k_m, return_indices=True)
    keep = np.isin(month_of(common), rg.MONTHS)
    common, i_s, i_m = common[keep], i_s[keep], i_m[keep]
    G, u, v = G_uv(m["msl"][i_m], m["lat"], m["lon"], False)
    out = dict(k=common, s=s[i_s], G=G, ug=u, vg=v)
    if "t2m" in z.files:
        dt, tt = daily_mean(t, z["t2m"][:, ri].astype(np.float64))
        kt = gs.ymd_of(dt)
        out["kt"], out["t2"] = kt, tt
    return out


def attach_t2(d, fn):
    z = np.load(fn)
    d["kt"], d["t2"] = gs.ymd_of(z["day"].astype("datetime64[D]")), z["t2"]
    return d


def load_tco95(run):
    fs = sorted(glob.glob(f"{ts.PARTS}/{run}_*.npz"))
    Z = [np.load(f) for f in fs]
    day = np.concatenate([z["day"] for z in Z])
    o = np.argsort(day)
    s10 = np.concatenate([z["s10"] for z in Z])[o]
    t2 = np.concatenate([z["t2"] for z in Z])[o]
    ms = np.concatenate([z["mslbox"] for z in Z])
    md = np.concatenate([z["mslday"] for z in Z])
    om = np.argsort(md)
    assert np.array_equal(md[om], day[o])
    G, u, v = G_uv(ms[om], xm.LATC, xm.LONC, False)
    k = gs.ymd_of(day[o].astype("datetime64[D]"))
    keep = np.isin(month_of(k), rg.MONTHS)
    return dict(k=k[keep], s=s10[keep], G=G[keep], ug=u[keep], vg=v[keep], kt=k[keep], t2=t2[keep])


def sub(d, w0, w1):
    w = gs.winter_of(d["k"])
    m = (w >= w0) & (w <= w1)
    return {k_: (v[m] if k_ in ("k", "s", "G", "ug", "vg") else v) for k_, v in d.items()}


def T_period(d, w0, w1):
    w = gs.winter_of(d["kt"])
    m = (w >= w0) & (w <= w1) & np.isin(month_of(d["kt"]), rg.MONTHS)
    return np.nanmean(d["t2"][m], axis=0)


def ours():
    rng = np.random.default_rng(3)
    spec = {
        "100km": dict(ctl=("tco95", "ctl1950g", 1953, 2099), cmp=(("hi1950g", 1951, 1980), ("hi1950g", 2071, 2099))),
        "31km": dict(ctl=("regions", "tco319_ctl1950d", 1953, 1969), cmp=(("tco319_ctl1950d", 1953, 1969), ("tco319_ssp585d", 2081, 2099))),
        "9km": dict(ctl=("regions", "dart9km_1950c", 1953, 1969), cmp=(("dart9km_1950c", 1953, 1969), ("dart9km_2080c", 2083, 2092))),
        "eerie_ctl": dict(ctl=("regions", "eerie_ctl", 1972, 2050), cmp=(("eerie_ctl", 1972, 2050), ("eerie_ssp245", 2041, 2070))),
    }
    cache = {}

    def get(run):
        if run in cache:
            return cache[run]
        if run.startswith("ctl1950g") or run.startswith("hi1950g"):
            d = load_tco95(run)
        elif run.startswith("tco319"):
            d = attach_t2(load_regions_run(run, "spd10"), f"{DO}/{run}_t2.npz")
        elif run.startswith("dart9km"):
            d = load_regions_run(run, "spd10", path=f"{rg.OUT_DIR}/transfer_checks_{run}_spd10.npz")
            d = attach_t2(d, f"{DO}/{run}_t2.npz")
        else:
            d = load_regions_run(run, "spd10")
        cache[run] = d
        return d

    res = {}
    for name, sp in spec.items():
        _, run, w0, w1 = sp["ctl"]
        d = get(run)
        winters = list(range(w0, w1 + 1))
        tab = winter_table(d["k"], d["s"], d["G"], d["ug"], d["vg"], winters, 150)
        T = winter_T_daily(d["kt"], d["t2"], winters, 150)
        out = dict(winters=[w0, w1], interannual=all_sens(tab, T, winters))
        (cn, c0, c1), (fn, f0, f1) = sp["cmp"]
        a, b = sub(get(cn), c0, c1), sub(get(fn), f0, f1)
        out["forced"] = forced(a["k"], a["s"], a["G"], b["k"], b["s"], b["G"], T_period(get(cn), c0, c1), T_period(get(fn), f0, f1), rng)
        out["forced"]["windows"] = [cn, c0, c1, fn, f0, f1]
        if name == "100km":      # also the 100km comparison of 13b/13d (ctl1950g 1953-2099 vs hi1950g 2081-2099)
            a, b = sub(get("ctl1950g"), 1953, 2099), sub(get("hi1950g"), 2081, 2099)
            out["forced_ctl_vs_hi"] = forced(a["k"], a["s"], a["G"], b["k"], b["s"], b["G"], T_period(get("ctl1950g"), 1953, 2099),
                                             T_period(get("hi1950g"), 2081, 2099), rng)
            # interannual in the transient run too (1951-2099, quadratic detrend)
            dh = get("hi1950g")
            wh = list(range(1951, 2100))
            th = winter_table(dh["k"], dh["s"], dh["G"], dh["ug"], dh["vg"], wh, 150)
            out["interannual_hi1950g_quadratic"] = sensitivity(th, winter_T_daily(dh["kt"], dh["t2"], wh, 150), wh, order=2)
        res[name] = out
        p, f = out["interannual"]["primary"], out["forced"]
        print(f"  {name:10s} interannual {p['eu']:+.2f} %/K [{p['ci90'][0]:+.2f}, {p['ci90'][1]:+.2f}] ({p['n_winters']} winters)  "
              f"forced within-G {f['within_G']:+.2f}% dT {f['dT_region']:.2f} K -> {f['within_per_K']:+.2f} %/K", flush=True)
    return res


# ------------------------------------------------------------------ CMIP6
def cmip6():
    rng = np.random.default_rng(2)
    inv = json.load(open(f"{DC}/inventory.json"))
    res = {}
    for fn in sorted(glob.glob(f"{DI}/*.npz")):
        name = os.path.basename(fn)[:-4]
        z = np.load(fn)
        kw = z["wind_ymd"] @ np.array([10000, 100, 1])
        kp = z["psl_ymd"] @ np.array([10000, 100, 1])
        common, i_w, i_p = np.intersect1d(kw, kp, return_indices=True)
        s, G, u, v = z["s"][i_w], z["G"][i_p], z["ug"][i_p], z["vg"][i_p]
        Tm = winter_T_monthly(z["historical|tas_ym"], z["historical|tas"], list(range(1951, 2015)))
        out = {}
        for tag, (w0, w1) in (("1951_2014", (1951, 2014)), ("1951_1980", (1951, 1980))):
            winters = list(range(w0, w1 + 1))
            tab = winter_table(common, s, G, u, v, winters, 150)
            T = Tm[w0 - 1951:w1 - 1951 + 1]
            out[f"interannual_{tag}"] = all_sens(tab, T, winters) if tag == "1951_2014" else dict(primary=sensitivity(tab, T, winters))
        out["interannual"] = out.pop("interannual_1951_2014")
        # forced: the transfer_stability_cmip6 inputs
        P = np.load(f"{DC}/psl_{name}.npz")
        Z = np.load(f"{DC}/daily_cmip6_{name}.npz")
        ri = [rg.REGIONS.index(r) for r in EU]
        S = {}
        for exp, (y0, y1) in (("historical", (1951, 1980)), ("ssp585", (2071, 2100))):
            a = P[f"{exp}|ymd"] @ np.array([10000, 100, 1])
            b = Z[f"{exp}|ymd"] @ np.array([10000, 100, 1])
            cm, i_a, i_b = np.intersect1d(a, b, return_indices=True)
            w = gs.winter_of(cm)
            keep = (w >= y0) & (w <= y1) & np.isin(month_of(cm), rg.MONTHS)
            cm, i_a, i_b = cm[keep], i_a[keep], i_b[keep]
            Gx = gs.region_G(gs.geostrophic(P[f"{exp}|psl"][i_a].astype(np.float64), P["lat"], P["lon"], True), P["lat"], P["lon"], EU)
            S[exp] = (cm, Z[f"{exp}|spd10"][i_b][:, ri].astype(np.float64), Gx)
        Tc = np.nanmean(Tm[:30], axis=0)
        Tf = np.nanmean(winter_T_monthly(z["ssp585|tas_ym"], z["ssp585|tas"], list(range(2071, 2101))), axis=0)
        (kc, sc, Gc), (kf, sf, Gf) = S["historical"], S["ssp585"]
        out["forced"] = forced(kc, sc, Gc, kf, sf, Gf, Tc, Tf, rng)
        # consistency of the 1951-2014 extraction with the 1951-1980 inputs of the forced test
        cm_, i1, i2 = np.intersect1d(common, kc, return_indices=True)
        out["check_max_abs_diff"] = dict(s=float(np.nanmax(np.abs(s[i1] - sc[i2]))), G=float(np.nanmax(np.abs(G[i1] - Gc[i2]))))
        out["wind_only"] = bool(inv.get(name, {}).get("wind_only"))
        res[name] = out
        p, f = out["interannual"]["primary"], out["forced"]
        print(f"  {name:18s} interannual {p['eu']:+.2f} %/K [{p['ci90'][0]:+.2f}, {p['ci90'][1]:+.2f}]  forced within-G {f['within_G']:+.2f}% "
              f"dT {f['dT_region']:.2f} K -> {f['within_per_K']:+.2f} %/K  check {out['check_max_abs_diff']}", flush=True)
    return res


# ------------------------------------------------------------------ observations
def stations_sel(z, var, winters, min_frac=0.8, min_winters=30):
    day = z["day"]
    k = gs.ymd_of(day)
    w = gs.winter_of(k)
    x = z[var].astype(np.float64)
    good = np.zeros((x.shape[0], len(winters)), bool)
    for j, y in enumerate(winters):
        m = w == y
        good[:, j] = np.isfinite(x[:, m]).mean(axis=1) >= min_frac
    keep = good.sum(axis=1) >= min_winters
    # blank the station's bad winters
    wi = np.searchsorted(winters, w)
    bad = ~good[:, np.clip(wi, 0, len(winters) - 1)]
    x = np.where(bad, np.nan, x)
    return x, keep, good


def obs():
    z = np.load(f"{D}/obs/hadisd_daily.npz")
    sel = json.load(open(f"{D}/obs/hadisd_selection.json"))
    ids = list(z["ids"])
    winters = list(range(1980, 2026))         # ERA5 msl of data/regimes starts in January 1979
    # ERA5 G on the station days
    m = np.load(f"{xm.OUT}/msl_era5_-.npz")
    km = gs.ymd_of(m["dates"].astype("datetime64[D]"))
    k_all = gs.ymd_of(z["day"])
    k, i_m, i_z = np.intersect1d(km, k_all, return_indices=True)
    keep = np.isin(gs.winter_of(k), winters)
    k, i_m, i_z = k[keep], i_m[keep], i_z[keep]
    z = dict(day=z["day"][i_z], ids=z["ids"], **{v_: z[v_][:, i_z] for v_ in ("ws", "ws3h", "ta")})
    mon = month_of(k)
    G, u, v = G_uv(m["msl"][i_m], m["lat"], m["lon"], False)
    e5 = np.load(f"{D}/regimes/interannual_era5_t.npz")
    T_e5 = winter_T_monthly(e5["ym"], e5["t2m"], winters)
    INV = winter_T_monthly(e5["ym"], e5["t850"], winters) - T_e5
    res = dict(source="HadISD v3.4.3.2025f (Dunn et al. 2012, 2014, 2016, 2019; Smith et al. 2011)", winters=[1980, 2025])
    variants = {}
    fe = f"{D}/regimes/interannual_era5_wind.npz"
    E = None
    if os.path.exists(fe):                    # ERA5 10 m wind: region means and at the station points
        E = np.load(fe)
        assert list(E["ids"]) == ids
        ke = gs.ymd_of(E["day"].astype("datetime64[D]"))
        _, ia, ib = np.intersect1d(k, ke, return_indices=True)
        assert len(ia) == len(k)
        E_st, E_reg = E["s_station"][ib].T.astype(np.float64), E["s_region"][ib].astype(np.float64)
    for wvar in ("ws", "ws3h"):
        ws, keep_w, _ = stations_sel(z, wvar, winters)
        ta, keep_t, _ = stations_sel(z, "ta", winters)
        S = np.full((len(k), len(EU)), np.nan)
        Tst = np.full((len(winters), len(EU)), np.nan)
        yB = np.full((len(winters), len(EU)), np.nan)
        nst = {}
        SE_m, SE_a = np.full_like(S, np.nan), np.full_like(S, np.nan)

        def index(x, n):
            """normalised speed per station (divided by its calendar-month mean), mean over the stations"""
            x = x.copy()
            for mo in rg.MONTHS:
                mm = mon == mo
                x[:, mm] = x[:, mm] / np.nanmean(x[:, mm], axis=1, keepdims=True)
            return x, np.where(np.isfinite(x).sum(axis=0) >= 0.5 * n, np.nanmean(x, axis=0), np.nan)

        for r, reg in enumerate(EU):
            idx = [ids.index(i) for i in sel["regions"][reg] if keep_w[ids.index(i)]]
            nst[reg] = len(idx)
            x, S[:, r] = index(ws[idx], len(idx))
            if E is not None and wvar == "ws":
                SE_m[:, r] = index(np.where(np.isfinite(ws[idx]), E_st[idx], np.nan), len(idx))[1]
                SE_a[:, r] = index(E_st[idx], len(idx))[1]
            # station temperature: mean of per-station winter anomalies
            it = [ids.index(i) for i in sel["regions"][reg] if keep_t[ids.index(i)]]
            Tw = winter_T_daily(k, ta[it].T, winters, 1)            # (W, n) winter means of valid days
            wcov = np.stack([np.isfinite(ta[it][:, gs.winter_of(k) == y]).mean(axis=1) for y in winters])
            Tw = np.where(wcov >= 0.8, Tw, np.nan)
            Ta = Tw - np.nanmean(Tw, axis=0, keepdims=True)
            Tst[:, r] = np.where(np.isfinite(Ta).sum(axis=1) >= 0.5 * len(it), np.nanmean(Ta, axis=1), np.nan)
            # variant B: per-station y_q5 (station days, region G), detrended per station, then the mean
            ys = []
            for i in range(len(idx)):
                tb = winter_table(k, x[i][:, None], G[:, [r]], None, None, winters, int(0.8 * 200))
                y = detrend(tb["y_q5"][:, 0], np.asarray(winters, float), 1)
                ys.append(y)
            ys = np.array(ys)
            yB[:, r] = np.where(np.isfinite(ys).sum(axis=0) >= 0.5 * len(idx), np.nanmean(ys, axis=0), np.nan)
        tab = winter_table(k, S, G, u, v, winters, 150)
        v_ = dict(n_stations=nst, n_stations_total=int(len({i for reg in EU for i in sel["regions"][reg] if keep_w[ids.index(i)]})))
        v_["T_station"] = all_sens(tab, Tst, winters)
        v_["T_era5"] = all_sens(tab, T_e5, winters)
        v_["inversion_era5"] = sensitivity(tab, Tst, winters, "y_q5", ("lnG",), 1, Xextra=prepare(tab, INV, winters, "y_q5", (), 1)[1][:, :, 0])
        # T and inversion together: y' = a + b T' + c lnG' + d INV'
        y_, X_ = prepare(tab, Tst, winters, "y_q5", ("lnG",), 1)
        Xi = prepare(tab, INV, winters, "y_q5", (), 1)[1][:, :, :1]
        XX = np.concatenate([X_, Xi], axis=2)
        v_["T_station_and_inversion"] = dict(b_T=float(np.nanmean([ols_b(y_[:, r], XX[:, r])[0] for r in range(len(EU))]) * 100),
                                             b_inv=float(np.nanmean([ols_b(y_[:, r], XX[:, r])[2] for r in range(len(EU))]) * 100))
        v_["corr_Tstation_Tera5"] = float(np.nanmean([np.corrcoef(*[detrend(a[:, r], np.asarray(winters, float), 1) for a in (Tst, T_e5)])[0, 1] for r in range(len(EU))]))
        v_["corr_T_inversion"] = float(np.nanmean([np.corrcoef(*[detrend(a[:, r], np.asarray(winters, float), 1) for a in (T_e5, INV)])[0, 1] for r in range(len(EU))]))
        tabB = dict(tab, y_q5=yB)
        v_["per_station_B"] = sensitivity(tabB, Tst, winters, "y_q5", ("lnG",), 1)
        v_["per_station_B_era5T"] = sensitivity(tabB, T_e5, winters, "y_q5", ("lnG",), 1)
        if E is not None and wvar == "ws":
            chk = {}
            for nm, SS in (("era5_at_stations_station_days", SE_m), ("era5_at_stations_all_days", SE_a), ("era5_region_mean", E_reg)):
                tb = winter_table(k, SS, G, u, v, winters, 150)
                chk[nm] = dict(T_era5=sensitivity(tb, T_e5, winters), T_station=sensitivity(tb, Tst, winters),
                               q5_dir_T_era5=sensitivity(tb, T_e5, winters, "y_q5", ("lnG", "ug", "vg"), 1),
                               first_diff_T_era5=sensitivity(tb, T_e5, winters, "y_q5", ("lnG",), "diff", nboot=200))
                yo, yv = tab["y_q5"], tb["y_q5"]
                chk[nm]["corr_y_with_stations"] = float(np.nanmean([np.corrcoef(*[detrend(a[:, r], np.asarray(winters, float), 1) for a in (yo, yv)])[0, 1] for r in range(len(EU))]))
                if nm == "era5_at_stations_station_days":       # paired: stations minus ERA5 at the same points and days
                    td = dict(tb, y_q5=yo - yv)
                    chk["stations_minus_era5_at_stations"] = dict(T_era5=sensitivity(td, T_e5, winters), T_station=sensitivity(td, Tst, winters))
                c_ = chk[nm]
                print(f"    ERA5 check {nm:32s} T_era5 {c_['T_era5']['eu']:+.2f} [{c_['T_era5']['ci90'][0]:+.2f}, {c_['T_era5']['ci90'][1]:+.2f}]  "
                      f"T_station {c_['T_station']['eu']:+.2f}  dir {c_['q5_dir_T_era5']['eu']:+.2f}  r(y, stations) {c_['corr_y_with_stations']:.2f}", flush=True)
            v_["era5_wind_check"] = chk
            d_ = chk["stations_minus_era5_at_stations"]["T_era5"]
            print(f"    stations minus ERA5 at stations (paired), T_era5: {d_['eu']:+.2f} [{d_['ci90'][0]:+.2f}, {d_['ci90'][1]:+.2f}]", flush=True)
        variants[wvar] = v_
        p = v_["T_station"]["primary"]
        print(f"  obs {wvar}: {v_['n_stations_total']} stations {nst}; station T {p['eu']:+.2f} %/K [{p['ci90'][0]:+.2f}, {p['ci90'][1]:+.2f}], "
              f"ERA5 T {v_['T_era5']['primary']['eu']:+.2f} [{v_['T_era5']['primary']['ci90'][0]:+.2f}, {v_['T_era5']['primary']['ci90'][1]:+.2f}], "
              f"B {v_['per_station_B']['eu']:+.2f} [{v_['per_station_B']['ci90'][0]:+.2f}, {v_['per_station_B']['ci90'][1]:+.2f}], "
              f"inv {v_['inversion_era5']['eu']:+.2f} %/K-of-inversion", flush=True)
    res.update(variants)
    return res


# ------------------------------------------------------------------ cross-model
def cross(res, key="primary"):
    names = list(res["cmip6"])
    x = np.array([res["cmip6"][n]["interannual"][key]["eu"] for n in names])
    y = np.array([res["cmip6"][n]["forced"]["within_per_K"] for n in names])
    se = np.array([res["cmip6"][n]["interannual"][key]["se"] for n in names])
    out = {}
    for tag, sel_ in (("cmip6_21", np.ones(len(names), bool)), ("cmip6_20_noECE", np.array([n != "EC-Earth3" for n in names]))):
        xs, ys = x[sel_], y[sel_]
        r = np.corrcoef(xs, ys)[0, 1]
        sl, ic = np.polyfit(xs, ys, 1)
        bs_r, bs_s = [], []
        for _ in range(2000):
            i = RNG.integers(0, len(xs), len(xs))
            if np.std(xs[i]) > 0:
                bs_r.append(np.corrcoef(xs[i], ys[i])[0, 1])
                bs_s.append(np.polyfit(xs[i], ys[i], 1)[0])
        perm = np.array([abs(np.corrcoef(RNG.permutation(xs), ys)[0, 1]) for _ in range(5000)])
        out[tag] = dict(n=int(sel_.sum()), r=float(r), slope=float(sl), intercept=float(ic),
                        r_ci90=[float(v) for v in np.percentile(bs_r, [5, 95])], slope_ci90=[float(v) for v in np.percentile(bs_s, [5, 95])],
                        p_perm=float((perm >= abs(r)).mean()),
                        spearman=float(np.corrcoef(np.argsort(np.argsort(xs)), np.argsort(np.argsort(ys)))[0, 1]),
                        x_spread_sd=float(np.std(xs)), mean_sampling_se=float(np.mean(se[sel_])))
    # pooled over models x regions (region forced change per K of that region's warming), and per region
    px = np.array([[res["cmip6"][n]["interannual"][key]["per_region"][r] for r in EU] for n in names])
    py = np.array([[res["cmip6"][n]["forced"]["per_region"][r]["within_G"] / res["cmip6"][n]["forced"]["per_region"][r]["dT"] for r in EU] for n in names])
    ok = np.isfinite(px) & np.isfinite(py)
    out["pooled_model_region"] = dict(n=int(ok.sum()), r=float(np.corrcoef(px[ok], py[ok])[0, 1]),
                                      per_region_r={r: float(np.corrcoef(px[:, j], py[:, j])[0, 1]) for j, r in enumerate(EU)})
    ox = [res["ours"][n]["interannual"][key]["eu"] for n in res["ours"]]
    oy = [res["ours"][n]["forced"]["within_per_K"] for n in res["ours"]]
    xa, ya = np.concatenate([x, ox]), np.concatenate([y, oy])
    out["all_25"] = dict(r=float(np.corrcoef(xa, ya)[0, 1]), slope=float(np.polyfit(xa, ya, 1)[0]), intercept=float(np.polyfit(xa, ya, 1)[1]))
    return out


def plot(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for a, key in zip(ax, ("primary", "q5_dir")):
        for n, v in res["cmip6"].items():
            p = v["interannual"][key]
            c = "#c0504d" if n == "EC-Earth3" else "#7f7f7f"
            a.errorbar(p["eu"], v["forced"]["within_per_K"], xerr=[[p["eu"] - p["ci90"][0]], [p["ci90"][1] - p["eu"]]], fmt="o", color=c, ms=4, lw=0.6)
            a.annotate(n, (p["eu"], v["forced"]["within_per_K"]), fontsize=6, color=c)
        for n, v in res["ours"].items():
            p = v["interannual"][key]
            a.errorbar(p["eu"], v["forced"]["within_per_K"], xerr=[[p["eu"] - p["ci90"][0]], [p["ci90"][1] - p["eu"]]], fmt="s", color="#3b6ea8", ms=6, lw=0.8)
            a.annotate(n, (p["eu"], v["forced"]["within_per_K"]), fontsize=7, color="#3b6ea8")
        o = res["obs"]["ws"]["T_station"][key]
        a.axvspan(o["ci90"][0], o["ci90"][1], color="#9bbb59", alpha=0.25, lw=0)
        a.axvline(o["eu"], color="#4f6228", lw=1.2, label=f"HadISD, station T: {o['eu']:+.2f} %/K")
        oe = res["obs"]["ws"]["T_era5"][key]
        a.axvline(oe["eu"], color="#4f6228", lw=1.0, ls="--", label=f"HadISD, ERA5 T: {oe['eu']:+.2f} %/K")
        ck = res["obs"]["ws"].get("era5_wind_check")
        if ck:
            kk = "T_era5" if key == "primary" else "q5_dir_T_era5"
            a.axvline(ck["era5_region_mean"][kk]["eu"], color="#8064a2", lw=1.0, ls=":", label=f"ERA5 wind, region mean: {ck['era5_region_mean'][kk]['eu']:+.2f} %/K")
            a.axvline(ck["era5_at_stations_station_days"][kk]["eu"], color="#8064a2", lw=1.0, ls="-.", label=f"ERA5 wind at the stations: {ck['era5_at_stations_station_days'][kk]['eu']:+.2f} %/K")
        cr = res["cross"][key]["cmip6_21"]
        xx = np.array(a.get_xlim())
        a.plot(xx, cr["intercept"] + cr["slope"] * xx, color="k", lw=0.8, label=f"CMIP6 fit r = {cr['r']:+.2f}")
        a.axhline(0, color="k", lw=0.4)
        a.axvline(0, color="k", lw=0.4)
        a.set_xlabel("interannual sensitivity of the transfer at fixed G (% per K of winter T)")
        a.set_ylabel("forced transfer change at fixed G (% per K of regional warming)")
        a.set_title({"primary": "G-quintile standardised, control ln G", "q5_dir": "plus winter-mean geostrophic direction"}[key], fontsize=10)
        a.legend(fontsize=7, loc="lower left")
    fig.suptitle("8 European regions, Sep-Mar")
    fig.savefig(f"{OUT}/interannual_transfer_scatter.png", dpi=120)
    plt.close(fig)


def main():
    res = dict(meta=dict(regions=EU, n_boot=NBOOT, primary="y_q5 on T and lnG, linear detrend, % per K",
                         note="forced = within-G-quintile transfer change per K of regional winter land warming"))
    print("Observations", flush=True)
    res["obs"] = obs()
    print("Our runs", flush=True)
    res["ours"] = ours()
    print("CMIP6", flush=True)
    res["cmip6"] = cmip6()
    res["cross"] = {k: cross(res, k) for k in ("primary", "q5_dir", "res", "dir_daily", "detrend_quadratic", "first_diff", "q5_noG")}
    ens = [v for v in res["cmip6"].values() if not v["wind_only"]]
    res["cmip6_ensemble"] = {k: dict(median=float(np.median([f(v) for v in ens])), p10_p90=[float(x) for x in np.percentile([f(v) for v in ens], [10, 90])])
                             for k, f in (("interannual_pct_per_K", lambda v: v["interannual"]["primary"]["eu"]),
                                          ("interannual_dir_pct_per_K", lambda v: v["interannual"]["q5_dir"]["eu"]),
                                          ("forced_within_G_pct", lambda v: v["forced"]["within_G"]),
                                          ("forced_within_per_K", lambda v: v["forced"]["within_per_K"]),
                                          ("dT_region", lambda v: v["forced"]["dT_region"]))}
    # implied forced change from the observed sensitivity
    imp = {}
    for key in ("primary", "q5_dir"):
        cr = res["cross"][key]["cmip6_21"]
        for tv in ("T_station", "T_era5"):
            o = res["obs"]["ws"][tv][key]
            imp[f"{key}|{tv}"] = dict(obs=o["eu"], obs_ci90=o["ci90"], implied_forced_per_K=cr["intercept"] + cr["slope"] * o["eu"],
                                      implied_range=sorted([cr["intercept"] + cr["slope"] * o["ci90"][0], cr["intercept"] + cr["slope"] * o["ci90"][1]]))
    res["implied"] = imp
    json.dump(res, open(f"{OUT}/interannual_transfer.json", "w"), indent=1)
    print("\nCross-model (CMIP6 21): " + "; ".join(f"{k} r {v['cmip6_21']['r']:+.2f} slope {v['cmip6_21']['slope']:+.2f} p {v['cmip6_21']['p_perm']:.3f}" for k, v in res["cross"].items()))
    print("CMIP6 ensemble: " + json.dumps(res["cmip6_ensemble"]))
    print("Implied: " + json.dumps(imp))
    plot(res)
    print(f"Saved {OUT}/interannual_transfer.json")


if __name__ == "__main__":
    main()
