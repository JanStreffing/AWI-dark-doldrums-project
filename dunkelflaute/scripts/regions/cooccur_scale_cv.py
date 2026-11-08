"""
Why does the per-K increase in low-wind events grow with aggregation scale in CMIP6 (methods_notes
15), and does it hold for Europe? Test of the variance explanation: aggregation lowers the
variability of the series relative to its mean, so the same relative decline of the mean moves more
of the distribution below a threshold fixed at the control's base rate.

Series: the 5 x 5 deg land boxes of REGION_SET=grid5 (wind only, onshore CF of the log-law 100 m
wind), combined as in grid_per_K.py into 10, 20 and 40 deg tiles (cos-lat weights, a tile needs
2/4/8 land boxes, a step half the tile weight) and into Europe (35-70 N, 10 W-30 E, boxes whose
centre lies inside). Two families of aggregates:
  NH       all 5 deg boxes, all 10/20/40 deg tiles (the scale series of grid_per_K.py)
  Europe   the 5/10/20 deg aggregates whose weighted centre lies in the Europe box, and Europe itself

Per aggregate (control and future, Sep-Mar winters of the comparison):
  CV       std / mean of the trailing 48 h mean (the series that is thresholded); also p5 / mean
  thr      matched to TARGET = 5 events per control winter (200 levels, 0.1-40th percentile of the
           control, as grid_per_K.box_change)
  ratio    future / control events at thr (all)
  mr       mean removed: future x mean_ctl / mean_fut (as grid_per_K)
  mo       mean only: control x mean_fut / mean_ctl, i.e. only the mean change applied to the control
  elast    d ln(events) / d ln(thr) at thr in the control, from thr x (1 -/+ 0.03): the sensitivity
           of the event count to a proportional shift, so ln(mo) ~ -elast x ln(mean_fut / mean_ctl)
Per scale: weighted means (weights = tile cos-lat land weight) over the valid aggregates (control
reaches 0.9 TARGET) of ln ratio, ln mr, ln mo, CV, p5/mean, elast and the mean change; per K =
100 (exp(mean ln ratio / dT) - 1).

Comparisons: our runs on the 6-hourly footing (grid_per_K.OURS), the 20 CMIP6 models (historical
1951-1980 against ssp585 2071-2100) and the EERIE / TCO95 bridges on the daily footing; ERA5
winters 1980-2026 (control statistics only) on both footings.

Usage: python3 cooccur_scale_cv.py      (module python3; reads ~20 GB, compute node, NPROC workers)
Output: dunkelflaute/results/regions/cooccur_scale_cv.json
"""
import glob
import json
import os
import sys
from multiprocessing import Pool

os.environ["REGION_SET"] = "grid5"
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("regions", "analysis", "cmip6"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/{p}")
import regions as rg                                          # noqa: E402
import region_events as ev                                    # noqa: E402
import region_per_K as rpk                                    # noqa: E402
import grid_per_K as gp                                       # noqa: E402
from robustness_mockert_events import trailing_mean           # noqa: E402

OUT = f"{_ROOT}/dunkelflaute/results/regions"
ERA5 = f"{_ROOT}/dunkelflaute/data/regions/cooccur_era5_wind.npz"
TARGET = 5.0
NGRID = 200
EPS = 0.03
EUROPE = (35.0, 70.0, -10.0, 30.0)
SCALES = ("5deg", "10deg", "20deg", "40deg")
ESCALES = ("5deg", "10deg", "20deg", "Europe")


def families():
    """{family: {scale: (labels, A (nAgg, nB))}}"""
    agg = gp.aggregations()
    la = np.array([rg.BOXES[b][0] for b in gp.B]) + 2.5
    lo = np.array([rg.BOXES[b][2] for b in gp.B]) + 2.5
    nh = {s: agg[s] for s in SCALES}
    eu = {}
    for s in ("5deg", "10deg", "20deg"):
        labs, A = agg[s]
        cla = (A @ la) / A.sum(1)
        clo = (A @ lo) / A.sum(1)
        m = (cla >= EUROPE[0]) & (cla <= EUROPE[1]) & (clo >= EUROPE[2]) & (clo <= EUROPE[3])
        eu[s] = ([l for l, k in zip(labs, m) if k], A[m])
    labs, A = agg["named"]
    eu["Europe"] = (["Europe"], A[[labs.index("Europe")]])
    return dict(NH=nh, Europe=eu)


def roll(ws, win):
    return np.concatenate([trailing_mean(w, win)[win - 1:] for w in ws])


def cnt(ws, thr, win, f=1.0):
    ev.WIN = win
    return np.array([ev.counts_grid(w * f, np.atleast_1d(thr)) for w in ws], float).mean(0)


def one(xc, xf, win):
    """Control / future winter lists of one aggregate -> dict of statistics (NaN if not usable)."""
    out = {}
    if any(not np.isfinite(w).all() for w in xc):
        return None
    rc = roll(xc, win)
    mc = float(np.mean(np.concatenate(xc)))
    out.update(cv_ctl=float(rc.std() / rc.mean()), p5_ctl=float(np.percentile(rc, 5) / rc.mean()), mean_ctl=mc)
    allc = np.concatenate(xc)
    grid = np.linspace(*np.percentile(allc, [0.1, 40]), NGRID)
    ev.WIN = win
    cc = np.stack([ev.counts_grid(w, grid) for w in xc]).mean(0)
    k = int(np.argmax(cc >= TARGET)) if (cc >= TARGET).any() else int(np.argmax(cc))
    thr = grid[k]
    a = cc[k]
    out.update(thr=float(thr), rate=float(a), valid=bool(a >= 0.9 * TARGET))
    lo, hi = cnt(xc, thr * (1 - EPS), win)[0], cnt(xc, thr * (1 + EPS), win)[0]
    out["elast"] = float(np.log(hi / lo) / np.log((1 + EPS) / (1 - EPS))) if lo > 0 and hi > 0 else np.nan
    if xf is None:
        return out
    if any(not np.isfinite(w).all() for w in xf):
        return None
    rf = roll(xf, win)
    mf = float(np.mean(np.concatenate(xf)))
    out.update(cv_fut=float(rf.std() / rf.mean()), p5_fut=float(np.percentile(rf, 5) / rf.mean()), mean_fut=mf,
               mean_pct=float(100 * (mf / mc - 1)))
    b = cnt(xf, thr, win)[0]
    bm = cnt(xf, thr, win, mc / mf)[0]
    am = cnt(xc, thr, win, mf / mc)[0]
    with np.errstate(divide="ignore"):
        out.update(lr=float(np.log(b / a)) if a > 0 else np.nan, lr_mr=float(np.log(bm / a)) if a > 0 else np.nan,
                   lr_mo=float(np.log(am / a)) if a > 0 else np.nan)
    return out


def scale_stats(Xc, Xf, mc, mf, win, dT, fam):
    """Per family and scale: weighted means over valid aggregates."""
    res = {}
    for fname, scales in fam.items():
        res[fname] = {}
        for s, (labs, A) in scales.items():
            Yc = gp.aggregate(Xc, A)
            Yf = gp.aggregate(Xf, A) if Xf is not None else None
            rows, wts = [], []
            for j in range(A.shape[0]):
                o = one([Yc[m, j] for m in mc], None if Yf is None else [Yf[m, j] for m in mf], win)
                if o is None or not o["valid"]:
                    continue
                if Yf is not None and not np.isfinite(o["lr"]):
                    continue
                rows.append(o)
                wts.append(A[j].sum())
            w = np.array(wts)
            keys = ["cv_ctl", "p5_ctl", "elast"] + (["cv_fut", "p5_fut", "mean_pct", "lr", "lr_mr", "lr_mo"] if Yf is not None else [])
            st = dict(n=len(rows))
            for k in keys:
                v = np.array([r[k] for r in rows])
                ok = np.isfinite(v)
                st[k] = float((v[ok] * w[ok]).sum() / w[ok].sum()) if ok.any() else np.nan
            if Yf is not None and rows:
                for k, t in (("lr", "perK"), ("lr_mr", "perK_mr"), ("lr_mo", "perK_mo")):
                    st[t] = float(100 * (np.exp(st[k] / dT) - 1))
                v = np.array([r["lr_mo"] for r in rows])
                pred = np.array([-r["elast"] * np.log(1 + r["mean_pct"] / 100) for r in rows])
                ok = np.isfinite(v) & np.isfinite(pred)
                st["lr_mo_pred"] = float((pred[ok] * w[ok]).sum() / w[ok].sum()) if ok.any() else np.nan
                st["corr_mo_pred"] = float(np.corrcoef(v[ok], pred[ok])[0, 1]) if ok.sum() > 2 else np.nan
            res[fname][s] = st
    return res


FAM = None


def job(spec):
    kind = spec[0]
    if kind == "6h":
        _, name = spec
        (cn, c0, c1), (fn, f0, f1) = ev.COMPARISONS[name]
        dc, df = ev.load(cn), ev.load(fn)
        mc = [m for _, m in ev.winters(dc, c0, c1)]
        mf = [m for _, m in ev.winters(df, f0, f1)]
        T = rpk.global_t()
        dT = rpk.mean_years(T, fn, f0, f1)[0] - rpk.mean_years(T, cn, c0, c1)[0]
        r = scale_stats(gp.metric6(dc, "wind"), gp.metric6(df, "wind"), mc, mf, 8, dT, FAM)
        ev._CACHE.clear()
        return f"6h|{name}", dict(dT=dT, n_winters=[len(mc), len(mf)], scales=r)
    if kind == "day":
        _, name, fc, ec, c0, c1, ff, ef, f0, f1, Tc, Tf = spec
        dc, df = gp.load_daily(fc, ec), gp.load_daily(ff, ef)
        if Tc is None:
            Tc = dict(zip(dc["tas_years"].tolist(), dc["tas"].tolist()))
            Tf = dict(zip(df["tas_years"].tolist(), df["tas"].tolist()))
        dT = float(np.mean([Tf[y] for y in range(f0, f1 + 1) if y in Tf]) - np.mean([Tc[y] for y in range(c0, c1 + 1) if y in Tc]))
        mc, mf = gp.dmasks(dc, c0, c1), gp.dmasks(df, f0, f1)
        r = scale_stats(gp.daily_metric(dc, "wind"), gp.daily_metric(df, "wind"), mc, mf, 2, dT, FAM)
        return f"day|{name}", dict(dT=dT, n_winters=[len(mc), len(mf)], scales=r)
    if kind == "era5":
        z = np.load(ERA5)
        names = list(z["regions"])
        ix = [names.index(b) for b in gp.B]
        out = {}
        for tag, tk, xk, win, nmin in (("6h", "time6", "wind6", 8, ev.MIN_STEPS), ("day", "day", "wind_day", 2, 200)):
            t = z[tk]
            if tag == "6h":
                wint = ev.winter_of(t.astype("datetime64[s]"))
            else:
                d = t.astype("datetime64[D]")
                m = d.astype("datetime64[M]").astype(int) % 12 + 1
                y = d.astype("datetime64[Y]").astype(int) + 1970
                wint = np.where(m >= 9, y + 1, y)
            X = z[xk][:, ix].astype(float)
            ms = [wint == yy for yy in range(1980, 2027) if (wint == yy).sum() >= nmin]
            out[tag] = scale_stats(X, None, ms, None, win, np.nan, FAM)
        return "era5", dict(scales=out)


def main():
    global FAM
    FAM = families()
    specs = [("6h", n) for n in gp.OURS]
    DC = gp.DC
    gt = {}
    for mod in ("eerie", "tco95"):
        gt.update({k: {int(y): t for y, t in v.items()} for k, v in json.load(open(f"{rg.BASE_DIR}/global_t2m_{mod}.json")).items()})
    for fn in sorted(glob.glob(f"{DC}/daily_cmip6_*.npz")):
        m = os.path.basename(fn)[len("daily_cmip6_"):-4]
        specs.append(("day", f"cmip6|{m}", fn, "historical", 1951, 1980, fn, "ssp585", 2071, 2100, None, None))
    E = lambda r: f"{DC}/daily_eerie_{r}.npz"          # noqa: E731
    T95 = lambda r: f"{DC}/daily_tco95_{r}.npz"        # noqa: E731
    specs.append(("day", "bridge|EERIE ctl->ssp245", E("ctl"), "run", 1972, 2050, E("ssp245"), "run", 2041, 2070, gt["eerie_ctl"], gt["eerie_ssp245"]))
    specs.append(("day", "bridge|EERIE hist->ssp245", E("hist"), "run", 1971, 2000, E("ssp245"), "run", 2041, 2070, gt["eerie_hist"], gt["eerie_ssp245"]))
    specs.append(("day", "bridge|TCO95 ctl->hi", T95("ctl1950g"), "run", 1953, 2099, T95("hi1950g"), "run", 2081, 2099, gt["tco95_ctl1950g"], gt["tco95_hi1950g"]))
    if os.path.exists(ERA5):
        specs.append(("era5",))
    # big 6-hourly jobs first
    with Pool(int(os.environ.get("NPROC", 16))) as pool:
        res = dict(pool.imap_unordered(job, specs))
    out = dict(target=TARGET, scales_NH=list(SCALES), scales_Europe=list(ESCALES), comparisons=res,
               n_aggregates={f: {s: len(v[0]) for s, v in sc.items()} for f, sc in FAM.items()})
    # CMIP6 ensemble per family and scale
    names = [n for n in res if n.startswith("day|cmip6|")]
    ens = {}
    for f, scl in (("NH", SCALES), ("Europe", ESCALES)):
        ens[f] = {}
        for s in scl:
            e = {}
            for k in ("perK", "perK_mr", "perK_mo", "cv_ctl", "cv_fut", "p5_ctl", "p5_fut", "elast", "mean_pct", "lr_mo", "lr_mo_pred"):
                v = np.array([res[n]["scales"][f][s].get(k, np.nan) for n in names])
                e[k] = dict(median=float(np.nanmedian(v)), p10=float(np.nanpercentile(v, 10)), p90=float(np.nanpercentile(v, 90)))
            ens[f][s] = e
    out["cmip6_ensemble"] = ens
    # across models and scales: does the mean-only part follow elast x mean decline?
    xs, ys = [], []
    for n in names:
        for f, scl in (("NH", SCALES),):
            for s in scl:
                st = res[n]["scales"][f][s]
                xs.append(st["lr_mo_pred"])
                ys.append(st["lr_mo"])
    out["cmip6_mo_vs_pred_corr"] = float(np.corrcoef(xs, ys)[0, 1])
    json.dump(out, open(f"{OUT}/cooccur_scale_cv.json", "w"), indent=1)

    for f, scl in (("NH", SCALES), ("Europe", ESCALES)):
        print(f"\n=== {f}: per K all [mean removed] {{mean only}}, CV ctl->fut, p5/mean ctl, elasticity, mean change")
        for n in sorted(res):
            if n == "era5" or n.startswith("day|cmip6|"):
                continue
            sc = res[n]["scales"][f]
            print(f"  {n:30s} " + "  ".join(f"{s} {sc[s]['perK']:+.1f} [{sc[s]['perK_mr']:+.1f}] {{{sc[s]['perK_mo']:+.1f}}} cv {sc[s]['cv_ctl']:.2f}->{sc[s]['cv_fut']:.2f} e {sc[s]['elast']:.1f} dm {sc[s]['mean_pct']:+.1f}"
                                             for s in scl if sc[s]["n"]))
        for s in scl:
            e = ens[f][s]
            print(f"  CMIP6 median {s:7s} perK {e['perK']['median']:+.1f} ({e['perK']['p10']:+.1f}..{e['perK']['p90']:+.1f}) mr {e['perK_mr']['median']:+.1f} mo {e['perK_mo']['median']:+.1f} "
                  f"cv {e['cv_ctl']['median']:.3f}->{e['cv_fut']['median']:.3f} p5/mean {e['p5_ctl']['median']:.3f}->{e['p5_fut']['median']:.3f} "
                  f"elast {e['elast']['median']:.2f} dmean {e['mean_pct']['median']:+.1f}% lr_mo {e['lr_mo']['median']:.3f} pred {e['lr_mo_pred']['median']:.3f}")
        if "era5" in res:
            for tag in ("6h", "day"):
                sc = res["era5"]["scales"][tag][f]
                print(f"  ERA5 {tag:4s} " + "  ".join(f"{s} cv {sc[s]['cv_ctl']:.3f} p5 {sc[s]['p5_ctl']:.3f} e {sc[s]['elast']:.2f} (n {sc[s]['n']})" for s in scl))
    print(f"\ncorr(mean-only ln ratio, -elast x ln(mean ratio)) across CMIP6 models x NH scales: {out['cmip6_mo_vs_pred_corr']:.2f}")
    print(f"Saved {OUT}/cooccur_scale_cv.json")


if __name__ == "__main__":
    main()
