"""
Gridded per-K maps of the change in winter low-output events: 5 x 5 deg land boxes over 20-70 N
(regions.py, REGION_SET=grid5), for our runs on their native 6-hourly footing and for CMIP6 (plus the
EERIE / TCO95 daily bridges) on the daily footing.

Per box, comparison and metric (wind = onshore CF of the log-law 100 m wind; comb = 0.44 solar +
0.56 wind): threshold matched so the comparison's control has TARGET events per winter; ratio of
events per winter future / control; per K = 100 (exp(ln ratio / dT) - 1); mean-removed ratio (future
series rescaled to the control mean, same threshold); mean change of the metric. A box is valid if
its control reaches 0.9 TARGET.

Chance range per box, from the TCO95 constant-1950 control (147 winters): all placements of two
disjoint blocks with the comparison's own winter counts (capped at 147 in total), 5-95 % of ln ratio.
For CMIP6: the share of models with ln ratio > 0, and the multi-model median of the per-K change.

Footings: 6-hourly (WIN 8) for data/regions_grid5/*.npz, daily (WIN 2) for data/cmip6_grid5/*.npz.

Usage: python3 grid_per_K.py     (module python3 or legacy; reads ~30 GB, run on a compute node)
Output: dunkelflaute/results/maps/grid_per_K.npz, grid_per_K_summary.json
"""
import glob
import json
import os
import sys

os.environ["REGION_SET"] = "grid5"
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("regions", "analysis", "cmip6"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/{p}")
import regions as rg            # noqa: E402
import region_events as ev      # noqa: E402
import region_per_K as rpk      # noqa: E402

OUT = f"{_ROOT}/dunkelflaute/results/maps"
DC = f"{_ROOT}/dunkelflaute/data/cmip6_grid5"
TARGET = 5.0
NGRID = 200
W_SOLAR = 0.44
KINDS = ("wind", "comb")
B = [r for r in rg.REGIONS if r != "DE"]
BI = [rg.REGIONS.index(r) for r in B]
OURS = ["9km", "31km", "100km", "eerie_ctl", "eerie_hist"]


def box_change(Xc, Xf, mc, mf, win):
    """Xc, Xf (T, nB); mc, mf lists of winter masks. Returns dict of (nB,) arrays."""
    ev.WIN = win
    nB = Xc.shape[1]
    out = {k: np.full(nB, np.nan) for k in ("thr", "rate", "ctl", "fut", "ratio", "ratio_mr", "mean_pct")}
    for j in range(nB):
        xc, xf = Xc[:, j], Xf[:, j]
        if not np.isfinite(xc).all() or not np.isfinite(xf).all():
            continue
        grid = np.linspace(*np.percentile(xc, [0.1, 40]), NGRID)
        cc = np.stack([ev.counts_grid(xc[m], grid) for m in mc]).mean(0)
        k = int(np.argmax(cc >= TARGET)) if (cc >= TARGET).any() else int(np.argmax(cc))
        thr = np.array([grid[k]])
        a = np.array([ev.counts_grid(xc[m], thr)[0] for m in mc], float)
        b = np.array([ev.counts_grid(xf[m], thr)[0] for m in mf], float)
        mcm = np.mean(np.concatenate([xc[m] for m in mc]))
        mfm = np.mean(np.concatenate([xf[m] for m in mf]))
        bm = np.array([ev.counts_grid(xf[m] * mcm / mfm, thr)[0] for m in mf], float)
        out["thr"][j], out["rate"][j] = grid[k], cc[k]
        out["ctl"][j], out["fut"][j] = a.mean(), b.mean()
        if a.mean() > 0:
            out["ratio"][j] = b.mean() / a.mean()
            out["ratio_mr"][j] = bm.mean() / a.mean()
        out["mean_pct"][j] = 100 * (mfm / mcm - 1)
    out["valid"] = out["rate"] >= 0.9 * TARGET
    return out


def metric6(d, kind):
    x = ev.metric(d, kind)
    return x[:, BI]


def null_ranges(counts, sizes):
    """counts (W, nB) at the control-matched threshold -> {(n_a, n_b): (p05, p95) (nB,)}, keyed by the
    comparison's own block sizes. When n_a + n_b leaves too few disjoint placements in W winters
    (the 100 km comparison uses all 147 control winters), n_a is cut to W - n_b - 40 (>= n_b): a
    shorter control block gives a WIDER range, so the test stays conservative."""
    out = {}
    W = len(counts)
    for na0, nb in sizes:
        na = na0 if na0 + nb <= W - 40 else max(W - nb - 40, nb)
        pairs = ev.block_pairs(W, na, nb)
        lr = ev.block_logratios(counts.astype(float), na, nb, pairs)        # (P, nB)
        out[(na0, nb)] = np.nanpercentile(lr, [5, 95], axis=0)
    return out


def aggregations():
    """name -> (labels, A (nAgg, nB) weights). Tiles of 10/20/40 deg (from 20 N, 180 W) built from the
    5 deg boxes (cos-lat weights, a tile needs >= 2/4/8 land boxes), the 24 named regions (boxes whose
    centre lies in the region box; DE: DE_BOX) and Europe (35-70 N, 10 W-30 E)."""
    la = np.array([rg.BOXES[b][0] for b in B]) + 2.5
    lo = np.array([rg.BOXES[b][2] for b in B]) + 2.5
    w = np.cos(np.deg2rad(la))
    out = {"5deg": (list(B), np.diag(w))}
    for s, nmin in ((10, 2), (20, 4), (40, 8)):
        ti = ((la - 20) // s).astype(int) * 1000 + ((lo + 180) // s).astype(int)
        labs, rows = [], []
        for t in np.unique(ti):
            m = ti == t
            if m.sum() >= nmin:
                labs.append(f"t{s}_{t}")
                rows.append(np.where(m, w, 0.0))
        out[f"{s}deg"] = (labs, np.array(rows))
    named = dict(rg.REGION_BOXES, DE=rg.DE_BOX, Europe=(35.0, 70.0, -10.0, 30.0))
    labs, rows = [], []
    for r, (a, b, c, d) in named.items():
        m = (la >= a) & (la <= b) & (lo >= c) & (lo <= d)
        if m.any():
            labs.append(r)
            rows.append(np.where(m, w, 0.0))
    out["named"] = (labs, np.array(rows))
    return out


AGG = None


def aggregate(X, A):
    ok = np.isfinite(X)
    num = np.where(ok, X, 0.0) @ A.T
    den = ok.astype(float) @ A.T
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den >= 0.5 * A.sum(1)[None], num / den, np.nan)


def scale_change(Xc, Xf, mc, mf, win, dT):
    """Area-mean per K at every aggregation scale (5 deg handled by the box maps)."""
    out = {}
    for name, (labs, A) in AGG.items():
        if name == "5deg":
            continue
        o = box_change(aggregate(Xc, A), aggregate(Xf, A), mc, mf, win)
        ok = o["valid"] & (o["ratio"] > 0)
        lr = np.where(ok, np.log(np.where(ok, o["ratio"], 1.0)), np.nan)
        lrm = np.where(ok & (o["ratio_mr"] > 0), np.log(np.where(ok & (o["ratio_mr"] > 0), o["ratio_mr"], 1.0)), np.nan)
        wt = A.sum(1)
        out[name] = dict(labels=labs, lnratio=lr.tolist(), n_valid=int(ok.sum()),
                         area_mean_perK=float(100 * (np.exp(np.nansum(lr * wt) / wt[ok].sum() / dT) - 1)) if ok.any() else np.nan,
                         area_mean_perK_meanremoved=float(100 * (np.exp(np.nansum(lrm * wt) / wt[np.isfinite(lrm)].sum() / dT) - 1)) if ok.any() else np.nan,
                         frac_up=float(np.mean(lr[ok] > 0)) if ok.any() else np.nan)
    return out


def ours(T):
    res = {}
    for name in OURS:
        (cn, c0, c1), (fn, f0, f1) = ev.COMPARISONS[name]
        if not all(os.path.exists(f"{rg.OUT_DIR}/{n}.npz") for n in (cn, fn)):
            print(f"skip {name}: missing input", flush=True)
            continue
        dc, df = ev.load(cn), ev.load(fn)
        mc = [m for _, m in ev.winters(dc, c0, c1)]
        mf = [m for _, m in ev.winters(df, f0, f1)]
        dT = rpk.mean_years(T, fn, f0, f1)[0] - rpk.mean_years(T, cn, c0, c1)[0]
        res[name] = dict(dT=dT, n_winters=(len(mc), len(mf)))
        for kind in KINDS:
            res[name][kind] = box_change(metric6(dc, kind), metric6(df, kind), mc, mf, 8)
            res[name][kind]["scales"] = scale_change(metric6(dc, kind), metric6(df, kind), mc, mf, 8, dT)
        print(f"  {name}: dT {dT:.2f}, winters {len(mc)}/{len(mf)}", flush=True)
        ev._CACHE.clear()
    return res


def null6(res):
    """Per-box null from the TCO95 6-hourly control with each comparison's block sizes."""
    d = ev.load("tco95_ctl1950g")
    ws = ev.winters(d, 1953, 2099)
    out = {}
    for kind in KINDS:
        X = metric6(d, kind)
        allm = np.concatenate([m[None] for _, m in ws]).any(0)
        thr = box_change(X, X, [m for _, m in ws], [m for _, m in ws], 8)["thr"]
        ev.WIN = 8
        C = np.full((len(ws), X.shape[1]), np.nan)
        for j in range(X.shape[1]):
            if np.isfinite(thr[j]):
                C[:, j] = [ev.counts_grid(X[m, j], np.array([thr[j]]))[0] for _, m in ws]
        sizes = sorted({tuple(v["n_winters"]) for v in res.values()})
        out[kind] = {f"{a}_{b}": v for (a, b), v in null_ranges(C, sizes).items()}
        del allm
    out["W"] = len(ws)
    ev._CACHE.clear()
    return out


def load_daily(fn, exp):
    z = np.load(fn)
    d = {k.split("|", 1)[1]: z[k] for k in z.files if k.startswith(exp + "|")}
    y, m = d["ymd"][:, 0], d["ymd"][:, 1]
    d["winter"] = np.where(m >= 9, y + 1, y)
    return d


def daily_metric(d, kind):
    x = d["wind"] if kind == "wind" else W_SOLAR * d["solar"] + (1 - W_SOLAR) * d["wind"]
    return x[:, BI]


def dmasks(d, y0, y1):
    return [d["winter"] == y for y in range(y0, y1 + 1) if (d["winter"] == y).sum() >= 200]


def daily():
    res = {}
    gt = {}
    for mod in ("eerie", "tco95"):
        f = f"{rg.BASE_DIR}/global_t2m_{mod}.json"
        gt.update({k: {int(y): t for y, t in v.items()} for k, v in json.load(open(f)).items()})
    jobs = []
    for fn in sorted(glob.glob(f"{DC}/daily_cmip6_*.npz")):
        m = os.path.basename(fn)[len("daily_cmip6_"):-4]
        jobs.append((f"cmip6|{m}|2071-2100", fn, "historical", 1951, 1980, fn, "ssp585", 2071, 2100, None, None))
        jobs.append((f"cmip6_mid|{m}|2041-2070", fn, "historical", 1951, 1980, fn, "ssp585", 2041, 2070, None, None))
    E = lambda r: f"{DC}/daily_eerie_{r}.npz"          # noqa: E731
    T95 = lambda r: f"{DC}/daily_tco95_{r}.npz"        # noqa: E731
    jobs.append(("bridge|EERIE ctl->ssp245", E("ctl"), "run", 1972, 2050, E("ssp245"), "run", 2041, 2070, gt["eerie_ctl"], gt["eerie_ssp245"]))
    jobs.append(("bridge|EERIE hist->ssp245", E("hist"), "run", 1971, 2000, E("ssp245"), "run", 2041, 2070, gt["eerie_hist"], gt["eerie_ssp245"]))
    jobs.append(("bridge|TCO95 ctl->hi", T95("ctl1950g"), "run", 1953, 2099, T95("hi1950g"), "run", 2081, 2099, gt["tco95_ctl1950g"], gt["tco95_hi1950g"]))
    for name, fc, ec, c0, c1, ff, ef, f0, f1, Tc, Tf in jobs:
        if not (os.path.exists(fc) and os.path.exists(ff)):
            print(f"skip {name}: missing input", flush=True)
            continue
        dc, df = load_daily(fc, ec), load_daily(ff, ef)
        if Tc is None:
            Tc = dict(zip(dc["tas_years"].tolist(), dc["tas"].tolist()))
            Tf = dict(zip(df["tas_years"].tolist(), df["tas"].tolist()))
        dT = np.mean([Tf[y] for y in range(f0, f1 + 1) if y in Tf]) - np.mean([Tc[y] for y in range(c0, c1 + 1) if y in Tc])
        mc, mf = dmasks(dc, c0, c1), dmasks(df, f0, f1)
        res[name] = dict(dT=float(dT), n_winters=(len(mc), len(mf)))
        for kind in KINDS:
            res[name][kind] = box_change(daily_metric(dc, kind), daily_metric(df, kind), mc, mf, 2)
            if not name.startswith("cmip6_mid"):
                res[name][kind]["scales"] = scale_change(daily_metric(dc, kind), daily_metric(df, kind), mc, mf, 2, dT)
        print(f"  {name}: dT {dT:.2f}", flush=True)
    return res


def per_K(r, dT):
    return 100 * (np.exp(np.log(r) / dT) - 1)


def wmean(x, w):
    ok = np.isfinite(x)
    return float((x[ok] * w[ok]).sum() / w[ok].sum()) if ok.any() else np.nan


def wcorr(a, b, w):
    ok = np.isfinite(a) & np.isfinite(b)
    a, b, w = a[ok], b[ok], w[ok]
    am, bm = (a * w).sum() / w.sum(), (b * w).sum() / w.sum()
    return float(((a - am) * (b - bm) * w).sum() / np.sqrt(((a - am) ** 2 * w).sum() * ((b - bm) ** 2 * w).sum()))


def main():
    os.makedirs(OUT, exist_ok=True)
    lat0 = np.array([rg.BOXES[b][0] for b in B])
    lon0 = np.array([rg.BOXES[b][2] for b in B])
    w = np.cos(np.deg2rad(lat0 + 2.5))
    global AGG
    AGG = aggregations()
    T = rpk.global_t()
    R6 = ours(T)
    N6 = null6(R6) if os.path.exists(f"{rg.OUT_DIR}/tco95_ctl1950g.npz") else {}
    RD = daily()

    arrays = dict(boxes=np.array(B), lat0=lat0, lon0=lon0)
    summ = dict(n_boxes=len(B), target=TARGET, ours={}, daily={}, cmip6={})
    for src, R in (("6h", R6), ("1d", RD)):
        for name, v in R.items():
            for kind in KINDS:
                o = v[kind]
                ok = o["valid"] & np.isfinite(o["ratio"]) & (o["ratio"] > 0)
                pk = np.where(ok, per_K(np.where(ok, o["ratio"], 1.0), v["dT"]), np.nan)
                pk_mr = np.where(ok & (o["ratio_mr"] > 0), per_K(np.where(ok & (o["ratio_mr"] > 0), o["ratio_mr"], 1.0), v["dT"]), np.nan)
                key = f"{src}|{name}|{kind}"
                for a in ("ratio", "ratio_mr", "mean_pct", "ctl", "fut"):
                    arrays[f"{key}|{a}"] = np.where(ok, o[a], np.nan)
                arrays[f"{key}|perK"] = pk
                arrays[f"{key}|perK_mr"] = pk_mr
                row = dict(dT=v["dT"], n_valid=int(ok.sum()), frac_up=float(np.mean(o["ratio"][ok] > 1)),
                           area_mean_lnratio_ratio=float(np.exp(wmean(np.log(np.where(ok, o["ratio"], np.nan)), w))),
                           area_mean_perK=wmean(pk, w), area_mean_perK_meanremoved=wmean(pk_mr, w))
                if src == "6h" and kind in N6:
                    nb = f"{v['n_winters'][0]}_{v['n_winters'][1]}"
                    if nb in N6[kind]:
                        p05, p95 = N6[kind][nb]
                        lr = np.log(np.where(ok, o["ratio"], np.nan))
                        arrays[f"{key}|null_p05"], arrays[f"{key}|null_p95"] = p05, p95
                        row["frac_above_null95"] = float(np.nanmean(lr[ok] > p95[ok]))
                        row["frac_below_null05"] = float(np.nanmean(lr[ok] < p05[ok]))
                (summ["ours"] if src == "6h" else summ["daily"])[f"{name}|{kind}"] = row

    # CMIP6 ensemble maps
    for grp in ("cmip6", "cmip6_mid"):
        names = [n for n in RD if n.startswith(grp + "|")]
        if not names:
            continue
        for kind in KINDS:
            P = np.stack([arrays[f"1d|{n}|{kind}|perK"] for n in names])
            L = np.stack([np.log(arrays[f"1d|{n}|{kind}|ratio"]) for n in names])
            nval = np.isfinite(P).sum(0)
            enough = nval >= 0.8 * len(names)
            med = np.where(enough, np.nanmedian(P, axis=0), np.nan)
            agree = np.where(enough, np.nanmean(L > 0, axis=0), np.nan)
            med_mr = np.where(enough, np.nanmedian(np.stack([arrays[f"1d|{n}|{kind}|perK_mr"] for n in names]), axis=0), np.nan)
            arrays[f"{grp}|{kind}|perK_median"] = med
            arrays[f"{grp}|{kind}|perK_mr_median"] = med_mr
            arrays[f"{grp}|{kind}|agree_up"] = agree
            arrays[f"{grp}|{kind}|perK_p10"] = np.where(enough, np.nanpercentile(P, 10, axis=0), np.nan)
            arrays[f"{grp}|{kind}|perK_p90"] = np.where(enough, np.nanpercentile(P, 90, axis=0), np.nan)
            row = dict(n_models=len(names), n_boxes=int(enough.sum()), area_mean_perK_median=wmean(med, w),
                       frac_boxes_median_up=float(np.nanmean(med[enough] > 0)),
                       frac_boxes_80pct_agree_up=float(np.nanmean(agree[enough] >= 0.8)),
                       frac_boxes_80pct_agree_down=float(np.nanmean(agree[enough] <= 0.2)),
                       area_mean_perK_meanremoved_median=wmean(med_mr, w))
            # pattern correlation with our runs (and between CMIP6 models, for scale)
            pc = {}
            for o in OURS:
                k6 = f"6h|{o}|{kind}|perK"
                if k6 in arrays:
                    pc[o] = wcorr(arrays[k6], med, w)
            for bname in [n for n in RD if n.startswith("bridge|")]:
                pc[bname] = wcorr(arrays[f"1d|{bname}|{kind}|perK"], med, w)
            loo = []
            for i, n in enumerate(names):
                rest = np.nanmedian(np.delete(P, i, axis=0), axis=0)
                loo.append(wcorr(P[i], rest, w))
            row["pattern_corr_with_median"] = pc
            row["pattern_corr_single_model_vs_rest_median"] = [float(np.median(loo)), float(np.min(loo)), float(np.max(loo))]
            # zonal means of the median per-K
            zon = {}
            for la in np.unique(lat0):
                s = lat0 == la
                zon[f"{la:.0f}-{la + 5:.0f}N"] = wmean(med[s], w[s])
            row["zonal_mean_perK_median"] = zon
            summ["cmip6"][f"{grp}|{kind}"] = row

    # per K against aggregation scale
    sc = {}
    for src, R in (("6h", R6), ("1d", RD)):
        for name, v in R.items():
            if name.startswith("cmip6_mid"):
                continue
            for kind in KINDS:
                s5 = summ["ours" if src == "6h" else "daily"][f"{name}|{kind}"]
                row = {"5deg": dict(area_mean_perK=s5["area_mean_perK"], area_mean_perK_meanremoved=s5["area_mean_perK_meanremoved"], frac_up=s5["frac_up"])}
                for scl, x in v[kind]["scales"].items():
                    row[scl] = {k: x[k] for k in ("area_mean_perK", "area_mean_perK_meanremoved", "frac_up", "n_valid")}
                sc[f"{src}|{name}|{kind}"] = row
    ens_sc = {}
    for kind in KINDS:
        names = [n for n in RD if n.startswith("cmip6|")]
        if not names:
            continue
        ens_sc[kind] = {}
        for scl in ("5deg", "10deg", "20deg", "40deg", "named"):
            vals = np.array([sc[f"1d|{n}|{kind}"][scl]["area_mean_perK"] for n in names])
            vmr = np.array([sc[f"1d|{n}|{kind}"][scl]["area_mean_perK_meanremoved"] for n in names])
            ens_sc[kind][scl] = dict(median=float(np.nanmedian(vals)), p10_p90=[float(x) for x in np.nanpercentile(vals, [10, 90])],
                                     median_meanremoved=float(np.nanmedian(vmr)))
        labs = RD[names[0]][kind]["scales"]["named"]["labels"]
        L = np.array([RD[n][kind]["scales"]["named"]["lnratio"] for n in names])
        dTs = np.array([RD[n]["dT"] for n in names])[:, None]
        ens_sc[kind]["named_regions_perK_median"] = dict(zip(labs, [float(x) for x in np.nanmedian(100 * (np.exp(L / dTs) - 1), axis=0)]))
    summ["scales"] = sc
    summ["cmip6_scales"] = ens_sc

    np.savez(f"{OUT}/grid_per_K.npz", **arrays)
    json.dump(summ, open(f"{OUT}/grid_per_K_summary.json", "w"), indent=1)
    for sec in ("ours", "daily"):
        print(f"\n{sec}:")
        for k, v in summ[sec].items():
            if sec == "daily" and k.startswith("cmip6"):
                continue
            print(f"  {k:34s} dT {v['dT']:.2f}  boxes {v['n_valid']:3d}  up {v['frac_up']:.0%}  area-mean ratio {v['area_mean_lnratio_ratio']:.2f}  "
                  f"perK {v['area_mean_perK']:+.1f} (mean removed {v['area_mean_perK_meanremoved']:+.1f})"
                  + (f"  >null95 {v['frac_above_null95']:.0%} <null05 {v['frac_below_null05']:.0%}" if "frac_above_null95" in v else ""))
    for k, v in summ["cmip6"].items():
        print(f"\n{k}: {v['n_models']} models, {v['n_boxes']} boxes, area-mean median {v['area_mean_perK_median']:+.1f} %/K "
              f"(mean removed {v['area_mean_perK_meanremoved_median']:+.1f}); median up in {v['frac_boxes_median_up']:.0%}; "
              f">=80% agree up {v['frac_boxes_80pct_agree_up']:.0%}, down {v['frac_boxes_80pct_agree_down']:.0%}")
        print("  pattern corr with CMIP6 median: " + ", ".join(f"{a} {b:.2f}" for a, b in v["pattern_corr_with_median"].items())
              + f"; single model vs rest: median {v['pattern_corr_single_model_vs_rest_median'][0]:.2f} "
              f"({v['pattern_corr_single_model_vs_rest_median'][1]:.2f}..{v['pattern_corr_single_model_vs_rest_median'][2]:.2f})")
        print("  zonal: " + ", ".join(f"{a} {b:+.1f}" for a, b in v["zonal_mean_perK_median"].items()))
    print("\nArea-mean per K (%/K) against aggregation scale [mean removed]:")
    for k, row in sc.items():
        if k.startswith("1d|cmip6|"):
            continue
        print(f"  {k:34s} " + "  ".join(f"{s} {x['area_mean_perK']:+.1f} [{x['area_mean_perK_meanremoved']:+.1f}]" for s, x in row.items()))
    for kind, row in ens_sc.items():
        print(f"  CMIP6 median {kind:5s}                 " + "  ".join(f"{s} {x['median']:+.1f} [{x['median_meanremoved']:+.1f}]" for s, x in row.items() if s != "named_regions_perK_median"))
        print("    named regions (CMIP6 median %/K): " + ", ".join(f"{a} {b:+.1f}" for a, b in row["named_regions_perK_median"].items()))
    print(f"Saved {OUT}/grid_per_K.npz")


if __name__ == "__main__":
    main()
