"""
Per-K change in winter low-wind and low wind-plus-solar events on a DAILY footing: CMIP6 (about 20
models), and the same analysis for EERIE and TCO95 aggregated to days, which bridges to the
6-hourly results of dunkelflaute/scripts/regions/.

Events: 48 h = 2 daily steps, Mockert window expansion, Sep-Mar winters labelled by their Jan-Mar
year, threshold matched so the comparison's control has TARGET events per winter (per region and
metric). Metrics: wind (onshore CF of the log-law 100 m daily-mean wind) and comb (0.44 solar +
0.56 wind).
Warming dT: global annual-mean 2 m temperature of the future years minus the control years
(CMIP6: Amon/tas of the same member; EERIE and TCO95: global_t2m_*.json).
Per region: ln(ratio)/dT -> %/K. Field mean over regions valid in every comparison (control reaches
90% of TARGET and the region has cells), so all models share one region set.
Null for the field mean: the TCO95 constant-1950 control on the daily footing (147 winters), all
disjoint placements of two 30-winter blocks, the same winters in every region.
Across models: spread of ln(ratio) against spread of ln(ratio)/dT, and the regression of ln(ratio)
on dT through the model ensemble.

Usage: python3 daily_events.py
Output: dunkelflaute/results/cmip6/daily_events_summary.json
"""
import glob
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/analysis")
import regions as rg                                          # noqa: E402
import region_events as ev                                    # noqa: E402
from robustness_mockert_events import trailing_mean           # noqa: E402

ev.WIN = 2
D = f"{_ROOT}/dunkelflaute/data/cmip6"
OUT = f"{_ROOT}/dunkelflaute/results/cmip6"
TARGET = 5.0
MIN_DAYS = 200
W_SOLAR = 0.44
R = rg.REGIONS


def winters_of(ymd):
    return np.where(ymd[:, 1] >= 9, ymd[:, 0] + 1, ymd[:, 0])


def load(fn, exp):
    z = np.load(fn)
    d = {k.split("|", 1)[1]: z[k] for k in z.files if k.startswith(exp + "|")}
    d["winter"] = winters_of(d["ymd"])
    return d


def metric(d, kind):
    return d["wind"] if kind == "wind" else W_SOLAR * d["solar"] + (1 - W_SOLAR) * d["wind"]


def winter_masks(d, y0, y1):
    out = []
    for y in range(y0, y1 + 1):
        m = d["winter"] == y
        if m.sum() >= MIN_DAYS:
            out.append(m)
    return out


def counts(d, kind, i, y0, y1, thr):
    x = metric(d, kind)[:, i]
    return np.stack([ev.counts_grid(x[m], thr) for m in winter_masks(d, y0, y1)])      # (W, K)


def matched(d, kind, i, y0, y1):
    x = metric(d, kind)[:, i]
    if not np.isfinite(x).any():
        return np.nan, 0.0
    grid = np.linspace(*np.nanpercentile(x, [0.1, 40]), 400)
    c = counts(d, kind, i, y0, y1, grid).mean(axis=0)
    k = int(np.argmax(c >= TARGET)) if (c >= TARGET).any() else int(np.argmax(c))
    return float(grid[k]), float(c[k])


def kinds_of(group):
    return ("wind",) if group.startswith("cmip6_windonly") else ("wind", "comb")


def mean_T(T, y0, y1):
    ys = [y for y in range(y0, y1 + 1) if y in T]
    return float(np.mean([T[y] for y in ys])) if ys else np.nan


def comparisons():
    """name -> (file, ctl_exp, c0, c1, file, fut_exp, f0, f1, T_ctl dict, T_fut dict, group)"""
    out = {}
    inv = json.load(open(f"{D}/inventory.json"))
    for fn in sorted(glob.glob(f"{D}/daily_cmip6_*.npz")):
        m = os.path.basename(fn)[len("daily_cmip6_"):-4]
        z = np.load(fn)
        T = {e: dict(zip(z[f"{e}|tas_years"].tolist(), z[f"{e}|tas"].tolist())) for e in ("historical", "ssp585")}
        g = "cmip6_windonly" if inv.get(m, {}).get("wind_only") else "cmip6"     # wind-only models stay out of the ensemble
        out[f"{m}|2071-2100"] = (fn, "historical", 1951, 1980, fn, "ssp585", 2071, 2100, T["historical"], T["ssp585"], g)
        out[f"{m}|2041-2070"] = (fn, "historical", 1951, 1980, fn, "ssp585", 2041, 2070, T["historical"], T["ssp585"], g + "_mid")
    gt = {}
    for mod in ("eerie", "tco95"):
        f = f"{rg.BASE_DIR}/global_t2m_{mod}.json"
        if os.path.exists(f):
            gt.update({k: {int(y): t for y, t in v.items()} for k, v in json.load(open(f)).items()})
    E = lambda r: f"{D}/daily_eerie_{r}.npz"     # noqa: E731
    T95 = lambda r: f"{D}/daily_tco95_{r}.npz"   # noqa: E731
    if os.path.exists(E("ctl")) and os.path.exists(E("ssp245")):
        out["EERIE|ctl->ssp245 2041-2070"] = (E("ctl"), "run", 1972, 2050, E("ssp245"), "run", 2041, 2070, gt["eerie_ctl"], gt["eerie_ssp245"], "bridge")
    if os.path.exists(E("hist")) and os.path.exists(E("ssp245")):
        out["EERIE|hist 1971-2000->ssp245 2041-2070"] = (E("hist"), "run", 1971, 2000, E("ssp245"), "run", 2041, 2070, gt["eerie_hist"], gt["eerie_ssp245"], "bridge")
    if os.path.exists(T95("ctl1950g")) and os.path.exists(T95("hi1950g")):
        out["TCO95|ctl->hi 2081-2099"] = (T95("ctl1950g"), "run", 1953, 2099, T95("hi1950g"), "run", 2081, 2099, gt["tco95_ctl1950g"], gt["tco95_hi1950g"], "bridge")
        out["TCO95|hi 1951-1980->2071-2100"] = (T95("hi1950g"), "run", 1953, 1980, T95("hi1950g"), "run", 2071, 2099, gt["tco95_hi1950g"], gt["tco95_hi1950g"], "bridge_cmipstyle")
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    C = comparisons()
    res = {}
    store = {}
    for name, (fc, ec, c0, c1, ff, ef, f0, f1, Tc, Tf, group) in C.items():
        dc, df = load(fc, ec), load(ff, ef)
        dT = mean_T(Tf, f0, f1) - mean_T(Tc, c0, c1)
        row = dict(group=group, dT=dT, metrics={})
        for kind in kinds_of(group):
            per = {}
            for i, r in enumerate(R):
                thr, rate = matched(dc, kind, i, c0, c1)
                if not np.isfinite(thr):
                    per[r] = dict(valid=False)
                    continue
                a = counts(dc, kind, i, c0, c1, np.array([thr]))[:, 0]
                b = counts(df, kind, i, f0, f1, np.array([thr]))[:, 0]
                xc, xf = metric(dc, kind)[:, i], metric(df, kind)[:, i]
                mc = np.nanmean(np.concatenate([xc[m] for m in winter_masks(dc, c0, c1)]))
                mf = np.nanmean(np.concatenate([xf[m] for m in winter_masks(df, f0, f1)]))
                b_adj = np.array([ev.counts_grid(xf[m] * mc / mf, np.array([thr]))[0] for m in winter_masks(df, f0, f1)])
                per[r] = dict(valid=bool(rate >= 0.9 * TARGET), threshold=thr, control=float(a.mean()), future=float(b.mean()),
                              ratio=float(b.mean() / a.mean()) if a.mean() else np.nan,
                              ratio_meanremoved=float(b_adj.mean() / a.mean()) if a.mean() else np.nan,
                              mean_change_pct=float(100 * (mf / mc - 1)))
                store[(name, kind, r)] = (a, b)
            row["metrics"][kind] = per
        res[name] = row
        print(f"  {name:38s} dT {dT:+.2f} K", flush=True)

    # common valid region set across all comparisons, per metric (wind-only models do not shrink it)
    common = {}
    for kind in ("wind", "comb"):
        common[kind] = [r for r in R if all(res[n]["metrics"][kind][r].get("valid") for n in res if not res[n]["group"].startswith("cmip6_windonly"))]

    # field statistics
    summary = {}
    for name, row in res.items():
        summary[name] = dict(group=row["group"], dT=row["dT"])
        for kind in kinds_of(row["group"]):
            V = common[kind]
            lr = np.array([np.log(row["metrics"][kind][r]["ratio"]) for r in V])
            la = np.array([np.log(row["metrics"][kind][r]["ratio_meanremoved"]) for r in V])
            ok = np.isfinite(lr)
            m = float(lr[ok].mean())
            summary[name][kind] = dict(n_regions=int(ok.sum()), mean_ratio=float(np.exp(m)),
                                       pct_per_K=float(100 * (np.exp(m / row["dT"]) - 1)),
                                       n_regions_up=int((lr[ok] > 0).sum()),
                                       share_from_mean=float(1 - np.nanmean(la[ok]) / m) if m else np.nan,
                                       mean_change_pct=float(np.mean([row["metrics"][kind][r]["mean_change_pct"] for r in V])))

    # null from the TCO95 daily control, 30 vs 30 winters, same winters in all regions
    null = {}
    t95 = f"{D}/daily_tco95_ctl1950g.npz"
    if os.path.exists(t95):
        d = load(t95, "run")
        for kind in ("wind", "comb"):
            V = common[kind]
            cols = []
            for r in V:
                i = R.index(r)
                thr, _ = matched(d, kind, i, 1953, 2099)
                cols.append(counts(d, kind, i, 1953, 2099, np.array([thr]))[:, 0])
            c = np.stack(cols, axis=1).astype(float)
            pairs = ev.block_pairs(len(c), 30, 30)
            lr = ev.block_logratios(c, 30, 30, pairs)
            fm = np.nanmean(lr, axis=1)
            null[kind] = dict(n_winters=len(c), n_pairs=len(pairs), field_ratio_p05_p95=[float(x) for x in np.exp(np.percentile(fm, [5, 95]))])

    # ensemble: does normalising per K shrink the spread?
    ens = {}
    for grp in ("cmip6", "cmip6_mid"):
        names = [n for n in summary if summary[n]["group"] == grp]
        if not names:
            continue
        ens[grp] = {}
        for kind in ("wind", "comb"):
            lnr = np.array([np.log(summary[n][kind]["mean_ratio"]) for n in names])
            dT = np.array([summary[n]["dT"] for n in names])
            pk = np.array([summary[n][kind]["pct_per_K"] for n in names])
            slope, icpt = np.polyfit(dT, lnr, 1)
            ens[grp][kind] = dict(n_models=len(names), mean_ratio_median=float(np.exp(np.median(lnr))),
                                  ratio_range=[float(np.exp(lnr.min())), float(np.exp(lnr.max()))],
                                  dT_range=[float(dT.min()), float(dT.max())],
                                  pct_per_K_median=float(np.median(pk)), pct_per_K_p10_p90=[float(x) for x in np.percentile(pk, [10, 90])],
                                  cv_lnratio=float(lnr.std() / abs(lnr.mean())), cv_lnratio_per_K=float((lnr / dT).std() / abs((lnr / dT).mean())),
                                  corr_lnratio_dT=float(np.corrcoef(dT, lnr)[0, 1]),
                                  regression_pct_per_K=float(100 * (np.exp(slope) - 1)), regression_intercept_ratio=float(np.exp(icpt)),
                                  n_models_up=int((lnr > 0).sum()))

    out = dict(target=TARGET, window_days=2, common_regions=common, comparisons=summary, null_tco95_daily=null,
               ensemble=ens, regions=res)
    json.dump(out, open(f"{OUT}/daily_events_summary.json", "w"), indent=1)

    for kind in ("wind", "comb"):
        print(f"\n=== {kind} (daily footing, {len(common[kind])} common regions): mean ratio, %/K, regions up, share from mean")
        for n, s in sorted(summary.items(), key=lambda x: (x[1]["group"], x[0])):
            if kind not in s:
                continue
            x = s[kind]
            print(f"  {n:38s} dT {s['dT']:+.2f}  ratio {x['mean_ratio']:.2f}  {x['pct_per_K']:+5.1f} %/K  up {x['n_regions_up']:2d}/{x['n_regions']}  mean {x['share_from_mean']:.0%}  dmeanCF {x['mean_change_pct']:+.1f}%")
        if kind in null:
            print(f"  null (TCO95 daily control, 30 vs 30 winters): field ratio 5-95% {null[kind]['field_ratio_p05_p95']}")
        for g, e in ens.items():
            x = e[kind]
            print(f"  ensemble {g}: {x['n_models']} models, ratio median {x['mean_ratio_median']:.2f} (range {x['ratio_range'][0]:.2f}-{x['ratio_range'][1]:.2f}), "
                  f"{x['pct_per_K_median']:+.1f} %/K (10-90% {x['pct_per_K_p10_p90'][0]:+.1f}..{x['pct_per_K_p10_p90'][1]:+.1f}); "
                  f"CV of ln ratio {x['cv_lnratio']:.2f} vs per K {x['cv_lnratio_per_K']:.2f}; r(dT, ln ratio) {x['corr_lnratio_dT']:.2f}; "
                  f"regression {x['regression_pct_per_K']:+.1f} %/K, intercept {x['regression_intercept_ratio']:.2f}; up in {x['n_models_up']}")
    print(f"Saved {OUT}/daily_events_summary.json")


if __name__ == "__main__":
    main()
