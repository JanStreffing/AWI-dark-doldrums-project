"""
Event ratios per K of global warming, and what part of the wind distribution drives them.

1. Per K: for every comparison in region_events.py, dT = global annual-mean 2t (global_t2m.py) of
   the future years minus that of the control years. Sensitivity per region = ln(ratio) / dT,
   reported as percent per K: 100 (exp(ln(ratio)/dT) - 1). Field mean over the valid regions with a
   winter bootstrap that resamples the SAME winters in every region (keeps spatial correlation;
   within-run, so it does not include decadal variability beyond the sampled winters).
2. Mean shift or tail change: the future CF series is rescaled so its Sep-Mar mean equals the
   control's (multiplicative, x * mean_ctl / mean_fut, per region), and events are recounted at the
   same matched threshold. ratio_meanremoved ~ 1 means the increase is a mean decline ("stilling");
   ratio_meanremoved ~ ratio means the low tail / persistence changed with the mean unchanged.
   Also the change in mean CF and in the 5th percentile of the 48 h mean.

Usage: python3 region_per_K.py
Output: dunkelflaute/results/regions/region_per_K_summary.json
"""
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/analysis")
import regions as rg                                                           # noqa: E402
import region_events as ev                                                     # noqa: E402
from robustness_mockert_events import trailing_mean                            # noqa: E402

RES = f"{_ROOT}/dunkelflaute/results/regions"
RNG = np.random.default_rng(2)
NBOOT = 5000


def global_t():
    out = {}
    for m in ("dart9km", "tco319", "tco95", "eerie"):
        fn = f"{rg.BASE_DIR}/global_t2m_{m}.json"
        if os.path.exists(fn):
            for k, v in json.load(open(fn)).items():
                out[k] = {int(y): t for y, t in v.items()}
    return out


def mean_years(T, run, y0, y1):
    ys = [y for y in range(y0, y1 + 1) if y in T.get(run, {})]
    return float(np.mean([T[run][y] for y in ys])) if ys else np.nan, len(ys)


def main():
    S = json.load(open(f"{RES}/region_events_summary.json"))
    C = np.load(f"{RES}/region_events_counts.npz")
    T = global_t()
    R = S["regions"]
    out = dict(dT={}, per_K={}, field_per_K={}, mean_vs_tail={})
    for name, ((cn, c0, c1), (fn, f0, f1)) in ev.COMPARISONS.items():
        if name not in S["comparisons"]:
            continue
        ta, na = mean_years(T, cn, c0, c1)
        tb, nb = mean_years(T, fn, f0, f1)
        dT = tb - ta
        out["dT"][name] = dict(dT=dT, control_C=ta - 273.15, future_C=tb - 273.15, n_years=[na, nb])
        if not np.isfinite(dT) or abs(dT) < 0.3:
            continue
        out["per_K"][name] = {}
        out["field_per_K"][name] = {}
        for kind in ("comb", "wind"):
            V = S["valid_regions"][kind]
            lr = {r: np.log(S["comparisons"][name][kind][r]["ratio"]) for r in R}
            out["per_K"][name][kind] = {r: float(100 * (np.exp(lr[r] / dT) - 1)) for r in R if np.isfinite(lr[r])}
            a = np.stack([C[f"{name}|{kind}|{r}|ctl"] for r in V], axis=1).astype(float)     # (Wc, R)
            b = np.stack([C[f"{name}|{kind}|{r}|fut"] for r in V], axis=1).astype(float)
            obs = float(np.mean(np.log(b.mean(0) / a.mean(0))))
            boots = []
            for _ in range(NBOOT):
                ia = RNG.integers(0, len(a), len(a))
                ib = RNG.integers(0, len(b), len(b))
                with np.errstate(divide="ignore"):
                    x = np.log(b[ib].mean(0) / a[ia].mean(0))
                boots.append(np.mean(x[np.isfinite(x)]))
            lo, hi = np.percentile(boots, [5, 95])
            pk = lambda x: float(100 * (np.exp(x / dT) - 1))   # noqa: E731
            out["field_per_K"][name][kind] = dict(n_regions=len(V), mean_ratio=float(np.exp(obs)), pct_per_K=pk(obs),
                                                  pct_per_K_p05_p95=[pk(lo), pk(hi)])

    # 2. mean shift vs tail/persistence
    for name, ((cn, c0, c1), (fn, f0, f1)) in ev.COMPARISONS.items():
        if name not in S["comparisons"] or name == "eerie_drift":
            continue
        dc, dfu = ev.load(cn), ev.load(fn)
        res = {}
        for kind in ("wind", "comb"):
            rows = {}
            for i, r in enumerate(R):
                if r not in S["valid_regions"][kind]:
                    continue
                thr = S["comparisons"][name][kind][r]["threshold"]
                xc = ev.metric(dc, kind)[:, i]
                xf = ev.metric(dfu, kind)[:, i]
                wc = ev.winters(dc, c0, c1)
                wf = ev.winters(dfu, f0, f1)
                mc = np.mean(np.concatenate([xc[m] for _, m in wc]))
                mf = np.mean(np.concatenate([xf[m] for _, m in wf]))
                rc = np.concatenate([trailing_mean(xc[m], ev.WIN)[ev.WIN - 1:] for _, m in wc])
                rf = np.concatenate([trailing_mean(xf[m], ev.WIN)[ev.WIN - 1:] for _, m in wf])
                n_c = np.mean([ev.counts_grid(xc[m], np.array([thr]))[0] for _, m in wc])
                n_f = np.mean([ev.counts_grid(xf[m], np.array([thr]))[0] for _, m in wf])
                n_f_adj = np.mean([ev.counts_grid(xf[m] * mc / mf, np.array([thr]))[0] for _, m in wf])
                rows[r] = dict(ratio=float(n_f / n_c), ratio_meanremoved=float(n_f_adj / n_c),
                               mean_change_pct=float(100 * (mf / mc - 1)),
                               p05_48h_change_pct=float(100 * (np.percentile(rf, 5) / np.percentile(rc, 5) - 1)),
                               cv_48h_ctl=float(rc.std() / rc.mean()), cv_48h_fut=float(rf.std() / rf.mean()))
            lr = np.log([v["ratio"] for v in rows.values()])
            la = np.log([v["ratio_meanremoved"] for v in rows.values()])
            res[kind] = dict(regions=rows, field=dict(mean_ratio=float(np.exp(lr.mean())),
                                                      mean_ratio_meanremoved=float(np.exp(la.mean())),
                                                      share_explained_by_mean=float(1 - la.mean() / lr.mean()),
                                                      mean_change_pct=float(np.mean([v["mean_change_pct"] for v in rows.values()])),
                                                      p05_48h_change_pct=float(np.mean([v["p05_48h_change_pct"] for v in rows.values()]))))
        out["mean_vs_tail"][name] = res

    json.dump(out, open(f"{RES}/region_per_K_summary.json", "w"), indent=1)

    print("Global warming between compared periods (annual-mean 2t):")
    for n, d in out["dT"].items():
        print(f"  {n:15s} {d['control_C']:.2f} -> {d['future_C']:.2f} C   dT = {d['dT']:+.2f} K  (years {d['n_years']})")
    for kind in ("wind", "comb"):
        print(f"\n=== {kind}: field-mean change in events per K of global warming (5-95% winter bootstrap)")
        for n, d in out["field_per_K"].items():
            x = d[kind]
            print(f"  {n:15s} ratio {x['mean_ratio']:.2f}  {x['pct_per_K']:+5.1f} %/K  [{x['pct_per_K_p05_p95'][0]:+5.1f}, {x['pct_per_K_p05_p95'][1]:+5.1f}]  ({x['n_regions']} regions)")
    for kind in ("wind", "comb"):
        print(f"\n=== {kind}: mean shift vs tail (field means)")
        for n, d in out["mean_vs_tail"].items():
            f = d[kind]["field"]
            print(f"  {n:15s} ratio {f['mean_ratio']:.2f} -> {f['mean_ratio_meanremoved']:.2f} with the mean change removed "
                  f"(share from mean {f['share_explained_by_mean']:.0%}); mean CF {f['mean_change_pct']:+.1f}%, "
                  f"5th pct of 48 h mean {f['p05_48h_change_pct']:+.1f}%")
    print(f"Saved {RES}/region_per_K_summary.json")


if __name__ == "__main__":
    main()
