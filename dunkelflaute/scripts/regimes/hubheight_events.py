"""
Low-wind events from interpolated 100 m wind against extrapolated 10 m wind, in the four CMIP6 models with daily
model-level wind (Gerrit Lohmann's review comment, 2026-11-17).

Inputs (all region means, 24 regions, daily, historical 1951-1980 against ssp585 2071-2100, Sep-Mar winters with at
least 200 days):
  cf    the paper's onshore CF of the log-law 100 m wind (data/cmip6/daily_cmip6_<model>.npz 'wind')
  s10   the daily-mean 10 m wind speed of the same file ('spd10'); the log law with fixed roughness is a constant
        factor per cell, so it ranks days like the 10 m wind
  v100  the 100 m wind interpolated linearly in ln z between model levels (hubheight_extract.py,
        data/regimes/hubheight_parts/hubheight_ml_<model>.npz 'V100ml_all'; speed of the daily-mean components)
Events: 48 h = 2 days, Mockert window expansion, threshold matched separately for each metric and region so that the
control has 5 events per winter (the low tail of each series), future counted at the same threshold. Matched
thresholds make the result independent of any monotone transform of a series, so s10 stands for the extrapolated
wind; cf differs from s10 only through the power curve and averaging over cells.
Per K: 100 (exp(mean over regions of ln ratio / dT) - 1), dT from results/cmip6/daily_events_summary.json.
Check: cf on the common days reproduces the stored per-region rates of daily_events_summary.json.
Caveat: v100 is the speed of daily-mean components, which is lower than the daily-mean speed on variable days.

Environment: levante, module load python3/2025.01-gcc-13.3.0
Usage: python3 hubheight_events.py
Output: dunkelflaute/results/regimes/hubheight_events.json
"""
import json
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/analysis")
from robustness_mockert_events import mockert_events, trailing_mean   # noqa: E402

MODELS = ["AWI-ESM-1-REcoM", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "NorESM2-LM"]
DC = f"{_ROOT}/dunkelflaute/data/cmip6"
HP = f"{_ROOT}/dunkelflaute/data/regimes/hubheight_parts"
OUT = f"{_ROOT}/dunkelflaute/results/regimes/hubheight_events.json"
EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
WIN, TARGET, MIN_DAYS = 2, 5.0, 200
C0, C1, F0, F1 = 1951, 1980, 2071, 2100


def winters_of(ymd):
    return np.where(ymd[:, 1] >= 9, ymd[:, 0] + 1, ymd[:, 0])


def key(ymd):
    return ymd @ np.array([10000, 100, 1])


def series(model, exp, regions):
    z, h = np.load(f"{DC}/daily_cmip6_{model}.npz"), np.load(f"{HP}/hubheight_ml_{model}.npz")
    assert list(z["regions"]) == list(h["regions"]) == regions
    kz, kh = key(z[f"{exp}|ymd"]), key(h[f"{exp}|ymd"])
    common, iz, ih = np.intersect1d(kz, kh, return_indices=True)
    ymd = z[f"{exp}|ymd"][iz]
    return ymd, dict(cf=z[f"{exp}|wind"][iz].astype(float), s10=z[f"{exp}|spd10"][iz].astype(float),
                     v100=h[f"{exp}|V100ml_all"][ih].astype(float)), len(kz) - len(common)


def winter_list(ymd, x, y0, y1):
    w = winters_of(ymd)
    return [x[w == y] for y in range(y0, y1 + 1) if (w == y).sum() >= MIN_DAYS]


def n_events(x, thr):
    return len(mockert_events(trailing_mean(x, WIN), WIN, thr)[1])


def matched(ws):
    allx = np.concatenate(ws)
    grid = np.linspace(*np.nanpercentile(allx, [0.1, 60]), 600)
    for g in grid:
        if np.mean([n_events(w, g) for w in ws]) >= TARGET:
            return float(g)
    return float(grid[-1])


def main():
    S = json.load(open(f"{_ROOT}/dunkelflaute/results/cmip6/daily_events_summary.json"))["regions"]
    regions = list(np.load(f"{DC}/daily_cmip6_{MODELS[0]}.npz")["regions"])
    out = dict(models={}, regions=regions, eu=EU)
    for m in MODELS:
        st = S[f"{m}|2071-2100"]
        dT = st["dT"]
        yc, Xc, dropc = series(m, "historical", regions)
        yf, Xf, dropf = series(m, "ssp585", regions)
        res = dict(dT=dT, days_dropped=[int(dropc), int(dropf)], per_region={})
        for i, r in enumerate(regions):
            row = {}
            for v in ("cf", "s10", "v100"):
                wc, wf = winter_list(yc, Xc[v][:, i], C0, C1), winter_list(yf, Xf[v][:, i], F0, F1)
                if not all(np.isfinite(np.concatenate(wc))) or not all(np.isfinite(np.concatenate(wf))):
                    row[v] = None
                    continue
                thr = matched(wc)
                a = np.mean([n_events(w, thr) for w in wc])
                b = np.mean([n_events(w, thr) for w in wf])
                row[v] = dict(threshold=thr, control=float(a), future=float(b), ratio=float(b / a))
                if v == "s10" or v == "v100":
                    row[v]["mean_change_pct"] = float(100 * (np.nanmean(np.concatenate(wf)) / np.nanmean(np.concatenate(wc)) - 1))
            stored = st["metrics"]["wind"][r]
            row["cf_stored"] = dict(control=stored["control"], future=stored["future"], ratio=stored["ratio"])
            res["per_region"][r] = row
        for tag, regs in (("all24", regions), ("eu8", EU)):
            agg = {}
            for v in ("cf", "s10", "v100"):
                L = [np.log(res["per_region"][r][v]["ratio"]) for r in regs if res["per_region"][r][v]]
                agg[v] = dict(n_regions=len(L), perK=float(100 * (np.exp(np.mean(L) / dT) - 1)))
                if v != "cf":
                    mc = [res["per_region"][r][v]["mean_change_pct"] for r in regs if res["per_region"][r][v]]
                    agg[v]["mean_change_pct"] = float(np.mean(mc))
            res[tag] = agg
        dev = max(abs(res["per_region"][r]["cf"]["future"] - res["per_region"][r]["cf_stored"]["future"])
                  for r in regions if res["per_region"][r]["cf"])
        res["cf_reproduce_max_abs_future_rate"] = float(dev)
        out["models"][m] = res
        print(f"{m:16s} dT {dT:.2f} dropped {dropc}/{dropf} days; cf repro max |d| {dev:.3f}")
        for tag in ("all24", "eu8"):
            a = res[tag]
            print(f"   {tag}: per K  cf {a['cf']['perK']:+.1f}  s10 {a['s10']['perK']:+.1f}  v100 {a['v100']['perK']:+.1f}"
                  f"   mean wind  s10 {a['s10']['mean_change_pct']:+.1f} %  v100 {a['v100']['mean_change_pct']:+.1f} %"
                  f"   (n {a['cf']['n_regions']}/{a['s10']['n_regions']}/{a['v100']['n_regions']})", flush=True)
    json.dump(out, open(OUT, "w"), indent=1)
    print(OUT)


if __name__ == "__main__":
    main()
