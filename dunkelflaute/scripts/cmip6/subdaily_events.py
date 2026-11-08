"""
Does the daily footing bias the CMIP6 per-K change? For the four CMIP6 models with 3-hourly wind
and solar, the per-K change in low wind and low wind-plus-solar events on three footings of the SAME
winters (extract_3hr.py): 6h (WIN 8, our runs' footing), d3h (daily from the 3-hourly data, WIN 2)
and day (the day/ files, WIN 2, as daily_events.py). Historical 1961-1990 against ssp585 2071-2100.
Thresholds matched per footing, region and metric (TARGET events per control winter); field means of
ln ratio over the regions valid on every footing in every model.

Usage: python3 subdaily_events.py
Output: dunkelflaute/results/cmip6/subdaily_summary.json
"""
import glob
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("regions", "analysis"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/{p}")
import regions as rg            # noqa: E402
import region_events as ev      # noqa: E402

D = f"{_ROOT}/dunkelflaute/data/cmip6"
OUT = f"{_ROOT}/dunkelflaute/results/cmip6"
TARGET = 5.0
W_SOLAR = 0.44
WINTERS = {"historical": (1961, 1990), "ssp585": (2071, 2100)}
FOOT = {"6h": 8, "d3h": 2, "day": 2}
MIN_FRAC = 0.9


def winter_masks(z, exp, foot):
    if foot == "6h":
        t = z[f"{exp}|time6"].astype("datetime64[s]")
        w = ev.winter_of(t)
        full = 4 * 200
    else:
        if foot == "d3h":
            d = z[f"{exp}|day_d3h"].astype("datetime64[D]")
            y = d.astype("datetime64[Y]").astype(int) + 1970
            m = d.astype("datetime64[M]").astype(int) % 12 + 1
        else:
            y, m = z[f"{exp}|ymd_day"][:, 0], z[f"{exp}|ymd_day"][:, 1]
        w = np.where(m >= 9, y + 1, y)
        full = 200
    y0, y1 = WINTERS[exp]
    return [w == y for y in range(y0, y1 + 1) if (w == y).sum() >= full]


def series(z, exp, foot, kind):
    sfx = {"6h": "6", "d3h": "_d3h", "day": "_day"}[foot]
    wd, so = z[f"{exp}|wind{sfx}"], z[f"{exp}|solar{sfx}"]
    return wd if kind == "wind" else W_SOLAR * so + (1 - W_SOLAR) * wd


def main():
    res = {}
    for fn in sorted(glob.glob(f"{D}/sub3h_*.npz")):
        model = os.path.basename(fn)[6:-4]
        z = np.load(fn)
        assert list(z["regions"]) == rg.REGIONS
        T = {e: dict(zip(z[f"{e}|tas_years"].tolist(), z[f"{e}|tas"].tolist())) for e in WINTERS}
        dT = np.mean([T["ssp585"][y] for y in range(2071, 2101) if y in T["ssp585"]]) - np.mean([T["historical"][y] for y in range(1961, 1991) if y in T["historical"]])
        res[model] = dict(dT=float(dT))
        for foot, win in FOOT.items():
            ev.WIN = win
            mc, mf = winter_masks(z, "historical", foot), winter_masks(z, "ssp585", foot)
            res[model][foot] = {"n_winters": [len(mc), len(mf)]}
            for kind in ("wind", "comb"):
                xc_all, xf_all = series(z, "historical", foot, kind), series(z, "ssp585", foot, kind)
                per = {}
                for i, r in enumerate(rg.REGIONS):
                    xc, xf = xc_all[:, i], xf_all[:, i]
                    if not np.isfinite(xc).all():
                        continue
                    grid = np.linspace(*np.percentile(xc, [0.1, 40]), 400)
                    cc = np.stack([ev.counts_grid(xc[m], grid) for m in mc]).mean(0)
                    k = int(np.argmax(cc >= TARGET)) if (cc >= TARGET).any() else int(np.argmax(cc))
                    thr = np.array([grid[k]])
                    a = np.array([ev.counts_grid(xc[m], thr)[0] for m in mc], float)
                    b = np.array([ev.counts_grid(xf[m], thr)[0] for m in mf], float)
                    per[r] = dict(valid=bool(cc[k] >= MIN_FRAC * TARGET), ratio=float(b.mean() / a.mean()) if a.mean() else np.nan)
                res[model][foot][kind] = per
        print(f"  {model}: dT {dT:.2f}", flush=True)

    common = {kind: [r for r in rg.REGIONS if all(res[m][f][kind].get(r, {}).get("valid") for m in res for f in FOOT)]
              for kind in ("wind", "comb")}
    summ = {}
    for m, v in res.items():
        summ[m] = dict(dT=v["dT"])
        for f in FOOT:
            for kind in ("wind", "comb"):
                lr = np.array([np.log(v[f][kind][r]["ratio"]) for r in common[kind]])
                summ[m][f"{f}|{kind}"] = dict(mean_ratio=float(np.exp(lr.mean())), pct_per_K=float(100 * (np.exp(lr.mean() / v["dT"]) - 1)),
                                              n_up=int((lr > 0).sum()), n=len(lr))
        for kind in ("wind", "comb"):
            l6 = np.array([np.log(v["6h"][kind][r]["ratio"]) for r in common[kind]])
            ld = np.array([np.log(v["day"][kind][r]["ratio"]) for r in common[kind]])
            l3 = np.array([np.log(v["d3h"][kind][r]["ratio"]) for r in common[kind]])
            summ[m][f"corr_regions_6h_day|{kind}"] = float(np.corrcoef(l6, ld)[0, 1])
            summ[m][f"corr_regions_6h_d3h|{kind}"] = float(np.corrcoef(l6, l3)[0, 1])
    out = dict(winters=WINTERS, target=TARGET, common_regions=common, models=summ, regions=res)
    json.dump(out, open(f"{OUT}/subdaily_summary.json", "w"), indent=1)
    for kind in ("wind", "comb"):
        print(f"\n=== {kind} ({len(common[kind])} common regions): %/K on 6h | daily from 3h | daily files")
        for m, s in summ.items():
            print(f"  {m:15s} dT {s['dT']:.2f}  " + "  ".join(f"{f} {s[f'{f}|{kind}']['pct_per_K']:+5.1f} (ratio {s[f'{f}|{kind}']['mean_ratio']:.2f}, up {s[f'{f}|{kind}']['n_up']}/{s[f'{f}|{kind}']['n']})" for f in FOOT)
                  + f"  r(6h,day) {s[f'corr_regions_6h_day|{kind}']:.2f}")
    print(f"Saved {OUT}/subdaily_summary.json")


if __name__ == "__main__":
    main()
