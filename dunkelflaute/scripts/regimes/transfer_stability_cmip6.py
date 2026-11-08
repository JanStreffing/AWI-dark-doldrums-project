"""
G-dependence part of the transfer change in CMIP6, for comparison with the TCO95 stability study
(transfer_stability.py). Per model and European region (regime_analysis.EU), daily sfcWind (s) and
geostrophic speed G from daily psl exactly as geostrophic.cmip6 (global grid, periodic, region
boxes; historical 1951-1980 against ssp585 2071-2100, Sep-Mar), split into control quintiles of G:
  between  T(f_fut, s_ctl, G_ctl) / T_ctl     the part a weaker G alone gives, s / G per G class fixed
  within   T_fut / T(f_fut, s_ctl, G_ctl)     the change of s / G at fixed G
If the G part is as large in CMIP6 as in TCO95, the model difference sits in the within-G part.

Usage: python3 transfer_stability_cmip6.py    (legacy python, login node, a few minutes)
Output: dunkelflaute/results/regimes/transfer_stability_cmip6.json
"""
import glob
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regimes")
os.environ.setdefault("NB", "200")
import transfer_stability as ts          # noqa: E402
gs, rg = ts.gs, ts.rg
DC = f"{_ROOT}/dunkelflaute/data/cmip6"
OUT = f"{_ROOT}/dunkelflaute/results/regimes/transfer_stability_cmip6.json"


def model(fn, rng):
    name = os.path.basename(fn)[4:-4]
    P = np.load(fn)
    Z = np.load(f"{DC}/daily_cmip6_{name}.npz")
    lat, lon = P["lat"], P["lon"]
    ri = [rg.REGIONS.index(r) for r in ts.EU]
    S = {}
    for exp, (y0, y1) in (("historical", (1951, 1980)), ("ssp585", (2071, 2100))):
        kp = P[f"{exp}|ymd"] @ np.array([10000, 100, 1])
        kz = Z[f"{exp}|ymd"] @ np.array([10000, 100, 1])
        common, i_p, i_z = np.intersect1d(kp, kz, return_indices=True)
        w = gs.winter_of(common)
        keep = (w >= y0) & (w <= y1) & np.isin(common // 100 % 100, rg.MONTHS)
        common, i_p, i_z = common[keep], i_p[keep], i_z[keep]
        G = gs.region_G(gs.geostrophic(P[f"{exp}|psl"][i_p].astype(np.float64), lat, lon, True), lat, lon, ts.EU)
        S[exp] = dict(G=G, s=Z[f"{exp}|spd10"][i_z][:, ri].astype(np.float64), winter=gs.winter_of(common))
    dc, df = S["historical"], S["ssp585"]
    selc, self_ = np.ones(len(dc["G"]), bool), np.ones(len(df["G"]), bool)
    cc, cf, _ = ts.quintile_classes(dc["G"], df["G"], selc)
    r = ts.run_decomp(dc, df, selc, self_, cc, cf, 5, "s", rng)
    tot = ts.tr_change(dc, df, selc, self_, "s")
    out = dict(total=float(tot.mean()), between_G=r["between"]["eu_mean"], within_G=r["within"]["eu_mean"],
               between_ci90=r["between"]["ci90"], within_ci90=r["within"]["ci90"],
               mean_G_change_pct=float(np.mean(100 * (df["G"].mean(0) / dc["G"].mean(0) - 1))),
               tr_per_G_class_ctl=r["classes"]["tr_c"], G_class_freq_fut=r["classes"]["freq_f"],
               per_region_total={k: float(v) for k, v in zip(ts.EU, tot)})
    print(f"  {name:18s} dG {out['mean_G_change_pct']:+5.2f}%  transfer {out['total']:+5.2f}%  between-G {out['between_G']:+5.2f}%  within-G {out['within_G']:+5.2f}%", flush=True)
    return name, out


def main():
    rng = np.random.default_rng(2)
    inv = json.load(open(f"{DC}/inventory.json"))
    res = {}
    for fn in sorted(glob.glob(f"{DC}/psl_*.npz")):
        n, o = model(fn, rng)
        o["wind_only"] = bool(inv.get(n, {}).get("wind_only"))
        res[n] = o
    ens = [v for v in res.values() if not v["wind_only"]]
    summ = {k: dict(median=float(np.median([v[k] for v in ens])), p10_p90=[float(x) for x in np.percentile([v[k] for v in ens], [10, 90])])
            for k in ("total", "between_G", "within_G", "mean_G_change_pct")}
    summ["n_models"] = len(ens)
    json.dump(dict(models=res, ensemble=summ), open(OUT, "w"), indent=1)
    print("Ensemble (" + str(len(ens)) + " models, EC-Earth3 excluded): " + "; ".join(f"{k} {v['median']:+.2f} [{v['p10_p90'][0]:+.2f}, {v['p10_p90'][1]:+.2f}]" for k, v in summ.items() if isinstance(v, dict)))
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
