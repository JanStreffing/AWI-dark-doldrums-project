"""
The G split of transfer_stability_cmip6.py (between-G and at-fixed-G parts of the transfer change,
control quintiles of G per region) on the 22 regions of geostrophic.CMIP_REGIONS, for comparison
with the monthly stability diagnostics of stability_cmip6_extract.py. Same inputs, periods and
classes as transfer_stability_cmip6.py; only the region list differs (the 8 European regions are
re-done as a check against transfer_stability_cmip6.json). Bootstrap reduced (NB, default 100).

Usage: python3 stability_cmip6_split22.py      (legacy python, login node, a few minutes)
Output: dunkelflaute/results/regimes/stability_cmip6_split22.json
"""
import glob
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regimes")
os.environ.setdefault("NB", "100")
import transfer_stability as ts          # noqa: E402
gs, rg = ts.gs, ts.rg
DC = f"{_ROOT}/dunkelflaute/data/cmip6"
OUT = f"{_ROOT}/dunkelflaute/results/regimes/stability_cmip6_split22.json"


def split(P, Z, regions, rng):
    lat, lon = P["lat"], P["lon"]
    ri = [rg.REGIONS.index(r) for r in regions]
    S = {}
    for exp, (y0, y1) in (("historical", (1951, 1980)), ("ssp585", (2071, 2100))):
        kp = P[f"{exp}|ymd"] @ np.array([10000, 100, 1])
        kz = Z[f"{exp}|ymd"] @ np.array([10000, 100, 1])
        common, i_p, i_z = np.intersect1d(kp, kz, return_indices=True)
        w = gs.winter_of(common)
        keep = (w >= y0) & (w <= y1) & np.isin(common // 100 % 100, rg.MONTHS)
        common, i_p, i_z = common[keep], i_p[keep], i_z[keep]
        G = gs.region_G(gs.geostrophic(P[f"{exp}|psl"][i_p].astype(np.float64), lat, lon, True), lat, lon, regions)
        S[exp] = dict(G=G, s=Z[f"{exp}|spd10"][i_z][:, ri].astype(np.float64), winter=gs.winter_of(common))
    dc, df = S["historical"], S["ssp585"]
    ts.EU = list(regions)                                  # winter_sums / tr_change loop over ts.EU
    selc, self_ = np.ones(len(dc["G"]), bool), np.ones(len(df["G"]), bool)
    cc, cf, _ = ts.quintile_classes(dc["G"], df["G"], selc)
    r = ts.run_decomp(dc, df, selc, self_, cc, cf, 5, "s", rng)
    tot = ts.tr_change(dc, df, selc, self_, "s")
    return dict(total=float(tot.mean()), between_G=r["between"]["eu_mean"], within_G=r["within"]["eu_mean"],
                within_ci90=r["within"]["ci90"],
                mean_G_change_pct=float(np.mean(100 * (df["G"].mean(0) / dc["G"].mean(0) - 1))),
                per_region_total={k: float(v) for k, v in zip(regions, tot)},
                per_region_within={k: float(v) for k, v in zip(regions, r["within"]["per_region"].values())}
                if isinstance(r["within"].get("per_region"), dict) else None)


def main():
    rng = np.random.default_rng(2)
    inv = json.load(open(f"{DC}/inventory.json"))
    res = {}
    for fn in sorted(glob.glob(f"{DC}/psl_*.npz")):
        name = os.path.basename(fn)[4:-4]
        P, Z = np.load(fn), np.load(f"{DC}/daily_cmip6_{name}.npz")
        o = dict(cmip22=split(P, Z, gs.CMIP_REGIONS, rng), eu8=split(P, Z, ts_EU8, rng),
                 wind_only=bool(inv.get(name, {}).get("wind_only")))
        res[name] = o
        a, b = o["cmip22"], o["eu8"]
        print(f"  {name:18s} 22: total {a['total']:+5.2f} between {a['between_G']:+5.2f} within {a['within_G']:+5.2f} | "
              f"8: total {b['total']:+5.2f} between {b['between_G']:+5.2f} within {b['within_G']:+5.2f}", flush=True)
    ens = {}
    for g in ("cmip22", "eu8"):
        for k in ("total", "between_G", "within_G", "mean_G_change_pct"):
            x = [v[g][k] for v in res.values() if not v["wind_only"]]
            ens[f"{g}_{k}"] = dict(median=float(np.median(x)), p10_p90=[float(q) for q in np.percentile(x, [10, 90])])
    json.dump(dict(regions22=gs.CMIP_REGIONS, models=res, ensemble=ens), open(OUT, "w"), indent=1)
    print("\n".join(f"  {k}: {v['median']:+.2f} [{v['p10_p90'][0]:+.2f}, {v['p10_p90'][1]:+.2f}]" for k, v in ens.items()))
    print(f"Saved {OUT}")


ts_EU8 = list(ts.EU)

if __name__ == "__main__":
    main()
