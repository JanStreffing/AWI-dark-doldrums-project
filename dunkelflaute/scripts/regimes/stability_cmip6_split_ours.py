"""
The G split of transfer_stability_cmip6.py (between-G and at-fixed-G parts of the transfer change,
control quintiles of G per region) for our runs, so they can be placed on the CMIP6 relations of
stability_cmip6_analysis.py. Comparisons of region_events.COMPARISONS used in geostrophic.py
(9km, 31km, 100km, eerie_ctl, eerie_hist). Speeds and G exactly as transfer_checks.py:
  spd10 | bilin_global   raw 10 m speed, G from msl bilinear on the global 2.5 deg grid (CMIP6 way);
                         8 European regions and the 22 CMIP_REGIONS
  spd100 | box           log-law 100 m speed, box-mean G (the 13b / transfer_stability.py footing), 8 EU
Bootstrap reduced (NB, default 100).

Usage: python3 stability_cmip6_split_ours.py      (legacy python, login node, a few minutes)
Output: dunkelflaute/results/regimes/stability_cmip6_split_ours.json
"""
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regimes")
os.environ.setdefault("NB", "100")
import transfer_stability as ts          # noqa: E402
import transfer_checks as tc             # noqa: E402
gs, rg, ev = tc.gs, tc.rg, tc.ev
OUT = f"{_ROOT}/dunkelflaute/results/regimes/stability_cmip6_split_ours.json"


def split(A, cn, fn, c0, c1, f0, f1, regions, rng):
    S = {}
    for run, (y0, y1) in ((cn, (c0, c1)), (fn, (f0, f1))):
        k, s, G = A[run]
        w = gs.winter_of(k)
        sel = (w >= y0) & (w <= y1) & np.isfinite(s).all(1) & np.isfinite(G).all(1)
        S[run] = dict(G=G[sel], s=s[sel], winter=w[sel])
    dc, df = S[cn], S[fn]
    ts.EU = list(regions)
    selc, self_ = np.ones(len(dc["G"]), bool), np.ones(len(df["G"]), bool)
    cc, cf, _ = ts.quintile_classes(dc["G"], df["G"], selc)
    r = ts.run_decomp(dc, df, selc, self_, cc, cf, 5, "s", rng)
    tot = ts.tr_change(dc, df, selc, self_, "s")
    return dict(total=float(tot.mean()), between_G=r["between"]["eu_mean"], within_G=r["within"]["eu_mean"],
                within_ci90=r["within"]["ci90"], between_ci90=r["between"]["ci90"],
                mean_G_change_pct=float(np.mean(100 * (df["G"].mean(0) / dc["G"].mean(0) - 1))),
                per_region_within=r["within"]["per_region"], per_region_total={k: float(v) for k, v in zip(regions, tot)})


def main():
    rng = np.random.default_rng(3)
    EU8 = list(gs.EU)
    res = {}
    for name in gs.OURS:
        (cn, c0, c1), (fn, f0, f1) = ev.COMPARISONS[name]
        o = {}
        for var, how, sets in (("spd10", "bilin_global", (("eu8", EU8), ("cmip22", gs.CMIP_REGIONS))),
                               ("spd100", "box", (("eu8", EU8),))):
            Sp = {run: tc.speeds(run, var) for run in (cn, fn)}
            if any(v is None for v in Sp.values()):
                continue
            for tag, regions in sets:
                GG = {run: tc.G_of(run, how, regions) for run in (cn, fn)}
                if any(v is None for v in GG.values()):
                    continue
                A = {}
                for run in (cn, fn):
                    (ks, s), (kg, G) = Sp[run], GG[run]
                    common, i_s, i_g = np.intersect1d(ks, kg, return_indices=True)
                    keep = np.isin(common // 100 % 100, rg.MONTHS)
                    A[run] = (common[keep], s[i_s[keep]][:, [rg.REGIONS.index(r) for r in regions]], G[i_g[keep]])
                o[f"{var}|{how}|{tag}"] = x = split(A, cn, fn, c0, c1, f0, f1, regions, rng)
                print(f"  {name:10s} {var:6s} {how:12s} {tag:6s} total {x['total']:+5.2f} between {x['between_G']:+5.2f} "
                      f"within {x['within_G']:+5.2f} [{x['within_ci90'][0]:+.2f}, {x['within_ci90'][1]:+.2f}] dG {x['mean_G_change_pct']:+5.2f}", flush=True)
        res[name] = o
    json.dump(dict(comparisons={k: ev.COMPARISONS[k] for k in gs.OURS}, runs=res), open(OUT, "w"), indent=1)
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
