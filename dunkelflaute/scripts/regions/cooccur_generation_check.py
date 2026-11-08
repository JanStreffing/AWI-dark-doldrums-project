"""
Quick check of low-wind co-occurrence in observed national wind generation against ERA5, for the
three of the 8 European regions that have generation downloads in this repository:
  DE  SMARD onshore CF (smard_validation/data/cf_smard_fullyear_2015_2026.csv, 3-hourly)
  PL  energy-charts onshore generation (poland_validation/data/energycharts_pl_production_*.csv)
  SE  energy-charts onshore generation (sweden_validation/data/energycharts_se_production_*.csv),
      compared with the SC box (56-64 N, 5-25 E, which also covers southern Norway and Finland)
PL and SE CF = generation / installed onshore capacity (energycharts_*_installed_power.csv, annual
values taken as year-end and interpolated linearly in time). All series sampled at 00/06/12/18 UTC
(the 6-hourly footing), Sep-Mar winters 2016-2026.

Events as cooccur_events.py: 48 h trailing mean (8 steps), window expansion, threshold matched to 5
events per winter in each series itself. Statistics: Jaccard index of the event masks per pair,
P(N >= 2) and P(N = 3) of the three, and the correlation of the 48 h means. ERA5 (cooccur_era5_wind.npz,
regions DE, PL, SC) over the same winters, and over all winters 1980-2026.

Usage: python3 cooccur_generation_check.py      (module python3, login node, < 1 min)
Output: dunkelflaute/results/regions/cooccur_generation_check.json
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
import cooccur_events as ce                                   # noqa: E402

DD = f"{_ROOT}/dunkelflaute"
OUT = f"{DD}/results/regions/cooccur_generation_check.json"
Y0, Y1 = 2016, 2026


def ec_cf(cc):
    fs = sorted(glob.glob(f"{DD}/{dict(pl='poland', se='sweden')[cc]}_validation/data/energycharts_{cc}_production_*.csv"))
    g = pd.concat([pd.read_csv(f, usecols=["time_utc", "Wind onshore"]) for f in fs])
    g["t"] = pd.to_datetime(g["time_utc"])
    g = g.drop_duplicates("t").set_index("t")["Wind onshore"].astype(float)
    cap = pd.read_csv(f"{DD}/{dict(pl='poland', se='sweden')[cc]}_validation/data/energycharts_{cc}_installed_power.csv")
    cap = cap[cap["time"].astype(str).str.fullmatch(r"\d{4}")]
    tc = pd.to_datetime(cap["time"].astype(int).astype(str) + "-12-31").values.astype("datetime64[s]").astype(float)
    tt = g.index.values.astype("datetime64[s]").astype(float)
    gw = np.interp(tt, tc, cap["Wind onshore"].astype(float).values)
    return pd.Series(g.values / (gw * 1000.0), index=g.index)


def de_cf():
    s = pd.read_csv(f"{DD}/smard_validation/data/cf_smard_fullyear_2015_2026.csv")
    t = pd.to_datetime(s["timestamp_utc"]).dt.tz_convert(None)
    return pd.Series(s["cf_onshore_smard"].values, index=t)


def on6h(series, t6):
    return np.column_stack([s.reindex(t6).values for s in series])


def stats(ws):
    thr, rate = ce.match(ws, 8)
    ms = ce.masks(ws, thr, 8)
    m = np.concatenate(ms)
    n = m.sum(1)
    J = ce.jaccard(ms)
    x = np.concatenate([np.column_stack([ce.trailing_mean(w[:, j], 8) for j in range(w.shape[1])])[7:] for w in ws])
    ok = np.isfinite(x).all(1)
    C = np.corrcoef(x[ok].T)
    pairs = {"DE-PL": (0, 1), "DE-SE": (0, 2), "PL-SE": (1, 2)}
    return dict(rate=rate.tolist(), jaccard={k: float(J[a, b]) for k, (a, b) in pairs.items()},
                corr48h={k: float(C[a, b]) for k, (a, b) in pairs.items()},
                P_ge2=float(np.mean(n >= 2)), P_eq3=float(np.mean(n == 3)), P_ge1=float(np.mean(n >= 1)))


def winters(X, t6, y0, y1):
    w = ce.ev.winter_of(t6)
    out = []
    for y in range(y0, y1 + 1):
        m = w == y
        if m.sum() >= ce.ev.MIN_STEPS and np.isfinite(X[m]).mean() > 0.97:
            x = X[m].copy()
            for j in range(x.shape[1]):                     # short gaps: linear in time
                bad = ~np.isfinite(x[:, j])
                if bad.any():
                    x[bad, j] = np.interp(np.where(bad)[0], np.where(~bad)[0], x[~bad, j])
            out.append((y, x))
    return out


def main():
    z = np.load(ce.ERA5)
    names = list(z["regions"])
    t6 = z["time6"].astype("datetime64[s]")
    E = z["wind6"][:, [names.index(r) for r in ("DE", "PL", "SC")]].astype(float)
    G = on6h([de_cf(), ec_cf("pl"), ec_cf("se")], pd.DatetimeIndex(t6))
    wg = winters(G, t6, Y0, Y1)
    ys = [y for y, _ in wg]
    we = [x for y, x in winters(E, t6, Y0, Y1) if y in ys]
    res = dict(winters=ys, generation=stats([x for _, x in wg]), era5_same_winters=stats(we),
               era5_1980_2026=stats([x for _, x in winters(E, t6, 1980, 2026)]))
    json.dump(res, open(OUT, "w"), indent=1)
    print(f"winters {ys}")
    for k, v in res.items():
        if k == "winters":
            continue
        print(f"  {k:20s} J " + " ".join(f"{a} {b:.3f}" for a, b in v["jaccard"].items())
              + "  r48h " + " ".join(f"{a} {b:.2f}" for a, b in v["corr48h"].items())
              + f"  P(N>=2) {v['P_ge2']:.4f} P(N=3) {v['P_eq3']:.4f} rate {np.round(v['rate'], 2)}")
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
