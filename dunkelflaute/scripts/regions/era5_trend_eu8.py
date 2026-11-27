"""ERA5 low-wind event trend for the 8 European regions, same estimator as era5_trend.py (methods_notes 15d).
Usage: python3 era5_trend_eu8.py (legacy python). Output: dunkelflaute/results/regions/era5_trend_eu8.json"""
import sys, json, numpy as np
sys.path.insert(0, "/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/scripts/regions")
import era5_trend as t
R = t._ROOT
EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
e = np.load(f"{R}/dunkelflaute/data/regions/cooccur_era5_wind.npz", allow_pickle=True)
ereg = [str(x) for x in e["regions"]]
day = e["day"].astype("datetime64[D]")
ymd = np.stack([day.astype("datetime64[Y]").astype(int) + 1970, day.astype("datetime64[M]").astype(int) % 12 + 1], 1)
ewin = t.winters_of(ymd)
T_era = t.winter_T_era5(t.era5_gmst_monthly()); years = [y for y in range(1980, 2027) if y in T_era]
y0, y1 = years[0], years[-1]
Er, Eu = t.field(e["wind_day"], ewin, [ereg.index(r) for r in EU], y0, y1)
lnE = np.log(Er.mean(0)); s_era = t.slope(lnE, [T_era[y] for y in years])
gt = json.load(open(f"{R}/dunkelflaute/data/regions/global_t2m_tco95.json"))
out = {"era5_per_K": s_era, "n_regions": len(Eu)}
for run in ("hi1950g", "ctl1950g"):
    z = np.load(f"{R}/dunkelflaute/data/cmip6/daily_tco95_{run}.npz"); mreg = [str(x) for x in z["regions"]]
    cols = [mreg.index(r) for r in EU]; win = t.winters_of(z["run|ymd"])
    if run == "hi1950g":
        T = t.winter_T_model({int(k): v for k, v in gt["tco95_hi1950g"].items()})
        Rr, _ = t.field(z["run|wind"], win, cols, y0, y1)
        out["hi_same_winters_per_K"] = t.slope(np.log(Rr.mean(0)), [T[y] for y in years])
    else:
        n = len(years); Tx = [T_era[y] for y in years]; null = []
        for a in range(1953, 2099 - n + 2):
            Rr, u = t.field(z["run|wind"], win, cols, a, a + n - 1)
            if len(u) == len(cols):
                null.append(t.slope(np.log(Rr.mean(0)), Tx))
        null = np.array(null)
        out["null_5_50_95"] = np.percentile(null, [5, 50, 95]).tolist(); out["null_sd"] = float(null.std())
        out["p_one_sided"] = float((null >= s_era).mean())
print(json.dumps(out, indent=1))
json.dump(out, open(f"{R}/dunkelflaute/results/regions/era5_trend_eu8.json", "w"), indent=1)
