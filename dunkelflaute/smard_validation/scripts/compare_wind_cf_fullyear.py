"""
Compare Germany wind CF three ways, full year, 2015-2026:
  recon: ERA5 log-law-reconstructed 100m wind -> wind_cf (dunkelflaute)
  real:  ERA5's own genuine 100m wind -> wind_cf (dunkelflaute)
  smard: real generation / year-specific installed capacity (ground truth)

recon-vs-smard and real-vs-smard both matter: real-vs-smard isolates how much
gap is inherent to "wind speed -> simplified power curve" (turbine downtime,
curtailment, wake losses, grid constraints -- things no wind-speed-only model
can capture), while recon-vs-smard additionally includes the log-law
reconstruction's own error on top of that.
"""
import json

import numpy as np
import pandas as pd

ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
ERA5_PATH = f"{ROOT}/dunkelflaute/data/germany_era5/wind_cf_era5_fullyear_2015_2026.npz"
SMARD_PATH = f"{ROOT}/dunkelflaute/smard_validation/data/cf_smard_fullyear_2015_2026.csv"
OUT_PATH = f"{ROOT}/dunkelflaute/smard_validation/data/wind_cf_comparison_fullyear_2015_2026.json"

era5 = np.load(ERA5_PATH, allow_pickle=True)
era5_time = pd.to_datetime([str(t) for t in era5["time"]]).tz_localize("UTC")
era5_df = pd.DataFrame({
    "cf_on_recon": era5["cf_on_recon"], "cf_off_recon": era5["cf_off_recon"],
    "cf_on_real": era5["cf_on_real"], "cf_off_real": era5["cf_off_real"],
}, index=era5_time)

smard = pd.read_csv(SMARD_PATH, parse_dates=["timestamp_utc"]).set_index("timestamp_utc")

df = era5_df.join(smard, how="inner")
print(f"joined n={len(df)}  {df.index.min()} .. {df.index.max()}")

# flag known capacity-step artifact years (offshore buildout) rather than silently drop
df["offshore_capacity_step_year"] = df.index.year.isin([2015, 2016, 2017])


def metrics(pred, true):
    e = pred - true
    rmse = float(np.sqrt(np.mean(e ** 2)))
    bias = float(np.mean(e))
    mae = float(np.mean(np.abs(e)))
    r = float(np.corrcoef(pred, true)[0, 1])
    sst = float(np.sum((true - true.mean()) ** 2))
    r2 = float(1 - np.sum(e ** 2) / sst) if sst > 0 else float("nan")
    return dict(rmse=rmse, bias=bias, mae=mae, r=r, r2=r2, mean_true=float(true.mean()), n=int(len(true)))


results = {"pooled": {}, "by_year": {}}
pairs = {
    "onshore_recon_vs_smard": ("cf_on_recon", "cf_onshore_smard", None),
    "onshore_real_vs_smard": ("cf_on_real", "cf_onshore_smard", None),
    "offshore_recon_vs_smard": ("cf_off_recon", "cf_offshore_smard", None),
    "offshore_real_vs_smard": ("cf_off_real", "cf_offshore_smard", None),
    "offshore_recon_vs_smard_excl_stepyears": ("cf_off_recon", "cf_offshore_smard", "excl_step"),
    "offshore_real_vs_smard_excl_stepyears": ("cf_off_real", "cf_offshore_smard", "excl_step"),
}

for name, (predcol, truecol, filt) in pairs.items():
    sub = df
    if filt == "excl_step":
        sub = df[~df["offshore_capacity_step_year"]]
    valid = sub[[predcol, truecol]].dropna()
    results["pooled"][name] = metrics(valid[predcol].values, valid[truecol].values)

for year, g in df.groupby(df.index.year):
    yr_res = {}
    for name, (predcol, truecol, filt) in pairs.items():
        if filt == "excl_step":
            continue
        valid = g[[predcol, truecol]].dropna()
        if len(valid) < 10:
            continue
        yr_res[name] = metrics(valid[predcol].values, valid[truecol].values)
    results["by_year"][int(year)] = yr_res

json.dump(results, open(OUT_PATH, "w"), indent=2)

print("\n=== POOLED (2015-01 .. 2026-08) ===")
for name, m in results["pooled"].items():
    print(f"{name:42s} RMSE {m['rmse']:.3f}  bias {m['bias']:+.3f}  r {m['r']:.3f}  R2 {m['r2']:.3f}  n={m['n']}")

print("\n=== BY YEAR: onshore recon vs smard ===")
for year, yr_res in sorted(results["by_year"].items()):
    m = yr_res.get("onshore_recon_vs_smard")
    if m:
        print(f"{year}  RMSE {m['rmse']:.3f}  bias {m['bias']:+.3f}  r {m['r']:.3f}")

print("\n=== BY YEAR: offshore recon vs smard (raw, includes capacity-step years) ===")
for year, yr_res in sorted(results["by_year"].items()):
    m = yr_res.get("offshore_recon_vs_smard")
    if m:
        print(f"{year}  RMSE {m['rmse']:.3f}  bias {m['bias']:+.3f}  r {m['r']:.3f}")

print(f"\nSaved {OUT_PATH}")
