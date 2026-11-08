"""
How robust is the 1950C -> 2080C change in Mockert Dunkelflaute events (Germany, Sep-Mar)?

Reads only the finished combined-CF series in dunkelflaute/results/ (1950c, 2080c,
validation_era5_2015_2026), so it runs in seconds without touching the TB-scale
reconstructions. Re-detects events with the same window-expansion rule as
dart/compute_germany_dunkelflaute_dart_mockert.py, then reports:
  - significance of the change: Welch t, Mann-Whitney, exact conditional Poisson test,
    and a year-block bootstrap of the 2080C/1950C ratio (whole winters resampled)
  - threshold (0.04..0.10) x window (24/48/72 h) sensitivity, all five sources
  - the same with Germany's 2024 weights (57.7/36.9/5.4) instead of Mockert's 2018 mix
  - a Li-style count: capacity-weighted wind CF < 0.2 and solar CF < 0.2, run > 24 h
  - low tails of the 48 h mean, component climatology, event onset month
  - recon-vs-real ERA5: per-winter count correlation and CF bias binned by real 48 h CF

Output: dunkelflaute/results/robustness/robustness_summary.json (and a printed summary).
"""
import json
import os

import numpy as np
from scipy import stats

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
os.chdir(_ROOT)

RESULTS = "dunkelflaute/results"
OUT = f"{RESULTS}/robustness/robustness_summary.json"
W_MOCKERT = dict(solar=0.44, on=0.50, off=0.06)   # Germany 2018 (IRENA), Mockert et al. 2023
W_2024 = dict(solar=0.577, on=0.369, off=0.054)    # Germany 2024 (BNetzA), proposal Sec. 5.4
DT_H = 3.0
THRESHOLDS = (0.04, 0.05, 0.06, 0.07, 0.08, 0.10)
WINDOWS_H = (24, 48, 72)
N_BOOT = 20000


def load():
    out = {}
    for run in ("1950c", "2080c"):
        d = np.load(f"{RESULTS}/{run}/combined_cf_{run}_sepmar_mockertweights.npz")
        out[run] = dict(winter=d["winter"], time=d["time"].astype("datetime64[s]"),
                        solar=d["solar"].astype(float), on=d["on"], off=d["off"])
    d = np.load(f"{RESULTS}/validation_era5_2015_2026/combined_cf_2015_2026_sepmar_mockertweights.npz",
                allow_pickle=True)
    for src, keys in (("era5_recon", ("cf_solar_recon", "cf_on_recon", "cf_off_recon")),
                      ("era5_real", ("cf_solar_real", "cf_on_real", "cf_off_real")),
                      ("smard", ("cf_solar_smard", "cf_onshore_smard", "cf_offshore_smard"))):
        out[src] = dict(winter=d["winter"], time=d["time"].astype("datetime64[s]"),
                        solar=d[keys[0]], on=d[keys[1]], off=d[keys[2]])
    return out


def combined(ds, w):
    return w["solar"] * ds["solar"] + w["on"] * ds["on"] + w["off"] * ds["off"]


def trailing_mean(x, win):
    m = np.full(len(x), np.nan)
    if len(x) >= win:
        c = np.cumsum(np.insert(x, 0, 0.0))
        m[win - 1:] = (c[win:] - c[:-win]) / win
    return m


def mockert_events(roll, win, thr):
    """Window expansion (Mockert Sec. 2.4): every below-threshold mean pulls its whole window in."""
    n = len(roll)
    flagged = np.zeros(n, bool)
    for i in np.where(np.nan_to_num(roll, nan=np.inf) < thr)[0]:
        flagged[max(0, i - win + 1):i + 1] = True
    events, start = [], None
    for i in range(n):
        if flagged[i] and start is None:
            start = i
        elif not flagged[i] and start is not None:
            events.append((start, i)); start = None
    if start is not None:
        events.append((start, n))
    return flagged, events


def per_winter(ds, win=16, thr=0.06, w=W_MOCKERT):
    comb = combined(ds, w)
    res = {}
    for y in np.unique(ds["winter"]):
        m = ds["winter"] == y
        roll = trailing_mean(comb[m], win)
        flagged, events = mockert_events(roll, win, thr)
        res[int(y)] = dict(n=len(events), dur_d=[(e - s) * DT_H / 24 for s, e in events],
                           flagged=float(flagged.mean()), roll=roll,
                           starts=[ds["time"][m][s] for s, _ in events])
    return res


def li_like(ds, thr=0.2, min_h=24.0):
    wind = (W_MOCKERT["on"] * ds["on"] + W_MOCKERT["off"] * ds["off"]) / (W_MOCKERT["on"] + W_MOCKERT["off"])
    counts = []
    for y in np.unique(ds["winter"]):
        m = ds["winter"] == y
        below = (wind[m] < thr) & (ds["solar"][m] < thr)
        n, run = 0, 0
        for b in np.append(below, False):
            if b:
                run += 1
            else:
                n += run * DT_H > min_h
                run = 0
        counts.append(n)
    return float(np.mean(counts))


def counts(res):
    return np.array([r["n"] for r in res.values()])


def main():
    data = load()
    base = {k: per_winter(v) for k, v in data.items()}
    out = {"weights_mockert": W_MOCKERT, "weights_2024": W_2024, "sources": {}}

    for k, r in base.items():
        n = counts(r)
        dur = np.concatenate([x["dur_d"] for x in r.values()])
        rolls = np.concatenate([x["roll"] for x in r.values()])
        rolls = rolls[np.isfinite(rolls)]
        st = np.concatenate([np.array(x["starts"], dtype="datetime64[s]") for x in r.values()])
        mon = st.astype("datetime64[M]").astype(int) % 12 + 1
        ds = data[k]
        out["sources"][k] = dict(
            n_winters=len(n), events_per_winter=float(n.mean()), sd=float(n.std(ddof=1)),
            per_winter={y: x["n"] for y, x in r.items()},
            flagged_pct=100 * float(np.mean([x["flagged"] for x in r.values()])),
            dur_median_d=float(np.median(dur)), dur_p90_d=float(np.percentile(dur, 90)), dur_max_d=float(dur.max()),
            roll48_p05=float(np.percentile(rolls, 5)), roll48_p10=float(np.percentile(rolls, 10)),
            roll48_sd=float(rolls.std()),
            mean_cf=dict(solar=float(np.nanmean(ds["solar"])), on=float(np.nanmean(ds["on"])),
                         off=float(np.nanmean(ds["off"])), combined=float(np.nanmean(combined(ds, W_MOCKERT)))),
            frac_onshore_cf_below_0p05=float(np.nanmean(ds["on"] < 0.05)),
            onset_month_frac={int(m): float(np.mean(mon == m)) for m in (9, 10, 11, 12, 1, 2, 3)},
            events_per_winter_2024_weights=float(counts(per_winter(ds, w=W_2024)).mean()),
            li_like_events_per_winter=li_like(ds),
        )

    a, b = counts(base["1950c"]), counts(base["2080c"])
    rng = np.random.default_rng(0)
    boot = np.array([rng.choice(b, len(b)).mean() / rng.choice(a, len(a)).mean() for _ in range(N_BOOT)])
    out["change_1950c_to_2080c"] = dict(
        ratio=float(b.mean() / a.mean()),
        welch_p=float(stats.ttest_ind(b, a, equal_var=False).pvalue),
        mannwhitney_p=float(stats.mannwhitneyu(b, a, alternative="two-sided").pvalue),
        poisson_exact_p=float(stats.binomtest(int(b.sum()), int(a.sum() + b.sum()), len(b) / (len(a) + len(b))).pvalue),
        dispersion_index=dict(c1950=float(a.var(ddof=1) / a.mean()), c2080=float(b.var(ddof=1) / b.mean())),
        bootstrap_ratio_p05_p50_p95=[float(x) for x in np.percentile(boot, [5, 50, 95])],
        bootstrap_p_ratio_le_1=float(np.mean(boot <= 1)),
        c1950_split_half=[float(a[:8].mean()), float(a[8:].mean())],
    )

    out["sensitivity_events_per_winter"] = {
        f"{win_h}h_{thr:.2f}": {k: float(counts(per_winter(v, win=win_h // 3, thr=thr)).mean()) for k, v in data.items()}
        for win_h in WINDOWS_H for thr in THRESHOLDS}

    # recon vs real ERA5 (same weather, different wind/solar fields)
    pr, pe = counts(base["era5_recon"]), counts(base["era5_real"])
    d = np.load(f"{RESULTS}/validation_era5_2015_2026/combined_cf_2015_2026_sepmar_mockertweights.npz", allow_pickle=True)
    rr, re = d["roll_recon"], d["roll_real"]
    ok = np.isfinite(rr) & np.isfinite(re)
    on_r, on_e = d["cf_on_recon"], d["cf_on_real"]
    out["recon_vs_real_era5"] = dict(
        per_winter_count_r=float(np.corrcoef(pr, pe)[0, 1]),
        roll48_ratio_by_real_bin={f"{lo}-{hi}": float(np.mean(rr[ok & (re >= lo) & (re < hi)]) /
                                                      np.mean(re[ok & (re >= lo) & (re < hi)]))
                                  for lo, hi in ((0, .06), (.06, .1), (.1, .2), (.2, 1))},
        onshore_cf_ratio_by_real_bin={f"{lo}-{hi}": float(np.mean(on_r[(on_e >= lo) & (on_e < hi)]) /
                                                          np.mean(on_e[(on_e >= lo) & (on_e < hi)]))
                                      for lo, hi in ((0, .05), (.05, .15), (.15, .4), (.4, 1.01))},
    )

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)

    c = out["change_1950c_to_2080c"]
    print(f"1950C {a.mean():.2f} -> 2080C {b.mean():.2f} events/winter, ratio {c['ratio']:.2f}, "
          f"Poisson p={c['poisson_exact_p']:.3f}, bootstrap 90% {c['bootstrap_ratio_p05_p50_p95'][0]:.2f}-"
          f"{c['bootstrap_ratio_p05_p50_p95'][2]:.2f}")
    print("ratio 2080C/1950C at 48 h:", {t: round(out["sensitivity_events_per_winter"][f"48h_{t:.2f}"]["2080c"] /
                                                  out["sensitivity_events_per_winter"][f"48h_{t:.2f}"]["1950c"], 2)
                                         for t in THRESHOLDS})
    print("Li-like:", {k: round(v["li_like_events_per_winter"], 2) for k, v in out["sources"].items()})
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
