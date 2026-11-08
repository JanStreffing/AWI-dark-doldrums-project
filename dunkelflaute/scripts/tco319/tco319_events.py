"""
Mockert events from the TCo319-DART (~31 km) Germany series, and the 100 / 31 / 9 km comparison
on a common 6-hourly footing.

Sources, all 6-hourly with the 48 h window = 8 steps, Mockert 2018 weights:
  31 km   compute_germany_tco319.py series: v7 solar + calibrated log law ("qm") and + uncalibrated ("ll")
  9 km    dunkelflaute/results/{1950c,2080c} 3 h series degraded to 6 h (wind sampled at 00/06/12/18,
          solar summed over 6 h), as in analysis/sampling_3h_vs_6h.py; v5 solar + uncalibrated log
          law, plus the calibrated wind from results/wind_calibration/dart_calibrated_*.npz
  100 km  TCO95 perfect-model series degraded the same way: native fields and the calibrated recipe
Periods: control winters 1953-1969 (the 9 km spin-up rule) and 1951-1969 (31 km, all); future
2083-2092 (the 9 km span) and 2081-2099 (31 and 100 km).

Usage: python3 tco319_events.py
Output: dunkelflaute/results/tco319/tco319_events_summary.json
"""
import json
import os
import sys

import numpy as np
from scipy import stats

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/analysis")
from robustness_mockert_events import W_MOCKERT, mockert_events, trailing_mean  # noqa: E402

os.chdir(_ROOT)
OUT = "dunkelflaute/results/tco319/tco319_events_summary.json"
THRESHOLDS = (0.04, 0.05, 0.06, 0.07, 0.08, 0.10)
WIN6 = 8
H3 = np.timedelta64(3, "h")


def winter_of(t):
    t1 = t - np.timedelta64(1, "s")
    m = (t1.astype("datetime64[M]").astype(int) % 12) + 1
    y = t1.astype("datetime64[Y]").astype(int) + 1970
    return np.where(m >= 9, y + 1, y)


def to_6h(t, comps):
    """3 h series -> 6 h: instantaneous wind at 00/06/12/18, solar mean of the two 3 h values."""
    hour = (t - t.astype("datetime64[D]")) // np.timedelta64(1, "h")
    pos = {v: i for i, v in enumerate(t)}
    end = np.array([i for i in np.where(hour % 6 == 0)[0] if t[i] - H3 in pos])
    prev = np.array([pos[t[i] - H3] for i in end])
    out = dict(time=t[end], on=comps["on"][end], off=comps["off"][end],
               solar=0.5 * (comps["solar"][end] + comps["solar"][prev]))
    return out


def counts(series, thr, w0, w1, min_steps=800):
    comb = W_MOCKERT["solar"] * series["solar"] + W_MOCKERT["on"] * series["on"] + W_MOCKERT["off"] * series["off"]
    w = winter_of(series["time"])
    out = {}
    for y in range(w0, w1 + 1):
        m = w == y
        if m.sum() < min_steps:
            continue
        _, ev = mockert_events(trailing_mean(comb[m], WIN6), WIN6, thr)
        out[y] = len(ev)
    return out


def change(ctl, fut):
    a, b = np.array(list(ctl.values())), np.array(list(fut.values()))
    p = stats.binomtest(int(b.sum()), int(a.sum() + b.sum()), len(b) / (len(a) + len(b))).pvalue if a.sum() + b.sum() else np.nan
    rng = np.random.default_rng(0)
    boot = np.array([rng.choice(b, len(b)).mean() / max(rng.choice(a, len(a)).mean(), 1e-9) for _ in range(10000)])
    return dict(control=float(a.mean()), future=float(b.mean()), ratio=float(b.mean() / a.mean()) if a.mean() else None,
                n_control=len(a), n_future=len(b), poisson_p=float(p),
                bootstrap_ratio_p05_p95=[float(x) for x in np.percentile(boot, [5, 95])])


def load_sources():
    src = {}
    for run in ("ctl1950d", "ssp585d"):
        d = np.load(f"dunkelflaute/data/tco319/{run}_germany_6h.npz")
        t = d["time"].astype("datetime64[s]")
        src[f"31km_qm_{run}"] = dict(time=t, solar=d["solar_v7"], on=d["on_qm"], off=d["off_qm"])
        src[f"31km_ll_{run}"] = dict(time=t, solar=d["solar_v7"], on=d["on_ll"], off=d["off_ll"])
    for run in ("1950c", "2080c"):
        d = np.load(f"dunkelflaute/results/{run}/combined_cf_{run}_sepmar_mockertweights.npz")
        t = d["time"].astype("datetime64[s]")
        src[f"9km_ll_{run}"] = to_6h(t, dict(solar=d["solar"].astype(float), on=d["on"], off=d["off"]))
        c = np.load(f"dunkelflaute/results/wind_calibration/dart_calibrated_{run}.npz")
        tc = c["time"].astype("datetime64[s]")
        assert np.array_equal(tc, t), "calibrated DART series must share the stored time axis"
        src[f"9km_qm_{run}"] = to_6h(t, dict(solar=d["solar"].astype(float), on=c["on_cal"], off=c["off_cal"]))
    for run in ("ctl1950g", "hi1950g"):
        d = np.load(f"dunkelflaute/data/tco95/{run}_germany_3h.npz")
        t = d["time"].astype("datetime64[s]")
        src[f"100km_native_{run}"] = to_6h(t, dict(solar=d["solar_native"], on=d["on_native"], off=d["off_native"]))
        src[f"100km_qm_{run}"] = to_6h(t, dict(solar=d["solar_v5"], on=d["on_qm"], off=d["off_qm"]))
    return src


def main():
    src = load_sources()
    pairs = {
        "9km_ll":     ("9km_ll_1950c", (1953, 1969), "9km_ll_2080c", (2083, 2092)),
        "9km_qm":     ("9km_qm_1950c", (1953, 1969), "9km_qm_2080c", (2083, 2092)),
        "31km_qm_9kmspans": ("31km_qm_ctl1950d", (1953, 1969), "31km_qm_ssp585d", (2083, 2092)),
        "31km_qm":    ("31km_qm_ctl1950d", (1953, 1969), "31km_qm_ssp585d", (2081, 2099)),
        "31km_qm_allctl": ("31km_qm_ctl1950d", (1951, 1969), "31km_qm_ssp585d", (2081, 2099)),
        "31km_ll":    ("31km_ll_ctl1950d", (1953, 1969), "31km_ll_ssp585d", (2081, 2099)),
        "100km_native": ("100km_native_ctl1950g", (1953, 2099), "100km_native_hi1950g", (2081, 2099)),
        "100km_qm":   ("100km_qm_ctl1950g", (1953, 2099), "100km_qm_hi1950g", (2081, 2099)),
        "100km_qm_9kmspans": ("100km_qm_ctl1950g", (1953, 1969), "100km_qm_hi1950g", (2083, 2092)),
    }
    out = dict(weights=W_MOCKERT, window_steps_6h=WIN6, comparisons={})
    for name, (c, (c0, c1), f, (f0, f1)) in pairs.items():
        out["comparisons"][name] = {f"{thr:.2f}": change(counts(src[c], thr, c0, c1), counts(src[f], thr, f0, f1))
                                    for thr in THRESHOLDS}
    # mean CF climatology per source (Sep-Mar), to see where resolutions differ
    clim = {}
    for k, s in src.items():
        clim[k] = {c: float(np.mean(s[c])) for c in ("solar", "on", "off")}
    out["mean_cf"] = clim
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(out, open(OUT, "w"), indent=1)

    print("Future / control events per winter, 6-hourly, by threshold (control -> future, ratio, Poisson p):")
    for name in pairs:
        row = out["comparisons"][name]
        print(f"  {name:20s} " + "  ".join(f"{t}: {r['control']:.2f}->{r['future']:.2f} ({r['ratio']:.2f}, p={r['poisson_p']:.3f})"
                                          for t, r in row.items() if t in ("0.04", "0.06", "0.08", "0.10")))
    print("\nMean Sep-Mar CF (solar / onshore / offshore):")
    for k, c in clim.items():
        print(f"  {k:24s} {c['solar']:.3f} / {c['on']:.3f} / {c['off']:.3f}")
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
