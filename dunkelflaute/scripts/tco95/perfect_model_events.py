"""
Perfect-model test, step 2: what the TCO95 series (perfect_model_series.py) say about the
9 km Dunkelflaute result.

Mockert events (48 h trailing mean of combined CF < threshold, window expansion, per Sep-Mar
winter, Mockert 2018 weights) for every source combination:
  native      native 100 m wind + native ssrd              (the model's truth)
  recon_ll    log law (ERA5 z0, Charnock) + v5 solar        (the uncalibrated DART recipe)
  recon_qm    calibrated log law + v5 solar                 (the calibrated DART recipe)
  wind_ll / wind_qm / solar_v5   one component reconstructed, the other native, to attribute errors
Periods: control = ctl1950g winters 1953-2099 (constant 1950, 147 winters), and its 1953-1969
subset (the 9 km 1950C span); future = hi1950g winters 2083-2092 (the 9 km 2080C span) and
2081-2099.

Answers:
  1. does the reconstruction error change with climate (mean CF bias, control vs future)?
  2. is the future/control event ratio the same from reconstructed fields as from native ones?
  3. natural variability: in the 147-winter control, how often do a 17-winter and a separate
     10-winter block of consecutive winters differ by a ratio as large as the 9 km result (1.58)?
  4. time evolution: native events per winter by decade along hi1950g, against the control.

Usage: python3 perfect_model_events.py
Output: dunkelflaute/results/tco95/perfect_model_summary.json
"""
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/analysis")
from robustness_mockert_events import W_MOCKERT, mockert_events, trailing_mean  # noqa: E402

os.chdir(_ROOT)
DATA = "dunkelflaute/data/tco95"
OUT = "dunkelflaute/results/tco95/perfect_model_summary.json"
THRESHOLDS = (0.04, 0.05, 0.06, 0.07, 0.08, 0.10)
WIN = 16
RATIO_9KM = 5.50 / 3.47
SOURCES = {
    "native":   ("solar_native", "on_native", "off_native"),
    "recon_ll": ("solar_v5", "on_ll", "off_ll"),
    "recon_qm": ("solar_v5", "on_qm", "off_qm"),
    "wind_ll":  ("solar_native", "on_ll", "off_ll"),
    "wind_qm":  ("solar_native", "on_qm", "off_qm"),
    "solar_v5": ("solar_v5", "on_native", "off_native"),
}
PERIODS = {"control_all": ("ctl1950g", 1953, 2099), "control_1953_1969": ("ctl1950g", 1953, 1969),
           "future_2083_2092": ("hi1950g", 2083, 2092), "future_2081_2099": ("hi1950g", 2081, 2099)}


def load(run):
    d = dict(np.load(f"{DATA}/{run}_germany_3h.npz"))
    t = d["time"].astype("datetime64[s]")
    m = ((t - np.timedelta64(1, "s")).astype("datetime64[M]").astype(int) % 12) + 1
    y = (t - np.timedelta64(1, "s")).astype("datetime64[Y]").astype(int) + 1970
    d["winter"] = np.where(m >= 9, y + 1, y)
    return d


def per_winter_counts(d, src, thr, w0, w1):
    s, on, off = SOURCES[src]
    comb = W_MOCKERT["solar"] * d[s] + W_MOCKERT["on"] * d[on] + W_MOCKERT["off"] * d[off]
    out = {}
    for w in range(w0, w1 + 1):
        m = d["winter"] == w
        if m.sum() < 1600:                       # incomplete Sep-Mar block
            continue
        _, ev = mockert_events(trailing_mean(comb[m], WIN), WIN, thr)
        out[w] = len(ev)
    return out


def mean_bias(d, w0, w1):
    m = (d["winter"] >= w0) & (d["winter"] <= w1)
    day = m & (d["solar_native"] > 0.005)
    r = lambda a, b: float(np.mean(d[a][m]) / np.mean(d[b][m]))   # noqa: E731
    return dict(on_ll_over_native=r("on_ll", "on_native"), on_qm_over_native=r("on_qm", "on_native"),
                on_llfsr_over_native=r("on_llfsr", "on_native"),
                off_ll_over_native=r("off_ll", "off_native"), off_qm_over_native=r("off_qm", "off_native"),
                solar_v5_over_native=r("solar_v5", "solar_native"),
                solar_v5_pre_over_native=float(np.mean(d["solar_v5_pre"][m]) / np.mean(d["solar_native"][m])),
                solar_v5_r=float(np.corrcoef(d["solar_v5"][day], d["solar_native"][day])[0, 1]),
                on_qm_r=float(np.corrcoef(d["on_qm"][m], d["on_native"][m])[0, 1]),
                off_qm_r=float(np.corrcoef(d["off_qm"][m], d["off_native"][m])[0, 1]),
                spd100_on_ll_over_native=float(np.mean(d["spd100_on_ll"][m]) / np.mean(d["spd100_on_native"][m])),
                spd100_off_ll_over_native=float(np.mean(d["spd100_off_ll"][m]) / np.mean(d["spd100_off_native"][m])))


def block_ratio_null(counts, n_a=17, n_b=10, thr_ratio=RATIO_9KM):
    """All placements of a block of n_a consecutive winters and a disjoint later or earlier block of n_b."""
    c = np.array([counts[w] for w in sorted(counts)], float)
    n = len(c)
    ratios = []
    for a in range(n - n_a + 1):
        ma = c[a:a + n_a].mean()
        for b in range(n - n_b + 1):
            if b + n_b <= a or b >= a + n_a:
                ratios.append(c[b:b + n_b].mean() / ma if ma else np.nan)
    ratios = np.array(ratios)
    ratios = ratios[np.isfinite(ratios)]
    return dict(n_pairs=int(len(ratios)), p05_p50_p95=[float(x) for x in np.percentile(ratios, [5, 50, 95])],
                frac_ge_9km_ratio=float(np.mean(ratios >= thr_ratio)), frac_le_inverse=float(np.mean(ratios <= 1 / thr_ratio)))


def main():
    runs = {r: load(r) for r in ("ctl1950g", "hi1950g")}
    out = dict(weights=W_MOCKERT, ratio_9km=RATIO_9KM, bias={}, events={}, ratios={})
    for p, (run, w0, w1) in PERIODS.items():
        out["bias"][p] = mean_bias(runs[run], w0, w1)
        out["events"][p] = {src: {f"{thr:.2f}": float(np.mean(list(per_winter_counts(runs[run], src, thr, w0, w1).values())))
                                  for thr in THRESHOLDS} for src in SOURCES}
    for fut in ("future_2083_2092", "future_2081_2099"):
        for ctl in ("control_all", "control_1953_1969"):
            out["ratios"][f"{fut}/{ctl}"] = {src: {t: out["events"][fut][src][t] / out["events"][ctl][src][t]
                                                    for t in out["events"][ctl][src]} for src in SOURCES}
    # per-winter agreement of reconstructed with native counts at 0.06
    agree = {}
    for p, (run, w0, w1) in PERIODS.items():
        nat = per_winter_counts(runs[run], "native", 0.06, w0, w1)
        agree[p] = {src: float(np.corrcoef([nat[w] for w in nat], [c for c in per_winter_counts(runs[run], src, 0.06, w0, w1).values()])[0, 1])
                    for src in SOURCES if src != "native"}
    out["per_winter_r_vs_native_0p06"] = agree
    ctl_counts = {src: per_winter_counts(runs["ctl1950g"], src, 0.06, 1953, 2099) for src in ("native", "recon_qm")}
    out["natural_variability_17_vs_10_winters_0p06"] = {src: block_ratio_null(c) for src, c in ctl_counts.items()}
    out["control_per_winter_native_0p06"] = ctl_counts["native"]
    # TCO95 has ~7 events/winter at 0.06, the 9 km 1950C 3.47: relative noise depends on the base
    # rate, so repeat the null at thresholds where the control has a 9 km-like base rate
    matched = {}
    for thr in (0.040, 0.045, 0.050, 0.055):
        c = per_winter_counts(runs["ctl1950g"], "native", thr, 1953, 2099)
        fut = per_winter_counts(runs["hi1950g"], "native", thr, 2081, 2099)
        matched[f"{thr:.3f}"] = dict(control_events_per_winter=float(np.mean(list(c.values()))),
                                     future_2081_2099_events_per_winter=float(np.mean(list(fut.values()))),
                                     **block_ratio_null(c))
    out["natural_variability_matched_base_rate"] = matched
    evo = {}
    for run in ("ctl1950g", "hi1950g"):
        c = per_winter_counts(runs[run], "native", 0.06, 1953, 2099)
        evo[run] = {f"{d}s": float(np.mean([c[w] for w in c if d <= w < d + 10])) for d in range(1960, 2100, 10)}
        evo[run]["1953-1959"] = float(np.mean([c[w] for w in c if w < 1960]))
    out["evolution_native_0p06_by_decade"] = evo

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(out, open(OUT, "w"), indent=1)

    print("Mean CF ratio reconstructed / native:")
    for p in PERIODS:
        b = out["bias"][p]
        print(f"  {p:20s} on ll {b['on_ll_over_native']:.3f} qm {b['on_qm_over_native']:.3f} llfsr {b['on_llfsr_over_native']:.3f} | "
              f"off ll {b['off_ll_over_native']:.3f} qm {b['off_qm_over_native']:.3f} | solar v5 {b['solar_v5_over_native']:.4f} "
              f"(pre {b['solar_v5_pre_over_native']:.3f}, r {b['solar_v5_r']:.4f})")
    print("\nEvents per winter at 0.06:")
    for p in PERIODS:
        print(f"  {p:20s} " + "  ".join(f"{s}={out['events'][p][s]['0.06']:.2f}" for s in SOURCES))
    print("\nFuture / control ratio, 2083-2092 vs 1953-1969 (the 9 km spans), by threshold:")
    rr = out["ratios"]["future_2083_2092/control_1953_1969"]
    for s in SOURCES:
        print(f"  {s:9s} " + " ".join(f"{t}:{rr[s][t]:.2f}" for t in rr[s]))
    rr = out["ratios"]["future_2081_2099/control_all"]
    print("Future 2081-2099 / control 1953-2099:")
    for s in SOURCES:
        print(f"  {s:9s} " + " ".join(f"{t}:{rr[s][t]:.2f}" for t in rr[s]))
    print("\nNatural variability (control, 17 vs 10 consecutive winters, 0.06):",
          json.dumps(out["natural_variability_17_vs_10_winters_0p06"]))
    print("Evolution native 0.06:", json.dumps(evo))
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
