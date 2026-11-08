"""
Multi-region test of the 9 km Dunkelflaute increase, and the EERIE IFS-FESOM2-SR comparison.

Inputs: dunkelflaute/data/regions/*.npz from the extract_*.py scripts. Everything on a 6-hourly
footing (3-hourly sources: wind at 00/06/12/18, solar mean of the two 3 h windows), Sep-Mar winters
labelled by the year of their Jan-Mar, 48 h trailing mean = 8 steps, Mockert window expansion.

Metrics per region (land cells only):
  comb   0.44 solar + 0.56 wind (Mockert 2018 solar share; onshore + offshore wind shares pooled)
  wind   wind CF alone (no solar method in the loop)
Wind is the uncalibrated log law in every model (the QM table was fitted over Germany only).

Thresholds: per comparison, region and metric, the threshold at which the comparison's own control
has TARGET events per winter (base-rate matching; removes mean-CF differences between models and
regions). Also fixed 0.06 for comb.

Tests:
  per region   Poisson exact test; block null from the TCO95 constant-1950 control (147 winters) and
               from the EERIE control (79 winters): all placements of disjoint 17- and 10-winter blocks
  field        statistic = mean over regions of log(ratio); for 9 km against 31 km, the mean of
               log(r9 / r31); null by drawing block placements from a control, the SAME winters in
               every region, so the spatial correlation of the real samples is kept

Usage: python3 region_events.py
Output: dunkelflaute/results/regions/region_events_summary.json (+ per-winter counts .npz)
"""
import json
import os
import sys

import numpy as np
from scipy import stats

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/analysis")
import regions as rg                                                     # noqa: E402
from robustness_mockert_events import mockert_events, trailing_mean     # noqa: E402

D = rg.OUT_DIR
OUT = f"{_ROOT}/dunkelflaute/results/regions"
WIN = 8
TARGET = 5.0
W_SOLAR = 0.44
H3 = np.timedelta64(3, "h")
MIN_STEPS = 800
RNG = np.random.default_rng(1)

COMPARISONS = {
    "9km":            (("dart9km_1950c", 1953, 1969), ("dart9km_2080c", 2083, 2092)),
    "31km_9kmspans":  (("tco319_ctl1950d", 1953, 1969), ("tco319_ssp585d", 2083, 2092)),
    "31km":           (("tco319_ctl1950d", 1953, 1969), ("tco319_ssp585d", 2081, 2099)),
    "100km_9kmspans": (("tco95_ctl1950g", 1953, 1969), ("tco95_hi1950g", 2083, 2092)),
    "100km":          (("tco95_ctl1950g", 1953, 2099), ("tco95_hi1950g", 2081, 2099)),
    "eerie_ctl":      (("eerie_ctl", 1972, 2050), ("eerie_ssp245", 2041, 2070)),
    "eerie_hist":     (("eerie_hist", 1971, 2000), ("eerie_ssp245", 2041, 2070)),
    "eerie_drift":    (("eerie_ctl", 1972, 2000), ("eerie_ctl", 2021, 2050)),
}
NULLS = {"tco95_ctl": ("tco95_ctl1950g", 1953, 2099), "eerie_ctl": ("eerie_ctl", 1972, 2050)}


def winter_of(t):
    t1 = t - np.timedelta64(1, "s")
    m = (t1.astype("datetime64[M]").astype(int) % 12) + 1
    y = t1.astype("datetime64[Y]").astype(int) + 1970
    return np.where(m >= 9, y + 1, y)


def to_6h(d):
    t = d["time"]
    hour = (t - t.astype("datetime64[D]")) // np.timedelta64(1, "h")
    pos = {v: i for i, v in enumerate(t)}
    end = np.array([i for i in np.where(hour % 6 == 0)[0] if t[i] - H3 in pos])
    prev = np.array([pos[t[i] - H3] for i in end])
    out = {k: (v[end] if k != "solar" else 0.5 * (v[end] + v[prev])) for k, v in d.items() if k != "time"}
    out["time"] = t[end]
    return out


_CACHE = {}


def load(name):
    if name in _CACHE:
        return _CACHE[name]
    z = np.load(f"{D}/{name}.npz")
    regions = list(z["regions"])
    assert regions == rg.REGIONS, f"{name}: region order differs"
    d = {k: z[k] for k in ("time", "solar", "wind", "wind_off_DE") if k in z.files}
    for k in ("wind_qm", "wind_off_DE_qm", "t2m"):
        if k in z.files:
            d[k] = z[k]
    d["time"] = d["time"].astype("datetime64[s]")
    dt_h = json.loads(str(z["meta"]))["dt_h"]
    if dt_h == 3:
        d = to_6h(d)
    d["winter"] = winter_of(d["time"])
    d["n_cells"] = z["n_cells"]
    _CACHE[name] = d
    return d


def metric(d, kind):
    if kind == "comb":
        return W_SOLAR * d["solar"] + (1 - W_SOLAR) * d["wind"]
    if kind == "wind":
        return d["wind"]
    raise ValueError(kind)


def winters(d, y0, y1):
    out = []
    for y in range(y0, y1 + 1):
        m = d["winter"] == y
        if m.sum() >= MIN_STEPS:
            out.append((y, m))
    return out


def counts_grid(x, thr):
    """Events per winter-series x (T,) for each threshold in thr (K,), Mockert window expansion,
    vectorised: flagged[i] = any(roll[i .. i+WIN-1] < thr)."""
    roll = trailing_mean(x, WIN)
    below = (np.nan_to_num(roll, nan=np.inf)[None, :] < thr[:, None]).astype(np.int32)       # (K, T)
    c = np.concatenate([np.zeros((len(thr), 1), np.int32), np.cumsum(below, axis=1)], axis=1)
    n = below.shape[1]
    hi = np.minimum(np.arange(n) + WIN, n)
    flagged = (c[:, hi] - c[:, :n]) > 0
    starts = flagged[:, 0].astype(int) + (flagged[:, 1:] & ~flagged[:, :-1]).sum(axis=1)
    return starts


def per_winter(d, kind, region, y0, y1, thr):
    x = metric(d, kind)[:, region]
    ws = winters(d, y0, y1)
    return np.array([y for y, _ in ws]), np.stack([counts_grid(x[m], thr) for _, m in ws])     # (W,), (W, K)


def matched_threshold(d, kind, region, y0, y1):
    """Lowest threshold on a fine grid at which the control mean reaches TARGET events per winter."""
    x = metric(d, kind)[:, region]
    qs = np.nanpercentile(x, [0.1, 40])
    grid = np.linspace(qs[0], qs[1], 400)
    _, c = per_winter(d, kind, region, y0, y1, grid)
    mean = c.mean(axis=0)
    k = int(np.argmax(mean >= TARGET)) if (mean >= TARGET).any() else int(np.argmax(mean))
    return float(grid[k])


def change(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    tot = a.sum() + b.sum()
    p = stats.binomtest(int(b.sum()), int(tot), len(b) / (len(a) + len(b))).pvalue if tot else np.nan
    return dict(control=float(a.mean()), future=float(b.mean()),
                ratio=float(b.mean() / a.mean()) if a.mean() else np.nan, poisson_p=float(p),
                n_control=len(a), n_future=len(b))


def block_pairs(n, n_a, n_b):
    pairs = []
    for a in range(n - n_a + 1):
        for b in range(n - n_b + 1):
            if b + n_b <= a or b >= a + n_a:
                pairs.append((a, b))
    return np.array(pairs)


def block_logratios(c, n_a, n_b, pairs):
    """c (W, R) per-winter counts -> (n_pairs, R) log ratios of block means."""
    cs = np.concatenate([np.zeros((1, c.shape[1])), np.cumsum(c, axis=0)])
    ma = (cs[pairs[:, 0] + n_a] - cs[pairs[:, 0]]) / n_a
    mb = (cs[pairs[:, 1] + n_b] - cs[pairs[:, 1]]) / n_b
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(mb / ma)


def main():
    os.makedirs(OUT, exist_ok=True)
    R = rg.REGIONS
    out = dict(regions=R, target_events_per_winter=TARGET, window_steps=WIN, solar_weight=W_SOLAR,
               comparisons={}, nulls={}, field={}, germany_mockert={})
    available = {n for c in COMPARISONS.values() for n, _, _ in c if os.path.exists(f"{D}/{n}.npz")}
    store = {}
    for name, ((cn, c0, c1), (fn, f0, f1)) in COMPARISONS.items():
        if cn not in available or fn not in available:
            print(f"skip {name}: missing input")
            continue
        dc, dfu = load(cn), load(fn)
        res = {}
        for kind in ("comb", "wind"):
            rows = {}
            for i, r in enumerate(R):
                thr = matched_threshold(dc, kind, i, c0, c1)
                _, a = per_winter(dc, kind, i, c0, c1, np.array([thr]))
                _, b = per_winter(dfu, kind, i, f0, f1, np.array([thr]))
                rows[r] = dict(threshold=thr, **change(a[:, 0], b[:, 0]))
                store[(name, kind, r)] = (a[:, 0], b[:, 0])
            res[kind] = rows
        rows = {}
        for i, r in enumerate(R):
            _, a = per_winter(dc, "comb", i, c0, c1, np.array([0.06]))
            _, b = per_winter(dfu, "comb", i, f0, f1, np.array([0.06]))
            rows[r] = change(a[:, 0], b[:, 0])
        res["comb_fixed006"] = rows
        out["comparisons"][name] = res

    # nulls from long controls, base-rate matched on the whole control
    null_logr = {}
    for nn, (src, y0, y1) in NULLS.items():
        if src not in available:
            continue
        d = load(src)
        for kind in ("comb", "wind"):
            cols = []
            for i in range(len(R)):
                thr = matched_threshold(d, kind, i, y0, y1)
                _, c = per_winter(d, kind, i, y0, y1, np.array([thr]))
                cols.append(c[:, 0])
            c = np.stack(cols, axis=1).astype(float)                    # (W, R)
            pairs = block_pairs(len(c), 17, 10)
            lr = block_logratios(c, 17, 10, pairs)                      # (P, R)
            null_logr[(nn, kind)] = (c, pairs, lr)
            out["nulls"][f"{nn}_{kind}"] = dict(n_winters=len(c), n_pairs=len(pairs),
                                                ratio_p05_p50_p95={r: [float(x) for x in np.exp(np.nanpercentile(lr[:, j], [5, 50, 95]))]
                                                                   for j, r in enumerate(R)})

    # a region enters the field tests for a metric only if every control reaches 90% of TARGET
    # (in sunny, low-wind regions the combined CF never gives TARGET separate events: lows merge)
    valid = {}
    for kind in ("comb", "wind"):
        v = np.ones(len(R), bool)
        for name in out["comparisons"]:
            v &= np.array([out["comparisons"][name][kind][r]["control"] >= 0.9 * TARGET for r in R])
        for nn in NULLS:
            if (nn, kind) in null_logr:
                v &= null_logr[(nn, kind)][0].mean(axis=0) >= 0.9 * TARGET
        valid[kind] = v
    out["valid_regions"] = {k: [r for r, x in zip(R, v) if x] for k, v in valid.items()}

    # per-region block-null p for the 9 km and the 9 km-span comparisons, and field tests
    for kind in ("comb", "wind"):
        fld = {}
        for nn in NULLS:
            if (nn, kind) not in null_logr:
                continue
            c, pairs, lr = null_logr[(nn, kind)]
            for name in ("9km", "31km_9kmspans", "100km_9kmspans", "eerie_ctl"):
                if name not in out["comparisons"]:
                    continue
                obs = np.array([np.log(out["comparisons"][name][kind][r]["ratio"]) for r in R])
                for j, r in enumerate(R):
                    out["comparisons"][name][kind][r][f"block_null_p_{nn}"] = float(np.nanmean(lr[:, j] >= obs[j]))
                ok = np.isfinite(obs) & valid[kind]
                stat = float(np.mean(obs[ok]))
                nullm = np.nanmean(lr[:, ok], axis=1)
                fld[f"{name}_mean_logratio_vs_{nn}"] = dict(observed=stat, observed_ratio=float(np.exp(stat)),
                                                            n_regions=int(ok.sum()),
                                                            p_one_sided=float(np.mean(nullm >= stat)),
                                                            null_p05_p95=[float(x) for x in np.exp(np.percentile(nullm, [5, 95]))],
                                                            n_regions_ratio_gt1=int((obs[ok] > 0).sum()))
            # 9 km against 31 km at the same block lengths: difference of two independent draws
            if "9km" in out["comparisons"] and "31km_9kmspans" in out["comparisons"]:
                o9 = np.array([np.log(out["comparisons"]["9km"][kind][r]["ratio"]) for r in R])
                o31 = np.array([np.log(out["comparisons"]["31km_9kmspans"][kind][r]["ratio"]) for r in R])
                ok = np.isfinite(o9) & np.isfinite(o31) & valid[kind]
                diff = o9[ok] - o31[ok]
                i1 = RNG.integers(0, len(pairs), 20000)
                i2 = RNG.integers(0, len(pairs), 20000)
                nd = np.nanmean(lr[i1][:, ok] - lr[i2][:, ok], axis=1)
                nsign = ((lr[i1][:, ok] - lr[i2][:, ok]) > 0).sum(axis=1)
                fld[f"9km_minus_31km_vs_{nn}"] = dict(observed_mean_logdiff=float(diff.mean()),
                                                      observed_ratio_of_ratios=float(np.exp(diff.mean())),
                                                      n_regions=int(ok.sum()), n_regions_9km_higher=int((diff > 0).sum()),
                                                      p_one_sided=float(np.mean(nd >= diff.mean())),
                                                      p_sign_count=float(np.mean(nsign >= (diff > 0).sum())),
                                                      pattern_corr_9km_31km=float(np.corrcoef(o9[ok], o31[ok])[0, 1]))
            if "9km" in out["comparisons"] and "eerie_ctl" in out["comparisons"]:
                o9 = np.array([np.log(out["comparisons"]["9km"][kind][r]["ratio"]) for r in R])
                oe = np.array([np.log(out["comparisons"]["eerie_ctl"][kind][r]["ratio"]) for r in R])
                ok = np.isfinite(o9) & np.isfinite(oe) & valid[kind]
                fld["pattern_corr_9km_eerie"] = float(np.corrcoef(o9[ok], oe[ok])[0, 1])
        out["field"][kind] = fld

    # Germany with offshore, Mockert 2018 weights, fixed 0.06 (the headline metric), where available
    for name, ((cn, c0, c1), (fn, f0, f1)) in COMPARISONS.items():
        if name not in out["comparisons"]:
            continue
        dc, dfu = load(cn), load(fn)
        row = {}
        for tag, wk, ok_ in (("ll", "wind", "wind_off_DE"), ("qm", "wind_qm", "wind_off_DE_qm")):
            if wk not in dc or ok_ not in dc:
                continue
            ab = []
            for d, y0, y1 in ((dc, c0, c1), (dfu, f0, f1)):
                x = 0.44 * d["solar"][:, 0] + 0.50 * d[wk][:, 0] + 0.06 * d[ok_]
                ab.append([len(mockert_events(trailing_mean(x[m], WIN), WIN, 0.06)[1]) for _, m in winters(d, y0, y1)])
            row[tag] = change(*ab)
        out["germany_mockert"][name] = row

    json.dump(out, open(f"{OUT}/region_events_summary.json", "w"), indent=1)
    np.savez(f"{OUT}/region_events_counts.npz",
             **{f"{n}|{k}|{r}|{p}": v[j] for (n, k, r), v in store.items() for j, p in enumerate(("ctl", "fut"))})

    # print
    for kind in ("comb", "wind"):
        print(f"\n=== {kind}: future/control ratio of events per winter (base rate {TARGET:.0f}/winter in control); * p<0.05 Poisson")
        names = [n for n in COMPARISONS if n in out["comparisons"]]
        print(f"{'region':7s}" + "".join(f"{n:>15s}" for n in names))
        for r in R:
            cells = []
            for n in names:
                x = out["comparisons"][n][kind][r]
                cells.append(f"{x['ratio']:.2f}{'*' if x['poisson_p'] < 0.05 else ' '}")
            print(f"{r:7s}" + "".join(f"{c:>15s}" for c in cells))
        print(f"  valid regions ({len(out['valid_regions'][kind])}): {' '.join(out['valid_regions'][kind])}")
        for k, v in out["field"].get(kind, {}).items():
            print(f"  {k}: {v}")
    print("\n=== Germany, Mockert mix with offshore, fixed 0.06")
    for n, row in out["germany_mockert"].items():
        print(f"  {n:15s} " + "  ".join(f"{t}: {x['control']:.2f}->{x['future']:.2f} ({x['ratio']:.2f}, p={x['poisson_p']:.3f})" for t, x in row.items()))
    print(f"Saved {OUT}/region_events_summary.json")


if __name__ == "__main__":
    main()
