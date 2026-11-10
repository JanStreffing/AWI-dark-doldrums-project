"""
Does the per-K change in low-wind events and Dunkelflauten depend on the base rate (TARGET events per
control winter)? Recomputes, for TARGET = 2, 5 and 10:
  (a) our runs on the 6-hourly footing (region_events.py + region_per_K.py logic): the four main
      comparisons (9km, 31km, 100km, eerie_ctl), wind and comb, field mean over the valid regions;
      validity = every control of every comparison in region_events.COMPARISONS and both null controls
      reaches 90 % of TARGET (as in region_events.py); through-origin fit of ln ratio on dT over the
      four main comparisons (as in paper/scripts/perK/plot_perK.py);
  (b) CMIP6 on the daily footing (daily_events.py logic): 20 models historical 1951-1980 against
      SSP5-8.5 2071-2100, common valid region set over every non-wind-only comparison of daily_events.py
      (CMIP6 end and mid century, EERIE and TCO95 bridges); EC-Earth3 wind only, outside the ensemble.
Everything else as in the source scripts (regions, winters, 48 h window, Mockert expansion, threshold
fixed between control and future). The threshold grid is that of the source scripts (400 points from the
0.1th to the 40th percentile of the full series); for TARGET != 5 it is extended with the same spacing
down to the series minimum and up to the 80th percentile, but only where the original grid does not
bracket TARGET (control rate at the lowest grid point already >= TARGET, or never reaching TARGET).
TARGET = 5 uses the original grid only and is checked against the stored results. Side check: our runs
at TARGET = 5 again with the extended grid (key ours_target5_extended_grid), because on the 6-hourly
footing the original grid's top (40th percentile of the full series, nights included) does not reach the
48 h-mean comb levels of the sunny regions, which is what made them invalid.

Does not import scripts/regimes/geostrophic.py (it sets region_events.WIN = 2).
Usage: python3 base_rate_sensitivity.py
Output: dunkelflaute/results/regions/base_rate_sensitivity.json, base_rate_sensitivity.png
"""
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/cmip6")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/analysis")
import regions as rg                                    # noqa: E402
import region_events as ev                              # noqa: E402
import region_per_K as rpk                              # noqa: E402

assert "geostrophic" not in sys.modules
assert ev.WIN == 8, ev.WIN

RES = f"{_ROOT}/dunkelflaute/results"
OUT_JSON = f"{RES}/regions/base_rate_sensitivity.json"
OUT_PNG = f"{RES}/regions/base_rate_sensitivity.png"
TARGETS = (2.0, 5.0, 10.0)
MAIN = ["9km", "31km", "100km", "eerie_ctl"]
KINDS = ("wind", "comb")
R = rg.REGIONS
pk = lambda s: float(100 * (np.exp(s) - 1))            # noqa: E731


def base_grid(x):
    return np.linspace(*np.nanpercentile(x, [0.1, 40]), 400)


def ext_grid(x, g):
    step = g[1] - g[0]
    lo, hi = np.nanmin(x), np.nanpercentile(x, 80)
    nb = int(np.ceil((g[0] - lo) / step)) if g[0] > lo else 0
    na = int(np.ceil((hi - g[-1]) / step)) if hi > g[-1] else 0
    return np.concatenate([g[0] - step * np.arange(nb, 0, -1), g, g[-1] + step * np.arange(1, na + 1)])


def pick(grid, mean, target):
    k = int(np.argmax(mean >= target)) if (mean >= target).any() else int(np.argmax(mean))
    return float(grid[k]), float(mean[k])


EXTEND_AT_5 = False     # True only for the side check of the 6-hourly comb validity at TARGET = 5


def matched(x, count_fn, target):
    """count_fn(grid) -> (W, K). Returns threshold, control rate there, grid tag."""
    g = base_grid(x)
    mean = count_fn(g).mean(axis=0)
    brackets = mean[0] < target and (mean >= target).any()
    if (target == 5.0 and not EXTEND_AT_5) or brackets:
        return (*pick(g, mean, target), "base" if brackets else "base_unbracketed")
    g2 = ext_grid(x, g)
    mean2 = count_fn(g2).mean(axis=0)
    tag = "extended" if (mean2[0] < target and (mean2 >= target).any()) else "extended_unbracketed"
    return (*pick(g2, mean2, target), tag)


# ------------------------------------------------------------------ (a) our runs, 6-hourly footing
def ours(target):
    assert ev.WIN == 8, ev.WIN
    T = rpk.global_t()
    comps, tags = {}, {}
    for name, ((cn, c0, c1), (fn, f0, f1)) in ev.COMPARISONS.items():
        dc, dfu = ev.load(cn), ev.load(fn)
        comps[name] = {}
        for kind in KINDS:
            rows = {}
            for i, r in enumerate(R):
                x = ev.metric(dc, kind)[:, i]
                thr, _, tag = matched(x, lambda g: ev.per_winter(dc, kind, i, c0, c1, g)[1], target)
                tags[tag] = tags.get(tag, 0) + 1
                _, a = ev.per_winter(dc, kind, i, c0, c1, np.array([thr]))
                _, b = ev.per_winter(dfu, kind, i, f0, f1, np.array([thr]))
                a, b = a[:, 0].astype(float), b[:, 0].astype(float)
                rows[r] = dict(threshold=thr, control=float(a.mean()), future=float(b.mean()),
                               ratio=float(b.mean() / a.mean()) if a.mean() else np.nan)
            comps[name][kind] = rows
    null_rate = {}
    for nn, (src, y0, y1) in ev.NULLS.items():
        d = ev.load(src)
        for kind in KINDS:
            v = []
            for i in range(len(R)):
                x = ev.metric(d, kind)[:, i]
                thr, _, _ = matched(x, lambda g: ev.per_winter(d, kind, i, y0, y1, g)[1], target)
                v.append(ev.per_winter(d, kind, i, y0, y1, np.array([thr]))[1][:, 0].mean())
            null_rate[(nn, kind)] = np.array(v)
    valid = {}
    for kind in KINDS:
        v = np.ones(len(R), bool)
        for name in comps:
            v &= np.array([comps[name][kind][r]["control"] >= 0.9 * target for r in R])
        for nn in ev.NULLS:
            v &= null_rate[(nn, kind)] >= 0.9 * target
        valid[kind] = [r for r, ok in zip(R, v) if ok]
    field = {}
    for name, ((cn, c0, c1), (fn, f0, f1)) in ev.COMPARISONS.items():
        dT = rpk.mean_years(T, fn, f0, f1)[0] - rpk.mean_years(T, cn, c0, c1)[0]
        if not np.isfinite(dT) or abs(dT) < 0.3:
            continue
        field[name] = dict(dT=float(dT))
        for kind in KINDS:
            lr = np.array([np.log(comps[name][kind][r]["ratio"]) for r in valid[kind]])
            m = float(np.mean(lr))
            field[name][kind] = dict(n_regions=len(lr), mean_ratio=float(np.exp(m)), pct_per_K=pk(m / dT),
                                     n_regions_up=int((lr > 0).sum()))
    fit = {}
    for kind in KINDS:
        x = np.array([field[n]["dT"] for n in MAIN])
        y = np.log([field[n][kind]["mean_ratio"] for n in MAIN])
        fit[kind] = pk((x * y).sum() / (x * x).sum())
    return dict(comparisons=comps, valid_regions=valid, n_valid={k: len(v) for k, v in valid.items()},
                field_per_K=field, origin_fit_main_pct_per_K=fit, grid_use=tags)


def ours_fixed_set(res, regions, target_res):
    """per-K for the main comparisons on a fixed region set (like-for-like across TARGETs)."""
    out = {}
    for kind in KINDS:
        V = regions[kind]
        out[kind] = {}
        for n in MAIN:
            dT = target_res["field_per_K"][n]["dT"]
            m = np.mean([np.log(target_res["comparisons"][n][kind][r]["ratio"]) for r in V])
            out[kind][n] = pk(m / dT)
        x = np.array([target_res["field_per_K"][n]["dT"] for n in MAIN])
        y = np.array([np.mean([np.log(target_res["comparisons"][n][kind][r]["ratio"]) for r in V]) for n in MAIN])
        out[kind]["origin_fit"] = pk((x * y).sum() / (x * x).sum())
        out[kind]["n_regions"] = len(V)
    return out


# ------------------------------------------------------------------ (b) CMIP6, daily footing
def cmip6(target, de):
    assert ev.WIN == 2, ev.WIN
    C = de.comparisons()
    res, tags = {}, {}
    for name, (fc, ec, c0, c1, ff, ef, f0, f1, Tc, Tf, group) in C.items():
        dc, df = de.load(fc, ec), de.load(ff, ef)
        dT = de.mean_T(Tf, f0, f1) - de.mean_T(Tc, c0, c1)
        row = dict(group=group, dT=float(dT), metrics={})
        for kind in de.kinds_of(group):
            per = {}
            for i, r in enumerate(R):
                x = de.metric(dc, kind)[:, i]
                if not np.isfinite(x).any():
                    per[r] = dict(valid=False)
                    continue
                thr, rate, tag = matched(x, lambda g: de.counts(dc, kind, i, c0, c1, g), target)
                tags[tag] = tags.get(tag, 0) + 1
                a = de.counts(dc, kind, i, c0, c1, np.array([thr]))[:, 0]
                b = de.counts(df, kind, i, f0, f1, np.array([thr]))[:, 0]
                per[r] = dict(valid=bool(rate >= 0.9 * target), threshold=thr, control=float(a.mean()),
                              future=float(b.mean()), ratio=float(b.mean() / a.mean()) if a.mean() else np.nan)
            row["metrics"][kind] = per
        res[name] = row
    common = {kind: [r for r in R if all(res[n]["metrics"][kind][r].get("valid") for n in res
                                         if not res[n]["group"].startswith("cmip6_windonly"))] for kind in KINDS}
    summ = {}
    for name, row in res.items():
        summ[name] = dict(group=row["group"], dT=row["dT"])
        for kind in de.kinds_of(row["group"]):
            lr = np.array([np.log(row["metrics"][kind][r]["ratio"]) for r in common[kind]])
            ok = np.isfinite(lr)
            m = float(lr[ok].mean())
            summ[name][kind] = dict(n_regions=int(ok.sum()), mean_ratio=float(np.exp(m)), pct_per_K=pk(m / row["dT"]),
                                    n_regions_up=int((lr[ok] > 0).sum()))
    names = [n for n in summ if summ[n]["group"] == "cmip6"]
    ens = {}
    for kind in KINDS:
        lnr = np.array([np.log(summ[n][kind]["mean_ratio"]) for n in names])
        dT = np.array([summ[n]["dT"] for n in names])
        p = np.array([summ[n][kind]["pct_per_K"] for n in names])
        ens[kind] = dict(n_models=len(names), pct_per_K_median=float(np.median(p)),
                         pct_per_K_p10_p90=[float(v) for v in np.percentile(p, [10, 90])],
                         n_models_up=int((lnr > 0).sum()), origin_fit_pct_per_K=pk((dT * lnr).sum() / (dT * dT).sum()),
                         per_model={n.split("|")[0]: float(v) for n, v in zip(names, p)})
    ece = {n.split("|")[0]: summ[n]["wind"]["pct_per_K"] for n in summ if summ[n]["group"] == "cmip6_windonly"}
    bridge = {n: {k: summ[n][k]["pct_per_K"] for k in KINDS} for n in summ if summ[n]["group"].startswith("bridge")}
    return dict(common_regions=common, n_common={k: len(v) for k, v in common.items()}, ensemble=ens,
                ecearth3_wind_pct_per_K=ece, bridges_pct_per_K=bridge, grid_use=tags), summ


# ------------------------------------------------------------------ reproduction check
def check_ours(o):
    S = json.load(open(f"{RES}/regions/region_events_summary.json"))
    P = json.load(open(f"{RES}/regions/region_per_K_summary.json"))
    md = 0.0
    for name in S["comparisons"]:
        for kind in KINDS:
            for r in R:
                s, x = S["comparisons"][name][kind][r], o["comparisons"][name][kind][r]
                for k in ("threshold", "control", "future", "ratio"):
                    if np.isfinite(s[k]) or np.isfinite(x[k]):
                        md = max(md, abs(s[k] - x[k]))
    same_valid = all(S["valid_regions"][k] == o["valid_regions"][k] for k in KINDS)
    dpk = max(abs(P["field_per_K"][n][k]["pct_per_K"] - o["field_per_K"][n][k]["pct_per_K"])
              for n in P["field_per_K"] for k in KINDS)
    return dict(max_abs_diff_region_values=md, same_valid_regions=same_valid, max_abs_diff_field_pct_per_K=dpk)


def check_cmip6(c, summ):
    S = json.load(open(f"{RES}/cmip6/daily_events_summary.json"))
    dpk = max(abs(S["comparisons"][n][k]["pct_per_K"] - summ[n][k]["pct_per_K"])
              for n in S["comparisons"] for k in KINDS if k in S["comparisons"][n])
    e = S["ensemble"]["cmip6"]
    dens = max(max(abs(e[k]["pct_per_K_median"] - c["ensemble"][k]["pct_per_K_median"]),
                   *[abs(a - b) for a, b in zip(e[k]["pct_per_K_p10_p90"], c["ensemble"][k]["pct_per_K_p10_p90"])])
               for k in KINDS)
    return dict(max_abs_diff_comparison_pct_per_K=dpk, max_abs_diff_ensemble=dens,
                same_common_regions=all(S["common_regions"][k] == c["common_regions"][k] for k in KINDS),
                same_n_models_up=all(e[k]["n_models_up"] == c["ensemble"][k]["n_models_up"] for k in KINDS))


# ------------------------------------------------------------------ figure
def figure(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {"9km": "#2a78d6", "31km": "#eb6834", "100km": "#1baf7a", "eerie_ctl": "#e87ba4"}
    lab = {"9km": "AWI-CM3 9 km", "31km": "AWI-CM3 31 km", "100km": "AWI-CM3 100 km", "eerie_ctl": "IFS-FESOM 9 km"}
    ts = [float(t) for t in TARGETS]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3), sharey=True, constrained_layout=True)
    for ax, kind in zip(axes, KINDS):
        A = [out["targets"][f"{t:g}"] for t in ts]
        med = [a["cmip6"]["ensemble"][kind]["pct_per_K_median"] for a in A]
        p10 = [a["cmip6"]["ensemble"][kind]["pct_per_K_p10_p90"][0] for a in A]
        p90 = [a["cmip6"]["ensemble"][kind]["pct_per_K_p10_p90"][1] for a in A]
        for m in A[0]["cmip6"]["ensemble"][kind]["per_model"]:
            ax.plot(ts, [a["cmip6"]["ensemble"][kind]["per_model"][m] for a in A], color="0.82", lw=0.6, zorder=1)
        ax.fill_between(ts, p10, p90, color="0.55", alpha=0.18, lw=0, zorder=0, label="CMIP6 10-90 %")
        ax.plot(ts, med, color="0.25", lw=2, marker="s", ms=5, zorder=3, label="CMIP6 median (daily)")
        for n in MAIN:
            ax.plot(ts, [a["ours"]["field_per_K"][n][kind]["pct_per_K"] for a in A], color=col[n], lw=1.4,
                    marker="o", ms=5, zorder=4, label=lab[n])
        ax.plot(ts, [a["ours"]["origin_fit_main_pct_per_K"][kind] for a in A], color="black", lw=2, ls="--",
                zorder=5, label="Fit, our runs (6-hourly)")
        e5 = out["ours_target5_extended_grid"]
        ax.scatter([5.0], [e5["origin_fit_main_pct_per_K"][kind]], s=40, facecolor="white", edgecolor="black", lw=1.2,
                   zorder=6, label="Fit, our runs, TARGET 5 on the extended grid")
        nreg_o = "/".join(str(a["ours"]["n_valid"][kind]) for a in A) + f" (5, extended grid: {e5['n_valid'][kind]})"
        nreg_c = "/".join(str(a["cmip6"]["n_common"][kind]) for a in A)
        ax.text(0.03, 0.97, f"regions ours {nreg_o}\nregions CMIP6 {nreg_c}", transform=ax.transAxes,
                va="top", fontsize=7, color="0.3")
        ax.set_xscale("log")
        ax.set_xticks(ts)
        ax.set_xticklabels([f"{t:g}" for t in ts])
        ax.minorticks_off()
        ax.set_xlabel("Base rate: events per control winter")
        ax.set_title("Wind" if kind == "wind" else "Wind + solar (0.56/0.44)", fontsize=9)
        ax.axhline(0, color="black", lw=0.6)
        ax.grid(axis="y", color="0.9", lw=0.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel("Change in events (% per K)")
    h, l_ = axes[0].get_legend_handles_labels()
    fig.legend(h, l_, loc="outside lower center", ncol=4, fontsize=7, frameon=False)
    fig.savefig(OUT_PNG, dpi=200)


def main():
    out = dict(targets={}, note=__doc__.strip().splitlines()[0:3])
    summs = {}
    for t in TARGETS:
        print(f"ours TARGET={t:g}", flush=True)
        out["targets"][f"{t:g}"] = dict(ours=ours(t))
    out["check_target5_ours"] = check_ours(out["targets"]["5"]["ours"])
    # side check: at TARGET = 5 on the 6-hourly footing, the comb regions TR, US_CA, CN_N, CN_E, KR, IN_N fail
    # validity on the original grid because its top (40th percentile of the full 6-hourly series, nights
    # included) lies below the 48 h-mean levels there; recount with the grid extended as for TARGET 2 and 10
    global EXTEND_AT_5
    EXTEND_AT_5 = True
    print("ours TARGET=5, extended grid", flush=True)
    o5 = ours(5.0)
    EXTEND_AT_5 = False
    out["ours_target5_extended_grid"] = dict(n_valid=o5["n_valid"], valid_regions=o5["valid_regions"], grid_use=o5["grid_use"],
                                             field_per_K=o5["field_per_K"], origin_fit_main_pct_per_K=o5["origin_fit_main_pct_per_K"])
    import daily_events as de                       # sets ev.WIN = 2 on import
    assert ev.WIN == 2, ev.WIN
    for t in TARGETS:
        print(f"cmip6 TARGET={t:g}", flush=True)
        c, summs[t] = cmip6(t, de)
        out["targets"][f"{t:g}"]["cmip6"] = c
    out["check_target5_cmip6"] = check_cmip6(out["targets"]["5"]["cmip6"], summs[5.0])
    # like-for-like: our runs on the regions valid at every TARGET
    inter = {k: [r for r in R if all(r in out["targets"][f"{t:g}"]["ours"]["valid_regions"][k] for t in TARGETS)] for k in KINDS}
    out["ours_on_regions_valid_at_all_targets"] = {f"{t:g}": ours_fixed_set(None, inter, out["targets"][f"{t:g}"]["ours"]) for t in TARGETS}
    # drop the bulky per-region values of our runs, keep thresholds/rates for the main comparisons
    for t in TARGETS:
        o = out["targets"][f"{t:g}"]["ours"]
        o["comparisons"] = {n: o["comparisons"][n] for n in MAIN}
    json.dump(out, open(OUT_JSON, "w"), indent=1)
    figure(out)

    print(json.dumps({k: out[k] for k in ("check_target5_ours", "check_target5_cmip6")}, indent=1))
    for t in TARGETS:
        a = out["targets"][f"{t:g}"]
        for kind in KINDS:
            o, e = a["ours"], a["cmip6"]["ensemble"][kind]
            print(f"T={t:<3g} {kind:4s} ours " + " ".join(f"{n} {o['field_per_K'][n][kind]['pct_per_K']:+.2f}" for n in MAIN)
                  + f" | fit {o['origin_fit_main_pct_per_K'][kind]:+.2f} ({o['n_valid'][kind]} reg)"
                  + f" | CMIP6 med {e['pct_per_K_median']:+.2f} [{e['pct_per_K_p10_p90'][0]:+.2f},{e['pct_per_K_p10_p90'][1]:+.2f}]"
                  + f" up {e['n_models_up']}/{e['n_models']} fit {e['origin_fit_pct_per_K']:+.2f} ({a['cmip6']['n_common'][kind]} reg)")
        print(f"   EC-Earth3 wind {a['cmip6']['ecearth3_wind_pct_per_K']}; grid ours {a['ours']['grid_use']} cmip6 {a['cmip6']['grid_use']}")
        print(f"   bridges {a['cmip6']['bridges_pct_per_K']}")
    for t, v in out["ours_on_regions_valid_at_all_targets"].items():
        print(f"fixed region set T={t}: " + "; ".join(f"{k} ({v[k]['n_regions']}) fit {v[k]['origin_fit']:+.2f} "
                                                      + " ".join(f"{n} {v[k][n]:+.2f}" for n in MAIN) for k in KINDS))
    e5 = out["ours_target5_extended_grid"]
    print(f"TARGET 5 extended grid: n_valid {e5['n_valid']}, grid {e5['grid_use']}, fit {e5['origin_fit_main_pct_per_K']}, "
          + "; ".join(f"{k} " + " ".join(f"{n} {e5['field_per_K'][n][k]['pct_per_K']:+.2f}" for n in MAIN) for k in KINDS))
    print(f"Saved {OUT_JSON}, {OUT_PNG}")


if __name__ == "__main__":
    main()
