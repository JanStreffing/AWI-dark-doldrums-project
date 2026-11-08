"""
Question C (methods notes 13d): what sets the at-fixed-G transfer change across CMIP6 models, and do
our runs sit on the CMIP6 relation?

Inputs
  monthly diagnostics   data/regimes/stability_cmip6_parts/cmip6_<model>.npz (stability_cmip6_extract.py)
                        and tco95_<run>.npz (stability_cmip6_tco95.py), optionally eerie_<run>.npz
  transfer terms        results/regimes/transfer_stability_cmip6.json (8 EU regions, at-fixed-G),
                        stability_cmip6_split22.json (same split, 22 regions), transfer_checks.json
                        (total transfer change, 22 regions, bilinear G), transfer_stability.json (TCO95)
Winter = Sep-Mar, equal weight per month, labelled by the year of its January. CMIP6 historical
1951-1980 against ssp585 2071-2100; TCO95 the 100km comparison (ctl1950g 1953-2099 against hi1950g
2081-2099) and hist_fut (hi1950g 1951-1980 against 2071-2099).
Region aggregates: plain mean over the 8 European regions (regime_analysis.EU) or the 22
geostrophic.CMIP_REGIONS, as the transfer terms are averaged.

Across models: Pearson r and Spearman rho of the at-fixed-G term with each diagnostic change (20
models; EC-Earth3 reported, also a 21-model variant), plus a pooled model x region regression with
region means removed (the per-region at-fixed-G terms of stability_cmip6_split22.json).

Usage: python3 stability_cmip6_analysis.py
Output: dunkelflaute/results/regimes/stability_cmip6.json, stability_cmip6_scatter_{eu8,cmip22}.png
"""
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("regions", "regimes"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/{p}")
import regions as rg                  # noqa: E402

RES = f"{_ROOT}/dunkelflaute/results/regimes"
PARTS = f"{_ROOT}/dunkelflaute/data/regimes/stability_cmip6_parts"
EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
BOX = dict(rg.BOXES, DE=rg.DE_BOX)
C22 = [r for r in rg.REGIONS if BOX[r][0] >= 29.0]
DIFF = ("th850_ths", "th850_tts", "dT850s", "tas", "tas_st", "ta850", "ts", "hfss", "hfss_st", "hfls", "snc", "snw",
        "lai", "treeFrac", "cropFrac", "grassFrac", "sfcWind", "ps")
TCO = {"100km": (("ctl1950g", 1953, 2099), ("hi1950g", 2081, 2099)),
       "hist_fut": (("hi1950g", 1951, 1980), ("hi1950g", 2071, 2099))}


def agg(x, regs, R):
    return float(np.nanmean([x[R.index(r)] for r in regs]))


def cmip6_model(name):
    z = np.load(f"{PARTS}/cmip6_{name}.npz")
    R = list(z["regions"])
    d = {}
    for v in DIFF:
        if f"historical|{v}" not in z.files:
            continue
        h, f = z[f"historical|{v}"].mean(0), z[f"ssp585|{v}"].mean(0)        # (R,)
        d[v] = (h, f)
    return R, d, z


def tco95_comp(comp):
    (cn, c0, c1), (fn, f0, f1) = TCO[comp]
    out = {}
    zc, zf = np.load(f"{PARTS}/tco95_{cn}.npz"), np.load(f"{PARTS}/tco95_{fn}.npz")
    R = list(zc["regions"])

    def wmean(z, y0, y1, v):
        yr = z["years"]
        x = z[v]                                   # (Y, 12, R)
        rows = []
        for w in range(y0, y1 + 1):
            a, b = np.where(yr == w - 1)[0], np.where(yr == w)[0]
            if len(a) and len(b):
                rows.append(np.concatenate([x[a[0], 8:12], x[b[0], 0:3]]).mean(0))
        return np.mean(rows, axis=0)
    for v in DIFF + ("blh", "fsr", "sd", "inv850_z", "skt"):
        key = {"ts": "skt"}.get(v, v)
        if key not in zc.files:
            continue
        out[v] = (wmean(zc, c0, c1, key), wmean(zf, f0, f1, key))
    return R, out


def summarise(R, d):
    """-> {region set: {var: change}}, plus control values."""
    res = {}
    for nm, regs in (("eu8", EU), ("cmip22", C22)):
        o = {}
        for v, (h, f) in d.items():
            if v in ("sfcWind", "ps", "lai", "blh", "fsr", "sd"):
                o[f"d_{v}_pct"] = agg(100 * (f / h - 1), regs, R)
            o[f"d_{v}"] = agg(f - h, regs, R)
            o[f"ctl_{v}"] = agg(h, regs, R)
        if "ta850" in d and "tas_st" in d:
            o["ratio_dtas_dta850"] = o["d_tas_st"] / o["d_ta850"]
        res[nm] = o
    return res


def corr(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5:
        return dict(n=int(m.sum()))
    rx, ry = np.argsort(np.argsort(x[m])), np.argsort(np.argsort(y[m]))
    r = float(np.corrcoef(x[m], y[m])[0, 1])
    b = np.polyfit(x[m], y[m], 1)
    # two-sided p from t with n-2 dof (normal approximation for the tail, scipy-free)
    n = int(m.sum())
    t = r * np.sqrt((n - 2) / max(1e-12, 1 - r * r))
    from math import erf, sqrt
    p = float(2 * (1 - 0.5 * (1 + erf(abs(t) / sqrt(2)))))
    return dict(n=n, pearson=r, spearman=float(np.corrcoef(rx, ry)[0, 1]), slope=float(b[0]), intercept=float(b[1]),
                p_normal_approx=p)


def main():
    inv = json.load(open(f"{_ROOT}/dunkelflaute/data/cmip6/inventory.json"))
    tsc = json.load(open(f"{RES}/transfer_stability_cmip6.json"))["models"]
    s22 = json.load(open(f"{RES}/stability_cmip6_split22.json"))["models"]
    tch = json.load(open(f"{RES}/transfer_checks.json"))
    tst = json.load(open(f"{RES}/transfer_stability.json"))
    models = {}
    pooled = []
    for name in inv:
        R, d, z = cmip6_model(name)
        s = summarise(R, d)
        info = json.loads(str(z["info"]))
        m = dict(wind_only=bool(inv[name].get("wind_only")), dT_global=tch["cmip6"]["models"][name]["dT"],
                 vars=info["vars"], grid_source=info["grid_source"])
        m["eu8"] = dict(within_G=tsc[name]["within_G"], between_G=tsc[name]["between_G"], total=tsc[name]["total"],
                        dG_pct=tsc[name]["mean_G_change_pct"], **s["eu8"])
        m["cmip22"] = dict(within_G=s22[name]["cmip22"]["within_G"], between_G=s22[name]["cmip22"]["between_G"],
                           total=tch["cmip6"]["models"][name]["cmip22"]["transfer_change_pct"],
                           dG_pct=tch["cmip6"]["models"][name]["cmip22"]["mean_G_change_pct"], **s["cmip22"])
        if "sftlf" in z.files:
            m["eu8"]["sftlf"] = agg(z["sftlf"], EU, R)
            m["cmip22"]["sftlf"] = agg(z["sftlf"], C22, R)
        if "n_cells_stab" in z.files:
            m["cells_stab_land"] = [int(z["n_cells_stab"].sum()), int(z["n_cells_land"].sum())]
        models[name] = m
        # pooled per region
        pw = s22[name]["cmip22"]["per_region_within"]
        for r in C22:
            i = R.index(r)
            row = dict(model=name, region=r, within=pw[r], wind_only=m["wind_only"])
            for v, (h, f) in d.items():
                row[f"d_{v}"] = float(f[i] - h[i])
            pooled.append(row)
    ours = {}
    for comp in TCO:
        R, d = tco95_comp(comp)
        s = summarise(R, d)
        g = tst[comp]["B2"]["G_q5"]
        s["eu8"].update(within_G=g["within"]["eu_mean"], between_G=g["between"]["eu_mean"], total=g["total"]["eu_mean"],
                        dG_pct=tst[comp]["B1"]["mean_G_change_pct"])
        if comp == "100km":
            a3 = tch["table"]["A3 bilinear G, 22 regions"]["100km"]
            s["cmip22"].update(total=a3["transfer_change_pct"], dG_pct=a3["mean_G_change_pct"])
        ours[f"tco95_{comp}"] = s
    eerie = f"{PARTS}/eerie_summary.json"
    if os.path.exists(eerie):
        ours.update(json.load(open(eerie)))
    # at-fixed-G split of our runs on the CMIP6 footing (10 m speed, bilinear global G), stability_cmip6_split_ours.json;
    # the TCO95 within_G of transfer_stability.json (log-law 100 m, box G) is kept as within_G_box
    so = json.load(open(f"{RES}/stability_cmip6_split_ours.json"))["runs"]
    alias = {"100km": "tco95_100km"}
    for k, v in so.items():
        o = ours.setdefault(alias.get(k, k), {"eu8": {}, "cmip22": {}})
        for g in ("eu8", "cmip22"):
            x = v[f"spd10|bilin_global|{g}"]
            o.setdefault(g, {})
            if "within_G" in o[g]:
                o[g]["within_G_box"] = o[g]["within_G"]
            o[g].update(within_G=x["within_G"], between_G=x["between_G"], total=x["total"], dG_pct=x["mean_G_change_pct"],
                        within_ci90=x["within_ci90"])
        o["eu8"]["within_G_spd100_box"] = v["spd100|box|eu8"]["within_G"]

    # ensemble stats and correlations
    ens = {n: m for n, m in models.items() if not m["wind_only"]}
    keys = ("within_G", "between_G", "total", "dG_pct", "d_th850_ths", "d_th850_tts", "d_dT850s", "d_tas", "d_tas_st",
            "d_ta850", "ratio_dtas_dta850", "d_hfss", "d_hfls", "d_snc", "d_snw", "d_lai", "d_lai_pct", "d_treeFrac",
            "d_cropFrac", "d_grassFrac", "d_sfcWind_pct", "ctl_th850_ths", "ctl_snc", "ctl_hfss", "sftlf")
    summary, cors = {}, {}
    for g in ("eu8", "cmip22"):
        summary[g] = {}
        for k in keys:
            x = [m[g].get(k, np.nan) for m in ens.values()]
            x = np.array(x, float)
            if np.isfinite(x).sum():
                summary[g][k] = dict(median=float(np.nanmedian(x)), p10_p90=[float(q) for q in np.nanpercentile(x, [10, 90])],
                                     n=int(np.isfinite(x).sum()), n_negative=int((x < 0).sum()))
        cors[g] = {}
        for tgt in ("within_G", "total"):
            for k in keys[3:]:
                a = [m[g].get(tgt, np.nan) for m in ens.values()]
                b = [m[g].get(k, np.nan) for m in ens.values()]
                c = corr(b, a)
                a21 = [m[g].get(tgt, np.nan) for m in models.values()]
                b21 = [m[g].get(k, np.nan) for m in models.values()]
                c21 = corr(b21, a21)
                cors[g][f"{tgt}~{k}"] = dict(cmip20=c, cmip21_with_ecearth3=c21)
            # per K of global warming
            for k in ("d_th850_ths", "d_hfss", "d_snc", "d_tas"):
                a = [m[g].get(tgt, np.nan) / m["dT_global"] for m in ens.values()]
                b = [m[g].get(k, np.nan) / m["dT_global"] for m in ens.values()]
                cors[g][f"{tgt}_perK~{k}_perK"] = dict(cmip20=corr(b, a))
    # pooled model x region, region means removed
    pr = {}
    P = [p for p in pooled if not p["wind_only"]]
    for v in ("th850_ths", "dT850s", "hfss", "snc", "snw", "tas", "lai", "treeFrac", "cropFrac", "grassFrac", "sfcWind"):
        rows = [p for p in P if np.isfinite(p.get(f"d_{v}", np.nan)) and np.isfinite(p["within"])]
        if len(rows) < 20:
            continue
        y = np.array([p["within"] for p in rows])
        x = np.array([p[f"d_{v}"] for p in rows])
        reg = np.array([p["region"] for p in rows])
        mod = np.array([p["model"] for p in rows])
        def demean(a, lab):
            a = a.copy()
            for u in np.unique(lab):
                a[lab == u] -= a[lab == u].mean()
            return a
        pr[v] = dict(raw=corr(x, y), region_demeaned=corr(demean(x, reg), demean(y, reg)),
                     model_demeaned=corr(demean(x, mod), demean(y, mod)), n_models=int(len(np.unique(mod))))
    # where our runs sit relative to the 20-model regressions: residual, in CMIP6 residual standard deviations
    resid = {}
    for g in ("eu8", "cmip22"):
        resid[g] = {}
        for k in ("d_th850_ths", "d_dT850s", "d_hfss", "d_snc", "ratio_dtas_dta850", "d_tas", "dG_pct",
                  "d_cropFrac", "d_treeFrac", "d_lai"):
            a = np.array([[m[g].get(k, np.nan), m[g]["within_G"]] for m in ens.values()], float)
            a = a[np.isfinite(a).all(1)]
            if len(a) < 5:
                continue
            b = np.polyfit(a[:, 0], a[:, 1], 1)
            sd = float(np.std(a[:, 1] - np.polyval(b, a[:, 0]), ddof=2))
            for n, o in ours.items():
                # IFS (OpenIFS 43r3; assumed for EERIE 48r1): climatological LAI and fixed land cover, so no change
                x = 0.0 if k in ("d_cropFrac", "d_treeFrac", "d_lai") else o.get(g, {}).get(k, np.nan)
                y = o.get(g, {}).get("within_G", np.nan)
                if np.isfinite(x) and np.isfinite(y):
                    r_ = float(y - np.polyval(b, x))
                    resid[g].setdefault(n, {})[k] = dict(x=float(x), within_G=float(y), predicted=float(np.polyval(b, x)),
                                                        residual=r_, residual_in_sd=r_ / sd)
    # partial correlations (20 models) of within_G with each land-cover change, controlling for d_tas and dG_pct
    partial = {}
    for g in ("eu8", "cmip22"):
        partial[g] = {}
        for k in ("d_cropFrac", "d_treeFrac", "d_lai", "d_th850_ths", "d_snc", "d_hfss"):
            rows = np.array([[m[g].get(k, np.nan), m[g]["within_G"], m[g]["d_tas"], m[g]["dG_pct"]] for m in ens.values()], float)
            rows = rows[np.isfinite(rows).all(1)]
            if len(rows) < 6:
                continue
            Z = np.column_stack([np.ones(len(rows)), rows[:, 2], rows[:, 3]])
            ex = rows[:, 0] - Z @ np.linalg.lstsq(Z, rows[:, 0], rcond=None)[0]
            ey = rows[:, 1] - Z @ np.linalg.lstsq(Z, rows[:, 1], rcond=None)[0]
            partial[g][k] = dict(n=len(rows), r_partial_given_dtas_dG=float(np.corrcoef(ex, ey)[0, 1]))
    out = dict(meta=dict(regions_eu8=EU, regions_cmip22=C22, cmip6_periods="historical 1951-1980 vs ssp585 2071-2100",
                         tco95=TCO, notes="d_x = future - control winter (Sep-Mar) mean; hfss positive upward (TCO95 converted from IFS sshf); "
                         "th850_ths = theta(850 hPa) - theta(tas, ps) over cells with control ps >= 900 hPa; snc in %; "
                         "ratio_dtas_dta850 uses tas on the same cells as ta850"),
               models=models, ours=ours, ensemble20=summary, correlations=cors, pooled_model_region=pr,
               ours_vs_cmip6_regression=resid, partial_correlations=partial)
    json.dump(out, open(f"{RES}/stability_cmip6.json", "w"), indent=1)
    table(models, ours, summary)
    print("\nCorrelations (20 models) of the at-fixed-G term:")
    for g in ("eu8", "cmip22"):
        for k, c in cors[g].items():
            c0 = c["cmip20"]
            if k.startswith("within_G") and "pearson" in c0:
                c1 = c.get("cmip21_with_ecearth3", {})
                print(f"  {g:6s} {k:34s} r {c0['pearson']:+.2f} rho {c0['spearman']:+.2f} n {c0['n']:2d} slope {c0['slope']:+.3f}"
                      + (f"   (21: r {c1['pearson']:+.2f})" if "pearson" in c1 else ""))
    print("\nPooled model x region (22 regions, 20 models):")
    for v, c in pr.items():
        print(f"  {v:10s} raw r {c['raw']['pearson']:+.2f}  region-demeaned r {c['region_demeaned']['pearson']:+.2f} "
              f"slope {c['region_demeaned']['slope']:+.3f}  model-demeaned r {c['model_demeaned']['pearson']:+.2f} (n {c['raw']['n']})")
    print("\nOur runs against the 20-model regressions (residual of within_G, in residual sd):")
    for g in resid:
        for n, d in resid[g].items():
            print(f"  {g:6s} {n:14s} " + "  ".join(f"{k} {v['residual']:+.2f} ({v['residual_in_sd']:+.1f} sd)" for k, v in d.items()))
    print("\nPartial r given d_tas and dG: " + json.dumps(partial))
    plot(models, ours)
    print(f"Saved {RES}/stability_cmip6.json")


def table(models, ours, summary):
    cols = ("within_G", "total", "d_th850_ths", "d_dT850s", "d_hfss", "d_snc", "d_tas", "d_ta850", "ratio_dtas_dta850", "d_lai", "d_treeFrac")
    for g in ("eu8", "cmip22"):
        print(f"\n{g}: " + " ".join(f"{c:>9s}" for c in cols))
        for n, m in sorted(models.items(), key=lambda kv: kv[1][g]["within_G"]):
            print(f"  {n:18s}" + " ".join(f"{m[g].get(c, np.nan):+9.2f}" for c in cols))
        print(f"  {'median':18s}" + " ".join(f"{summary[g].get(c, {}).get('median', np.nan):+9.2f}" for c in cols))
        for n, o in ours.items():
            if g in o:
                print(f"  {n:18s}" + " ".join(f"{o[g].get(c, np.nan):+9.2f}" for c in cols))


def plot(models, ours):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xs = (("d_th850_ths", "d theta(850) - theta(2 m) (K)"), ("d_hfss", "d hfss (W m-2, up +)"), ("d_snc", "d snow cover (%)"),
          ("d_lai", "d LAI (m2 m-2)"), ("d_treeFrac", "d tree fraction (%)"), ("d_cropFrac", "d crop fraction (%)"))
    for g in ("eu8", "cmip22"):
        fig, axs = plt.subplots(2, 3, figsize=(13, 8))
        for ax, (k, lab) in zip(axs.ravel(), xs):
            for n, m in models.items():
                x, y = m[g].get(k, np.nan), m[g]["within_G"]
                if not np.isfinite(x):
                    continue
                c = "tab:red" if m["wind_only"] else "0.3"
                ax.plot(x, y, "o", color=c, ms=5)
                ax.annotate(n.replace("HadGEM3-GC31-", "Had").replace("-ESM", "E").replace("1-2-", ""), (x, y), fontsize=6, color=c)
            veg = k in ("d_lai", "d_treeFrac", "d_cropFrac")
            for n, o in ours.items():
                x = 0.0 if veg else o.get(g, {}).get(k, np.nan)      # IFS: fixed land cover, climatological LAI
                y = o.get(g, {}).get("within_G", np.nan)
                if np.isfinite(x) and np.isfinite(y):
                    ax.plot(x, y, "*", color="tab:blue", ms=11, mfc="none" if veg else "tab:blue")
                    ax.annotate(n.replace("tco95_", "100km_"), (x, y), fontsize=6, color="tab:blue")
            a = np.array([[m[g].get(k, np.nan), m[g]["within_G"]] for m in models.values() if not m["wind_only"]])
            a = a[np.isfinite(a).all(1)]
            if len(a) > 4:
                r = np.corrcoef(a.T)[0, 1]
                b = np.polyfit(a[:, 0], a[:, 1], 1)
                xx = np.linspace(a[:, 0].min(), a[:, 0].max(), 10)
                ax.plot(xx, np.polyval(b, xx), "-", color="0.6", lw=1)
                ax.set_title(f"r = {r:+.2f} (n = {len(a)})", fontsize=9)
            ax.axhline(0, color="0.8", lw=0.5)
            ax.set_xlabel(lab)
            ax.set_ylabel("transfer change at fixed G (%)")
        fig.suptitle(f"At-fixed-G transfer change against winter diagnostics, {g} (CMIP6 grey, EC-Earth3 red, our runs blue; open stars placed at zero vegetation change)")
        fig.tight_layout()
        fig.savefig(f"{RES}/stability_cmip6_scatter_{g}.png", dpi=110)
        plt.close(fig)


if __name__ == "__main__":
    main()
