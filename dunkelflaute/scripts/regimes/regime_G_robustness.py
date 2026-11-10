"""
Is the within-regime weakening of the geostrophic speed G robust against sampling and natural variability?

Definitions are those of geostrophic.py and regime_analysis.py, imported rather than copied:
  G      region mean of |grad p| / (rho f) from daily msl on the 2.5 deg grid, for the 8 European regions
  days   full days of the region speed series that have msl, Sep-Mar (as geostrophic.ours())
  regime projection of msl anomalies (control seasonal cycle, both periods) on the ERA5 EOFs, nearest ERA5 centroid
  change 100 (mean_fut G / mean_ctl G - 1) per region over all regime days of the period, then the mean over regions
         ("all" = the same over all days, = field mean_G_change_pct)
The script first reproduces geostrophic_summary.json (ours -> field -> europe -> regimes -> mean_G_change_pct).

Tests (per comparison, per regime and "all"):
  1 winter bootstrap: control and future winters resampled with replacement, independently (NB draws)
  2 natural-variability null from the 147-winter TCO95 1950 control (tco95_ctl1950g 1953-2099, its own labels from the
    100 km comparison): all ordered placements of two disjoint contiguous blocks of winters with the comparison's
    winter counts; "change" = second block against first. 100 km: a 19-winter block against the rest of the control
    (closest to 147 vs 19), and 30 vs 19.
  3 interannual sd of the winter-mean within-regime G (Europe mean), signal / sd and signal / standard error
  4 assignment artefacts: future labelled with its own seasonal cycle; change in the core (nearest tercile of the
    control distance to the centroid) and with future days reweighted to the control distance terciles
  5 CMIP6 (20 models, historical 1951-1980 vs ssp585 2071-2100): same bootstrap per model

Usage: python3 regime_G_robustness.py
Output: dunkelflaute/results/regimes/regime_G_robustness.json, regime_G_robustness/*.png
"""
import glob
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regimes")
import geostrophic as gs                   # noqa: E402  (also puts regions/analysis/cmip6 on the path)
import regime_analysis as ra               # noqa: E402
import regions as rg                       # noqa: E402
import region_events as ev                 # noqa: E402
import region_per_K as rpk                 # noqa: E402

import matplotlib                          # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt            # noqa: E402

OUT = f"{_ROOT}/dunkelflaute/results/regimes"
FIG = f"{OUT}/regime_G_robustness"
RUNS = ["9km", "31km", "100km", "eerie_ctl"]
NULL_RUN = "100km"            # its control is the 147-winter constant-forcing null
NB = 5000
EU = ra.EU
MIN_WDAYS = 5                 # regime days a winter needs for a winter-mean regime G (test 3)


def per_k(pct, dT):
    return float(100 * (np.exp(np.log1p(pct / 100) / dT) - 1))


def run_days(run):
    """ymd, G (D, 8), dates, msl for the days geostrophic.ours() uses."""
    z = np.load(f"{rg.OUT_DIR}/{run}.npz")
    t = z["time"].astype("datetime64[s]")
    day = (t - np.timedelta64(1, "s")).astype("datetime64[D]")
    u, cnt = np.unique(day, return_counts=True)
    u = u[cnt == cnt.max()]
    m = np.load(f"{ra.D}/msl_{run}.npz")
    md = m["dates"].astype("datetime64[D]")
    common, _, i_m = np.intersect1d(gs.ymd_of(u), gs.ymd_of(md), return_indices=True)
    keep = np.isin(common // 100 % 100, rg.MONTHS)
    common, i_m = common[keep], i_m[keep]
    msl = m["msl"][i_m].astype(np.float64)
    G = gs.region_G(gs.geostrophic(msl, m["lat"], m["lon"], False), m["lat"], m["lon"], EU)
    return common, G, md[i_m], msl


def classify(pat, dates, msl, coef):
    """ra.assign plus the distance to the assigned centroid in PC space."""
    pcs = (ra.anomalies(dates, msl, coef) * pat["wts"].ravel()) @ pat["eofs"].T
    d2 = ((pcs[:, None, :] - pat["cent"][None]) ** 2).sum(-1)
    return d2.argmin(1), np.sqrt(d2.min(1))


def winter_sums(win, lab, G, years):
    """Per winter: G sums (W, 5, R) and day counts (W, 5) for the 4 regimes and all days."""
    S = np.zeros((len(years), 5, G.shape[1]))
    N = np.zeros((len(years), 5))
    idx = np.searchsorted(years, win)
    assert np.all(years[idx] == win)
    for j in range(4):
        m = lab == j
        np.add.at(S[:, j], idx[m], G[m])
        np.add.at(N[:, j], idx[m], 1)
    np.add.at(S[:, 4], idx, G)
    np.add.at(N[:, 4], idx, 1)
    return S, N


def change(wc, Sc, Nc, wf, Sf, Nf):
    """Weights over winters (B, Wc), (B, Wf) -> change in % (B, 5), mean over regions of the per-region change."""
    gc = np.einsum("bw,wjr->bjr", wc, Sc) / (wc @ Nc)[:, :, None]
    gf = np.einsum("bw,wjr->bjr", wf, Sf) / (wf @ Nf)[:, :, None]
    return 100 * (gf / gc - 1).mean(-1)


def boot_weights(rng, n, B):
    W = np.zeros((B, n))
    np.add.at(W, (np.repeat(np.arange(B), n), rng.integers(0, n, B * n)), 1)
    return W


def bootstrap(P, rng):
    (Sc, Nc), (Sf, Nf) = P["ctl"]["SN"], P["fut"]["SN"]
    return change(boot_weights(rng, len(Sc), NB), Sc, Nc, boot_weights(rng, len(Sf), NB), Sf, Nf)


def block_null(S, N, n1, n2, rest=False):
    """All ordered placements of disjoint contiguous blocks (n1 control, n2 'future'); rest: control = all other winters."""
    W = len(S)
    pairs = []
    for j in range(W - n2 + 1):
        if rest:
            a = np.ones(W)
            a[j:j + n2] = 0
            pairs.append((a, j))
            continue
        for i in range(W - n1 + 1):
            if j + n2 <= i or j >= i + n1:
                a = np.zeros(W)
                a[i:i + n1] = 1
                pairs.append((a, j))
    wc = np.stack([a for a, _ in pairs])
    wf = np.zeros_like(wc)
    for k, (_, j) in enumerate(pairs):
        wf[k, j:j + n2] = 1
    return change(wc, S, N, wf, S, N)


def winter_means(S, N):
    """Winter-mean regime G, Europe mean of the region means (W, 5); NaN where fewer than MIN_WDAYS days."""
    with np.errstate(invalid="ignore", divide="ignore"):
        x = (S / N[:, :, None]).mean(-1)
    x[N < MIN_WDAYS] = np.nan
    return x


def pctl(x):
    return [float(v) for v in np.percentile(x, [5, 95])]


def comparison(name, pat):
    (cn, c0, c1), (fn, f0, f1) = ev.COMPARISONS[name]
    P = {}
    for tag, run, y0, y1 in (("ctl", cn, c0, c1), ("fut", fn, f0, f1)):
        k, G, d, msl = run_days(run)
        sel = np.isin(gs.winter_of(k), range(y0, y1 + 1))
        assert np.isfinite(G[sel]).all()
        P[tag] = dict(k=k[sel], G=G[sel], d=d[sel], msl=msl[sel], win=gs.winter_of(k[sel]))
    coef = ra.seasonal_fit(P["ctl"]["d"], P["ctl"]["msl"])
    for tag in P:
        P[tag]["lab"], P[tag]["dist"] = classify(pat, P[tag]["d"], P[tag]["msl"], coef)
    # labels identical to geostrophic.regime_labels
    lc, lf = gs.regime_labels(pat, P["ctl"]["d"], P["ctl"]["msl"], P["fut"]["d"], P["fut"]["msl"])
    assert (lc == P["ctl"]["lab"]).all() and (lf == P["fut"]["lab"]).all()
    lf_own, _ = classify(pat, P["fut"]["d"], P["fut"]["msl"], ra.seasonal_fit(P["fut"]["d"], P["fut"]["msl"]))
    for tag in P:
        P[tag]["years"] = np.unique(P[tag]["win"])
        P[tag]["SN"] = winter_sums(P[tag]["win"], P[tag]["lab"], P[tag]["G"], P[tag]["years"])
    P["fut"]["SN_own"] = winter_sums(P["fut"]["win"], lf_own, P["fut"]["G"], P["fut"]["years"])
    T = rpk.global_t()
    P["dT"] = rpk.mean_years(T, fn, f0, f1)[0] - rpk.mean_years(T, cn, c0, c1)[0]
    return P


def obs_change(P, fut_key="SN"):
    (Sc, Nc), (Sf, Nf) = P["ctl"]["SN"], P["fut"][fut_key]
    return change(np.ones((1, len(Sc))), Sc, Nc, np.ones((1, len(Sf))), Sf, Nf)[0]


def distance_checks(P):
    """Per regime: future share of days in the control distance terciles; change in the core tercile; change with the
    future reweighted to equal tercile shares."""
    out = {}
    c, f = P["ctl"], P["fut"]
    for j, nm in enumerate(ra.NAMES_ORDER):
        mc, mf = c["lab"] == j, f["lab"] == j
        e = np.percentile(c["dist"][mc], [100 / 3, 200 / 3])
        tc, tf = np.digitize(c["dist"][mc], e), np.digitize(f["dist"][mf], e)
        Gc = np.stack([c["G"][mc][tc == t].mean(0) for t in range(3)])     # (3, R)
        Gf = np.stack([f["G"][mf][tf == t].mean(0) for t in range(3)])
        out[nm] = dict(fut_tercile_share=[float(np.mean(tf == t)) for t in range(3)],
                       ctl_G_by_tercile=[float(v) for v in Gc.mean(1)],
                       core_change_pct=float(100 * (Gf[0] / Gc[0] - 1).mean()),
                       outer_change_pct=float(100 * (Gf[2] / Gc[2] - 1).mean()),
                       reweighted_change_pct=float(100 * (Gf.mean(0) / Gc.mean(0) - 1).mean()),
                       mean_dist_change_pct=float(100 * (f["dist"][mf].mean() / c["dist"][mc].mean() - 1)))
    return out


def interannual(P, chg):
    (Sc, Nc), (Sf, Nf) = P["ctl"]["SN"], P["fut"]["SN"]
    xc, xf = winter_means(Sc, Nc), winter_means(Sf, Nf)
    out = {}
    for j, nm in enumerate(ra.NAMES_ORDER + ["all"]):
        a, b = xc[:, j][np.isfinite(xc[:, j])], xf[:, j][np.isfinite(xf[:, j])]
        sc, sf = 100 * a.std(ddof=1) / a.mean(), 100 * b.std(ddof=1) / a.mean()
        se = np.sqrt(sc ** 2 / len(a) + sf ** 2 / len(b))
        out[nm] = dict(sd_ctl_pct=float(sc), sd_fut_pct=float(sf), n_ctl=len(a), n_fut=len(b), se_pct=float(se),
                       change_over_sd=float(chg[j] / sc), change_over_se=float(chg[j] / se))
    return out


def cmip6_model(path, pat):
    """Same footing as geostrophic.cmip6(), European regions only."""
    model = os.path.basename(path)[4:-4]
    c0, c1, f0, f1 = 1951, 1980, 2071, 2100
    Pz, Z = np.load(path), np.load(f"{gs.DC}/daily_cmip6_{model}.npz")
    lat, lon = Pz["lat"], Pz["lon"]
    ii = np.where(lat >= 20.0)[0]
    jj = np.where((lon >= -80.0) & (lon <= 40.0))[0]
    P = {}
    for tag, exp, y0, y1 in (("ctl", "historical", c0, c1), ("fut", "ssp585", f0, f1)):
        kp = Pz[f"{exp}|ymd"] @ np.array([10000, 100, 1])
        kz = Z[f"{exp}|ymd"] @ np.array([10000, 100, 1])
        common, i_p, _ = np.intersect1d(kp, kz, return_indices=True)
        w = gs.winter_of(common)
        keep = (w >= y0) & (w <= y1)
        common, i_p = common[keep], i_p[keep]
        psl = Pz[f"{exp}|psl"][i_p].astype(np.float64)
        G = gs.region_G(gs.geostrophic(psl, lat, lon, True), lat, lon, EU)
        P[tag] = dict(k=common, G=G, msl=psl[:, ii][:, :, jj], win=gs.winter_of(common), d=gs.dates_of(common))
    lc, lf = gs.regime_labels(pat, P["ctl"]["d"], P["ctl"]["msl"], P["fut"]["d"], P["fut"]["msl"])
    for tag, lab in (("ctl", lc), ("fut", lf)):
        P[tag]["years"] = np.unique(P[tag]["win"])
        P[tag]["SN"] = winter_sums(P[tag]["win"], lab, P[tag]["G"], P[tag]["years"])
        del P[tag]["msl"]
    T = {e: dict(zip(Z[f"{e}|tas_years"].tolist(), Z[f"{e}|tas"].tolist())) for e in ("historical", "ssp585")}
    P["dT"] = float(np.mean([T["ssp585"][y] for y in range(f0, f1 + 1) if y in T["ssp585"]])
                    - np.mean([T["historical"][y] for y in range(c0, c1 + 1) if y in T["historical"]]))
    return model, P


def main():
    os.makedirs(FIG, exist_ok=True)
    pat = np.load(f"{ra.OUT}/regime_patterns.npz")
    ra.NAMES_ORDER = [str(x) for x in pat["names"]]
    COLS = ra.NAMES_ORDER + ["all"]
    stored = json.load(open(f"{OUT}/geostrophic_summary.json"))
    rng = np.random.default_rng(11)
    res = dict(columns=COLS, n_boot=NB, runs={}, cmip6={}, null_source="tco95_ctl1950g 1953-2099 (100 km comparison labels)")

    comps = {}
    for name in RUNS:
        P = comparison(name, pat)
        comps[name] = P
        obs = obs_change(P)
        st = stored["ours"][name]["field"]
        ref = [st["europe"]["regimes"][nm]["mean_G_change_pct"] for nm in ra.NAMES_ORDER] + [st["mean_G_change_pct"]]
        dmax = float(np.max(np.abs(obs - np.array(ref))))
        print(f"{name}: {len(P['ctl']['years'])} vs {len(P['fut']['years'])} winters, dT {P['dT']:.3f} (stored {stored['ours'][name]['dT']:.3f}); "
              f"max |obs - stored| {dmax:.2e} %", flush=True)
        bs = bootstrap(P, rng)
        own = obs_change(P, "SN_own")
        res["runs"][name] = dict(
            n_winters=[len(P["ctl"]["years"]), len(P["fut"]["years"])], dT=P["dT"], reproduce_max_abs_diff_pct=dmax,
            change_pct=dict(zip(COLS, obs.tolist())), stored_change_pct=dict(zip(COLS, ref)),
            change_per_K=dict(zip(COLS, [per_k(v, P["dT"]) for v in obs])),
            boot_p05_p95=dict(zip(COLS, [pctl(bs[:, j]) for j in range(5)])),
            boot_same_sign=dict(zip(COLS, [float(np.mean(np.sign(bs[:, j]) == np.sign(obs[j]))) for j in range(5)])),
            interannual=interannual(P, obs),
            ownclim_change_pct=dict(zip(COLS, own.tolist())),
            distance=distance_checks(P),
            regime_days=dict(ctl=P["ctl"]["SN"][1].sum(0).tolist(), fut=P["fut"]["SN"][1].sum(0).tolist()))
        res["runs"][name]["_boot"] = bs

    # natural-variability null from the TCO95 control
    S0, N0 = comps[NULL_RUN]["ctl"]["SN"]
    yrs0 = comps[NULL_RUN]["ctl"]["years"]
    x0 = winter_means(S0, N0)
    tr = np.polyfit(yrs0, 100 * x0[:, 4] / x0[:, 4].mean(), 1)[0] * 100
    res["null_control"] = dict(n_winters=len(yrs0), all_days_trend_pct_per_century=float(tr),
                               interannual_sd_pct={nm: float(100 * np.nanstd(x0[:, j], ddof=1) / np.nanmean(x0[:, j])) for j, nm in enumerate(COLS)})
    print(f"null control: {len(yrs0)} winters, all-day G trend {tr:+.2f} % per century", flush=True)
    for name in RUNS:
        nc, nf = res["runs"][name]["n_winters"]
        obs = np.array([res["runs"][name]["change_pct"][c] for c in COLS])
        variants = {"rest_vs_19": block_null(S0, N0, None, nf, rest=True), "30_vs_19": block_null(S0, N0, 30, nf)} \
            if name == NULL_RUN else {f"{nc}_vs_{nf}": block_null(S0, N0, nc, nf)}
        R = {}
        for vn, nl in variants.items():
            R[vn] = dict(n_placements=len(nl),
                         p05_p95=dict(zip(COLS, [pctl(nl[:, j]) for j in range(5)])),
                         sd=dict(zip(COLS, [float(nl[:, j].std()) for j in range(5)])),
                         p_one_sided=dict(zip(COLS, [float(np.mean(nl[:, j] <= obs[j]) if obs[j] < 0 else np.mean(nl[:, j] >= obs[j])) for j in range(5)])),
                         p_two_sided=dict(zip(COLS, [float(np.mean(np.abs(nl[:, j]) >= abs(obs[j]))) for j in range(5)])))
            R[vn]["_null"] = nl
        res["runs"][name]["block_null"] = R

    # CMIP6
    for path in sorted(glob.glob(f"{gs.DC}/psl_*.npz")):
        model = os.path.basename(path)[4:-4]
        if model not in stored["cmip6"]:
            continue
        model, P = cmip6_model(path, pat)
        obs = obs_change(P)
        st = stored["cmip6"][model]["field"]["europe"]["regimes"]
        ref = [st[nm]["mean_G_change_pct"] for nm in ra.NAMES_ORDER]
        dmax = float(np.max(np.abs(obs[:4] - np.array(ref))))
        bs = bootstrap(P, rng)
        res["cmip6"][model] = dict(
            dT=P["dT"], reproduce_max_abs_diff_pct=dmax, change_pct=dict(zip(COLS, obs.tolist())),
            change_per_K=dict(zip(COLS, [per_k(v, P["dT"]) for v in obs])),
            boot_p05_p95=dict(zip(COLS, [pctl(bs[:, j]) for j in range(5)])),
            boot_same_sign=dict(zip(COLS, [float(np.mean(np.sign(bs[:, j]) == np.sign(obs[j]))) for j in range(5)])),
            interannual=interannual(P, obs))
        print(f"  {model:18s} ScBL {obs[2]:+6.2f} % [{np.percentile(bs[:, 2], 5):+6.2f}, {np.percentile(bs[:, 2], 95):+6.2f}]  "
              f"all {obs[4]:+6.2f} %  max|d stored| {dmax:.1e}", flush=True)
    ens = {}
    for j, nm in enumerate(COLS):
        ch = np.array([v["change_pct"][nm] for v in res["cmip6"].values()])
        hi = np.array([v["boot_p05_p95"][nm][1] for v in res["cmip6"].values()])
        lo = np.array([v["boot_p05_p95"][nm][0] for v in res["cmip6"].values()])
        n = len(ch)
        nneg = int((ch < 0).sum())
        from math import comb
        ens[nm] = dict(n_models=n, n_negative=nneg, n_boot90_below_zero=int((hi < 0).sum()), n_boot90_above_zero=int((lo > 0).sum()),
                       median_change_pct=float(np.median(ch)),
                       sign_test_p_two_sided=float(min(1.0, 2 * sum(comb(n, k) for k in range(max(nneg, n - nneg), n + 1)) / 2 ** n)))
    res["cmip6_ensemble"] = ens

    plots(res, comps)
    for name in RUNS:
        res["runs"][name].pop("_boot")
        for v in res["runs"][name]["block_null"].values():
            v.pop("_null")
    json.dump(res, open(f"{OUT}/regime_G_robustness.json", "w"), indent=1)
    report(res)


def plots(res, comps):
    COLS = res["columns"]
    lab = {"9km": "9 km", "31km": "31 km", "100km": "100 km", "eerie_ctl": "IFS-FESOM"}
    # 1: change with bootstrap and null ranges
    fig, axs = plt.subplots(1, len(RUNS), figsize=(13, 3.6), sharey=True)
    for ax, name in zip(axs, RUNS):
        r = res["runs"][name]
        nv = list(r["block_null"])[0]
        for j, c in enumerate(COLS):
            nl = r["block_null"][nv]["p05_p95"][c]
            ax.add_patch(plt.Rectangle((j - 0.3, nl[0]), 0.6, nl[1] - nl[0], color="0.85", lw=0))
            b = r["boot_p05_p95"][c]
            ax.plot([j, j], b, color="C0", lw=2)
            ax.plot(j, r["change_pct"][c], "o", color="C0")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xticks(range(len(COLS)), COLS)
        ax.set_title(f"{lab[name]} ({r['n_winters'][0]} vs {r['n_winters'][1]} winters)", fontsize=9)
    axs[0].set_ylabel("change in mean G (%)")
    fig.text(0.5, 0.005, "dot: change; blue bar: winter bootstrap 5-95 %; grey box: TCO95-control block null 5-95 %"
             f" ({', '.join(list(res['runs'][n]['block_null'])[0] for n in RUNS)})", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(f"{FIG}/change_boot_null.png", dpi=150)
    plt.close(fig)

    # 2: null histograms, ScBL and all days
    fig, axs = plt.subplots(2, len(RUNS), figsize=(13, 5.5))
    for i, name in enumerate(RUNS):
        r = res["runs"][name]
        nv = list(r["block_null"])[0]
        for row, c in enumerate(("ScBL", "all")):
            j = COLS.index(c)
            ax = axs[row, i]
            ax.hist(r["block_null"][nv]["_null"][:, j], bins=40, color="0.6")
            ax.hist(r["_boot"][:, j], bins=40, color="C0", alpha=0.45)
            ax.axvline(r["change_pct"][c], color="C3")
            ax.axvline(0, color="k", lw=0.6)
            ax.set_title(f"{lab[name]} {c}: p={r['block_null'][nv]['p_one_sided'][c]:.3f}", fontsize=9)
            ax.set_yticks([])
    axs[1, 0].set_xlabel("change in mean G (%)")
    fig.text(0.5, 0.005, "grey: TCO95-control block null; blue: winter bootstrap; red: observed", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(f"{FIG}/null_hist_ScBL_all.png", dpi=150)
    plt.close(fig)

    # 3: CMIP6 per model, all regimes
    ms = sorted(res["cmip6"], key=lambda m: res["cmip6"][m]["change_pct"]["ScBL"])
    fig, axs = plt.subplots(1, len(COLS), figsize=(14, 5), sharey=True)
    for ax, c in zip(axs, COLS):
        for i, m in enumerate(ms):
            v = res["cmip6"][m]
            b = v["boot_p05_p95"][c]
            col = "C3" if b[1] < 0 else ("C0" if b[0] > 0 else "0.4")
            ax.plot(b, [i, i], color=col, lw=1.5)
            ax.plot(v["change_pct"][c], i, "o", color=col, ms=4)
        ax.axvline(0, color="k", lw=0.6)
        ax.set_title(f"{c}: {res['cmip6_ensemble'][c]['n_boot90_below_zero']}/{len(ms)} below 0", fontsize=9)
        ax.set_xlabel("change in mean G (%)")
    axs[0].set_yticks(range(len(ms)), ms, fontsize=7)
    fig.text(0.5, 0.005, "CMIP6 1951-80 vs 2071-2100 ssp585, winter bootstrap 5-95 %; red: range below zero", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(f"{FIG}/cmip6_boot.png", dpi=150)
    plt.close(fig)

    # 4: winter-mean ScBL and all-day G in the TCO95 control, with the comparison periods
    P = comps[NULL_RUN]
    x = winter_means(*P["ctl"]["SN"])
    yrs = P["ctl"]["years"]
    fig, ax = plt.subplots(figsize=(10, 3.3))
    for j, c in ((2, "ScBL"), (4, "all")):
        ax.plot(yrs, x[:, j], label=f"{c} (sd {res['null_control']['interannual_sd_pct'][c]:.1f} %)")
    ax.set_ylabel("winter-mean G, Europe (m/s)")
    ax.set_title("TCO95 1950 control: winter-mean within-regime G", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{FIG}/tco95_control_winter_G.png", dpi=150)
    plt.close(fig)


def report(res):
    COLS = res["columns"]
    print("\nrun        col    chg%   %/K  boot5-95        same  null5-95        p1    p2    sd%   chg/sd chg/se  own%   core%  rewt%")
    for name, r in res["runs"].items():
        nv = list(r["block_null"])[0]
        for c in COLS:
            b, nl = r["boot_p05_p95"][c], r["block_null"][nv]["p05_p95"][c]
            ia = r["interannual"][c]
            dd = r["distance"].get(c, {})
            print(f"{name:10s} {c:5s} {r['change_pct'][c]:+6.2f} {r['change_per_K'][c]:+5.2f} [{b[0]:+6.2f},{b[1]:+6.2f}] {r['boot_same_sign'][c]:.2f}"
                  f"  [{nl[0]:+6.2f},{nl[1]:+6.2f}] {r['block_null'][nv]['p_one_sided'][c]:.3f} {r['block_null'][nv]['p_two_sided'][c]:.3f}"
                  f" {ia['sd_ctl_pct']:5.2f} {ia['change_over_sd']:+6.2f} {ia['change_over_se']:+6.2f} {r['ownclim_change_pct'][c]:+6.2f}"
                  + (f" {dd['core_change_pct']:+6.2f} {dd['reweighted_change_pct']:+6.2f}  fut terciles {['%.2f' % v for v in dd['fut_tercile_share']]}" if dd else ""))
        if len(r["block_null"]) > 1:
            for vn, v in list(r["block_null"].items())[1:]:
                print(f"{'':10s} null {vn}: " + "  ".join(f"{c} [{v['p05_p95'][c][0]:+.2f},{v['p05_p95'][c][1]:+.2f}] p {v['p_one_sided'][c]:.3f}" for c in COLS))
    print(f"\nnull control: {res['null_control']}")
    print("\nCMIP6 ensemble:")
    for c, v in res["cmip6_ensemble"].items():
        print(f"  {c:5s} {v}")
    print(f"Saved {OUT}/regime_G_robustness.json and {FIG}/*.png")


if __name__ == "__main__":
    main()
