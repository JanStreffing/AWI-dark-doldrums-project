"""
How robust is the split of the change in low-wind events into a regime-frequency term and a within-regime term
(fig:mechanism c of the paper)?

Definitions are those of regime_analysis.py (imported, not copied): ERA5 EOFs and centroids from
results/regimes/regime_patterns.npz; model msl anomalies from the seasonal cycle of the comparison's control,
applied to both periods; nearest ERA5 centroid; wind events (48 h, window expansion, control-matched threshold from
results/regions/region_events_summary.json) on the 8 European regions; an event's regime = most frequent regime over
its calendar days (ties go to the lowest regime index, numpy argmax).
Per region, with f_k regime days per winter and lam_k events per regime day,
  total = sum_k lam_k^ctl (f_k^fut - f_k^ctl)  [frequency term]  +  sum_k f_k^fut (lam_k^fut - lam_k^ctl)  [within]
then the mean over the 8 regions, divided by the comparison's dT (results/regions/region_per_K_summary.json), as
paper/scripts/mechanism/reduce_mechanism.py.

Steps (comparisons 9km, 31km, 100km, eerie_ctl of region_events.COMPARISONS):
  0 per winter, region and regime: event counts e[w, r, k] and regime days n[w, k]; any set of winter weights then
    gives f_k = sum_w W_w n[w, k] / sum_w W_w and e_k likewise, so bootstrap and null are weight matrices
  1 reproduce regime_summary.json (frequency_term, within_regime_term, total_change per region) with unit weights
  2 winter bootstrap: control and future winters resampled with replacement, independently, NB draws, the same
    draw for all 8 regions; 8-region means per K of the two terms, the total and the per-regime parts
  3 block null from the 147-winter TCO95 constant-1950 control (tco95_ctl1950g 1953-2099, labels and thresholds of the
    100 km comparison): all ordered placements of two disjoint contiguous blocks with the comparison's winter counts
    (as regime_G_robustness.block_null); 100 km: a 19-winter block against the rest, and 30 vs 19
  4 sensitivities: symmetric (midpoint) decomposition; ties in the event regime split equally between the tied regimes
  5 side check for fig:mechanism d: per-regime change in mean G with bilinear G (transfer_checks_msl_<run>.npz,
    the CMIP6 footing) against box-mean G (the paper's footing for our runs), labels as regime_G_robustness

Environment: levante, module load python3/2025.01-gcc-13.3.0
Usage: python3 regime_decomposition_robustness.py
Output: dunkelflaute/results/regimes/regime_decomposition_robustness.json and .npz (bootstrap and null draws)
"""
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regimes")
import regime_analysis as ra               # noqa: E402  (puts regions/ and analysis/ on the path)
import region_events as ev                 # noqa: E402
WIN_6H = ev.WIN                            # 8 six-hourly steps = 48 h; geostrophic.py sets ev.WIN = 2 on import
import geostrophic as gs                   # noqa: E402
import regime_G_robustness as rgr          # noqa: E402
ev.WIN = WIN_6H
assert ev.WIN == 8

OUT = f"{_ROOT}/dunkelflaute/results/regimes"
RUNS = ["9km", "31km", "100km", "eerie_ctl"]
NULL_RUN = "100km"
EU = ra.EU
K = ra.K
NB = 5000
TERMS = ["total", "frequency", "within"]


def labels(name, pat):
    """Regime labels exactly as regime_analysis.main (all msl days of the winter range, control seasonal cycle)."""
    (cn, c0, c1), (fn, f0, f1) = ev.COMPARISONS[name]
    wts = pat["wts"].ravel()
    out = {}
    dc, xc, *_ = ra.load_msl(cn)
    df_, xf, *_ = ra.load_msl(fn)
    wc, wf = ra.winter_of_day(dc), ra.winter_of_day(df_)
    mc, mf = (wc >= c0) & (wc <= c1), (wf >= f0) & (wf <= f1)
    dc, xc, wc, df_, xf, wf = dc[mc], xc[mc], wc[mc], df_[mf], xf[mf], wf[mf]
    cc = ra.seasonal_fit(dc, xc)
    out["ctl"] = (dc, wc, ra.assign(ra.anomalies(dc, xc, cc) * wts, pat["eofs"], pat["cent"]))
    out["fut"] = (df_, wf, ra.assign(ra.anomalies(df_, xf, cc) * wts, pat["eofs"], pat["cent"]))
    return out


def winter_arrays(name, pat, S):
    """Per period: years (W,), e (W, R, K) events by regime, e_split (ties shared), n (W, K) regime days,
    plus counts of tied and unlabelled events."""
    (cn, c0, c1), (fn, f0, f1) = ev.COMPARISONS[name]
    L = labels(name, pat)
    R = S["regions"]
    out = {}
    for per, run, y0, y1 in (("ctl", cn, c0, c1), ("fut", fn, f0, f1)):
        dates, wins, lab = L[per]
        rd = dict(zip(dates, lab))
        src = ev.load(run)
        years = None
        e = e_split = None
        ties = unlab = tot = 0
        for ir, r in enumerate(EU):
            thr = S["comparisons"][name]["wind"][r]["threshold"]
            E = ra.events_with_days(src, "wind", R.index(r), y0, y1, thr)
            if years is None:
                years = np.array(sorted(E))
                e = np.zeros((len(years), len(EU), K))
                e_split = np.zeros_like(e)
            assert list(years) == sorted(E)
            for iw, y in enumerate(years):
                for days in E[y]:
                    tot += 1
                    rr = [rd[dd] for dd in days if dd in rd]
                    if not rr:
                        unlab += 1
                        continue
                    bc = np.bincount(rr, minlength=K)
                    e[iw, ir, bc.argmax()] += 1
                    top = bc == bc.max()
                    ties += int(top.sum() > 1)
                    e_split[iw, ir] += top / top.sum()
        idx = np.searchsorted(years, wins)
        ok = (idx < len(years)) & (years[np.minimum(idx, len(years) - 1)] == wins)
        n = np.zeros((len(years), K))
        np.add.at(n, (idx[ok], lab[ok]), 1)
        out[per] = dict(years=years, e=e, e_split=e_split, n=n, n_events=tot, n_tied=ties, n_unlabelled=unlab)
    return out


def terms(wc, ec, nc, wf, ef, nf, symmetric=False):
    """Weights (B, Wc), (B, Wf) -> per region and regime (B, R, K): frequency part, within part."""
    fc = (wc @ nc) / wc.sum(1, keepdims=True)                              # (B, K) regime days per winter
    ff = (wf @ nf) / wf.sum(1, keepdims=True)
    Ec = np.einsum("bw,wrk->brk", wc, ec) / wc.sum(1)[:, None, None]       # (B, R, K) events per winter
    Ef = np.einsum("bw,wrk->brk", wf, ef) / wf.sum(1)[:, None, None]
    fc3, ff3 = fc[:, None, :], ff[:, None, :]
    lc = np.divide(Ec, fc3, out=np.zeros_like(Ec), where=fc3 > 0)
    lf = np.divide(Ef, ff3, out=np.zeros_like(Ef), where=ff3 > 0)
    if symmetric:
        return 0.5 * (lc + lf) * (ff3 - fc3), 0.5 * (fc3 + ff3) * (lf - lc)
    return lc * (ff3 - fc3), ff3 * (lf - lc)


def summarise(fq, wi, dT):
    """(B, R, K) parts -> dict of (B,) 8-region means per K for total, frequency, within, and per-regime parts."""
    F, Wt = fq.mean(1) / dT, wi.mean(1) / dT                               # (B, K)
    return dict(total=F.sum(1) + Wt.sum(1), frequency=F.sum(1), within=Wt.sum(1), freq_k=F, within_k=Wt)


def boot_weights(rng, n, B):
    return rgr.boot_weights(rng, n, B)


def placements(W, n1, n2, rest=False):
    """Ordered disjoint contiguous block pairs (as regime_G_robustness.block_null) -> weights (P, W), (P, W)."""
    wc, wf = [], []
    for j in range(W - n2 + 1):
        b = np.zeros(W)
        b[j:j + n2] = 1
        if rest:
            wc.append(1 - b)
            wf.append(b)
            continue
        for i in range(W - n1 + 1):
            if j + n2 <= i or j >= i + n1:
                a = np.zeros(W)
                a[i:i + n1] = 1
                wc.append(a)
                wf.append(b)
    return np.stack(wc), np.stack(wf)


def pctl(x, q=(5, 95)):
    return [float(v) for v in np.percentile(x, q)]


def p_one(null, obs):
    """One-sided p in the direction of the observed sign."""
    return float(np.mean(null >= obs) if obs > 0 else np.mean(null <= obs))


def bilinear_regime_G(name, pat):
    """Per-regime change in mean G (%), box-mean vs bilinear G, days and labels as regime_G_robustness.comparison."""
    P = rgr.comparison(name, pat)
    (cn, _, _), (fn, _, _) = ev.COMPARISONS[name]
    res = {}
    for tag, run in (("ctl", cn), ("fut", fn)):
        m = np.load(f"{ra.D}/transfer_checks_msl_{run}.npz")
        k = gs.ymd_of(m["dates"].astype("datetime64[D]"))
        common, i_p, i_m = np.intersect1d(P[tag]["k"], k, return_indices=True)
        assert len(common) == len(P[tag]["k"]), (name, tag, len(common), len(P[tag]["k"]))
        Gb = gs.region_G(gs.geostrophic(m["msl"][i_m].astype(np.float64), m["lat"], m["lon"], True), m["lat"], m["lon"], EU)
        res[tag] = (P[tag]["lab"][i_p], P[tag]["G"][i_p], Gb)
    out = {}
    for how, j in (("box", 1), ("bilinear", 2)):
        row = {}
        for kk, nm in enumerate(ra.NAMES_ORDER + ["all"]):
            mc = np.ones(len(res["ctl"][0]), bool) if nm == "all" else res["ctl"][0] == kk
            mf = np.ones(len(res["fut"][0]), bool) if nm == "all" else res["fut"][0] == kk
            row[nm] = float(np.mean(100 * (res["fut"][j][mf].mean(0) / res["ctl"][j][mc].mean(0) - 1)))
        out[how] = row
    out["dT"] = float(P["dT"])
    out["per_K"] = {how: {nm: rgr.per_k(v, P["dT"]) for nm, v in out[how].items()} for how in ("box", "bilinear")}
    return out


def main():
    pat = np.load(f"{OUT}/regime_patterns.npz")
    ra.NAMES_ORDER = [str(x) for x in pat["names"]]
    names = ra.NAMES_ORDER
    assert names == ["NAO+", "NAO-", "ScBL", "AR"]
    S = json.load(open(f"{_ROOT}/dunkelflaute/results/regions/region_events_summary.json"))
    stored = json.load(open(f"{OUT}/regime_summary.json"))
    pk = json.load(open(f"{_ROOT}/dunkelflaute/results/regions/region_per_K_summary.json"))
    gstored = json.load(open(f"{OUT}/geostrophic_summary.json"))
    rng = np.random.default_rng(2026)
    res = dict(regimes=names, regions=EU, n_boot=NB, units="events per winter per K (8-region mean) unless noted",
               null_source="tco95_ctl1950g 1953-2099, labels and thresholds of the 100 km comparison", runs={})
    draws = {}
    A = {}
    for name in RUNS:
        A[name] = X = winter_arrays(name, pat, S)
        c, f = X["ctl"], X["fut"]
        dT = pk["dT"][name]["dT"]
        one_c, one_f = np.ones((1, len(c["years"]))), np.ones((1, len(f["years"])))
        fq, wi = terms(one_c, c["e"], c["n"], one_f, f["e"], f["n"])
        tot = (fq + wi).sum(-1)[0]
        # 1 reproduction of regime_summary.json, per region
        rep = {}
        for ir, r in enumerate(EU):
            st = stored["events"][name][f"wind|{r}"]
            rep[r] = dict(frequency_term=[float(fq[0, ir].sum()), st["frequency_term"]],
                          within_regime_term=[float(wi[0, ir].sum()), st["within_regime_term"]],
                          total_change=[float(tot[ir]), st["total_change"]],
                          lam_ctl=[np.divide(c["e"].mean(0)[ir], c["n"].mean(0)).tolist(),
                                   st["events_per_regime_day_ctl"]])
        dmax = max(abs(v[0] - v[1]) for r in rep.values() for kname, v in r.items() if kname != "lam_ctl")
        dlam = max(float(np.max(np.abs(np.array(r["lam_ctl"][0]) - np.array(r["lam_ctl"][1])))) for r in rep.values())
        assert len(c["years"]) == stored["events"][name]["wind|DE"]["ctl"]["n_winters"]
        assert len(f["years"]) == stored["events"][name]["wind|DE"]["fut"]["n_winters"]
        obs = summarise(fq, wi, dT)
        obs = {k: (v[0].tolist() if v.ndim > 1 else float(v[0])) for k, v in obs.items()}
        # 2 bootstrap
        wc, wf = boot_weights(rng, len(c["years"]), NB), boot_weights(rng, len(f["years"]), NB)
        bs = summarise(*terms(wc, c["e"], c["n"], wf, f["e"], f["n"]), dT)
        draws[f"{name}|boot"] = np.column_stack([bs["total"], bs["frequency"], bs["within"], bs["freq_k"], bs["within_k"]])
        # 4 sensitivities
        sym = summarise(*terms(one_c, c["e"], c["n"], one_f, f["e"], f["n"], symmetric=True), dT)
        tie = summarise(*terms(one_c, c["e_split"], c["n"], one_f, f["e_split"], f["n"]), dT)
        # example table: per region, per regime f, lam (per 100 regime days), parts per K
        ctl_f, fut_f = c["n"].mean(0), f["n"].mean(0)
        reg = {}
        for ir, r in enumerate(EU):
            reg[r] = dict(f_ctl=ctl_f.tolist(), f_fut=fut_f.tolist(),
                          e_ctl=c["e"].mean(0)[ir].tolist(), e_fut=f["e"].mean(0)[ir].tolist(),
                          lam100_ctl=(100 * c["e"].mean(0)[ir] / ctl_f).tolist(),
                          lam100_fut=(100 * f["e"].mean(0)[ir] / fut_f).tolist(),
                          freq_k_perK=(fq[0, ir] / dT).tolist(), within_k_perK=(wi[0, ir] / dT).tolist(),
                          frequency_perK=float(fq[0, ir].sum() / dT), within_perK=float(wi[0, ir].sum() / dT),
                          total_perK=float(tot[ir] / dT))
        res["runs"][name] = dict(
            n_winters=[len(c["years"]), len(f["years"])], dT=dT,
            reproduce=dict(max_abs_diff_terms_events_per_winter=float(dmax), max_abs_diff_lam_ctl=dlam,
                           per_region=rep),
            observed=obs,
            boot=dict(p05_p95={t: pctl(bs[t]) for t in TERMS},
                      freq_k_p05_p95=[pctl(bs["freq_k"][:, k]) for k in range(K)],
                      within_k_p05_p95=[pctl(bs["within_k"][:, k]) for k in range(K)],
                      share_positive={t: float(np.mean(bs[t] > 0)) for t in TERMS},
                      within_ScBL_share_of_within_median=float(np.median(bs["within_k"][:, 2] / bs["within"]))),
            symmetric={k: (v[0].tolist() if v.ndim > 1 else float(v[0])) for k, v in sym.items()},
            ties_split={k: (v[0].tolist() if v.ndim > 1 else float(v[0])) for k, v in tie.items()},
            events=dict(ctl=dict(n=c["n_events"], tied=c["n_tied"], unlabelled=c["n_unlabelled"]),
                        fut=dict(n=f["n_events"], tied=f["n_tied"], unlabelled=f["n_unlabelled"])),
            within_positive_regions=int(sum(v["within_perK"] > 0 for v in reg.values())),
            regions=reg)
        print(f"{name}: {len(c['years'])} vs {len(f['years'])} winters, dT {dT:.3f}; reproduce max|d| {dmax:.2e} "
              f"(lam {dlam:.1e}); per K total {obs['total']:+.3f} freq {obs['frequency']:+.3f} within {obs['within']:+.3f}; "
              f"boot within {pctl(bs['within'])}", flush=True)

    # 3 block null
    c0 = A[NULL_RUN]["ctl"]
    W0 = len(c0["years"])
    res["null_control"] = dict(n_winters=W0, years=[int(c0["years"][0]), int(c0["years"][-1])],
                               events_per_winter=float(c0["e"].sum((1, 2)).mean() / len(EU)))
    for name in RUNS:
        nc, nf = res["runs"][name]["n_winters"]
        dT = res["runs"][name]["dT"]
        obs = res["runs"][name]["observed"]
        variants = {"rest_vs_19": (None, nf, True), "30_vs_19": (30, nf, False)} if name == NULL_RUN \
            else {f"{nc}_vs_{nf}": (nc, nf, False)}
        R = {}
        for vn, (n1, n2, rest) in variants.items():
            wc, wf = placements(W0, n1, n2, rest)
            nl = summarise(*terms(wc, c0["e"], c0["n"], wf, c0["e"], c0["n"]), dT)
            draws[f"{name}|null|{vn}"] = np.column_stack([nl["total"], nl["frequency"], nl["within"], nl["freq_k"], nl["within_k"]])
            R[vn] = dict(n_placements=len(wc),
                         p05_p95={t: pctl(nl[t]) for t in TERMS},
                         p_one_sided={t: p_one(nl[t], obs[t]) for t in TERMS},
                         p_two_sided={t: float(np.mean(np.abs(nl[t]) >= abs(obs[t]))) for t in TERMS},
                         within_k_p05_p95=[pctl(nl["within_k"][:, k]) for k in range(K)],
                         within_k_p_one_sided=[p_one(nl["within_k"][:, k], obs["within_k"][k]) for k in range(K)],
                         freq_k_p05_p95=[pctl(nl["freq_k"][:, k]) for k in range(K)])
        res["runs"][name]["block_null"] = R
        v = list(R.values())[0]
        print(f"  null {name} {list(R)[0]} ({v['n_placements']}): " + "  ".join(
            f"{t} {obs[t]:+.3f} [{v['p05_p95'][t][0]:+.3f},{v['p05_p95'][t][1]:+.3f}] p {v['p_one_sided'][t]:.4f}" for t in TERMS), flush=True)

    # 5 bilinear vs box-mean G per regime
    res["G_footing"] = {}
    for name in RUNS:
        g = bilinear_regime_G(name, pat)
        ref = gstored["ours"][name]["field"]["europe"]["regimes"]
        g["box_matches_geostrophic_summary_max_abs_diff_pct"] = float(max(abs(g["box"][nm] - ref[nm]["mean_G_change_pct"]) for nm in names))
        res["G_footing"][name] = g
        print(f"  G {name}: box " + " ".join(f"{nm} {g['box'][nm]:+.2f}" for nm in g["box"]) + " | bilinear "
              + " ".join(f"{nm} {g['bilinear'][nm]:+.2f}" for nm in g["bilinear"])
              + f" | box vs stored {g['box_matches_geostrophic_summary_max_abs_diff_pct']:.1e}", flush=True)

    json.dump(res, open(f"{OUT}/regime_decomposition_robustness.json", "w"), indent=1)
    np.savez(f"{OUT}/regime_decomposition_robustness.npz",
             columns=np.array(["total", "frequency", "within"] + [f"freq_{n}" for n in names] + [f"within_{n}" for n in names]),
             **{k.replace("|", "__"): v for k, v in draws.items()})
    print(f"Saved {OUT}/regime_decomposition_robustness.json and .npz")


if __name__ == "__main__":
    main()
