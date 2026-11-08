"""
Field significance of the gridded per-K maps (5 x 5 deg land boxes, 20-70 N; methods 15) for our runs.

Null: the TCO95 constant-1950 control (ctl1950g, 147 winters 1953-2099), on the 6-hourly footing, with
the threshold of every box matched to 5 events per control winter exactly as in grid_per_K.py. For each
comparison, all placements of two disjoint blocks of whole consecutive winters with the comparison's
own winter counts (the control block cut to W - n_b - 40 where too few placements remain, as in
grid_per_K.null_ranges) give a null ln(ratio) per box. The SAME placement is used in every box, so the
null keeps the spatial correlation of the real samples, as the region field test of methods 10c/14
(region_events.block_pairs / block_logratios, reused here).

Per box: two-sided p = min(1, 2 min(p_up, p_down)), p_up = (1 + #{null >= obs}) / (1 + P), and
likewise p_down (empirical; P placements). A Gaussian variant (z from the mean and standard deviation
of the null placements) is computed as a check on the empirical tails.

Multiple testing: false discovery rate control (Benjamini and Hochberg 1995), as recommended for
gridded fields by Wilks (2016, BAMS 97, 2263-2273) with alpha_FDR = 0.10: sort the N valid p-values,
find the largest k with p_(k) <= alpha_FDR k / N, and mark every box with p <= p_(k). The field is
significant (Wilks' global test) if at least one box is marked.

Field significance of the area mean: statistic = cos-lat weighted mean of ln(ratio) over the boxes
valid in the comparison (all null placements are finite there); null = the same weighted mean for
every placement; one-sided p = (1 + #{null >= obs}) / (1 + P).

CMIP6 is not tested this way (it has no common control); the map marks model agreement instead.
This script also lists every CMIP6 model's own area-mean per K, from grid_per_K.npz.

Usage: sbatch --job-name=fieldsig dunkelflaute/jobs/regions_py.sbatch regions/fieldsig_grid.py
       (loads the 6 GB TCO95 grid5 control; the null counts are cached in results/maps/fieldsig_null_counts.npz,
       later runs only redo the tests, a few seconds on the login node)
Output: dunkelflaute/results/maps/fieldsig_grid.npz, fieldsig_grid_summary.json
"""
import json
import math
import os
import sys

os.environ["REGION_SET"] = "grid5"
import numpy as np                                                     # noqa: E402

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
import grid_per_K as gpk                                               # noqa: E402

ev, rg = gpk.ev, gpk.rg
OUT = f"{_ROOT}/dunkelflaute/results/maps"
NULL_CACHE = f"{OUT}/fieldsig_null_counts.npz"
ALPHA_FDR = 0.10
KINDS = ("wind", "comb")
# winters (control, future) per comparison, as counted by grid_per_K.py (job 28763196 log)
SIZES = {"9km": (17, 10), "31km": (17, 19), "100km": (147, 19), "eerie_ctl": (79, 30), "eerie_hist": (30, 30)}


def null_counts():
    """Per-winter event counts (W, nB) in the TCO95 control at the per-box matched threshold."""
    if os.path.exists(NULL_CACHE):
        z = np.load(NULL_CACHE)
        return {k: z[f"{k}|counts"] for k in KINDS}, {k: z[f"{k}|thr"] for k in KINDS}
    d = ev.load("tco95_ctl1950g")
    ws = [m for _, m in ev.winters(d, 1953, 2099)]
    C, TH = {}, {}
    for kind in KINDS:
        X = gpk.metric6(d, kind)
        thr = gpk.box_change(X, X, ws, ws, 8)["thr"]
        ev.WIN = 8
        c = np.full((len(ws), X.shape[1]), np.nan)
        for j in range(X.shape[1]):
            if np.isfinite(thr[j]):
                c[:, j] = [ev.counts_grid(X[m, j], np.array([thr[j]]))[0] for m in ws]
        C[kind], TH[kind] = c, thr
        print(f"  null counts {kind}: {len(ws)} winters, {np.isfinite(thr).sum()} boxes, "
              f"mean {np.nanmean(c):.2f} events/winter", flush=True)
    np.savez(NULL_CACHE, **{f"{k}|counts": C[k] for k in KINDS}, **{f"{k}|thr": TH[k] for k in KINDS})
    return C, TH


def null_logratios(c, na0, nb):
    W = len(c)
    na = na0 if na0 + nb <= W - 40 else max(W - nb - 40, nb)
    pairs = ev.block_pairs(W, na, nb)
    return ev.block_logratios(c.astype(float), na, nb, pairs), na


def fdr(p, alpha):
    """Benjamini-Hochberg on the finite entries of p; returns (mask, p threshold or nan)."""
    ok = np.isfinite(p)
    ps = np.sort(p[ok])
    n = len(ps)
    below = ps <= alpha * np.arange(1, n + 1) / n
    if not below.any():
        return np.zeros_like(ok), np.nan
    pk = ps[np.nonzero(below)[0].max()]
    return ok & (p <= pk), float(pk)


def gauss_p(obs, lr):
    mu, sd = np.nanmean(lr, axis=0), np.nanstd(lr, axis=0, ddof=1)
    z = np.abs(obs - mu) / sd
    return np.array([math.erfc(x / math.sqrt(2)) if np.isfinite(x) else np.nan for x in z])


def main():
    G = np.load(f"{OUT}/grid_per_K.npz")
    lat0 = G["lat0"]
    w = np.cos(np.deg2rad(lat0 + 2.5))
    C, TH = null_counts()
    arrays, summ = {}, dict(alpha_fdr=ALPHA_FDR, null="TCO95 ctl1950g 1953-2099, 6-hourly, block placements", ours={}, cmip6={})
    for kind in KINDS:
        for name, (na0, nb) in SIZES.items():
            key = f"6h|{name}|{kind}"
            if f"{key}|ratio" not in G.files:
                continue
            obs = np.log(G[f"{key}|ratio"])
            lr, na = null_logratios(C[kind], na0, nb)
            P = len(lr)
            fin = np.isfinite(lr).all(0)
            ok = np.isfinite(obs) & fin
            obs_ok = np.where(ok, obs, np.nan)
            n_up = (lr >= obs_ok[None]).sum(0)
            n_dn = (lr <= obs_ok[None]).sum(0)
            p_up = np.where(ok, (1 + n_up) / (1 + P), np.nan)
            p_dn = np.where(ok, (1 + n_dn) / (1 + P), np.nan)
            p = np.minimum(1.0, 2 * np.minimum(p_up, p_dn))
            sig, pk = fdr(p, ALPHA_FDR)
            pg = np.where(ok, gauss_p(obs_ok, lr), np.nan)
            sig_g, pk_g = fdr(pg, ALPHA_FDR)
            # consistency with grid_per_K.npz: the same 5-95 % chance range
            p05, p95 = np.nanpercentile(lr, [5, 95], axis=0)
            dev = np.nanmax(np.abs(np.r_[p05 - G[f"{key}|null_p05"], p95 - G[f"{key}|null_p95"]]))
            # area mean
            stat = float((obs[ok] * w[ok]).sum() / w[ok].sum())
            nullm = lr[:, ok] @ w[ok] / w[ok].sum()
            p_field = float((1 + (nullm >= stat).sum()) / (1 + P))
            dT = json.load(open(f"{OUT}/grid_per_K_summary.json"))["ours"][f"{name}|{kind}"]["dT"]
            perK = lambda x: float(100 * (np.exp(x / dT) - 1))        # noqa: E731
            up = sig & (obs > 0)
            dn = sig & (obs < 0)
            row = dict(n_control=na0, n_future=nb, n_control_null=na, n_placements=P, n_boxes=int(ok.sum()),
                       fdr_p_threshold=pk, n_sig=int(sig.sum()), n_sig_up=int(up.sum()), n_sig_down=int(dn.sum()),
                       frac_sig=float(sig.sum() / ok.sum()), frac_area_sig_up=float(w[up].sum() / w[ok].sum()),
                       frac_area_sig_down=float(w[dn].sum() / w[ok].sum()),
                       field_significant_wilks=bool(sig.any()),
                       n_sig_gauss=int(sig_g.sum()), n_sig_up_gauss=int((sig_g & (obs > 0)).sum()),
                       n_sig_down_gauss=int((sig_g & (obs < 0)).sum()),
                       n_p_lt_005=int((p[ok] < 0.05).sum()), expected_false_at_005=float(0.05 * ok.sum()),
                       area_mean_lnratio=stat, area_mean_ratio=float(np.exp(stat)), area_mean_perK_of_mean_ln=perK(stat),
                       null_area_mean_ratio_p05_p50_p95=[float(np.exp(x)) for x in np.percentile(nullm, [5, 50, 95])],
                       null_area_mean_perK_p95=perK(float(np.percentile(nullm, 95))),
                       p_field_area_mean=p_field, n_null_area_mean_ge_obs=int((nullm >= stat).sum()),
                       max_dev_from_grid_per_K_null_range=float(dev))
            summ["ours"][f"{name}|{kind}"] = row
            arrays[f"{key}|p"] = p
            arrays[f"{key}|p_gauss"] = pg
            arrays[f"{key}|fdr_sig"] = sig
            arrays[f"{key}|fdr_sig_gauss"] = sig_g
            arrays[f"{key}|null_area_mean"] = nullm
            print(f"{key:24s} ctl/fut {na0}/{nb} (null ctl {na}), {P} placements, {ok.sum()} boxes: "
                  f"FDR sig {sig.sum()} ({up.sum()} up, {dn.sum()} down; p_FDR {pk:.4f}); gauss {sig_g.sum()}; "
                  f"p<0.05 {(p[ok] < 0.05).sum()} (chance {0.05 * ok.sum():.0f}); area-mean ratio {np.exp(stat):.3f} "
                  f"({perK(stat):+.2f} %/K) vs null 95 % {np.exp(np.percentile(nullm, 95)):.3f}, p_field {p_field:.5f}; "
                  f"null-range check {dev:.1e}", flush=True)

    # CMIP6: every model's own area mean (cos-lat weighted mean of per-box per K) and its agreement map
    for kind in KINDS:
        names = sorted({k.split("|")[2] for k in G.files if k.startswith("1d|cmip6|") and k.endswith(f"|{kind}|perK")})
        am = {}
        for m in names:
            v = G[f"1d|cmip6|{m}|2071-2100|{kind}|perK"]
            ok = np.isfinite(v)
            am[m] = float((v[ok] * w[ok]).sum() / w[ok].sum())
        a = G[f"cmip6|{kind}|agree_up"]
        okc = np.isfinite(a)
        summ["cmip6"][kind] = dict(n_models=len(names), area_mean_perK_per_model=am,
                                   n_models_area_mean_up=int(sum(v > 0 for v in am.values())),
                                   frac_boxes_agree_up=float(np.mean(a[okc] >= 0.8)),
                                   frac_boxes_agree_down=float(np.mean(a[okc] <= 0.2)),
                                   frac_area_agree_up=float(w[okc & (a >= 0.8)].sum() / w[okc].sum()),
                                   frac_area_agree_down=float(w[okc & (a <= 0.2)].sum() / w[okc].sum()))
        print(f"CMIP6 {kind}: area mean up in {summ['cmip6'][kind]['n_models_area_mean_up']}/{len(names)} models; "
              f">= 80 % agree up in {summ['cmip6'][kind]['frac_boxes_agree_up']:.0%} of boxes "
              f"({summ['cmip6'][kind]['frac_area_agree_up']:.0%} of area), down {summ['cmip6'][kind]['frac_boxes_agree_down']:.0%}")
    np.savez(f"{OUT}/fieldsig_grid.npz", lat0=lat0, lon0=G["lon0"], **arrays)
    json.dump(summ, open(f"{OUT}/fieldsig_grid_summary.json", "w"), indent=1)
    print(f"Saved {OUT}/fieldsig_grid.npz, fieldsig_grid_summary.json")


if __name__ == "__main__":
    main()
