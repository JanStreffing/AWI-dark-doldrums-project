"""
What raises the surface-wind transfer s / G in the TCO95 (~100 km) rerun under warming?
Input: per region-day diagnostics from transfer_stability_extract.py (8 European regions, Sep-Mar).

  G          region mean daily geostrophic speed from the 2.5 deg daily msl (geostrophic.py functions)
  s          s100ll (log-law 100 m from 10 m, ERA5 z0; the series geostrophic.py uses), also s10, s100n
  transfer   mean(s) / mean(G) over the days of a period; change 100 ((sf/Gf)/(sc/Gc) - 1)

B1  baseline transfer change, s100ll / s10 / s100n
B2  stability classes per region-day, quintiles of the pooled control days of that region:
      ri      bulk Richardson number level 91 -> 85 (~10 -> ~305 m) from region-day means,
              g / T91 * dtheta * dz / |dV|^2
      ri_st   fraction of cell-steps with Ri(91 -> 85) > 0.25
      sshf    daily-mean sensible heat flux, W m-2, IFS sign (positive downward = stable)
      blh     daily-mean boundary-layer height
      t2skt   2t - skt
      inv     theta(L85) - theta(2 m) ~ T85 - T2 + g (h85 + 8 m) / cp
      G       G itself (does the transfer depend on G, so a weaker G raises s / G?)
      G x ri  5 x 5 joint classes
    and fixed physical classes (ri < 0, 0-0.25, >= 0.25; sshf < 0, >= 0).
    Transfer T = sum_k f_k s_k / sum_k f_k G_k (f class frequency, s_k, G_k class means).
      between  T(f_fut, s_ctl, G_ctl) / T_ctl          class frequencies change
      within   T_fut / T(f_fut, s_ctl, G_ctl)          per-class transfer changes
    (and the reverse order: within at control frequencies, then between). Regions averaged in %.
B3  per region and calendar month: transfer change against the change of snow cover, sd, fsr, sshf,
    2t - skt, inversion, Ri stable fraction; Pearson / Spearman across region-months.
B4  vertical chain s10/G = (V85/G)(V91/V85)(s10/V91), and s100n/s10.
B5  hi1950g 1951-1980 against 1991-2020 (and the control drift 1953-1980 against 2071-2099).
Uncertainty: block bootstrap over winters (control and future resampled separately), 1000 draws.

Usage: python3 transfer_stability.py
Output: dunkelflaute/results/regimes/transfer_stability.json, transfer_stability_*.png
"""
import glob
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("regions", "regimes"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/{p}")
import regions as rg                       # noqa: E402
import extract_msl as xm                   # noqa: E402
import geostrophic as gs                   # noqa: E402

EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
PARTS = f"{_ROOT}/dunkelflaute/data/regimes/transfer_stability_parts"
OUT = f"{_ROOT}/dunkelflaute/results/regimes"
G0, CP = 9.80665, 1004.7
NB = int(os.environ.get("NB", 1000))
COMPS = {
    "100km":       (("ctl1950g", 1953, 2099), ("hi1950g", 2081, 2099)),
    "hist_fut":    (("hi1950g", 1951, 1980), ("hi1950g", 2071, 2099)),
    "hist_trend":  (("hi1950g", 1951, 1980), ("hi1950g", 1991, 2020)),
    "ctl_drift":   (("ctl1950g", 1953, 1980), ("ctl1950g", 2071, 2099)),
}
LEV = [85, 86, 87, 88, 89, 90, 91]


def load(run):
    fs = sorted(glob.glob(f"{PARTS}/{run}_*.npz"))
    Z = [np.load(f) for f in fs]
    d = {}
    for k in Z[0].files:
        if k.startswith("m_"):
            continue
        d[k] = np.concatenate([z[k] for z in Z])
    o = np.argsort(d["day"])
    d = {k: v[o] for k, v in d.items() if k not in ("mslbox", "mslday")}
    ms = np.concatenate([z["mslbox"] for z in Z])
    md = np.concatenate([z["mslday"] for z in Z])
    om = np.argsort(md)
    ms, md = ms[om], md[om]
    assert np.array_equal(md, d["day"]), f"{run}: msl days differ"
    d["G"] = gs.region_G(gs.geostrophic(ms.astype(np.float64), xm.LATC, xm.LONC, False), xm.LATC, xm.LONC, EU)
    d["winter"] = gs.winter_of(gs.ymd_of(d["day"]))
    d["month"] = d["day"].astype("datetime64[M]").astype(int) % 12 + 1
    years = np.array([int(os.path.basename(f).split("_")[1][:4]) for f in fs])
    m = {k[2:]: np.stack([z[k] for z in Z]) for k in Z[0].files if k.startswith("m_")}   # (Y, 12, R)
    # stability indices (D, R)
    d["ri"] = G0 / d["T"][:, 6] * d["dth"][:, 0] * d["h"][:, 0] / d["sh2"][:, 0]
    d["ri_st"] = d["ri_st"][:, 0]
    d["t2skt"] = d["t2"] - d["skt"]
    d["inv"] = d["T"][:, 0] - d["t2"] + G0 * (d["h"][:, 0] + 8.0) / CP
    d["V85"], d["V91"] = d["V"][:, 0], d["V"][:, 6]
    # V40: wind at the 40 m blending height of the IFS 10 m diagnostic over land (z0m > 0.03 m),
    # linear in height between levels 90 and 89 (heights above ground: h + 10 m for level 91 at ~10 m)
    h90, h89 = d["h"][:, 5] + 10.0, d["h"][:, 4] + 10.0
    wgt = (40.0 - h90) / (h89 - h90)
    d["V40"] = d["V"][:, 5] + wgt * (d["V"][:, 4] - d["V"][:, 5])
    print(f"{run}: {len(fs)} years {years.min()}-{years.max()}, {len(d['day'])} days", flush=True)
    return d, years, m


# ------------------------------------------------------------------ class decomposition
def winter_sums(d, cls, K, s, g=None):
    """Per winter, region, class: sums of s and G and counts. cls (D, R) in 0..K-1."""
    ws = np.unique(d["winter"])
    g = d["G"] if g is None else g
    S = np.zeros((len(ws), len(EU), K))
    Gs, N = np.zeros_like(S), np.zeros_like(S)
    wi = np.searchsorted(ws, d["winter"])
    for r in range(len(EU)):
        np.add.at(S[:, r], (wi, cls[:, r]), s[:, r])
        np.add.at(Gs[:, r], (wi, cls[:, r]), g[:, r])
        np.add.at(N[:, r], (wi, cls[:, r]), 1.0)
    return ws, S, Gs, N


def decomp(Sc, Gc, Nc, Sf, Gf, Nf):
    """Totals (R, K) -> dict of per-region arrays (R,) and per-class stats."""
    fc, ff = Nc / Nc.sum(1, keepdims=True), Nf / Nf.sum(1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        sc, gc, sf, gf = Sc / Nc, Gc / Nc, Sf / Nf, Gf / Nf
    sc, gc, sf, gf = (np.nan_to_num(x) for x in (sc, gc, sf, gf))
    Tc = Sc.sum(1) / Gc.sum(1)
    Tf = Sf.sum(1) / Gf.sum(1)
    Tb = (ff * sc).sum(1) / (ff * gc).sum(1)          # future frequencies, control per-class means
    Tw = (fc * sf).sum(1) / (fc * gf).sum(1)          # control frequencies, future per-class means
    return dict(total=100 * (Tf / Tc - 1), between=100 * (Tb / Tc - 1), within=100 * (Tf / Tb - 1),
                within_ctlfreq=100 * (Tw / Tc - 1), between_futmeans=100 * (Tf / Tw - 1),
                freq_c=fc, freq_f=ff, tr_c=sc / np.where(gc > 0, gc, np.nan), tr_f=sf / np.where(gf > 0, gf, np.nan),
                G_c=gc, G_f=gf, s_c=sc, s_f=sf)


def run_decomp(dc, df, selc, self_, cls_c, cls_f, K, svar, rng, gvar="G"):
    wc, Sc, Gc, Nc = winter_sums({k: v[selc] for k, v in dc.items()}, cls_c[selc], K, dc[svar][selc], dc[gvar][selc])
    wf, Sf, Gf, Nf = winter_sums({k: v[self_] for k, v in df.items()}, cls_f[self_], K, df[svar][self_], df[gvar][self_])
    full = decomp(Sc.sum(0), Gc.sum(0), Nc.sum(0), Sf.sum(0), Gf.sum(0), Nf.sum(0))
    keys = ("total", "between", "within", "within_ctlfreq", "between_futmeans")
    boot = {k: [] for k in keys}
    for _ in range(NB):
        ic = rng.integers(0, len(wc), len(wc))
        jf = rng.integers(0, len(wf), len(wf))
        b = decomp(Sc[ic].sum(0), Gc[ic].sum(0), Nc[ic].sum(0), Sf[jf].sum(0), Gf[jf].sum(0), Nf[jf].sum(0))
        for k in keys:
            boot[k].append(np.nanmean(b[k]))
    res = {k: dict(eu_mean=float(np.nanmean(full[k])), ci90=[float(x) for x in np.nanpercentile(boot[k], [5, 95])],
                   per_region={r: float(full[k][i]) for i, r in enumerate(EU)}) for k in keys}
    res["classes"] = {k: [float(x) for x in np.nanmean(full[k], axis=0)] for k in
                      ("freq_c", "freq_f", "tr_c", "tr_f", "G_c", "G_f", "s_c", "s_f")}
    return res


def quintile_classes(xc, xf, selc, K=5):
    """Per region control-quantile edges (from the selected control days)."""
    cc, cf = np.zeros(xc.shape, int), np.zeros(xf.shape, int)
    edges = []
    for r in range(xc.shape[1]):
        e = np.nanpercentile(xc[selc, r], np.linspace(0, 100, K + 1)[1:-1])
        cc[:, r], cf[:, r] = np.digitize(xc[:, r], e), np.digitize(xf[:, r], e)
        edges.append(e)
    return cc, cf, np.array(edges)


def fixed_classes(x, edges):
    return np.digitize(x, edges)


# ------------------------------------------------------------------ helpers
def tr_change(dc, df, selc, self_, s):
    return np.array([100 * ((df[s][self_, r].mean() / df["G"][self_, r].mean()) /
                            (dc[s][selc, r].mean() / dc["G"][selc, r].mean()) - 1) for r in range(len(EU))])


def mean_change(dc, df, selc, self_, v, rel=True):
    a, b = dc[v][selc].mean(0), df[v][self_].mean(0)
    return 100 * (b / a - 1) if rel else b - a


def spearman(a, b):
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


def main():
    rng = np.random.default_rng(1)
    D = {r: load(r) for r in ("ctl1950g", "hi1950g")}
    res = dict(meta=dict(regions=EU, levels=LEV, n_boot=NB,
                         note="s100ll = log-law 100 m from 10 m speed with fixed ERA5 z0 (geostrophic.py series); "
                              "sshf IFS sign positive downward; Ri from region-day means level 91 -> 85"))
    figs = {}
    for name, ((cn, c0, c1), (fn, f0, f1)) in COMPS.items():
        dc, yc, mc = D[cn]
        df, yf, mf = D[fn]
        selc = (dc["winter"] >= c0) & (dc["winter"] <= c1)
        self_ = (df["winter"] >= f0) & (df["winter"] <= f1)
        out = dict(ctl=[cn, c0, c1], fut=[fn, f0, f1], n_winters=[int(len(np.unique(dc["winter"][selc]))), int(len(np.unique(df["winter"][self_])))])
        # B1
        b1 = {}
        for s in ("s100ll", "s10", "s100n", "V91", "V40", "V85"):
            x = tr_change(dc, df, selc, self_, s)
            b1[s] = dict(eu_mean=float(x.mean()), per_region={r: float(v) for r, v in zip(EU, x)})
        b1["mean_G_change_pct"] = float(mean_change(dc, df, selc, self_, "G").mean())
        b1["mean_s10_change_pct"] = float(mean_change(dc, df, selc, self_, "s10").mean())
        out["B1"] = b1
        # index means
        out["index_means"] = {v: dict(ctl=float(np.nanmean(dc[v][selc])), fut=float(np.nanmean(df[v][self_])))
                              for v in ("ri", "ri_st", "sshf", "sshf_pos", "blh", "t2skt", "inv", "t2", "skt", "G", "s10", "V91", "V40", "V85", "s100n")}
        out["level_heights_m"] = [float(x) for x in dc["h"][selc].mean((0, 2))] + [0.0]
        # B2
        b2 = {}
        for v in ("ri", "ri_st", "sshf", "blh", "t2skt", "inv", "G"):
            cc, cf, e = quintile_classes(dc[v], df[v], selc)
            b2[f"{v}_q5"] = run_decomp(dc, df, selc, self_, cc, cf, 5, "s100ll", rng)
            b2[f"{v}_q5"]["edges_eu_mean"] = [float(x) for x in e.mean(0)]
        cc, cf, _ = quintile_classes(dc["G"], df["G"], selc)
        rc, rf, _ = quintile_classes(dc["ri"], df["ri"], selc)
        b2["G_x_ri_q5x5"] = run_decomp(dc, df, selc, self_, cc * 5 + rc, cf * 5 + rf, 25, "s100ll", rng)
        # stratification at fixed G (ri is circular at fixed G: its denominator is the shear, ~ V85)
        for v in ("sshf", "inv", "t2skt", "blh"):
            rc, rf, _ = quintile_classes(dc[v], df[v], selc)
            b2[f"G_x_{v}_q5x5"] = run_decomp(dc, df, selc, self_, cc * 5 + rc, cf * 5 + rf, 25, "s100ll", rng)
        for v, e in (("ri", [0.0, 0.25]), ("sshf", [0.0])):
            b2[f"{v}_fixed"] = run_decomp(dc, df, selc, self_, fixed_classes(dc[v], e), fixed_classes(df[v], e), len(e) + 1, "s100ll", rng)
            b2[f"{v}_fixed"]["edges"] = e
        # chain links split by stability class: 10 m diagnostic (s10 / V40), resolved profile (V40 / V85),
        # boundary-layer top (V85 / G)
        for v in ("sshf", "inv", "blh", "G"):
            cc, cf, _ = quintile_classes(dc[v], df[v], selc)
            for num, den in (("s10", "V40"), ("V40", "V85"), ("V85", "G")):
                b2[f"{num}/{den}|{v}_q5"] = run_decomp(dc, df, selc, self_, cc, cf, 5, num, rng, gvar=den)
        # chain links split by joint G x stratification classes
        for v in ("inv", "blh"):
            rc, rf, _ = quintile_classes(dc[v], df[v], selc)
            for num, den in (("s10", "V40"), ("V40", "V85"), ("V85", "G")):
                b2[f"{num}/{den}|G_x_{v}_q5x5"] = run_decomp(dc, df, selc, self_, cc * 5 + rc, cf * 5 + rf, 25, num, rng, gvar=den)
        # stratification at fixed G and fixed s10 (25 joint quintile classes, control weights): does the
        # boundary layer destabilise beyond what the changed wind itself does?
        sc_, sf_, _ = quintile_classes(dc["s10"], df["s10"], selc)
        kc, kf = cc * 5 + sc_, cf * 5 + sf_
        fixed = {}
        for v in ("inv", "t2skt", "blh", "sshf", "ri_st"):
            per = []
            for i in range(len(EU)):
                a, b = kc[selc, i], kf[self_, i]
                xa, xb = dc[v][selc, i], df[v][self_, i]
                na = np.bincount(a, minlength=25).astype(float)
                ma = np.bincount(a, xa, 25) / np.maximum(na, 1)
                nb = np.bincount(b, minlength=25).astype(float)
                mb = np.bincount(b, xb, 25) / np.maximum(nb, 1)
                ok = (na >= 20) & (nb >= 20)
                per.append(float((na[ok] * (mb[ok] - ma[ok])).sum() / na[ok].sum()))
            fixed[v] = dict(change_fixed_G_s10=float(np.mean(per)),
                            change_raw=float(np.mean(df[v][self_].mean(0) - dc[v][selc].mean(0))),
                            per_region={r: x for r, x in zip(EU, per)})
        out["stratification_at_fixed_G_s10"] = fixed
        mc_ = np.where(dc["month"] >= 9, dc["month"] - 9, dc["month"] + 3)
        mf_ = np.where(df["month"] >= 9, df["month"] - 9, df["month"] + 3)
        b2["month"] = run_decomp(dc, df, selc, self_, np.repeat(mc_[:, None], 8, 1), np.repeat(mf_[:, None], 8, 1), 7, "s100ll", rng)
        out["B2"] = b2
        # B3 per region-month
        rows = []
        ycs = (yc >= c0 - 1) & (yc <= c1)
        yfs = (yf >= f0 - 1) & (yf <= f1)
        for i, r in enumerate(EU):
            for mo in rg.MONTHS:
                a = selc & (dc["month"] == mo)
                b = self_ & (df["month"] == mo)
                tr = 100 * ((df["s100ll"][b, i].mean() / df["G"][b, i].mean()) / (dc["s100ll"][a, i].mean() / dc["G"][a, i].mean()) - 1)
                # monthly files: winter months of the selected winters (Sep-Dec of y0-1 .. y1-1, Jan-Mar of y0 .. y1)
                def mmean(m, ys, y0, y1, key):
                    yy = ys[(ys >= (y0 - 1 if mo >= 9 else y0)) & (ys <= (y1 - 1 if mo >= 9 else y1))]
                    idx = np.isin(ys, yy)
                    return float(m[key][idx, mo - 1, i].mean())
                row = dict(region=r, month=mo, transfer_change_pct=float(tr),
                           dG_pct=float(100 * (df["G"][b, i].mean() / dc["G"][a, i].mean() - 1)),
                           ds10_pct=float(100 * (df["s10"][b, i].mean() / dc["s10"][a, i].mean() - 1)))
                for key in ("scf_depth", "scf_swe", "sd", "fsr", "sshf", "skt", "2t", "asn"):
                    x0, x1 = mmean(mc, yc, c0, c1, key), mmean(mf, yf, f0, f1, key)
                    row[f"{key}_ctl"], row[f"{key}_fut"] = x0, x1
                    row[f"d_{key}"] = x1 - x0
                row["d_sshf"] /= 10800.0
                row["sshf_ctl"] /= 10800.0
                row["sshf_fut"] /= 10800.0
                for v in ("t2skt", "inv", "ri_st", "blh", "sshf"):
                    row[f"dd_{v}"] = float(np.nanmean(df[v][b, i]) - np.nanmean(dc[v][a, i]))
                rows.append(row)
        tr = np.array([x["transfer_change_pct"] for x in rows])
        cor = {}
        for key in ("d_scf_depth", "d_scf_swe", "d_sd", "d_fsr", "d_sshf", "dd_t2skt", "dd_inv", "dd_ri_st", "dd_blh", "dG_pct", "d_2t"):
            x = np.array([r_[key] for r_ in rows])
            cor[key] = dict(pearson=float(np.corrcoef(x, tr)[0, 1]), spearman=spearman(x, tr),
                            mean=float(x.mean()), slope=float(np.polyfit(x, tr, 1)[0]) if x.std() > 0 else None)
        snowy = np.array([x["scf_depth_ctl"] > 0.1 for x in rows])
        X = np.column_stack([np.ones(len(rows)), [x["dG_pct"] for x in rows], [x["d_scf_depth"] for x in rows]])
        coef, *_ = np.linalg.lstsq(X, tr, rcond=None)
        resid_G = tr - X[:, :2] @ np.linalg.lstsq(X[:, :2], tr, rcond=None)[0]
        resid_S = X[:, 2] - X[:, :2] @ np.linalg.lstsq(X[:, :2], X[:, 2], rcond=None)[0]
        cor["regression_tr_on_dG_dscf"] = dict(intercept=float(coef[0]), per_pct_dG=float(coef[1]), per_unit_dscf=float(coef[2]),
                                               partial_r_dscf_given_dG=float(np.corrcoef(resid_G, resid_S)[0, 1]),
                                               r2=float(1 - ((tr - X @ coef) ** 2).sum() / ((tr - tr.mean()) ** 2).sum()))
        out["B3"] = dict(rows=rows, corr_with_transfer_change=cor,
                         transfer_change_snowy_ctl_scf_gt_0p1=float(tr[snowy].mean()) if snowy.any() else None,
                         transfer_change_not_snowy=float(tr[~snowy].mean()), n_snowy=int(snowy.sum()), n=len(rows),
                         fsr_changes_relative_pct=float(np.mean([100 * x["d_fsr"] / x["fsr_ctl"] for x in rows])))
        # B4
        def rat(dd, sel, a, b):
            return dd[a][sel].mean(0) / dd[b][sel].mean(0)
        b4 = {}
        for a, b in (("V85", "G"), ("V91", "V85"), ("s10", "V91"), ("s10", "G"), ("s100n", "s10"), ("s100ll", "s10"), ("V85", "V91"), ("s100n", "V91"),
                     ("V40", "V85"), ("s10", "V40"), ("s100n", "G"), ("V40", "G"), ("s100n", "V85")):
            rc_, rf_ = rat(dc, selc, a, b), rat(df, self_, a, b)
            x = 100 * (rf_ / rc_ - 1)
            b4[f"{a}/{b}"] = dict(ctl=float(rc_.mean()), fut=float(rf_.mean()), change_pct=float(x.mean()),
                                   per_region={r: float(v) for r, v in zip(EU, x)})
        for k in range(6):
            a_, b_ = dc["V"][selc][:, k].mean(0) / dc["V"][selc][:, 6].mean(0), df["V"][self_][:, k].mean(0) / df["V"][self_][:, 6].mean(0)
            b4[f"V{LEV[k]}/V91"] = dict(ctl=float(a_.mean()), fut=float(b_.mean()), change_pct=float((100 * (b_ / a_ - 1)).mean()))
        out["B4"] = b4
        res[name] = out
        figs[name] = out
        print(f"\n=== {name}: {cn} {c0}-{c1} vs {fn} {f0}-{f1}  winters {out['n_winters']}")
        print("  B1 transfer change: " + ", ".join(f"{k} {v['eu_mean']:+.2f}%" for k, v in b1.items() if isinstance(v, dict)) +
              f"; dG {b1['mean_G_change_pct']:+.2f}%, ds10 {b1['mean_s10_change_pct']:+.2f}%")
        print("  index means ctl -> fut: " + "; ".join(f"{k} {v['ctl']:.3g}->{v['fut']:.3g}" for k, v in out["index_means"].items()))
        for k, v in b2.items():
            print(f"  B2 {k:14s} total {v['total']['eu_mean']:+.2f}  between {v['between']['eu_mean']:+.2f} [{v['between']['ci90'][0]:+.2f},{v['between']['ci90'][1]:+.2f}]"
                  f"  within {v['within']['eu_mean']:+.2f} [{v['within']['ci90'][0]:+.2f},{v['within']['ci90'][1]:+.2f}]"
                  f"  | within@ctlfreq {v['within_ctlfreq']['eu_mean']:+.2f} between@futmeans {v['between_futmeans']['eu_mean']:+.2f}")
            if k.endswith("q5") or k.endswith("fixed"):
                c = v["classes"]
                print("       freq c/f " + " ".join(f"{a:.3f}/{b:.3f}" for a, b in zip(c["freq_c"], c["freq_f"])) +
                      " | tr c/f " + " ".join(f"{a:.3f}/{b:.3f}" for a, b in zip(c["tr_c"], c["tr_f"])))
        print("  B3 corr with transfer change: " + "; ".join(f"{k} r={v['pearson']:+.2f} rho={v['spearman']:+.2f} mean={v['mean']:.3g}" for k, v in cor.items() if "pearson" in v))
        print("  B3 regression: " + json.dumps(cor["regression_tr_on_dG_dscf"]))
        print("  stratification change at fixed G and s10 (raw): " + "; ".join(f"{k} {v['change_fixed_G_s10']:+.3f} ({v['change_raw']:+.3f})" for k, v in fixed.items()))
        print(f"  B3 transfer change snowy (ctl scf>0.1, n={out['B3']['n_snowy']}) {out['B3']['transfer_change_snowy_ctl_scf_gt_0p1']}, other {out['B3']['transfer_change_not_snowy']:+.2f}, fsr rel change {out['B3']['fsr_changes_relative_pct']:+.3f}%")
        print("  B4 " + "; ".join(f"{k} {v['ctl']:.3f}->{v['fut']:.3f} ({v['change_pct']:+.2f}%)" for k, v in b4.items()))
    json.dump(res, open(f"{OUT}/transfer_stability.json", "w"), indent=1)
    print(f"Saved {OUT}/transfer_stability.json")
    plot(res)


def plot(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    for name in ("100km", "hist_fut"):
        r = res[name]["B2"]
        keys = ["ri_q5", "sshf_q5", "blh_q5", "inv_q5", "t2skt_q5", "G_q5"]
        fig, ax = plt.subplots(2, len(keys), figsize=(3.0 * len(keys), 5.5), constrained_layout=True)
        for j, k in enumerate(keys):
            c = r[k]["classes"]
            x = np.arange(1, 6)
            ax[0, j].plot(x, c["tr_c"], "o-", color="#3b6ea8", label="control")
            ax[0, j].plot(x, c["tr_f"], "o-", color="#c0504d", label="future")
            ax[0, j].set_title(f"{k}\nbetween {r[k]['between']['eu_mean']:+.1f}%, within {r[k]['within']['eu_mean']:+.1f}%", fontsize=9)
            ax[1, j].bar(x - 0.2, c["freq_c"], 0.4, color="#3b6ea8")
            ax[1, j].bar(x + 0.2, c["freq_f"], 0.4, color="#c0504d")
            ax[1, j].set_xlabel("control quintile")
        ax[0, 0].set_ylabel("transfer s/G")
        ax[1, 0].set_ylabel("frequency")
        ax[0, 0].legend(fontsize=8)
        fig.suptitle(f"TCO95 {name}: transfer per stability class (8 EU regions)")
        fig.savefig(f"{OUT}/transfer_stability_classes_{name}.png", dpi=110)
        plt.close(fig)
        rows = res[name]["B3"]["rows"]
        fig, ax = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
        tr = [x["transfer_change_pct"] for x in rows]
        for a, key, lab in ((ax[0], "d_scf_depth", "change in snow cover fraction"), (ax[1], "dd_inv", "change in theta(L85)-theta(2m) (K)"),
                            (ax[2], "dd_ri_st", "change in Ri>0.25 fraction")):
            x = [q[key] for q in rows]
            sc = a.scatter(x, tr, c=[q["month"] if q["month"] >= 9 else q["month"] + 12 for q in rows], cmap="viridis", s=18)
            a.set_xlabel(lab)
            a.set_title(f"r = {res[name]['B3']['corr_with_transfer_change'][key]['pearson']:+.2f}")
            a.axhline(0, color="k", lw=0.5)
        ax[0].set_ylabel("transfer change (%)")
        fig.colorbar(sc, ax=ax[2], label="month (13=Jan)")
        fig.savefig(f"{OUT}/transfer_stability_regionmonth_{name}.png", dpi=110)
        plt.close(fig)


if __name__ == "__main__":
    main()
