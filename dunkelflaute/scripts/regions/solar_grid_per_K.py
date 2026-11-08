"""
Gridded per-K maps of winter low-SOLAR events on the 5 x 5 deg land boxes over 20-70 N (methods 15),
for our runs (6-hourly) and 20 CMIP6 models (daily), with the field significance of methods 15a.

Same machinery as grid_per_K.py and fieldsig_grid.py; only the metric differs:
  solar      solar CF (absolute; the series the combined metric uses: 9 km v5, 31 km v7, 100 km and
             IFS-FESOM2-SR native ssrd; CMIP6 daily-mean rsds / 1000)
  solar_ds   solar CF / control calendar-day climatology (solar_common.py: 31-day running mean, floored
             at 5 W m-2). A box whose climatology falls below the floor on any day (polar night) is not
             used for solar_ds.
Per box: threshold matched to 5 events per control winter on 200 levels between the 0.1st and 40th
percentiles of the control's 48 h means; valid if the control reaches 0.9 x 5. Ratio, per K, mean-removed
ratio, mean change of the CF, and the share of control event time in Nov-Jan.
Significance, our runs: block null from the TCO95 constant-1950 control (147 winters) with each
comparison's winter counts (control block cut to W - n_b - 40 where needed), the same placement in
every box; two-sided empirical p per box; Benjamini-Hochberg FDR at alpha 0.10 (Wilks 2016); area mean of
ln ratio against the same area mean of every placement.
CMIP6: median per-K map (boxes valid in >= 80 % of models), share of models with ln ratio > 0, every
model's area mean; pattern correlation (cos-lat weighted) of our runs with the CMIP6 median, and of each
model with the median of the others.

Usage: sbatch --job-name=solar_grid dunkelflaute/jobs/regions_py.sbatch regions/solar_grid_per_K.py
Output: dunkelflaute/results/regions/solar_grid_per_K.npz, solar_grid_per_K_summary.json
"""
import glob
import json
import os
import sys
from multiprocessing import Pool

os.environ["REGION_SET"] = "grid5"
import numpy as np                                                     # noqa: E402

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("regions", "analysis", "cmip6"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/{p}")
import regions as rg                                                   # noqa: E402
import region_events as ev                                             # noqa: E402
import region_per_K as rpk                                             # noqa: E402
import solar_common as sc                                              # noqa: E402
from fieldsig_grid import fdr                                          # noqa: E402

OUT = f"{_ROOT}/dunkelflaute/results/regions"
DC = f"{_ROOT}/dunkelflaute/data/cmip6_grid5"
TARGET = sc.TARGET
NGRID = 200
ALPHA_FDR = 0.10
KINDS = ("solar", "solar_ds")
B = [r for r in rg.REGIONS if r != "DE"]
BI = [rg.REGIONS.index(r) for r in B]
OURS = ["9km", "31km", "100km", "eerie_ctl", "eerie_hist"]


def load6(name):
    z = np.load(f"{rg.OUT_DIR}/{name}.npz")
    assert list(z["regions"]) == rg.REGIONS
    d = dict(time=z["time"].astype("datetime64[s]"), solar=z["solar"][:, BI].astype(float))
    if json.loads(str(z["meta"]))["dt_h"] == 3:
        d = ev.to_6h(d)
    d["winter"] = ev.winter_of(d["time"])
    d["di"], _ = sc.day_index(d["time"])
    t1 = d["time"] - np.timedelta64(1, "s")
    d["month"] = (t1.astype("datetime64[M]").astype(int) % 12) + 1
    return d


def wmasks6(d, y0, y1):
    return [d["winter"] == y for y in range(y0, y1 + 1) if (d["winter"] == y).sum() >= ev.MIN_STEPS]


def load1(fn, exp):
    z = np.load(fn)
    d = {k.split("|", 1)[1]: z[k] for k in z.files if k.startswith(exp + "|")}
    d["solar"] = d["solar"][:, BI].astype(float)
    y, m = d["ymd"][:, 0], d["ymd"][:, 1]
    d["winter"] = np.where(m >= 9, y + 1, y)
    d["di"] = sc.md_index(m, d["ymd"][:, 2])
    d["month"] = m
    return d


def wmasks1(d, y0, y1):
    return [d["winter"] == y for y in range(y0, y1 + 1) if (d["winter"] == y).sum() >= 200]


def box_change(Xc, Xf, mc, mf, win, months_c):
    """As grid_per_K.box_change for solar series, plus the Nov-Jan share of control event time."""
    ev.WIN = win
    nB = Xc.shape[1]
    out = {k: np.full(nB, np.nan) for k in ("thr", "rate", "ctl", "fut", "ratio", "ratio_mr", "mean_pct", "ndj")}
    mo = np.concatenate([months_c[m] for m in mc])
    for j in range(nB):
        xc, xf = Xc[:, j], Xf[:, j]
        if not np.isfinite(xc).all() or not np.isfinite(xf).all():
            continue
        roll = np.concatenate([sc.trailing(xc[m], win)[win - 1:] for m in mc])
        grid = np.linspace(*np.percentile(roll, [0.1, 40]), NGRID)
        cc = np.stack([ev.counts_grid(xc[m], grid) for m in mc]).mean(0)
        k = int(np.argmax(cc >= TARGET)) if (cc >= TARGET).any() else int(np.argmax(cc))
        thr = np.array([grid[k]])
        a = np.array([ev.counts_grid(xc[m], thr)[0] for m in mc], float)
        b = np.array([ev.counts_grid(xf[m], thr)[0] for m in mf], float)
        mcm = np.mean(np.concatenate([xc[m] for m in mc]))
        mfm = np.mean(np.concatenate([xf[m] for m in mf]))
        bm = np.array([ev.counts_grid(xf[m] * mcm / mfm, thr)[0] for m in mf], float)
        f = np.concatenate([sc.flagged(xc[m], grid[k], win) for m in mc])
        out["thr"][j], out["rate"][j] = grid[k], cc[k]
        out["ctl"][j], out["fut"][j] = a.mean(), b.mean()
        if a.mean() > 0:
            out["ratio"][j] = b.mean() / a.mean()
            out["ratio_mr"][j] = bm.mean() / a.mean()
        out["mean_pct"][j] = 100 * (mfm / mcm - 1)
        out["ndj"][j] = f[np.isin(mo, (11, 12, 1))].sum() / max(f.sum(), 1)
    out["valid"] = out["rate"] >= 0.9 * TARGET
    return out


def both_kinds(dc, df, mc, mf, win):
    clim = sc.climatology(dc["solar"], dc["di"], np.any(np.stack(mc), axis=0))
    ok_ds = np.nanmin(clim, axis=0) >= sc.CLIM_FLOOR
    res = {"solar": box_change(dc["solar"], df["solar"], mc, mf, win, dc["month"])}
    Xc, Xf = sc.deseason(dc["solar"], dc["di"], clim), sc.deseason(df["solar"], df["di"], clim)
    Xc[:, ~ok_ds] = np.nan
    Xf[:, ~ok_ds] = np.nan
    res["solar_ds"] = box_change(Xc, Xf, mc, mf, win, dc["month"])
    res["solar_ds"]["mean_pct"] = res["solar"]["mean_pct"].copy()          # the CF mean change, same in both
    res["clim_min"] = np.nanmin(clim, axis=0)
    return res


def run_ours(name):
    (cn, c0, c1), (fn, f0, f1) = ev.COMPARISONS[name]
    T = rpk.global_t()
    dc, df = load6(cn), load6(fn)
    mc, mf = wmasks6(dc, c0, c1), wmasks6(df, f0, f1)
    dT = rpk.mean_years(T, fn, f0, f1)[0] - rpk.mean_years(T, cn, c0, c1)[0]
    r = both_kinds(dc, df, mc, mf, 8)
    r.update(dT=dT, n_winters=(len(mc), len(mf)))
    print(f"  {name}: dT {dT:.2f}, winters {len(mc)}/{len(mf)}", flush=True)
    return name, r


def run_cmip6(fn):
    m = os.path.basename(fn)[len("daily_cmip6_"):-4]
    dc, df = load1(fn, "historical"), load1(fn, "ssp585")
    mc, mf = wmasks1(dc, 1951, 1980), wmasks1(df, 2071, 2100)
    Tc = dict(zip(dc["tas_years"].tolist(), dc["tas"].tolist()))
    Tf = dict(zip(df["tas_years"].tolist(), df["tas"].tolist()))
    dT = float(np.mean([Tf[y] for y in range(2071, 2101) if y in Tf]) - np.mean([Tc[y] for y in range(1951, 1981) if y in Tc]))
    r = both_kinds(dc, df, mc, mf, 2)
    r.update(dT=dT, n_winters=(len(mc), len(mf)))
    print(f"  cmip6 {m}: dT {dT:.2f}", flush=True)
    return m, r


def null_counts():
    """TCO95 control: per-winter counts (W, nB) per kind at the control-matched threshold."""
    d = load6("tco95_ctl1950g")
    ms = wmasks6(d, 1953, 2099)
    r = both_kinds(d, d, ms, ms, 8)
    ev.WIN = 8
    clim = sc.climatology(d["solar"], d["di"], np.any(np.stack(ms), axis=0))
    X = {"solar": d["solar"], "solar_ds": sc.deseason(d["solar"], d["di"], clim)}
    out = {}
    for kind in KINDS:
        thr = r[kind]["thr"]
        C = np.full((len(ms), len(B)), np.nan)
        for j in range(len(B)):
            if np.isfinite(thr[j]):
                C[:, j] = [ev.counts_grid(X[kind][m, j], np.array([thr[j]]))[0] for m in ms]
        out[kind] = C
        out[f"{kind}_valid"] = r[kind]["valid"]
    print(f"  null: {len(ms)} winters", flush=True)
    return out


def per_K(r, dT):
    return 100 * (np.exp(np.log(r) / dT) - 1)


def wmean(x, w):
    ok = np.isfinite(x)
    return float((x[ok] * w[ok]).sum() / w[ok].sum()) if ok.any() else np.nan


def wcorr(a, b, w):
    ok = np.isfinite(a) & np.isfinite(b)
    a, b, w = a[ok], b[ok], w[ok]
    am, bm = (a * w).sum() / w.sum(), (b * w).sum() / w.sum()
    return float(((a - am) * (b - bm) * w).sum() / np.sqrt(((a - am) ** 2 * w).sum() * ((b - bm) ** 2 * w).sum()))


def main():
    os.makedirs(OUT, exist_ok=True)
    lat0 = np.array([rg.BOXES[b][0] for b in B])
    lon0 = np.array([rg.BOXES[b][2] for b in B])
    w = np.cos(np.deg2rad(lat0 + 2.5))
    files = sorted(f for f in glob.glob(f"{DC}/daily_cmip6_*.npz") if "EC-Earth3" not in f)
    with Pool(12) as pool:
        jobs_c = pool.map_async(run_cmip6, files)
        jobs_o = pool.map_async(run_ours, OURS)
        jobs_n = pool.apply_async(null_counts)
        RC, RO, N = dict(jobs_c.get()), dict(jobs_o.get()), jobs_n.get()

    arrays = dict(boxes=np.array(B), lat0=lat0, lon0=lon0)
    summ = dict(n_boxes=len(B), target=TARGET, alpha_fdr=ALPHA_FDR, ours={}, cmip6_models={}, cmip6={})
    for src, RR in (("6h", RO), ("1d", RC)):
        for name, v in RR.items():
            arrays[f"{src}|{name}|clim_min"] = v["clim_min"]
            for kind in KINDS:
                o = v[kind]
                ok = o["valid"] & np.isfinite(o["ratio"]) & (o["ratio"] > 0)
                pk = np.where(ok, per_K(np.where(ok, o["ratio"], 1.0), v["dT"]), np.nan)
                okm = ok & (o["ratio_mr"] > 0)
                pk_mr = np.where(okm, per_K(np.where(okm, o["ratio_mr"], 1.0), v["dT"]), np.nan)
                key = f"{src}|{name}|{kind}"
                for a in ("ratio", "ratio_mr", "ctl", "fut", "ndj"):
                    arrays[f"{key}|{a}"] = np.where(ok, o[a], np.nan)
                arrays[f"{key}|mean_pct_perK"] = o["mean_pct"] / v["dT"]
                arrays[f"{key}|perK"] = pk
                arrays[f"{key}|perK_mr"] = pk_mr
                row = dict(dT=v["dT"], n_winters=list(v["n_winters"]), n_valid=int(ok.sum()), frac_up=float(np.mean(o["ratio"][ok] > 1)),
                           area_mean_perK=wmean(pk, w), area_mean_perK_meanremoved=wmean(pk_mr, w),
                           area_mean_lnratio=wmean(np.log(np.where(ok, o["ratio"], np.nan)), w),
                           area_mean_mean_change_perK_all=wmean(o["mean_pct"] / v["dT"], w),
                           area_mean_mean_change_perK_valid=wmean(np.where(ok, o["mean_pct"], np.nan) / v["dT"], w),
                           ndj_share_area_mean=wmean(np.where(ok, o["ndj"], np.nan), w))
                if src == "6h":
                    # field significance against the TCO95 control
                    na0, nb = v["n_winters"]
                    W = len(N[kind])
                    na = na0 if na0 + nb <= W - 40 else max(W - nb - 40, nb)
                    pairs = ev.block_pairs(W, na, nb)
                    lr = ev.block_logratios(N[kind], na, nb, pairs)
                    obs = np.log(np.where(ok, o["ratio"], np.nan))
                    fin = np.isfinite(lr).all(0) & N[f"{kind}_valid"]
                    okk = np.isfinite(obs) & fin
                    obs_ok = np.where(okk, obs, np.nan)
                    P = len(lr)
                    p_up = np.where(okk, (1 + (lr >= obs_ok[None]).sum(0)) / (1 + P), np.nan)
                    p_dn = np.where(okk, (1 + (lr <= obs_ok[None]).sum(0)) / (1 + P), np.nan)
                    p = np.minimum(1.0, 2 * np.minimum(p_up, p_dn))
                    sig, pthr = fdr(p, ALPHA_FDR)
                    stat = float((obs[okk] * w[okk]).sum() / w[okk].sum())
                    nullm = lr[:, okk] @ w[okk] / w[okk].sum()
                    arrays[f"{key}|fdr_sig"] = sig
                    arrays[f"{key}|p"] = p
                    p05, p95 = np.nanpercentile(lr, [5, 95], axis=0)
                    row.update(n_tested=int(okk.sum()), n_placements=P, n_control_null=na, fdr_p_threshold=pthr,
                               n_sig=int(sig.sum()), n_sig_up=int((sig & (obs > 0)).sum()), n_sig_down=int((sig & (obs < 0)).sum()),
                               n_p_lt_005=int((p[okk] < 0.05).sum()), expected_false_at_005=float(0.05 * okk.sum()),
                               frac_above_null95=float(np.mean(obs[okk] > p95[okk])), frac_below_null05=float(np.mean(obs[okk] < p05[okk])),
                               area_mean_ratio_tested=float(np.exp(stat)),
                               null_area_mean_ratio_p05_p50_p95=[float(np.exp(x)) for x in np.percentile(nullm, [5, 50, 95])],
                               p_field_up=float((1 + (nullm >= stat).sum()) / (1 + P)),
                               p_field_down=float((1 + (nullm <= stat).sum()) / (1 + P)))
                    summ["ours"][f"{name}|{kind}"] = row
                else:
                    summ["cmip6_models"][f"{name}|{kind}"] = row

    names = sorted(RC)
    for kind in KINDS:
        P = np.stack([arrays[f"1d|{n}|{kind}|perK"] for n in names])
        L = np.stack([np.log(arrays[f"1d|{n}|{kind}|ratio"]) for n in names])
        M = np.stack([arrays[f"1d|{n}|{kind}|mean_pct_perK"] for n in names])
        nval = np.isfinite(P).sum(0)
        enough = nval >= 0.8 * len(names)
        med = np.where(enough, np.nanmedian(P, axis=0), np.nan)
        agree = np.where(enough, np.nanmean(L > 0, axis=0), np.nan)
        arrays[f"cmip6|{kind}|perK_median"] = med
        arrays[f"cmip6|{kind}|agree_up"] = agree
        arrays[f"cmip6|{kind}|perK_mr_median"] = np.where(enough, np.nanmedian(np.stack([arrays[f"1d|{n}|{kind}|perK_mr"] for n in names]), axis=0), np.nan)
        arrays[f"cmip6|{kind}|mean_pct_perK_median"] = np.nanmedian(M, axis=0)
        am = {n: summ["cmip6_models"][f"{n}|{kind}"]["area_mean_perK"] for n in names}
        row = dict(n_models=len(names), n_boxes=int(enough.sum()), area_mean_perK_median_map=wmean(med, w),
                   area_mean_perK_models_median=float(np.median(list(am.values()))),
                   area_mean_perK_models_p10_p90=[float(x) for x in np.percentile(list(am.values()), [10, 90])],
                   n_models_area_mean_up=int(sum(v > 0 for v in am.values())),
                   frac_boxes_median_up=float(np.nanmean(med[enough] > 0)),
                   frac_boxes_80pct_agree_up=float(np.nanmean(agree[enough] >= 0.8)),
                   frac_boxes_80pct_agree_down=float(np.nanmean(agree[enough] <= 0.2)),
                   area_mean_perK_meanremoved_median_map=wmean(arrays[f"cmip6|{kind}|perK_mr_median"], w),
                   area_mean_mean_change_perK_median_map=wmean(arrays[f"cmip6|{kind}|mean_pct_perK_median"], w),
                   ndj_share_models_median=float(np.median([summ["cmip6_models"][f"{n}|{kind}"]["ndj_share_area_mean"] for n in names])))
        pc = {o: wcorr(arrays[f"6h|{o}|{kind}|perK"], med, w) for o in OURS}
        loo = [wcorr(P[i], np.nanmedian(np.delete(P, i, axis=0), axis=0), w) for i in range(len(names))]
        row["pattern_corr_with_median"] = pc
        row["pattern_corr_single_model_vs_rest_median"] = [float(np.median(loo)), float(np.min(loo)), float(np.max(loo))]
        # how much of the per-box event change is the mean change: correlation across boxes in the median maps
        row["corr_boxes_perK_vs_meanchange_median_maps"] = wcorr(med, arrays[f"cmip6|{kind}|mean_pct_perK_median"], w)
        row["corr_boxes_perK_vs_meanchange_ours"] = {o: wcorr(arrays[f"6h|{o}|{kind}|perK"], arrays[f"6h|{o}|{kind}|mean_pct_perK"], w) for o in OURS}
        zon = {}
        for la in np.unique(lat0):
            s_ = lat0 == la
            zon[f"{la:.0f}-{la + 5:.0f}N"] = wmean(med[s_], w[s_])
        row["zonal_mean_perK_median"] = zon
        summ["cmip6"][kind] = row

    np.savez(f"{OUT}/solar_grid_per_K.npz", **arrays)
    json.dump(summ, open(f"{OUT}/solar_grid_per_K_summary.json", "w"), indent=1)
    for k, v in summ["ours"].items():
        print(f"  {k:22s} dT {v['dT']:.2f} boxes {v['n_valid']:3d} up {v['frac_up']:.0%} perK {v['area_mean_perK']:+.2f} (mr {v['area_mean_perK_meanremoved']:+.2f}) "
              f"meanCF {v['area_mean_mean_change_perK_all']:+.2f} %/K NDJ {v['ndj_share_area_mean']:.2f} | FDR {v['n_sig']} ({v['n_sig_up']}+/{v['n_sig_down']}-) of {v['n_tested']}, "
              f"p<.05 {v['n_p_lt_005']} (chance {v['expected_false_at_005']:.0f}); area ratio {v['area_mean_ratio_tested']:.3f} null {v['null_area_mean_ratio_p05_p50_p95']} "
              f"p_up {v['p_field_up']:.4f} p_dn {v['p_field_down']:.4f}")
    for kind, v in summ["cmip6"].items():
        print(f"  CMIP6 {kind}: {v}")
    print(f"Saved {OUT}/solar_grid_per_K.npz")


if __name__ == "__main__":
    main()
