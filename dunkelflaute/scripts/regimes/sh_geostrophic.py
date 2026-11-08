"""
Southern-Hemisphere test of the geostrophic mechanism (methods notes 13b): if NH winter low-wind
events increase because the large-scale pressure gradients weaken, then SH winter (April-September)
events should change with the sign of the change in the geostrophic speed G, including where G
strengthens.

Per comparison and SH region (sh_regions.py; daily series from sh_extract.py), on the daily
footing of geostrophic.py and daily_events.py:
  G    region mean of |grad p| / (rho |f|) from daily-mean msl on the 2.5 deg SH grid (centred
       differences, periodic in longitude, rho = 1.25 kg m-3), grid points inside the region box
  s    daily-mean region-mean speed (CMIP6: sfcWind; our runs: log-law 100 m speed, as geostrophic.py)
  a_m  control transfer sum(s) / sum(G) per calendar month; s* = a_m G
  R_s, R_G   future / control events of s and s* (48 h = 2 days, Mockert window expansion, each
             series with its own threshold matched to 5 events per control winter)
  R_cf       the same for the onshore wind CF (the wind-only events of daily_events.py), with the
             change of the mean CF and the ratio after removing the mean change
  per K      100 (exp(ln R / dT) - 1); for G: 100 ((G_fut / G_ctl)^(1 / dT) - 1)
Winters are April-September of one year. Statistics:
  * our runs: winter bootstrap (5000, the same winters in every region) of R_s and R_cf per K
  * null: the TCO95 constant-1950 control (147 SH winters), all disjoint block placements of the
    comparison's control and future lengths (30/30 for CMIP6), the same winters in every region
Cross-hemisphere test: per CMIP6 model and region, event change per K (R_s) against G change per K,
NH from geostrophic_summary.json (22 regions) and SH here (7 named regions): correlation, OLS fit,
sign agreement, and the SH prediction from the NH fit alone.

Usage: python3 sh_geostrophic.py
Output: dunkelflaute/results/sh/sh_geostrophic_summary.json
"""
import glob
import json
import os
import sys

import numpy as np
from scipy import stats

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("regions", "analysis", "regimes", "cmip6"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/{p}")
import sh_regions as sh                    # noqa: E402
import region_events as ev                 # noqa: E402
import region_per_K as rpk                 # noqa: E402
import geostrophic as gs                   # noqa: E402

ev.WIN = 2
D = sh.OUT_DIR
OUT = f"{_ROOT}/dunkelflaute/results/sh"
TARGET = 5.0
MIN_DAYS = 170
NBOOT = 5000
RNG = np.random.default_rng(7)
BOXES = sh.sh_boxes()
NAMED = list(sh.NAMED)
OURS = {   # name: (file, control span), (file, future span), global-t keys
    "9km":        (("dart9km_1950c", 1953, 1969), ("dart9km_2080c", 2083, 2092)),
    "31km":       (("tco319_ctl1950d", 1953, 1969), ("tco319_ssp585d", 2081, 2099)),
    "100km":      (("tco95_ctl1950g", 1953, 2099), ("tco95_hi1950g", 2081, 2099)),
    "eerie_ctl":  (("eerie_ctl", 1972, 2050), ("eerie_ssp245", 2041, 2070)),
    "eerie_hist": (("eerie_hist", 1971, 2000), ("eerie_ssp245", 2041, 2070)),
}


def pct_per_K(lnr, dT):
    return float(100 * (np.exp(lnr / dT) - 1))


def load(fn, exp, speed):
    z = np.load(fn)
    d = {k.split("|", 1)[1]: z[k] for k in z.files if k.startswith(exp + "|")}
    d["regions"] = [str(x) for x in z["regions"]]
    assert d["regions"] == list(BOXES)
    d["k"] = d["ymd"] @ np.array([10000, 100, 1])
    d["year"], d["month"] = d["ymd"][:, 0], d["ymd"][:, 1]
    d["s"] = d[speed]
    vg = np.abs(gs.geostrophic(d["psl"].astype(np.float64), z["plat"], z["plon"], True))
    LA, LO = np.meshgrid(z["plat"], z["plon"], indexing="ij")
    d["G"] = np.stack([np.nanmean(vg[:, (LA >= a) & (LA <= b) & (LO >= c) & (LO <= e)], axis=1)
                       for a, b, c, e in BOXES.values()], axis=1)
    return d


def masks(d, y0, y1):
    return [d["year"] == y for y in range(y0, y1 + 1) if (d["year"] == y).sum() >= MIN_DAYS]


def matched(xc, mc):
    grid = np.linspace(*np.nanpercentile(xc, [0.1, 40]), 400)
    c = np.stack([ev.counts_grid(xc[m], grid) for m in mc]).mean(0)
    j = int(np.argmax(c >= TARGET)) if (c >= TARGET).any() else int(np.argmax(c))
    return float(grid[j]), float(c[j])


def per_winter(x, ms, thr):
    return np.array([ev.counts_grid(x[m], np.array([thr]))[0] for m in ms], float)


def analyse(dc, c0, c1, df, f0, f1, dT):
    mc, mf = masks(dc, c0, c1), masks(df, f0, f1)
    sel_c, sel_f = np.any(mc, axis=0), np.any(mf, axis=0)
    per, cnt = {}, {}
    for i, r in enumerate(BOXES):
        sc, sf, Gc, Gf = dc["s"][:, i], df["s"][:, i], dc["G"][:, i], df["G"][:, i]
        if not (np.isfinite(sc[sel_c]).all() and np.isfinite(Gc[sel_c]).all()):
            continue
        a = {m: np.sum(sc[sel_c & (dc["month"] == m)]) / np.sum(Gc[sel_c & (dc["month"] == m)]) for m in sh.MONTHS}
        xc = np.array([a[m] for m in dc["month"]]) * Gc
        xf = np.array([a[m] for m in df["month"]]) * Gf
        row = {}
        for key, (yc, yf) in dict(s=(sc, sf), G=(xc, xf), cf=(dc["wind"][:, i], df["wind"][:, i])).items():
            thr, rate = matched(yc, mc)
            A, B = per_winter(yc, mc, thr), per_winter(yf, mf, thr)
            row[f"ratio_{key}"] = float(B.mean() / A.mean())
            row[f"events_ctl_{key}"] = float(A.mean())
            row[f"events_fut_{key}"] = float(B.mean())
            row[f"valid_{key}"] = bool(rate >= 0.9 * TARGET)
            cnt[(r, key)] = (A, B)
            if key == "cf":
                mcf, mff = np.mean(yc[sel_c]), np.mean(yf[sel_f])
                Bm = per_winter(yf * mcf / mff, mf, thr)
                row["ratio_cf_meanremoved"] = float(Bm.mean() / A.mean())
                row["mean_cf_change_pct"] = float(100 * (mff / mcf - 1))
        ms_c, ms_f, mg_c, mg_f = sc[sel_c].mean(), sf[sel_f].mean(), Gc[sel_c].mean(), Gf[sel_f].mean()
        row.update(mean_s_change_pct=float(100 * (ms_f / ms_c - 1)), mean_G_change_pct=float(100 * (mg_f / mg_c - 1)),
                   transfer_ctl=float(ms_c / mg_c),
                   transfer_change_pct=float(100 * ((ms_f / mg_f) / (ms_c / mg_c) - 1)),
                   corr_s_G_ctl=float(np.corrcoef(sc[sel_c], Gc[sel_c])[0, 1]))
        row["pct_per_K_s"] = pct_per_K(np.log(row["ratio_s"]), dT)
        row["pct_per_K_G_events"] = pct_per_K(np.log(row["ratio_G"]), dT)
        row["pct_per_K_cf"] = pct_per_K(np.log(row["ratio_cf"]), dT)
        row["G_pct_per_K"] = float(100 * ((mg_f / mg_c) ** (1 / dT) - 1))
        row["share_circulation"] = float(np.log(row["ratio_G"]) / np.log(row["ratio_s"])) if row["ratio_s"] != 1 else np.nan
        per[r] = row
    return per, cnt, (len(mc), len(mf))


def field(per, regs, dT):
    V = [r for r in regs if r in per]
    ls = np.array([np.log(per[r]["ratio_s"]) for r in V])
    lg = np.array([np.log(per[r]["ratio_G"]) for r in V])
    lc = np.array([np.log(per[r]["ratio_cf"]) for r in V])
    dG = np.array([per[r]["mean_G_change_pct"] for r in V])
    return dict(n_regions=len(V), mean_ratio_s=float(np.exp(ls.mean())), mean_ratio_G=float(np.exp(lg.mean())),
                mean_ratio_cf=float(np.exp(lc.mean())), pct_per_K_s=pct_per_K(ls.mean(), dT),
                pct_per_K_G_events=pct_per_K(lg.mean(), dT), pct_per_K_cf=pct_per_K(lc.mean(), dT),
                share_circulation=float(lg.mean() / ls.mean()) if ls.mean() else np.nan,
                slope_lnRs_on_lnRG=float(np.sum(ls * lg) / np.sum(lg * lg)),
                r_lnRs_lnRG=float(np.corrcoef(ls, lg)[0, 1]) if len(V) > 2 else np.nan,
                mean_G_change_pct=float(dG.mean()),
                transfer_change_pct=float(np.mean([per[r]["transfer_change_pct"] for r in V])),
                n_G_up=int((dG > 0).sum()), n_s_events_up=int((ls > 0).sum()),
                n_sign_opposite=int((np.sign(ls) == -np.sign(dG)).sum()))


def bootstrap(cnt, regs, key, n_c, n_f, dT):
    V = [r for r in regs if (r, key) in cnt]
    A = np.stack([cnt[(r, key)][0] for r in V], axis=1)
    B = np.stack([cnt[(r, key)][1] for r in V], axis=1)
    ia = RNG.integers(0, n_c, (NBOOT, n_c))
    ib = RNG.integers(0, n_f, (NBOOT, n_f))
    with np.errstate(divide="ignore", invalid="ignore"):
        lr = np.log(B[ib].mean(1) / A[ia].mean(1))           # (NBOOT, R)
    pk = 100 * (np.exp(lr / dT) - 1)
    fm = 100 * (np.exp(np.nanmean(lr, 1) / dT) - 1)
    return dict(field_pct_per_K_p05_p95=[float(x) for x in np.nanpercentile(fm, [5, 95])],
                region_pct_per_K_p05_p95={r: [float(x) for x in np.nanpercentile(pk[:, j], [5, 95])] for j, r in enumerate(V)})


def null_counts(d, c0, c1):
    """per-winter counts (W, R) of s at the control-matched thresholds, TCO95 control."""
    mc = masks(d, c0, c1)
    out = {}
    for i, r in enumerate(BOXES):
        x = d["s"][:, i]
        thr, _ = matched(x, mc)
        out[r] = per_winter(x, mc, thr)
    return out


def null_dist(nc_, regs, n_a, n_b):
    V = [r for r in regs if r in nc_]
    c = np.stack([nc_[r] for r in V], axis=1)
    pairs = ev.block_pairs(len(c), n_a, n_b)
    lr = ev.block_logratios(c, n_a, n_b, pairs)
    return V, lr


def cmip_models():
    inv = json.load(open(f"{_ROOT}/dunkelflaute/data/cmip6/inventory.json"))
    out = {}
    for fn in sorted(glob.glob(f"{D}/sh_cmip6_*.npz")):
        m = os.path.basename(fn)[len("sh_cmip6_"):-4]
        out[m] = (fn, "cmip6_windonly" if inv.get(m, {}).get("wind_only") else "cmip6")
    return out


def fit(x, y):
    x, y = np.asarray(x), np.asarray(y)
    b, a = np.polyfit(x, y, 1)
    return dict(n=len(x), r=float(np.corrcoef(x, y)[0, 1]), slope=float(b), intercept=float(a),
                sign_agree=float(np.mean(np.sign(y) == -np.sign(x))))


def main():
    os.makedirs(OUT, exist_ok=True)
    res = dict(regions=list(BOXES), named=NAMED, long_names=sh.LONG, boxes=BOXES, cmip6={}, ours={})
    T = rpk.global_t()

    # null from the TCO95 constant-1950 control (SH winters 1953-2099)
    t95c = load(f"{D}/sh_tco95_ctl1950g.npz", "run", "spd100")
    NC = null_counts(t95c, 1953, 2099)
    grid = [r for r in BOXES if r not in NAMED]

    def nullinfo(n_a, n_b):
        o = {}
        for nm, regs in (("named", NAMED), ("grid", grid)):
            V, lr = null_dist(NC, regs, n_a, n_b)
            fm = np.nanmean(lr, 1)
            o[nm] = dict(n_pairs=len(lr), field_ratio_p05_p95=[float(x) for x in np.exp(np.percentile(fm, [5, 95]))])
            if nm == "named":
                o["region_ratio_p05_p95"] = {r: [float(x) for x in np.exp(np.nanpercentile(lr[:, j], [5, 95]))] for j, r in enumerate(V)}
        return o
    res["null_tco95"] = {"30/30": nullinfo(30, 30)}

    for m, (fn, grp) in cmip_models().items():
        z = np.load(fn)
        Tc = dict(zip(z["historical|tas_years"].tolist(), z["historical|tas"].tolist()))
        Tf = dict(zip(z["ssp585|tas_years"].tolist(), z["ssp585|tas"].tolist()))
        dT = float(np.mean([Tf[y] for y in range(2071, 2101) if y in Tf]) - np.mean([Tc[y] for y in range(1951, 1981) if y in Tc]))
        dc, df = load(fn, "historical", "spd10"), load(fn, "ssp585", "spd10")
        per, cnt, (n_c, n_f) = analyse(dc, 1951, 1980, df, 2071, 2100, dT)
        res["cmip6"][m] = dict(group=grp, dT=dT, n_winters=[n_c, n_f], field_named=field(per, NAMED, dT),
                               field_grid=field(per, grid, dT), regions=per)
        f = res["cmip6"][m]["field_named"]
        print(f"  {m:18s} dT {dT:.2f}  named: R_s {f['mean_ratio_s']:.2f} R_G {f['mean_ratio_G']:.2f} R_cf {f['mean_ratio_cf']:.2f} "
              f"dG {f['mean_G_change_pct']:+.1f}%  transfer {f['transfer_change_pct']:+.1f}%  G up {f['n_G_up']}/7", flush=True)

    for name, ((cn, c0, c1), (fn_, f0, f1)) in OURS.items():
        src_c, run_c = cn.split("_", 1)
        src_f, run_f = fn_.split("_", 1)
        fc, ff = f"{D}/sh_{src_c}_{run_c}.npz", f"{D}/sh_{src_f}_{run_f}.npz"
        if not (os.path.exists(fc) and os.path.exists(ff)):
            print(f"  {name}: missing input")
            continue
        dT = rpk.mean_years(T, fn_, f0, f1)[0] - rpk.mean_years(T, cn, c0, c1)[0]
        dc = t95c if fc.endswith("sh_tco95_ctl1950g.npz") else load(fc, "run", "spd100")
        df = load(ff, "run", "spd100")
        per, cnt, (n_c, n_f) = analyse(dc, c0, c1, df, f0, f1, dT)
        row = dict(dT=dT, n_winters=[n_c, n_f], field_named=field(per, NAMED, dT), field_grid=field(per, grid, dT), regions=per)
        for nm, regs in (("named", NAMED), ("grid", grid)):
            row[f"boot_{nm}_s"] = bootstrap(cnt, regs, "s", n_c, n_f, dT)
            row[f"boot_{nm}_cf"] = bootstrap(cnt, regs, "cf", n_c, n_f, dT)
        if fc.endswith("sh_tco95_ctl1950g.npz") and n_c > n_f * 2:
            row["null"] = "the control itself (147 winters)"
        key = f"{n_c}/{n_f}"
        if key not in res["null_tco95"] and n_c + n_f <= 147:
            res["null_tco95"][key] = nullinfo(n_c, n_f)
        res["ours"][name] = row
        f = row["field_named"]
        print(f"  {name:11s} dT {dT:.2f}  named: R_s {f['mean_ratio_s']:.2f} ({f['pct_per_K_s']:+.1f} %/K, boot {row['boot_named_s']['field_pct_per_K_p05_p95']}) "
              f"R_G {f['mean_ratio_G']:.2f} R_cf {f['mean_ratio_cf']:.2f} dG {f['mean_G_change_pct']:+.1f}% transfer {f['transfer_change_pct']:+.1f}%", flush=True)

    # CMIP6 ensemble per region
    C = {m: v for m, v in res["cmip6"].items() if v["group"] == "cmip6"}
    ens = {}
    for r in BOXES:
        rows = [(v["regions"][r], v["dT"]) for v in C.values() if r in v["regions"]]
        if len(rows) < 10:
            continue
        e = {}
        for k in ("pct_per_K_s", "pct_per_K_G_events", "pct_per_K_cf", "G_pct_per_K", "mean_G_change_pct",
                  "transfer_change_pct", "mean_s_change_pct", "mean_cf_change_pct", "corr_s_G_ctl", "transfer_ctl"):
            x = np.array([a[k] for a, _ in rows])
            e[k] = dict(median=float(np.median(x)), p10_p90=[float(v) for v in np.percentile(x, [10, 90])])
        e["n_models"] = len(rows)
        e["n_models_G_up"] = int(sum(a["mean_G_change_pct"] > 0 for a, _ in rows))
        e["n_models_events_down"] = int(sum(a["ratio_s"] < 1 for a, _ in rows))
        e["n_models_cf_events_down"] = int(sum(a["ratio_cf"] < 1 for a, _ in rows))
        ens[r] = e
    res["cmip6_ensemble_regions"] = ens
    fe = {}
    for nm in ("field_named", "field_grid"):
        fe[nm] = {}
        for k in ("pct_per_K_s", "pct_per_K_G_events", "pct_per_K_cf", "mean_G_change_pct", "transfer_change_pct", "share_circulation",
                  "slope_lnRs_on_lnRG", "r_lnRs_lnRG"):
            x = np.array([v[nm][k] for v in C.values()])
            fe[nm][k] = dict(median=float(np.nanmedian(x)), p10_p90=[float(v) for v in np.nanpercentile(x, [10, 90])])
    res["cmip6_ensemble_field"] = fe

    # cross-hemisphere test (CMIP6, per model and region): events per K against G per K, for the
    # speed-space events of 13b (R_s) and for the wind-CF events of daily_events.py
    g = json.load(open(f"{_ROOT}/dunkelflaute/results/regimes/geostrophic_summary.json"))
    de = json.load(open(f"{_ROOT}/dunkelflaute/results/cmip6/daily_events_summary.json"))
    res["cross_hemisphere"], res["points"] = {}, {}
    for kind in ("s", "cf"):
        nh = []
        for m, v in g["cmip6"].items():
            for r, a in v["regions"].items():
                gpk = 100 * ((1 + a["mean_G_change_pct"] / 100) ** (1 / v["dT"]) - 1)
                if kind == "s":
                    nh.append((m, r, gpk, pct_per_K(np.log(a["ratio_s"]), v["dT"])))
                else:
                    e = de["regions"][f"{m}|2071-2100"]
                    nh.append((m, r, gpk, pct_per_K(np.log(e["metrics"]["wind"][r]["ratio"]), e["dT"])))
        key = f"pct_per_K_{kind}"
        shp = [(m, r, v["regions"][r]["G_pct_per_K"], v["regions"][r][key]) for m, v in C.items() for r in NAMED if r in v["regions"]]
        shg = [(m, r, v["regions"][r]["G_pct_per_K"], v["regions"][r][key]) for m, v in C.items() for r in grid
               if r in v["regions"] and np.isfinite(v["regions"][r][key])]
        xn, yn = np.array([a[2] for a in nh]), np.array([a[3] for a in nh])
        xs, ys = np.array([a[2] for a in shp]), np.array([a[3] for a in shp])
        xg, yg = np.array([a[2] for a in shg]), np.array([a[3] for a in shg])
        test = dict(nh=fit(xn, yn), sh=fit(xs, ys), sh_grid=fit(xg, yg), pooled=fit(np.r_[xn, xs], np.r_[yn, ys]))
        fn_nh = test["nh"]
        pred = fn_nh["intercept"] + fn_nh["slope"] * xs
        test["sh_from_nh_fit"] = dict(r=float(np.corrcoef(pred, ys)[0, 1]), rmse=float(np.sqrt(np.mean((pred - ys) ** 2))),
                                      bias=float(np.mean(ys - pred)), sign_hit=float(np.mean(np.sign(pred) == np.sign(ys))),
                                      sd_obs=float(ys.std()))
        up = xs > 0
        test["sh_where_G_up"] = dict(n=int(up.sum()), frac_events_down=float(np.mean(ys[up] < 0)), median_events=float(np.median(ys[up])))
        test["sh_where_G_down"] = dict(n=int((~up).sum()), frac_events_up=float(np.mean(ys[~up] > 0)), median_events=float(np.median(ys[~up])))
        # one point per region: the CMIP6 median
        med_nh = {r: (np.median([a[2] for a in nh if a[1] == r]), np.median([a[3] for a in nh if a[1] == r])) for r in g["cmip6_regions"]}
        med_sh = {r: (np.median([a[2] for a in shp if a[1] == r]), np.median([a[3] for a in shp if a[1] == r])) for r in NAMED}
        med_g = {r: (np.median([a[2] for a in shg if a[1] == r]), np.median([a[3] for a in shg if a[1] == r]))
                 for r in grid if sum(a[1] == r for a in shg) >= 10}
        xm = np.array([v[0] for v in list(med_nh.values()) + list(med_sh.values())])
        ym = np.array([v[1] for v in list(med_nh.values()) + list(med_sh.values())])
        test["region_medians_pooled"] = fit(xm, ym)
        test["region_medians_pooled"]["p"] = float(stats.pearsonr(xm, ym)[1])
        test["region_medians_nh"] = fit([v[0] for v in med_nh.values()], [v[1] for v in med_nh.values()])
        test["region_medians_sh"] = fit([v[0] for v in med_sh.values()], [v[1] for v in med_sh.values()])
        test["grid_medians_sh"] = fit([v[0] for v in med_g.values()], [v[1] for v in med_g.values()])
        test["grid_medians_sh"]["p"] = float(stats.pearsonr([v[0] for v in med_g.values()], [v[1] for v in med_g.values()])[1])
        test["region_medians"] = dict(nh={r: [float(a), float(b)] for r, (a, b) in med_nh.items()},
                                      sh={r: [float(a), float(b)] for r, (a, b) in med_sh.items()},
                                      sh_grid={r: [float(a), float(b)] for r, (a, b) in med_g.items()})
        # within each model, across its 29 regions (NH 22 + SH 7)
        rm = []
        for m in C:
            a = [(x, y) for mm, _, x, y in nh + shp if mm == m]
            rm.append(np.corrcoef(*np.array(a).T)[0, 1])
        test["per_model_r"] = dict(median=float(np.median(rm)), p10_p90=[float(v) for v in np.percentile(rm, [10, 90])],
                                   n_negative=int((np.array(rm) < 0).sum()), n_models=len(rm))
        # our runs, 7 named SH regions each
        o = [(n, r, v["regions"][r]["G_pct_per_K"], v["regions"][r][key]) for n, v in res["ours"].items() for r in NAMED]
        test["ours_sh"] = fit([a[2] for a in o], [a[3] for a in o])
        res["cross_hemisphere"][kind] = test
        res["points"][kind] = dict(nh=[list(a) for a in nh], sh=[list(a) for a in shp], sh_grid=[list(a) for a in shg], ours_sh=[list(a) for a in o])
    test = res["cross_hemisphere"]
    json.dump(res, open(f"{OUT}/sh_geostrophic_summary.json", "w"), indent=1, default=float)

    print("\nCMIP6 ensemble per named SH region: events %/K (s, CF), G %/K, transfer change, models G up / events down")
    for r in NAMED:
        if r in ens:
            e = ens[r]
            print(f"  {r:6s} s {e['pct_per_K_s']['median']:+5.1f} [{e['pct_per_K_s']['p10_p90'][0]:+.1f},{e['pct_per_K_s']['p10_p90'][1]:+.1f}]  "
                  f"s* {e['pct_per_K_G_events']['median']:+5.1f}  cf {e['pct_per_K_cf']['median']:+5.1f}  G {e['G_pct_per_K']['median']:+5.2f} %/K "
                  f"(dG {e['mean_G_change_pct']['median']:+.1f}%)  transfer {e['transfer_change_pct']['median']:+.1f}%  "
                  f"G up {e['n_models_G_up']}/{e['n_models']}  events down {e['n_models_events_down']}  cf down {e['n_models_cf_events_down']}")
    print("CMIP6 field (named / grid):", json.dumps({k: {kk: round(vv['median'], 2) for kk, vv in v.items()} for k, v in fe.items()}))
    for kind in ("s", "cf"):
        t = test[kind]
        print(f"Cross-hemisphere [{kind}]: " + "; ".join(f"{k} r {t[k]['r']:+.2f} slope {t[k]['slope']:+.2f} icpt {t[k]['intercept']:+.2f} sign {t[k]['sign_agree']:.2f} n {t[k]['n']}"
              for k in ("nh", "sh", "sh_grid", "pooled", "region_medians_pooled", "region_medians_nh", "region_medians_sh", "grid_medians_sh", "ours_sh")))
        print(f"   p(region medians) {t['region_medians_pooled']['p']:.1e}, p(grid medians) {t['grid_medians_sh']['p']:.1e}; SH from NH fit {t['sh_from_nh_fit']}; "
              f"G up {t['sh_where_G_up']}; G down {t['sh_where_G_down']}; per-model r {t['per_model_r']}")
    print("Null TCO95:", json.dumps(res["null_tco95"])[:1500])
    print(f"Saved {OUT}/sh_geostrophic_summary.json")


if __name__ == "__main__":
    main()
