"""
Solar-only winter low-output events on the DAILY footing: 20 CMIP6 models (day/rsds), and the EERIE and
TCO95 daily bridges; plus the attribution of the combined (0.44 solar + 0.56 wind) change to wind and
solar. Companion of dunkelflaute/scripts/regions/solar_region_events.py (6-hourly, our runs), with the
machinery of daily_events.py (methods 14) unchanged: 48 h = 2 daily steps, Mockert window expansion,
Sep-Mar winters (>= 200 days), threshold matched to 5 events per control winter per region and metric,
common region set over all comparisons, null for the field mean from the TCO95 daily constant-1950
control (all disjoint 30/30-winter placements, same winters in every region).

Metrics: wind and comb (recomputed, checked against daily_events_summary.json), solar (daily-mean rsds /
1000, absolute) and solar_ds (divided by the control calendar-day climatology, 31-day running mean,
floored at 5 W m-2; future divided by the control climatology). For solar the threshold levels span the
0.1st-40th percentiles of the control's 2-day means.
Per model: field-mean ratio, per K, regions up, mean change of the CF (all 24 regions and the valid set),
mean-removed per K, 5th percentile of the 48 h mean, monthly mean change; 5-95 % from 5000 winter
resamples with the same winters in every region (as fieldsig_cmip6_bootstrap.py).
Ensemble: median per K, 10-90 %, models up, regression of ln ratio on dT, and the correlation of solar
per K with the mean solar change per K across models.
Attribution (comb): factorial recombination of wind and solar winters as in solar_region_events.py
(winters aligned by day of the season; 29 Feb dropped except in 360-day calendars), and the climatology
variants.

Usage: python3 solar_daily_events.py     (module python3/2025.01; about 15 min)
Output: dunkelflaute/results/cmip6/solar_daily_events_summary.json
"""
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/cmip6")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/analysis")
import daily_events as de                                     # noqa: E402
import solar_common as sc                                     # noqa: E402
from robustness_mockert_events import trailing_mean           # noqa: E402

ev = de.ev
ev.WIN = 2
WIN = 2
OUT = f"{_ROOT}/dunkelflaute/results/cmip6"
TARGET = sc.TARGET
W_SOLAR = 0.44
KINDS = ("wind", "comb", "solar", "solar_ds")
R = de.R
NBOOT = 5000
RNG = np.random.default_rng(12)


def prep(d):
    if "di" not in d:
        d["di"] = sc.md_index(d["ymd"][:, 1], d["ymd"][:, 2])
        d["feb29"] = (d["ymd"][:, 1] == 2) & (d["ymd"][:, 2] == 29)
        d["is360"] = bool(((d["ymd"][:, 1] == 2) & (d["ymd"][:, 2] == 30)).any())
    return d


def series(d, kind, clim):
    if kind in ("wind", "comb"):
        return de.metric(d, kind)
    if kind == "solar":
        return d["solar"]
    return sc.deseason(d["solar"], d["di"], clim)


def matched(x, ms, kind):
    if not np.isfinite(x).any():
        return np.nan, 0.0
    if kind in ("wind", "comb"):
        grid = np.linspace(*np.nanpercentile(x, [0.1, 40]), 400)
    else:
        roll = np.concatenate([trailing_mean(x[m], WIN)[WIN - 1:] for m in ms])
        grid = np.linspace(*np.nanpercentile(roll, [0.1, 40]), 400)
    c = np.stack([ev.counts_grid(x[m], grid) for m in ms]).mean(0)
    k = int(np.argmax(c >= TARGET)) if (c >= TARGET).any() else int(np.argmax(c))
    return float(grid[k]), float(c[k])


def counts(x, ms, thr):
    return np.array([ev.counts_grid(x[m], np.array([thr]))[0] for m in ms], float)


def aligned(d, x, ms, L):
    keep = np.ones(len(x), bool) if d["is360"] else ~d["feb29"]
    return np.stack([x[m & keep][:L] for m in ms])


def slot_clim(A):
    raw = np.nanmean(A, axis=0)
    return np.array([np.nanmean(raw[max(0, k - sc.CLIM_HALF):k + sc.CLIM_HALF + 1]) for k in range(len(raw))])


def attribution(Wc, Sc_, Wf, Sf):
    nc, nf = len(Wc), len(Wf)
    ws = 1 - W_SOLAR
    pcc, pfc = sc.pairings(nc, nc, True), sc.pairings(nf, nc, False)
    pcf, pff = sc.pairings(nc, nf, False), sc.pairings(nf, nf, True)
    base = ws * Wc[pcc[:, 0]] + W_SOLAR * Sc_[pcc[:, 1]]
    thr, rate = sc.matched_rows(base, WIN)
    E = lambda X: float(sc.counts_rows(X, thr, WIN).mean())                         # noqa: E731
    e_cc, e_fc = E(base), E(ws * Wf[pfc[:, 0]] + W_SOLAR * Sc_[pfc[:, 1]])
    e_cf, e_ff = E(ws * Wc[pcf[:, 0]] + W_SOLAR * Sf[pcf[:, 1]]), E(ws * Wf[pff[:, 0]] + W_SOLAR * Sf[pff[:, 1]])
    cc = ws * Wc + W_SOLAR * Sc_
    thr_c, rate_c = sc.matched_rows(cc, WIN)
    c_c, c_f = float(sc.counts_rows(cc, thr_c, WIN).mean()), float(sc.counts_rows(ws * Wf + W_SOLAR * Sf, thr_c, WIN).mean())
    scl, wcl = slot_clim(Sc_), slot_clim(Wc)
    xs_c, xs_f = ws * Wc + W_SOLAR * scl[None], ws * Wf + W_SOLAR * scl[None]
    t1, r1 = sc.matched_rows(xs_c, WIN)
    xw_c, xw_f = ws * wcl[None] + W_SOLAR * Sc_, ws * wcl[None] + W_SOLAR * Sf
    t2, r2 = sc.matched_rows(xw_c, WIN)
    lg = lambda a, b: float(np.log(b / a)) if a > 0 and b > 0 else np.nan            # noqa: E731
    return dict(base_rate=rate, ln_wind=lg(e_cc, e_fc), ln_solar=lg(e_cc, e_cf), ln_both=lg(e_cc, e_ff),
                coupled_rate=rate_c, ln_coupled=lg(c_c, c_f),
                clim_solarfixed_rate=r1, ln_clim_solarfixed=lg(sc.counts_rows(xs_c, t1, WIN).mean(), sc.counts_rows(xs_f, t1, WIN).mean()),
                clim_windfixed_rate=r2, ln_clim_windfixed=lg(sc.counts_rows(xw_c, t2, WIN).mean(), sc.counts_rows(xw_f, t2, WIN).mean()))


def comb_event_diagnostics(dc, df, mc, mf, i, thr):
    """As solar_region_events.comb_event_diagnostics, on the daily footing."""
    rw_c, rs_c, fl, mo = [], [], [], []
    for m in mc:
        w, s = dc["wind"][m, i], dc["solar"][m, i]
        rw_c.append(trailing_mean(w, WIN))
        rs_c.append(trailing_mean(s, WIN))
        fl.append(np.nan_to_num(trailing_mean((1 - W_SOLAR) * w + W_SOLAR * s, WIN), nan=np.inf) < thr)
        mo.append(dc["ymd"][m, 1])
    rw_c, rs_c, fl, mo = map(np.concatenate, (rw_c, rs_c, fl, mo))
    rw_f = np.concatenate([trailing_mean(df["wind"][m, i], WIN) for m in mf])
    okc, okf = np.isfinite(rw_c), np.isfinite(rw_f)
    if not fl.any() or not okc.any():
        return {}
    srt = np.sort(rw_c[okc])
    q = float(np.median(np.searchsorted(srt, rw_c[fl]) / len(srt)))
    qc = lambda p: np.quantile(rw_c[okc], p)                    # noqa: E731
    qf = lambda p: np.quantile(rw_f[okf], p)                    # noqa: E731
    wsh = (1 - W_SOLAR) * rw_c[fl] / ((1 - W_SOLAR) * rw_c[fl] + W_SOLAR * rs_c[fl])
    return dict(ctl_event_time_NDJ=float(np.isin(mo[fl], (11, 12, 1)).mean()), wind_share_of_comb_at_events=float(np.nanmean(wsh)),
                wind_rank_at_events_median=q, wind_q_change_pct_at_event_rank=float(100 * (qf(q) / qc(q) - 1)),
                wind_q05_change_pct=float(100 * (qf(0.05) / qc(0.05) - 1)),
                wind_mean_change_pct=float(100 * (rw_f[okf].mean() / rw_c[okc].mean() - 1)))


def main():
    S0 = json.load(open(f"{OUT}/daily_events_summary.json"))
    C = de.comparisons()
    res, store, attr = {}, {}, {}
    for name, (fc, ec, c0, c1, ff, ef, f0, f1, Tc, Tf, group) in C.items():
        if group.startswith("cmip6_windonly"):
            continue
        dc, df = prep(de.load(fc, ec)), prep(de.load(ff, ef))
        mc, mf = de.winter_masks(dc, c0, c1), de.winter_masks(df, f0, f1)
        clim = sc.climatology(dc["solar"], dc["di"], np.any(np.stack(mc), axis=0))
        dT = de.mean_T(Tf, f0, f1) - de.mean_T(Tc, c0, c1)
        row = dict(group=group, dT=dT, n_winters=[len(mc), len(mf)], metrics={})
        for kind in KINDS:
            Xc, Xf = series(dc, kind, clim), series(df, kind, clim)
            per = {}
            for i, r in enumerate(R):
                thr, rate = matched(Xc[:, i], mc, kind)
                if not np.isfinite(thr):
                    per[r] = dict(valid=False)
                    continue
                a, b = counts(Xc[:, i], mc, thr), counts(Xf[:, i], mf, thr)
                xc = np.concatenate([Xc[m, i] for m in mc])
                xf = np.concatenate([Xf[m, i] for m in mf])
                mcm, mfm = np.nanmean(xc), np.nanmean(xf)
                badj = counts(Xf[:, i] * mcm / mfm, mf, thr)
                rc = np.concatenate([trailing_mean(Xc[m, i], WIN)[WIN - 1:] for m in mc])
                rf = np.concatenate([trailing_mean(Xf[m, i], WIN)[WIN - 1:] for m in mf])
                p = dict(valid=bool(rate >= 0.9 * TARGET), threshold=thr, rate=rate, control=float(a.mean()), future=float(b.mean()),
                         ratio=float(b.mean() / a.mean()) if a.mean() else np.nan,
                         ratio_meanremoved=float(badj.mean() / a.mean()) if a.mean() else np.nan,
                         mean_change_pct=float(100 * (mfm / mcm - 1)),
                         p05_48h_change_pct=float(100 * (np.nanpercentile(rf, 5) / np.nanpercentile(rc, 5) - 1)))
                if kind in ("wind", "comb"):
                    p["check_ratio"] = float(p["ratio"] - S0["regions"][name]["metrics"][kind][r]["ratio"])
                if kind == "solar":
                    moc = np.concatenate([dc["ymd"][m, 1] for m in mc])
                    mof = np.concatenate([df["ymd"][m, 1] for m in mf])
                    p["monthly_mean_change_pct"] = {int(k): float(100 * (np.nanmean(xf[mof == k]) / np.nanmean(xc[moc == k]) - 1))
                                                    for k in (9, 10, 11, 12, 1, 2, 3)}
                    f = np.concatenate([sc.flagged(Xc[m, i], thr, WIN) for m in mc])
                    p["ctl_event_time_NDJ"] = float(f[np.isin(moc, (11, 12, 1))].sum() / max(f.sum(), 1))
                if kind == "solar_ds":
                    moc = np.concatenate([dc["ymd"][m, 1] for m in mc])
                    f = np.concatenate([sc.flagged(Xc[m, i], thr, WIN) for m in mc])
                    p["ctl_event_time_NDJ"] = float(f[np.isin(moc, (11, 12, 1))].sum() / max(f.sum(), 1))
                per[r] = p
                store[(name, kind, r)] = (a, b)
            row["metrics"][kind] = per
        for i, r in enumerate(R):
            if "threshold" in row["metrics"]["comb"][r]:
                row["metrics"]["comb"][r].update(comb_event_diagnostics(dc, df, mc, mf, i, row["metrics"]["comb"][r]["threshold"]))
        # attribution
        L = min(min(int((m & (np.ones(len(m), bool) if d["is360"] else ~d["feb29"])).sum()) for m in ms)
                for d, ms in ((dc, mc), (df, mf)))
        A = {}
        for i, r in enumerate(R):
            if not np.isfinite(dc["solar"][:, i]).any():
                continue
            A[r] = attribution(aligned(dc, dc["wind"][:, i], mc, L), aligned(dc, dc["solar"][:, i], mc, L),
                               aligned(df, df["wind"][:, i], mf, L), aligned(df, df["solar"][:, i], mf, L))
        attr[name] = A
        res[name] = row
        print(f"  {name:40s} dT {dT:+.2f}", flush=True)

    common = {k: [r for r in R if all(res[n]["metrics"][k][r].get("valid") for n in res)] for k in KINDS}

    # null from the TCO95 daily control (30/30 winters), per kind over the common set
    null = {}
    d = prep(de.load(f"{de.D}/daily_tco95_ctl1950g.npz", "run"))
    ms = de.winter_masks(d, 1953, 2099)
    clim = sc.climatology(d["solar"], d["di"], np.any(np.stack(ms), axis=0))
    pairs = ev.block_pairs(len(ms), 30, 30)
    for kind in KINDS:
        X = series(d, kind, clim)
        cols = []
        for r in common[kind]:
            i = R.index(r)
            thr, _ = matched(X[:, i], ms, kind)
            cols.append(counts(X[:, i], ms, thr))
        c = np.stack(cols, axis=1)
        fm = np.nanmean(ev.block_logratios(c, 30, 30, pairs), axis=1)
        null[kind] = dict(n_winters=len(c), n_pairs=len(pairs), field_ratio_p05_p95=[float(x) for x in np.exp(np.percentile(fm, [5, 95]))],
                          fm=fm)

    summary = {}
    for name, row in res.items():
        dT = row["dT"]
        pk = lambda x: float(100 * (np.exp(x / dT) - 1))                   # noqa: E731
        s = dict(group=row["group"], dT=dT)
        for kind in KINDS:
            V = common[kind]
            a = np.stack([store[(name, kind, r)][0] for r in V], axis=1)
            b = np.stack([store[(name, kind, r)][1] for r in V], axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                lr = np.log(b.mean(0) / a.mean(0))
            ok = np.isfinite(lr)
            obs = float(lr[ok].mean())
            boots = np.empty(NBOOT)
            for k in range(NBOOT):
                ia, ib = RNG.integers(0, len(a), len(a)), RNG.integers(0, len(b), len(b))
                with np.errstate(divide="ignore", invalid="ignore"):
                    x = np.log(b[ib].mean(0) / a[ia].mean(0))
                boots[k] = np.mean(x[np.isfinite(x)])
            lo, hi = np.percentile(boots, [5, 95])
            la = np.array([np.log(row["metrics"][kind][r]["ratio_meanremoved"]) for r in V])
            all24 = [row["metrics"][kind][r]["mean_change_pct"] for r in R if "mean_change_pct" in row["metrics"][kind][r]]
            s[kind] = dict(n_regions=int(ok.sum()), mean_ratio=float(np.exp(obs)), pct_per_K=pk(obs), pct_per_K_p05_p95=[pk(lo), pk(hi)],
                           ratio_p05_p95=[float(np.exp(lo)), float(np.exp(hi))], n_regions_up=int((lr[ok] > 0).sum()),
                           pct_per_K_meanremoved=pk(float(np.nanmean(la))), share_from_mean=float(1 - np.nanmean(la) / obs) if obs else np.nan,
                           mean_change_pct=float(np.mean([row["metrics"][kind][r]["mean_change_pct"] for r in V])),
                           mean_change_pct_all24=float(np.mean(all24)), mean_change_per_K_all24=float(np.mean(all24) / dT),
                           p05_48h_change_pct=float(np.mean([row["metrics"][kind][r]["p05_48h_change_pct"] for r in V])),
                           p_null_up=float(np.mean(null[kind]["fm"] >= obs)), p_null_down=float(np.mean(null[kind]["fm"] <= obs)))
            if kind in ("wind", "comb"):
                s[kind]["check_vs_daily_events"] = float(np.exp(obs) - S0["comparisons"][name][kind]["mean_ratio"]) \
                    if common[kind] == S0["common_regions"][kind] else None
            if kind == "solar":
                s[kind]["monthly_mean_change_pct_all24"] = {m: float(np.mean([row["metrics"]["solar"][r]["monthly_mean_change_pct"][m]
                                                                               for r in R if "monthly_mean_change_pct" in row["metrics"]["solar"][r]]))
                                                            for m in (9, 10, 11, 12, 1, 2, 3)}
            if kind in ("solar", "solar_ds"):
                s[kind]["ctl_event_time_NDJ"] = float(np.mean([row["metrics"][kind][r]["ctl_event_time_NDJ"] for r in V]))
        # attribution field means over the comb common set where the decoupled base reaches the rate
        A = attr[name]
        V = [r for r in common["comb"] if r in A and A[r]["base_rate"] >= 0.9 * TARGET]
        f = {k: float(np.nanmean([A[r][k] for r in V])) for k in ("ln_wind", "ln_solar", "ln_both", "ln_coupled",
                                                                    "ln_clim_solarfixed", "ln_clim_windfixed")}
        s["attribution"] = dict(n_regions=len(V), pct_per_K={k[3:]: pk(v) for k, v in f.items()},
                                pct_per_K_interaction=pk(f["ln_both"] - f["ln_wind"] - f["ln_solar"]),
                                share_wind=f["ln_wind"] / f["ln_both"] if f["ln_both"] else np.nan,
                                share_solar=f["ln_solar"] / f["ln_both"] if f["ln_both"] else np.nan,
                                pct_per_K_wind_alone_same_regions=pk(float(np.mean([np.log(row["metrics"]["wind"][r]["ratio"]) for r in V]))),
                                comb_events={k: float(np.nanmean([row["metrics"]["comb"][r].get(k, np.nan) for r in V]))
                                             for k in ("ctl_event_time_NDJ", "wind_share_of_comb_at_events", "wind_rank_at_events_median",
                                                       "wind_q_change_pct_at_event_rank", "wind_q05_change_pct", "wind_mean_change_pct")},
                                per_region=A)
        summary[name] = s

    ens = {}
    for grp in ("cmip6", "cmip6_mid"):
        names = [n for n in summary if summary[n]["group"] == grp]
        ens[grp] = {}
        for kind in KINDS:
            lnr = np.array([np.log(summary[n][kind]["mean_ratio"]) for n in names])
            dT = np.array([summary[n]["dT"] for n in names])
            pkv = np.array([summary[n][kind]["pct_per_K"] for n in names])
            ens[grp][kind] = dict(n_models=len(names), pct_per_K_median=float(np.median(pkv)),
                                  pct_per_K_p10_p90=[float(x) for x in np.percentile(pkv, [10, 90])],
                                  n_models_up=int((lnr > 0).sum()), corr_lnratio_dT=float(np.corrcoef(dT, lnr)[0, 1]),
                                  origin_fit_pct_per_K=float(100 * (np.exp((dT * lnr).sum() / (dT * dT).sum()) - 1)),
                                  n_models_beyond_null_up=int(sum(summary[n][kind]["p_null_up"] < 0.05 for n in names)),
                                  n_models_beyond_null_down=int(sum(summary[n][kind]["p_null_down"] < 0.05 for n in names)),
                                  mean_change_per_K_all24_median=float(np.median([summary[n][kind]["mean_change_per_K_all24"] for n in names])))
        sm = np.array([summary[n]["solar"]["mean_change_per_K_all24"] for n in names])
        for kind in ("solar", "solar_ds"):
            pkv = np.array([summary[n][kind]["pct_per_K"] for n in names])
            ens[grp][kind]["corr_with_mean_solar_change_per_K"] = float(np.corrcoef(sm, pkv)[0, 1])
            sl, ic = np.polyfit(sm, pkv, 1)
            ens[grp][kind]["fit_vs_mean_solar_change"] = [float(ic), float(sl)]
        at = [summary[n]["attribution"] for n in names]
        ens[grp]["attribution"] = {k: float(np.median([a["pct_per_K"][k] for a in at])) for k in at[0]["pct_per_K"]}
        ens[grp]["attribution"]["share_wind_median"] = float(np.median([a["share_wind"] for a in at]))
        ens[grp]["attribution"]["wind_alone_same_regions"] = float(np.median([a["pct_per_K_wind_alone_same_regions"] for a in at]))
        ens[grp]["attribution"]["comb_events"] = {k: float(np.median([a["comb_events"][k] for a in at])) for k in at[0]["comb_events"]}
    for k in null:
        null[k].pop("fm")
    out = dict(target=TARGET, window_days=WIN, clim_floor=sc.CLIM_FLOOR, common_regions=common, comparisons=summary,
               null_tco95_daily=null, ensemble=ens,
               regions={n: dict(group=r["group"], dT=r["dT"], metrics=r["metrics"]) for n, r in res.items()})
    json.dump(out, open(f"{OUT}/solar_daily_events_summary.json", "w"), indent=1)

    for kind in KINDS:
        print(f"\n=== {kind} ({len(common[kind])} common regions: {' '.join(common[kind])})")
        for n, s in sorted(summary.items(), key=lambda x: (x[1]["group"], x[0])):
            x = s[kind]
            print(f"  {n:40s} dT {s['dT']:+.2f} ratio {x['mean_ratio']:.3f} {x['pct_per_K']:+5.1f} %/K [{x['pct_per_K_p05_p95'][0]:+.1f},{x['pct_per_K_p05_p95'][1]:+.1f}] "
                  f"up {x['n_regions_up']:2d}/{x['n_regions']} meanCF {x['mean_change_pct_all24']:+5.1f}% ({x['mean_change_per_K_all24']:+.2f}/K) "
                  f"mr {x['pct_per_K_meanremoved']:+.1f} p05 {x['p05_48h_change_pct']:+.1f}% pnull {x['p_null_up']:.3f}/{x['p_null_down']:.3f}"
                  + (f" chk {x['check_vs_daily_events']:+.1e}" if x.get("check_vs_daily_events") is not None else ""))
        print(f"  null {null[kind]['field_ratio_p05_p95']}")
        for g, e in ens.items():
            x = e[kind]
            print(f"  ensemble {g}: median {x['pct_per_K_median']:+.1f} (10-90 {x['pct_per_K_p10_p90'][0]:+.1f}..{x['pct_per_K_p10_p90'][1]:+.1f}), "
                  f"up {x['n_models_up']}/{x['n_models']}, beyond null up {x['n_models_beyond_null_up']} down {x['n_models_beyond_null_down']}, "
                  f"origin fit {x['origin_fit_pct_per_K']:+.1f}, r(dT) {x['corr_lnratio_dT']:.2f}, mean change {x['mean_change_per_K_all24_median']:+.2f} %/K"
                  + (f", r(mean solar change) {x['corr_with_mean_solar_change_per_K']:.2f}, fit {x['fit_vs_mean_solar_change']}" if "corr_with_mean_solar_change_per_K" in x else ""))
    print("\n=== attribution of comb (%/K)")
    for n, s in sorted(summary.items(), key=lambda x: (x[1]["group"], x[0])):
        a = s["attribution"]
        p = a["pct_per_K"]
        print(f"  {n:40s} {a['n_regions']:2d}: coupled {p['coupled']:+.1f} both {p['both']:+.1f} = wind {p['wind']:+.1f} + solar {p['solar']:+.1f} "
              f"+ int {a['pct_per_K_interaction']:+.1f}; share wind {a['share_wind']:.0%}; clim solar fixed {p['clim_solarfixed']:+.1f}, wind fixed {p['clim_windfixed']:+.1f}; "
              f"wind alone {a['pct_per_K_wind_alone_same_regions']:+.1f}; comb events " + ", ".join(f"{k} {v:.2f}" for k, v in a["comb_events"].items()))
    for g, e in ens.items():
        print(f"  ensemble {g} medians: {e['attribution']}")
    print(f"Saved {OUT}/solar_daily_events_summary.json")


if __name__ == "__main__":
    main()
