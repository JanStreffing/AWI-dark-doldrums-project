"""
Solar-only winter low-output events in our runs (24 regions, 6-hourly footing), and the attribution of
the combined (0.44 solar + 0.56 wind) change to wind and solar.

Same machinery as region_events.py / region_per_K.py (methods 10c, 10d): the solar CF series that the
combined metric already uses (9 km v5 with exact monthly rescale, 31 km v7, 100 km native ssrd,
IFS-FESOM2-SR native ssrd; 3-hourly sources degraded to 6 h by region_events.to_6h), Sep-Mar winters,
48 h trailing mean = 8 steps, Mockert window expansion, threshold matched so the comparison's control
has 5 events per winter (lowest of 400 levels, as region_events.matched_threshold; for solar the
levels span the 0.1st-40th percentiles of the control's 48 h means, because half the raw 6-hourly
values are night zeros).

Metrics
  wind, comb   recomputed here as a check against region_events_summary.json (must agree exactly)
  solar        solar CF, absolute (consistent with the wind definition)
  solar_ds     solar CF divided by its control calendar-day climatology (31-day running mean, floored
               at 5 W m-2), i.e. the fraction of the normal irradiance for the date; the future is
               divided by the CONTROL climatology, so mean changes count. Its 48 h mean is a 2-day mean
               of whole days (8 steps of 6 h), so the diurnal cycle does not enter.
Per region: ratio of events per winter, Poisson p, mean change of the CF (Sep-Mar), mean-removed ratio
(future x mean_ctl / mean_fut, same threshold), change of the 5th percentile of the 48 h mean, monthly
mean changes, and in the control the share of event time per calendar month (does the absolute
definition only count December?).
Field means per K (region_per_K.py): ln(ratio)/dT, field mean over valid regions (control reaches 90 %
of 5 events per winter in every comparison and both long controls), 5-95 % from 5000 winter resamples
with the same winters in every region; null for the field mean from the TCO95 constant-1950 control
(147 winters), all placements of two disjoint blocks with the comparison's own winter counts (control
block cut to W - n_b - 40 where too few placements remain, as grid_per_K.null_ranges), same winters in
every region.

Attribution of comb (regions valid for comb in region_events_summary.json):
  factorial    wind and solar series aligned by calendar slot (29 Feb dropped) and recombined across
               winters: control wind of winter i with control solar of another winter j (base), future
               wind with control solar (wind change only), control wind with future solar (solar change
               only), future wind with future solar of another winter (both). 16 solar winters per wind
               winter (solar_common.pairings). Threshold matched in the decoupled base. ln R_both =
               ln R_wind + ln R_solar + interaction. The decoupling removes the wind-solar covariance in
               all four cells alike; the coupled ratio on the same aligned series is given alongside.
  climatology  the solar (or wind) part replaced in BOTH periods by the control calendar-slot
               climatology (31-day running mean per 6 h slot), threshold re-matched in that control.

Usage: python3 solar_region_events.py      (module python3/2025.01; about 10 min on a login node)
Output: dunkelflaute/results/regions/solar_region_events_summary.json, solar_region_events_counts.npz
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
import region_events as ev                                               # noqa: E402
import region_per_K as rpk                                               # noqa: E402
import solar_common as sc                                                # noqa: E402
from robustness_mockert_events import trailing_mean                      # noqa: E402

OUT = f"{_ROOT}/dunkelflaute/results/regions"
WIN = 8
TARGET = sc.TARGET
W_SOLAR = 0.44
NBOOT = 5000
RNG = np.random.default_rng(11)
COMP = {k: v for k, v in ev.COMPARISONS.items() if k != "eerie_drift"}
MAIN = ["9km", "31km", "100km", "eerie_ctl", "eerie_hist"]
KINDS = ("wind", "comb", "solar", "solar_ds")
R = rg.REGIONS
NULL = ("tco95_ctl1950g", 1953, 2099)
NULL2 = ("eerie_ctl", 1972, 2050)


def load(name):
    d = ev.load(name)
    if "di" not in d:
        d["di"], d["feb29"] = sc.day_index(d["time"])
        t1 = d["time"] - np.timedelta64(1, "s")
        d["slot"] = ((t1 - t1.astype("datetime64[D]")) // np.timedelta64(1, "h")).astype(int) // 6
        d["month"] = (t1.astype("datetime64[M]").astype(int) % 12) + 1
    return d


def wmasks(d, y0, y1):
    return [m for _, m in ev.winters(d, y0, y1)]


def series(d, kind, clim):
    if kind in ("wind", "comb"):
        return ev.metric(d, kind)
    if kind == "solar":
        return d["solar"]
    return sc.deseason(d["solar"], d["di"], clim)


def matched(x, ms):
    """As region_events.matched_threshold (400 levels, lowest reaching TARGET), for a solar series of the
    control run. The level grid spans the 0.1st-40th percentiles of the 48 h means of the control winters,
    not of the raw series: half of the raw 6-hourly solar values are night zeros, so its 40th percentile
    is 0."""
    roll = np.concatenate([trailing_mean(x[m], WIN)[WIN - 1:] for m in ms])
    grid = np.linspace(*np.nanpercentile(roll, [0.1, 40]), 400)
    c = np.stack([ev.counts_grid(x[m], grid) for m in ms]).mean(0)
    k = int(np.argmax(c >= TARGET)) if (c >= TARGET).any() else int(np.argmax(c))
    return float(grid[k]), float(c[k])


def counts(x, ms, thr):
    return np.array([ev.counts_grid(x[m], np.array([thr]))[0] for m in ms], float)


def poisson(a, b):
    tot = a.sum() + b.sum()
    return float(stats.binomtest(int(b.sum()), int(tot), len(b) / (len(a) + len(b))).pvalue) if tot else np.nan


def month_share(x, ms, months, thr):
    """Share of flagged (event) time per calendar month (Sep..Mar) in the given winters."""
    f = np.concatenate([sc.flagged(x[m], thr, WIN) for m in ms])
    mo = np.concatenate([months[m] for m in ms])
    tot = max(f.sum(), 1)
    return {int(k): float(f[mo == k].sum() / tot) for k in (9, 10, 11, 12, 1, 2, 3)}


def aligned(d, x, ms):
    """(W, NDAY*4) series aligned by calendar slot; 29 Feb dropped; missing slots NaN."""
    out = np.full((len(ms), sc.NDAY * 4), np.nan)
    for k, m in enumerate(ms):
        mm = m & ~d["feb29"] & (d["di"] >= 0)
        out[k, d["di"][mm] * 4 + d["slot"][mm]] = x[mm]
    return out


def slot_clim(A):
    """Calendar-slot climatology of aligned (W, L) series, 31-day running mean per 6 h slot."""
    raw = np.nanmean(A, axis=0).reshape(sc.NDAY, 4)
    sm = np.empty_like(raw)
    for k in range(sc.NDAY):
        lo, hi = max(0, k - sc.CLIM_HALF), min(sc.NDAY, k + sc.CLIM_HALF + 1)
        sm[k] = np.nanmean(raw[lo:hi], axis=0)
    return sm.reshape(-1)


def attribution(Wc, Sc_, Wf, Sf):
    """Factorial and climatology attribution of the combined change for one region (aligned arrays)."""
    nc, nf = len(Wc), len(Wf)
    ws = 1 - W_SOLAR
    pcc, pfc = sc.pairings(nc, nc, True), sc.pairings(nf, nc, False)
    pcf, pff = sc.pairings(nc, nf, False), sc.pairings(nf, nf, True)
    base = ws * Wc[pcc[:, 0]] + W_SOLAR * Sc_[pcc[:, 1]]
    thr, rate = sc.matched_rows(base, WIN)
    E = lambda X: float(sc.counts_rows(X, thr, WIN).mean())                         # noqa: E731
    e_cc = E(base)
    e_fc = E(ws * Wf[pfc[:, 0]] + W_SOLAR * Sc_[pfc[:, 1]])
    e_cf = E(ws * Wc[pcf[:, 0]] + W_SOLAR * Sf[pcf[:, 1]])
    e_ff = E(ws * Wf[pff[:, 0]] + W_SOLAR * Sf[pff[:, 1]])
    # coupled, same aligned series
    cc = ws * Wc + W_SOLAR * Sc_
    thr_c, rate_c = sc.matched_rows(cc, WIN)
    c_c = float(sc.counts_rows(cc, thr_c, WIN).mean())
    c_f = float(sc.counts_rows(ws * Wf + W_SOLAR * Sf, thr_c, WIN).mean())
    # climatology variants
    scl, wcl = slot_clim(Sc_), slot_clim(Wc)
    xs_c, xs_f = ws * Wc + W_SOLAR * scl[None], ws * Wf + W_SOLAR * scl[None]
    t1, r1 = sc.matched_rows(xs_c, WIN)
    xw_c, xw_f = ws * wcl[None] + W_SOLAR * Sc_, ws * wcl[None] + W_SOLAR * Sf
    t2, r2 = sc.matched_rows(xw_c, WIN)
    lg = lambda a, b: float(np.log(b / a)) if a > 0 and b > 0 else np.nan            # noqa: E731
    return dict(base_rate=rate, thr=thr, ln_wind=lg(e_cc, e_fc), ln_solar=lg(e_cc, e_cf), ln_both=lg(e_cc, e_ff),
                coupled_rate=rate_c, ln_coupled=lg(c_c, c_f),
                clim_solarfixed_rate=r1, ln_clim_solarfixed=lg(sc.counts_rows(xs_c, t1, WIN).mean(), sc.counts_rows(xs_f, t1, WIN).mean()),
                clim_windfixed_rate=r2, ln_clim_windfixed=lg(sc.counts_rows(xw_c, t2, WIN).mean(), sc.counts_rows(xw_f, t2, WIN).mean()))


def comb_event_diagnostics(dc, dfu, mc, mf, i, thr):
    """For the control's combined events (48 h mean of comb below thr): share of event time in Nov-Jan,
    mean wind share of the combined 48 h mean (0.56 W / comb), the median percentile rank of the wind 48 h
    mean within the control's wind 48 h means, and the change of the wind 48 h mean quantile at that rank
    and at the 5th percentile (future against control)."""
    rw_c, rs_c, fl, mo = [], [], [], []
    for m in mc:
        w, s = dc["wind"][m, i], dc["solar"][m, i]
        rw_c.append(trailing_mean(w, WIN))
        rs_c.append(trailing_mean(s, WIN))
        x = (1 - W_SOLAR) * w + W_SOLAR * s
        roll = trailing_mean(x, WIN)
        fl.append(np.nan_to_num(roll, nan=np.inf) < thr)          # the 48 h means below threshold
        mo.append(dc["month"][m])
    rw_c, rs_c, fl, mo = map(np.concatenate, (rw_c, rs_c, fl, mo))
    rw_f = np.concatenate([trailing_mean(dfu["wind"][m, i], WIN) for m in mf])
    okc, okf = np.isfinite(rw_c), np.isfinite(rw_f)
    srt = np.sort(rw_c[okc])
    rank = np.searchsorted(srt, rw_c[fl]) / len(srt)
    q = float(np.median(rank)) if fl.any() else np.nan
    qc = lambda p: np.quantile(rw_c[okc], p)                    # noqa: E731
    qf = lambda p: np.quantile(rw_f[okf], p)                    # noqa: E731
    wsh = (1 - W_SOLAR) * rw_c[fl] / ((1 - W_SOLAR) * rw_c[fl] + W_SOLAR * rs_c[fl])
    return dict(ctl_event_time_NDJ=float(np.isin(mo[fl], (11, 12, 1)).mean()) if fl.any() else np.nan,
                wind_share_of_comb_at_events=float(np.mean(wsh)) if fl.any() else np.nan,
                wind_rank_at_events_median=q,
                wind_q_change_pct_at_event_rank=float(100 * (qf(q) / qc(q) - 1)) if np.isfinite(q) else np.nan,
                wind_q05_change_pct=float(100 * (qf(0.05) / qc(0.05) - 1)),
                wind_mean_change_pct=float(100 * (rw_f[okf].mean() / rw_c[okc].mean() - 1)))


def null_counts(src, y0, y1):
    d = load(src)
    ms = wmasks(d, y0, y1)
    allm = np.any(np.stack(ms), axis=0)
    clim = sc.climatology(d["solar"], d["di"], allm)
    out = {}
    for kind in KINDS:
        X = series(d, kind, clim)
        cols = []
        for i in range(len(R)):
            thr, _ = matched(X[:, i], ms) if kind not in ("wind", "comb") else (ev.matched_threshold(d, kind, i, y0, y1), None)
            cols.append(counts(X[:, i], ms, thr))
        out[kind] = np.stack(cols, axis=1)
    return out


def null_field(c, na0, nb, valid):
    W = len(c)
    na = na0 if na0 + nb <= W - 40 else max(W - nb - 40, nb)
    pairs = ev.block_pairs(W, na, nb)
    lr = ev.block_logratios(c.astype(float), na, nb, pairs)
    return np.nanmean(lr[:, valid], axis=1), lr, na


def main():
    S0 = json.load(open(f"{OUT}/region_events_summary.json"))
    T = rpk.global_t()
    out = dict(regions=R, target=TARGET, window_steps=WIN, kinds=KINDS, clim_floor=sc.CLIM_FLOOR, comparisons={},
               valid_regions={}, field_per_K={}, attribution={}, nulls={})
    store, attr = {}, {}
    for name, ((cn, c0, c1), (fn, f0, f1)) in COMP.items():
        dc, dfu = load(cn), load(fn)
        mc, mf = wmasks(dc, c0, c1), wmasks(dfu, f0, f1)
        allc = np.any(np.stack(mc), axis=0)
        clim = sc.climatology(dc["solar"], dc["di"], allc)
        dT = rpk.mean_years(T, fn, f0, f1)[0] - rpk.mean_years(T, cn, c0, c1)[0]
        res = dict(dT=dT, n_winters=[len(mc), len(mf)], clim_min=dict(zip(R, np.nanmin(clim, axis=0).tolist())))
        for kind in KINDS:
            Xc, Xf = series(dc, kind, clim), series(dfu, kind, clim)
            rows = {}
            for i, r in enumerate(R):
                if kind in ("wind", "comb"):
                    thr = ev.matched_threshold(dc, kind, i, c0, c1)
                    rate = None
                else:
                    thr, rate = matched(Xc[:, i], mc)
                a, b = counts(Xc[:, i], mc, thr), counts(Xf[:, i], mf, thr)
                xc = np.concatenate([Xc[m, i] for m in mc])
                xf = np.concatenate([Xf[m, i] for m in mf])
                mcm, mfm = xc.mean(), xf.mean()
                badj = counts(Xf[:, i] * mcm / mfm, mf, thr)
                rc = np.concatenate([trailing_mean(Xc[m, i], WIN)[WIN - 1:] for m in mc])
                rf = np.concatenate([trailing_mean(Xf[m, i], WIN)[WIN - 1:] for m in mf])
                row = dict(threshold=thr, control=float(a.mean()), future=float(b.mean()),
                           ratio=float(b.mean() / a.mean()) if a.mean() else np.nan, poisson_p=poisson(a, b),
                           ratio_meanremoved=float(badj.mean() / a.mean()) if a.mean() else np.nan,
                           mean_change_pct=float(100 * (mfm / mcm - 1)),
                           p05_48h_change_pct=float(100 * (np.percentile(rf, 5) / np.percentile(rc, 5) - 1)))
                if kind in ("wind", "comb"):
                    ref = S0["comparisons"][name][kind][r]
                    row["check_ratio_vs_region_events"] = float(row["ratio"] - ref["ratio"])
                if kind in ("solar", "solar_ds", "wind"):
                    row["ctl_event_time_by_month"] = month_share(Xc[:, i], mc, dc["month"], thr)
                if kind == "solar":
                    mo_c = np.concatenate([dc["month"][m] for m in mc])
                    mo_f = np.concatenate([dfu["month"][m] for m in mf])
                    row["monthly_mean_change_pct"] = {int(k): float(100 * (xf[mo_f == k].mean() / xc[mo_c == k].mean() - 1))
                                                      for k in (9, 10, 11, 12, 1, 2, 3)}
                    row["mean_ctl"] = float(mcm)
                rows[r] = row
                store[(name, kind, r)] = (a, b)
            res[kind] = rows
        # where combined events sit in the wind and solar distributions (control), and how wind moves there
        for i, r in enumerate(R):
            res["comb"][r].update(comb_event_diagnostics(dc, dfu, mc, mf, i, res["comb"][r]["threshold"]))
        # attribution of the combined change
        if name in MAIN + ["31km_9kmspans", "100km_9kmspans"]:
            A = {}
            for i, r in enumerate(R):
                Wc_, Sc_ = aligned(dc, dc["wind"][:, i], mc), aligned(dc, dc["solar"][:, i], mc)
                Wf_, Sf_ = aligned(dfu, dfu["wind"][:, i], mf), aligned(dfu, dfu["solar"][:, i], mf)
                A[r] = attribution(Wc_, Sc_, Wf_, Sf_)
            attr[name] = A
        out["comparisons"][name] = res
        print(f"  {name}: dT {dT:.2f}, winters {len(mc)}/{len(mf)}", flush=True)
        ev._CACHE.clear()

    # long controls: validity and nulls
    NC = {"tco95_ctl": null_counts(*NULL)}
    ev._CACHE.clear()
    NC["eerie_ctl"] = null_counts(*NULL2)
    ev._CACHE.clear()
    for kind in KINDS:
        v = np.ones(len(R), bool)
        for name in COMP:
            v &= np.array([out["comparisons"][name][kind][r]["control"] >= 0.9 * TARGET for r in R])
        for nn in NC:
            v &= NC[nn][kind].mean(0) >= 0.9 * TARGET
        out["valid_regions"][kind] = [r for r, x in zip(R, v) if x]

    # field means per K with bootstrap and block null
    for name in COMP:
        res = out["comparisons"][name]
        dT = res["dT"]
        pk = lambda x: float(100 * (np.exp(x / dT) - 1))                    # noqa: E731
        out["field_per_K"][name] = {}
        for kind in KINDS:
            V = out["valid_regions"][kind]
            vm = np.array([r in V for r in R])
            a = np.stack([store[(name, kind, r)][0] for r in V], axis=1)
            b = np.stack([store[(name, kind, r)][1] for r in V], axis=1)
            with np.errstate(divide="ignore"):
                lr = np.log(b.mean(0) / a.mean(0))
            obs = float(np.mean(lr[np.isfinite(lr)]))
            boots = np.empty(NBOOT)
            for k in range(NBOOT):
                ia, ib = RNG.integers(0, len(a), len(a)), RNG.integers(0, len(b), len(b))
                with np.errstate(divide="ignore", invalid="ignore"):
                    x = np.log(b[ib].mean(0) / a[ia].mean(0))
                boots[k] = np.mean(x[np.isfinite(x)])
            lo, hi = np.percentile(boots, [5, 95])
            nm, lrn, na = null_field(NC["tco95_ctl"][kind], len(a), len(b), vm)
            # per-region two-sided block-null p
            obs_r = np.array([np.log(res[kind][r]["ratio"]) for r in R])
            p_up = (1 + (lrn >= obs_r[None]).sum(0)) / (1 + len(lrn))
            p_dn = (1 + (lrn <= obs_r[None]).sum(0)) / (1 + len(lrn))
            p2 = np.minimum(1, 2 * np.minimum(p_up, p_dn))
            for j, r in enumerate(R):
                res[kind][r]["block_null_p2_tco95"] = float(p2[j])
            mean_ch = float(np.mean([res[kind][r]["mean_change_pct"] for r in V]))
            lam = np.array([np.log(res[kind][r]["ratio_meanremoved"]) for r in V])
            out["field_per_K"][name][kind] = dict(
                n_regions=len(V), mean_ratio=float(np.exp(obs)), pct_per_K=pk(obs), pct_per_K_p05_p95=[pk(lo), pk(hi)],
                ratio_p05_p95=[float(np.exp(lo)), float(np.exp(hi))], n_regions_up=int((lr > 0).sum()),
                n_regions_sig_up=int(sum(res[kind][r]["block_null_p2_tco95"] < 0.05 and res[kind][r]["ratio"] > 1 for r in V)),
                n_regions_sig_down=int(sum(res[kind][r]["block_null_p2_tco95"] < 0.05 and res[kind][r]["ratio"] < 1 for r in V)),
                null_field_ratio_p05_p95=[float(np.exp(x)) for x in np.percentile(nm, [5, 95])], null_n_control=na,
                p_field_up=float((1 + (nm >= obs).sum()) / (1 + len(nm))), p_field_down=float((1 + (nm <= obs).sum()) / (1 + len(nm))),
                mean_change_pct=mean_ch, mean_change_per_K=mean_ch / dT,
                p05_48h_change_pct=float(np.mean([res[kind][r]["p05_48h_change_pct"] for r in V])),
                pct_per_K_meanremoved=pk(float(lam[np.isfinite(lam)].mean())),
                share_from_mean=float(1 - lam[np.isfinite(lam)].mean() / obs) if obs else np.nan)
    out["nulls"] = {nn: {k: float(v.mean()) for k, v in c.items()} for nn, c in NC.items()}

    # attribution field means (comb-valid regions with a decoupled base that reaches the rate)
    for name, A in attr.items():
        dT = out["comparisons"][name]["dT"]
        pk = lambda x: float(100 * (np.exp(x / dT) - 1))                    # noqa: E731
        V = [r for r in S0["valid_regions"]["comb"] if A[r]["base_rate"] >= 0.9 * TARGET]
        f = {k: float(np.nanmean([A[r][k] for r in V])) for k in ("ln_wind", "ln_solar", "ln_both", "ln_coupled",
                                                                    "ln_clim_solarfixed", "ln_clim_windfixed")}
        wind_alone = float(np.mean([np.log(out["comparisons"][name]["wind"][r]["ratio"]) for r in V]))
        solar_alone = float(np.mean([np.log(out["comparisons"][name]["solar"][r]["ratio"]) for r in V]))
        comb_real = float(np.mean([np.log(out["comparisons"][name]["comb"][r]["ratio"]) for r in V]))
        out["attribution"][name] = dict(
            n_regions=len(V), regions=V, per_region=A,
            pct_per_K={k[3:]: pk(v) for k, v in f.items()},
            pct_per_K_interaction=pk(f["ln_both"] - f["ln_wind"] - f["ln_solar"]),
            share_wind=f["ln_wind"] / f["ln_both"], share_solar=f["ln_solar"] / f["ln_both"],
            pct_per_K_same_regions=dict(wind_alone=pk(wind_alone), solar_alone=pk(solar_alone), comb=pk(comb_real)),
            comb_events={k: float(np.nanmean([out["comparisons"][name]["comb"][r][k] for r in V]))
                         for k in ("ctl_event_time_NDJ", "wind_share_of_comb_at_events", "wind_rank_at_events_median",
                                   "wind_q_change_pct_at_event_rank", "wind_q05_change_pct", "wind_mean_change_pct")})

    json.dump(out, open(f"{OUT}/solar_region_events_summary.json", "w"), indent=1)
    np.savez(f"{OUT}/solar_region_events_counts.npz",
             **{f"{n}|{k}|{r}|{p}": v[j] for (n, k, r), v in store.items() for j, p in enumerate(("ctl", "fut"))})

    # print
    for kind in KINDS:
        print(f"\n=== {kind}: valid {len(out['valid_regions'][kind])}: {' '.join(out['valid_regions'][kind])}")
        for n, f in out["field_per_K"].items():
            x = f[kind]
            print(f"  {n:15s} ratio {x['mean_ratio']:.3f} {x['pct_per_K']:+5.1f} %/K [{x['pct_per_K_p05_p95'][0]:+5.1f},{x['pct_per_K_p05_p95'][1]:+5.1f}] "
                  f"up {x['n_regions_up']:2d}/{x['n_regions']} sig+ {x['n_regions_sig_up']} sig- {x['n_regions_sig_down']}  "
                  f"null {x['null_field_ratio_p05_p95'][0]:.3f}-{x['null_field_ratio_p05_p95'][1]:.3f} p_up {x['p_field_up']:.4f} p_dn {x['p_field_down']:.4f}  "
                  f"mean {x['mean_change_pct']:+.1f}% ({x['mean_change_per_K']:+.2f} %/K) p05 {x['p05_48h_change_pct']:+.1f}%  "
                  f"meanremoved {x['pct_per_K_meanremoved']:+.1f} share {x['share_from_mean']:.0%}")
    chk = [abs(out["comparisons"][n][k][r]["check_ratio_vs_region_events"]) for n in COMP for k in ("wind", "comb") for r in R
           if np.isfinite(out["comparisons"][n][k][r]["check_ratio_vs_region_events"])]
    print(f"\ncheck wind/comb ratios against region_events_summary.json: max |diff| {max(chk):.2e}")
    print("\n=== attribution of comb (%/K)")
    for n, a in out["attribution"].items():
        p = a["pct_per_K"]
        print(f"  {n:15s} {a['n_regions']} regions: coupled {p['coupled']:+.1f}, decoupled both {p['both']:+.1f} = wind {p['wind']:+.1f} + solar "
              f"{p['solar']:+.1f} + int {a['pct_per_K_interaction']:+.1f}; share wind {a['share_wind']:.0%}; clim: solar fixed {p['clim_solarfixed']:+.1f}, "
              f"wind fixed {p['clim_windfixed']:+.1f}; same regions: wind alone {a['pct_per_K_same_regions']['wind_alone']:+.1f}, "
              f"solar alone {a['pct_per_K_same_regions']['solar_alone']:+.1f}, comb {a['pct_per_K_same_regions']['comb']:+.1f}")
        print(f"      comb events: NDJ {a['comb_events']['ctl_event_time_NDJ']:.2f}, wind share {a['comb_events']['wind_share_of_comb_at_events']:.2f}, "
              f"wind rank {a['comb_events']['wind_rank_at_events_median']:.2f}, wind quantile change there {a['comb_events']['wind_q_change_pct_at_event_rank']:+.1f}% "
              f"(at 5 %: {a['comb_events']['wind_q05_change_pct']:+.1f}%, mean {a['comb_events']['wind_mean_change_pct']:+.1f}%)")
    print("\n=== control event time in Nov-Jan (share), field mean; uniform would be 0.43")
    for kind in ("wind", "solar", "solar_ds"):
        for n in MAIN:
            sh = [sum(out["comparisons"][n][kind][r]["ctl_event_time_by_month"][m] for m in (11, 12, 1)) for r in out["valid_regions"][kind]]
            dj = [sum(out["comparisons"][n][kind][r]["ctl_event_time_by_month"][m] for m in (12, 1)) for r in out["valid_regions"][kind]]
            print(f"  {kind:8s} {n:10s} NDJ {np.mean(sh):.2f} (min {np.min(sh):.2f}, max {np.max(sh):.2f})  DJ {np.mean(dj):.2f}")
    print(f"Saved {OUT}/solar_region_events_summary.json")


if __name__ == "__main__":
    main()
