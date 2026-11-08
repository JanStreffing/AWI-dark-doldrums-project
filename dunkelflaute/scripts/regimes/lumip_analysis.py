"""
LUMIP land-use test of the CMIP6 at-fixed-G transfer term (methods notes 13e): the winter 10 m wind
difference between experiments that differ only in land use, per unit land-cover change, and the
part of the CMIP6 -2.0 % (8 European regions) it explains.

Pairs, land-use effect = A minus B:
  hist    historical  - hist-noLu              late 1985-2014, early 1951-1980 (and the change
                                                early -> late in each, difference in differences)
  ssp370  ssp370      - ssp370-ssp126Lu        SSP3 minus SSP1 land use, late 2071-2100, early 2016-2035
  ssp126  ssp126-ssp370Lu - ssp126             SSP3 minus SSP1 land use, same windows
  esm585  esm-ssp585  - esm-ssp585-ssp126Lu    SSP5 minus SSP1 land use (emission driven), same windows
Inputs: lumip_parts/lumip_mon_*.npz (lumip_extract_monthly.py), lumip_day_*.npz (lumip_extract_daily.py),
results/regimes/stability_cmip6.json (CMIP6 at-fixed-G term and ssp585 land-cover changes).

Winter = Sep-Mar, equal month weights, labelled by its January. Region aggregates are plain means
over the 8 European regions (eu8) or the 22 regions with a southern edge at or north of 29 N (c22),
as in stability_cmip6_analysis.py. Wind and G differences in %, land cover in percentage points of
cell area (lai in m2 m-2), tas in K.

Uncertainty of a difference of ensemble means (n_A, n_B members):
  se_member   sigma sqrt(1/n_A + 1/n_B), sigma = spread across members of the single-member window
              mean (field mean, %), pooled over both experiments (ddof 1); 90 % range +/- 1.645 se
  boot        1000 draws, winters resampled with replacement in each experiment (same draw for all
              members of an experiment); 5-95 %
  same_exp    differences between disjoint member groups of the same experiment (group sizes of the
              pair where possible, else 1 vs 1): what internal variability alone gives
Circulation control (monthly G of monthly-mean psl, only as a control): per region, the slope b of
monthly ln s on ln G anomalies (member x calendar month climatology removed, both experiments
pooled); at-fixed-circulation part = 100 (ln s_A/s_B - b ln G_A/G_B). Also s/G (b = 1).

Daily (where day/sfcWind exists in both experiments): mean wind CF change and the Mockert 48 h event
frequency (daily_events.py footing: WIN = 2, threshold matched to 5 events per winter in experiment
B, members pooled), ratio A/B, field mean of ln ratio; mean-removed ratio as in daily_events.py.

Usage: python3 lumip_analysis.py       (legacy python, login node, about a minute)
Output: dunkelflaute/results/regimes/lumip_landuse.json, lumip_landuse.png
"""
import glob
import itertools
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("regions", "analysis", "cmip6"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/{p}")
import regions as rg                     # noqa: E402
import region_events as ev               # noqa: E402

ev.WIN = 2
PARTS = f"{_ROOT}/dunkelflaute/data/regimes/lumip_parts"
RES = f"{_ROOT}/dunkelflaute/results/regimes"
R = rg.REGIONS
EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
BOX = dict(rg.BOXES, DE=rg.DE_BOX)
C22 = [r for r in R if BOX[r][0] >= 29.0]
SETS = {"eu8": [R.index(r) for r in EU], "c22": [R.index(r) for r in C22]}
PAIRS = {"hist": ("historical", "hist-noLu", (1985, 2014), (1951, 1980)),
         "ssp370": ("ssp370", "ssp370-ssp126Lu", (2071, 2100), (2016, 2035)),
         "ssp126": ("ssp126-ssp370Lu", "ssp126", (2071, 2100), (2016, 2035)),
         "esm585": ("esm-ssp585", "esm-ssp585-ssp126Lu", (2071, 2100), (2016, 2035))}
LC = ("cropFrac", "treeFrac", "grassFrac", "pastureFrac", "lai")
NB = 1000
TARGET = 5.0
RNG = np.random.default_rng(7)


def load_mon(model, exp):
    f = f"{PARTS}/lumip_mon_{model}_{exp}.npz"
    if not os.path.exists(f):
        return None
    z = np.load(f)
    mems = [str(m) for m in z["members"]]
    out = dict(winters=z["winters"], members=mems)
    for v in ("sfcWind", "G", "tas") + LC:
        if all(f"{m}|{v}" in z.files for m in mems):
            out[v] = np.stack([z[f"{m}|{v}"].astype(np.float64) for m in mems])     # (M, W, 7, R) or (M, W, R)
        else:
            have = [m for m in mems if f"{m}|{v}" in z.files]
            if have:
                out[v] = np.stack([z[f"{m}|{v}"].astype(np.float64) for m in have])
                out[f"{v}_members"] = have
    return out


def winter_mean(x):
    return x.mean(axis=2) if x.ndim == 4 else x                 # (M, W, R), NaN if any month missing


def window(d, v, w):
    """(M, Wsel, R) winter means in window w, winters complete in every member."""
    x = winter_mean(d[v])
    sel = (d["winters"] >= w[0]) & (d["winters"] <= w[1])
    x = x[:, sel]
    ok = np.isfinite(x[:, :, SETS["eu8"]]).all(axis=(0, 2))
    return x[:, ok]


def fm(pct, s):
    return float(np.nanmean(pct[SETS[s]]))


def pct_diff(a, b):
    return 100 * (a / b - 1)


def diff_stats(xa, xb, rel=True):
    """xa (Ma, Wa, R), xb (Mb, Wb, R) -> dict per set: d, se_member, ci90, boot, same_exp; plus per-region d."""
    A, B = np.nanmean(xa, axis=(0, 1)), np.nanmean(xb, axis=(0, 1))
    d = pct_diff(A, B) if rel else A - B
    out = {"per_region": {r: float(d[i]) for i, r in enumerate(R)}}
    ma, mb = xa.shape[0], xb.shape[0]
    ref = 0.5 * (A + B)
    # single-member window means as field-mean anomalies
    mem_a = np.nanmean(xa, axis=1)
    mem_b = np.nanmean(xb, axis=1)
    # bootstrap over winters
    ia = RNG.integers(0, xa.shape[1], (NB, xa.shape[1]))
    ib = RNG.integers(0, xb.shape[1], (NB, xb.shape[1]))
    Aw, Bw = np.nanmean(xa, axis=0), np.nanmean(xb, axis=0)           # (W, R) ensemble means per winter
    bA = Aw[ia].mean(axis=1)
    bB = Bw[ib].mean(axis=1)
    bd = pct_diff(bA, bB) if rel else bA - bB
    for s in SETS:
        def f_(v):
            return pct_diff(v, ref) if rel else v - ref
        fa = np.array([fm(f_(m), s) for m in mem_a])
        fb = np.array([fm(f_(m), s) for m in mem_b])
        dev = np.concatenate([fa - fa.mean(), fb - fb.mean()]) if (ma > 1 or mb > 1) else np.array([])
        dof = (ma - 1) + (mb - 1)
        sigma = float(np.sqrt((dev ** 2).sum() / dof)) if dof > 0 else np.nan
        se = sigma * np.sqrt(1 / ma + 1 / mb) if np.isfinite(sigma) else np.nan
        # interannual fallback: sd of field-mean winter anomalies of the ensemble mean / sqrt(W)
        wa_ = np.array([fm(f_(x), s) for x in Aw])
        wb_ = np.array([fm(f_(x), s) for x in Bw])
        se_w = float(np.sqrt(wa_.var(ddof=1) / len(wa_) + wb_.var(ddof=1) / len(wb_)))
        se_use = se if np.isfinite(se) else se_w
        o = dict(d=fm(d, s), sigma_member=sigma, se_member=float(se), se_winters=se_w, se=float(se_use),
                 se_source="member" if np.isfinite(se) else "winters",
                 ci90_member=[fm(d, s) - 1.645 * se_use, fm(d, s) + 1.645 * se_use],
                 boot_p05_p95=[float(q) for q in np.percentile(np.nanmean(bd[:, SETS[s]], axis=1), [5, 95])],
                 n_A=ma, n_B=mb, n_winters=[int(xa.shape[1]), int(xb.shape[1])])
        # same-experiment differences: disjoint member groups inside one experiment, sizes as close to the pair's as fit
        same = []
        nmin, nmax = min(ma, mb), max(ma, mb)
        for fx, mx in ((fa, ma), (fb, mb)):
            if mx < 2:
                continue
            k1 = max(1, min(nmin, mx // 2))
            k2 = min(nmax, mx - k1)
            for _ in range(300):
                p = RNG.permutation(mx)
                same.append(fx[p[:k1]].mean() - fx[p[k1:k1 + k2]].mean())
        if same:
            same = np.array(same)
            o["same_exp"] = dict(n=len(same), sd=float(same.std()), p05_p95=[float(q) for q in np.percentile(same, [5, 95])])
        if ma >= 2:
            o["two_member_diff_A"] = float(fa[0] - fa[1])
        if mb >= 2:
            o["two_member_diff_B"] = float(fb[0] - fb[1])
        out[s] = o
    return out


def circulation(da, db, wa, wb):
    """Monthly ln s on ln G slope per region (anomalies), and the at-fixed-circulation part."""
    def anoms(d, w):
        sel = (d["winters"] >= w[0]) & (d["winters"] <= w[1])
        ls, lg = np.log(d["sfcWind"][:, sel]), np.log(d["G"][:, sel])              # (M, W, 7, R)
        return ls - np.nanmean(ls, axis=1, keepdims=True), lg - np.nanmean(lg, axis=1, keepdims=True)
    sa, ga = anoms(da, wa)
    sb, gb = anoms(db, wb)
    s = np.concatenate([sa.reshape(-1, len(R)), sb.reshape(-1, len(R))])
    g = np.concatenate([ga.reshape(-1, len(R)), gb.reshape(-1, len(R))])
    ok = np.isfinite(s) & np.isfinite(g)
    b = np.array([np.sum(s[ok[:, i], i] * g[ok[:, i], i]) / np.sum(g[ok[:, i], i] ** 2) for i in range(len(R))])
    r = np.array([np.corrcoef(s[ok[:, i], i], g[ok[:, i], i])[0, 1] for i in range(len(R))])
    return b, r


def pair_monthly(model, pair):
    ea, eb, wl, we = PAIRS[pair]
    da, db = load_mon(model, ea), load_mon(model, eb)
    if da is None or db is None or "sfcWind" not in da or "sfcWind" not in db:
        return None
    out = dict(A=ea, B=eb, members_A=da["members"], members_B=db["members"], windows=dict(late=wl, early=we))
    for wn, w in (("late", wl), ("early", we)):
        xa, xb = window(da, "sfcWind", w), window(db, "sfcWind", w)
        if xa.shape[1] < 10 or xb.shape[1] < 10:
            continue
        o = {"sfcWind": diff_stats(xa, xb)}
        if "G" in da and "G" in db and "G_members" not in da and "G_members" not in db:
            o["G"] = diff_stats(window(da, "G", w), window(db, "G", w))
            b, rr = circulation(da, db, w, w)
            ds = np.array([o["sfcWind"]["per_region"][r] for r in R])
            dg = np.array([o["G"]["per_region"][r] for r in R])
            fixed = 100 * (np.exp(np.log1p(ds / 100) - b * np.log1p(dg / 100)) - 1)
            transfer = 100 * ((1 + ds / 100) / (1 + dg / 100) - 1)
            o["circulation"] = {s: dict(slope_b_mean=float(np.mean(b[SETS[s]])), corr_mean=float(np.mean(rr[SETS[s]])),
                                        d_fixed_circ=fm(fixed, s), d_transfer_sG=fm(transfer, s),
                                        circ_part=fm(ds, s) - fm(fixed, s)) for s in SETS}
            o["circulation"]["per_region_b"] = {r: float(b[i]) for i, r in enumerate(R)}
        if "tas" in da and "tas" in db and "tas_members" not in da and "tas_members" not in db:
            o["tas"] = diff_stats(window(da, "tas", w), window(db, "tas", w), rel=False)
        for v in LC:
            if v in da and v in db:
                xa_, xb_ = window(da, v, w) if f"{v}_members" not in da else None, window(db, v, w) if f"{v}_members" not in db else None
                if xa_ is None or xb_ is None:
                    xa_ = winter_mean(da[v])[:, (da["winters"] >= w[0]) & (da["winters"] <= w[1])]
                    xb_ = winter_mean(db[v])[:, (db["winters"] >= w[0]) & (db["winters"] <= w[1])]
                A, B = np.nanmean(xa_, axis=(0, 1)), np.nanmean(xb_, axis=(0, 1))
                o[v] = {s: fm(A - B, s) for s in SETS}
                o[v]["per_region"] = {r: float((A - B)[i]) for i, r in enumerate(R)}
                o[v]["ctl_B"] = {s: fm(B, s) for s in SETS}
                if v == "lai":
                    o[v].update({f"{s}_pct": fm(pct_diff(A, B), s) for s in SETS})
        out[wn] = o
    # change early -> late in each experiment (hist pair): difference in differences
    if pair == "hist" and "early" in out:
        dd = {}
        for nm, d in (("A", da), ("B", db)):
            e, l_ = window(d, "sfcWind", we), window(d, "sfcWind", wl)
            ch = pct_diff(np.nanmean(l_, axis=(0, 1)), np.nanmean(e, axis=(0, 1)))
            dd[nm] = {s: fm(ch, s) for s in SETS}
            if "G" in d and "G_members" not in d:
                chg = pct_diff(np.nanmean(window(d, "G", wl), axis=(0, 1)), np.nanmean(window(d, "G", we), axis=(0, 1)))
                dd[nm + "_G"] = {s: fm(chg, s) for s in SETS}
        dd["A_minus_B"] = {s: dd["A"][s] - dd["B"][s] for s in SETS}
        out["change_early_late"] = dd
    return out


# ---------------- daily: CF and events
def load_day(model, exp):
    f = f"{PARTS}/lumip_day_{model}_{exp}.npz"
    if not os.path.exists(f):
        return None
    z = np.load(f)
    return {str(m): dict(ymd=z[f"{m}|ymd"], wind=z[f"{m}|wind"].astype(np.float64), spd10=z[f"{m}|spd10"].astype(np.float64))
            for m in z["members"]}


def member_winters(d):
    """-> list over members of list of (W) arrays per region: [(T_w, R)]"""
    out = []
    for m in d.values():
        w = np.where(m["ymd"][:, 1] >= 9, m["ymd"][:, 0] + 1, m["ymd"][:, 0])
        out.append([m["wind"][w == y] for y in np.unique(w) if (w == y).sum() >= 200])
    return out


def events_pair(da, db):
    A, B = member_winters(da), member_winters(db)
    per = {}
    for i, r in enumerate(R):
        xb = np.concatenate([x[:, i] for mem in B for x in mem])
        if not np.isfinite(xb).any():
            continue
        grid = np.linspace(*np.nanpercentile(xb, [0.1, 40]), 400)
        cb = np.stack([ev.counts_grid(x[:, i], grid) for mem in B for x in mem]).mean(0)
        k = int(np.argmax(cb >= TARGET)) if (cb >= TARGET).any() else int(np.argmax(cb))
        thr = np.array([grid[k]])
        ca = [[ev.counts_grid(x[:, i], thr)[0] for x in mem] for mem in A]
        cbm = [[ev.counts_grid(x[:, i], thr)[0] for x in mem] for mem in B]
        mA, mB = np.nanmean(np.concatenate([x[:, i] for mem in A for x in mem])), np.nanmean(xb)
        cadj = np.array([ev.counts_grid(x[:, i] * mB / mA, thr)[0] for mem in A for x in mem])
        ea = np.mean(np.concatenate([np.array(c, float) for c in ca]))
        eb = np.mean(np.concatenate([np.array(c, float) for c in cbm]))
        per[r] = dict(threshold=float(thr[0]), events_B=float(eb), events_A=float(ea), ratio=float(ea / eb),
                      ratio_meanremoved=float(cadj.mean() / eb), dCF_pct=float(100 * (mA / mB - 1)), valid=bool(eb >= 0.9 * TARGET),
                      memb_events_A=[float(np.mean(c)) for c in ca], memb_events_B=[float(np.mean(c)) for c in cbm])
    out = {"per_region": per}
    for s, idx in SETS.items():
        V = [R[i] for i in idx if R[i] in per and per[R[i]]["valid"]]
        lr = np.array([np.log(per[r]["ratio"]) for r in V])
        la = np.array([np.log(per[r]["ratio_meanremoved"]) for r in V])
        # member-level field-mean ln(events / B ensemble events) -> spread for se
        def memf(key):
            M = len(per[V[0]][key])
            return np.array([np.mean([np.log(max(per[r][key][j], 0.05) / per[r]["events_B"]) for r in V]) for j in range(M)])
        fa, fb = memf("memb_events_A"), memf("memb_events_B")
        dev = np.concatenate([fa - fa.mean(), fb - fb.mean()])
        dof = len(fa) + len(fb) - 2
        sig = float(np.sqrt((dev ** 2).sum() / dof)) if dof > 0 else np.nan
        se = sig * np.sqrt(1 / len(fa) + 1 / len(fb))
        m = float(lr.mean())
        out[s] = dict(n_regions=len(V), mean_ratio=float(np.exp(m)), ratio_ci90_member=[float(np.exp(m - 1.645 * se)), float(np.exp(m + 1.645 * se))],
                      mean_ratio_meanremoved=float(np.exp(la.mean())),
                      dCF_pct=float(np.mean([per[r]["dCF_pct"] for r in V])), n_A=len(fa), n_B=len(fb))
    return out


def daily_part():
    res = {}
    models = sorted({os.path.basename(f)[len("lumip_day_"):-4].split("_")[0] for f in glob.glob(f"{PARTS}/lumip_day_*.npz")})
    for model in models:
        for pair, (ea, eb, wl, we) in PAIRS.items():
            da, db = load_day(model, ea), load_day(model, eb)
            if da is None or db is None:
                continue
            o = events_pair(da, db)
            # mean CF change with member se (same footing as the monthly diff_stats)
            def wm(d, key):
                rows = []
                for m in d.values():
                    w = np.where(m["ymd"][:, 1] >= 9, m["ymd"][:, 0] + 1, m["ymd"][:, 0])
                    rows.append(np.stack([m[key][w == y].mean(0) for y in np.unique(w) if (w == y).sum() >= 200]))
                n = min(len(x) for x in rows)
                return np.stack([x[:n] for x in rows])
            o["CF"] = diff_stats(wm(da, "wind"), wm(db, "wind"))
            o["spd10_daily"] = diff_stats(wm(da, "spd10"), wm(db, "spd10"))
            res[f"{model}|{pair}"] = o
            print(f"  daily {model:14s} {pair:7s} n {len(da)}/{len(db)}  dspd10 {o['spd10_daily']['eu8']['d']:+.2f}%  dCF {o['CF']['eu8']['d']:+.2f}% "
                  f"[{o['CF']['eu8']['ci90_member'][0]:+.2f}, {o['CF']['eu8']['ci90_member'][1]:+.2f}]  events ratio {o['eu8']['mean_ratio']:.3f} "
                  f"[{o['eu8']['ratio_ci90_member'][0]:.3f}, {o['eu8']['ratio_ci90_member'][1]:.3f}]  mean-removed {o['eu8']['mean_ratio_meanremoved']:.3f}"
                  f" | c22 dCF {o['CF']['c22']['d']:+.2f}% ratio {o['c22']['mean_ratio']:.3f}", flush=True)
    return res


# ---------------- per unit land cover and the CMIP6 attribution
def wls_origin(x, y, se):
    """Weighted least squares through the origin: slope, se (from weights), and the scatter-inflated se."""
    x, y, se = map(np.asarray, (x, y, se))
    w = 1 / se ** 2
    b = np.sum(w * x * y) / np.sum(w * x * x)
    se_b = 1 / np.sqrt(np.sum(w * x * x))
    chi2 = np.sum(w * (y - b * x) ** 2) / max(1, len(x) - 1)
    return float(b), float(se_b), float(se_b * np.sqrt(max(1.0, chi2))), float(chi2)


def wls(X, y, se):
    """WLS without intercept: beta, covariance scaled by max(1, reduced chi2), reduced chi2."""
    w = 1 / se ** 2
    XtW = X.T * w
    A = XtW @ X
    beta = np.linalg.solve(A, XtW @ y)
    dof = max(1, len(y) - X.shape[1])
    chi2 = float(np.sum(w * (y - X @ beta) ** 2) / dof)
    return beta, np.linalg.inv(A) * max(1.0, chi2), chi2


def lumip_points(mon, s, preds, windows, ykind="fixed"):
    pts = []
    for k, o in mon.items():
        model, pair = k.split("|")
        for wn in windows:
            if wn not in o or (wn == "early" and pair != "hist"):
                continue
            L = o[wn]
            if not all(v in L for v in preds):
                continue
            key = {"fixed": "d_fixed_circ", "sG": "d_transfer_sG", "raw": "none"}[ykind]
            y = L.get("circulation", {}).get(s, {}).get(key, L["sfcWind"][s]["d"])
            pts.append(dict(key=f"{k}|{wn}", model=model, x=[L[v][s] for v in preds], y=y, se=L["sfcWind"][s]["se"]))
    return pts


def multi_fits(mon, cm):
    """Per-unit land-cover sensitivity from the LUMIP pairs (WLS through the origin, at fixed circulation),
    with a cluster bootstrap over LUMIP models, and the land-use part of each CMIP6 model's at-fixed-G term."""
    out = {}
    ens = {n: m for n, m in cm["models"].items() if not m["wind_only"]}
    for s in SETS:
        g = "eu8" if s == "eu8" else "cmip22"
        wg20 = np.array([m[g]["within_G"] for m in ens.values()])
        out[s] = {}
        for preds, wset, ykind in [(p, w, "fixed") for p in (("cropFrac",), ("treeFrac",), ("cropFrac", "treeFrac")) for w in (("late",), ("late", "early"))] + \
                [(p, ("late",), y) for p in (("cropFrac",), ("treeFrac",), ("cropFrac", "treeFrac")) for y in ("raw", "sG")]:
            if True:
                pts = lumip_points(mon, s, preds, wset, ykind)
                if len(pts) < len(preds) + 2:
                    continue
                X = np.array([p["x"] for p in pts])
                y = np.array([p["y"] for p in pts])
                se = np.array([p["se"] for p in pts])
                beta, cov, chi2 = wls(X, y, se)
                mods = sorted({p["model"] for p in pts})
                pm = np.array([p["model"] for p in pts])
                boot = []
                for _ in range(2000):
                    pick = RNG.choice(mods, len(mods), replace=True)
                    idx = np.concatenate([np.where(pm == mm)[0] for mm in pick])
                    try:
                        boot.append(wls(X[idx], y[idx], se[idx])[0])
                    except np.linalg.LinAlgError:
                        continue
                boot = np.array(boot)
                rows = {n: ([m[g].get(f"d_{v}", np.nan) for v in preds], m[g]["within_G"]) for n, m in ens.items()}
                rows = {n: r for n, r in rows.items() if np.all(np.isfinite(r[0]))}
                DX = np.array([r[0] for r in rows.values()])
                WG = np.array([r[1] for r in rows.values()])
                pred = DX @ beta
                bp = boot @ DX.T                                                   # (B, n_models)
                gp = (beta + RNG.multivariate_normal(np.zeros(len(beta)), cov, 4000)) @ DX.T
                name = "+".join(preds) + "|" + "+".join(wset) + ("" if ykind == "fixed" else f"|{ykind}")
                out[s][name] = dict(n_points=len(pts), lumip_models=mods, beta=dict(zip(preds, beta.tolist())),
                                    beta_se=dict(zip(preds, np.sqrt(np.diag(cov)).tolist())), reduced_chi2=chi2,
                                    beta_boot_p05_p95={v: [float(q) for q in np.percentile(boot[:, i], [5, 95])] for i, v in enumerate(preds)},
                                    n_cmip6=len(rows), cmip6_within_G_subset_mean=float(WG.mean()), cmip6_within_G_subset_median=float(np.median(WG)),
                                    cmip6_within_G_20_median=float(np.median(wg20)),
                                    lu_part_mean=float(pred.mean()), lu_part_median=float(np.median(pred)),
                                    lu_part_mean_p05_p95_boot=[float(q) for q in np.percentile(bp.mean(1), [5, 95])],
                                    lu_part_mean_p05_p95_gauss=[float(q) for q in np.percentile(gp.mean(1), [5, 95])],
                                    share_of_subset_mean=float(pred.mean() / WG.mean()),
                                    share_p05_p95_boot=[float(q) for q in np.percentile(bp.mean(1) / WG.mean(), [5, 95])],
                                    residual_mean=float((WG - pred).mean()), residual_median=float(np.median(WG - pred)),
                                    points={p["key"]: dict(x=p["x"], y=p["y"], se=p["se"]) for p in pts},
                                    per_model={n: dict(d=r[0], within_G=r[1], lu_part=float(np.dot(r[0], beta))) for n, r in rows.items()})
                o = out[s][name]
                print(f"fit {s} {name:28s} n {len(pts):2d} ({len(mods)} models) beta " + ", ".join(f"{v} {b:+.3f} +/- {e:.3f} [boot {o['beta_boot_p05_p95'][v][0]:+.2f}, {o['beta_boot_p05_p95'][v][1]:+.2f}]" for v, b, e in zip(preds, beta, np.sqrt(np.diag(cov))))
                      + f" chi2 {chi2:.1f} | CMIP6 {len(rows)} models: within_G mean {WG.mean():+.2f}, LU part mean {pred.mean():+.2f} [boot {o['lu_part_mean_p05_p95_boot'][0]:+.2f}, {o['lu_part_mean_p05_p95_boot'][1]:+.2f}] "
                        f"median {np.median(pred):+.2f}, share {o['share_of_subset_mean']:.0%} [{o['share_p05_p95_boot'][0]:.0%}, {o['share_p05_p95_boot'][1]:.0%}], residual mean {o['residual_mean']:+.2f}", flush=True)
    return out


def main():
    models =sorted({os.path.basename(f)[len("lumip_mon_"):-4].split("_")[0] for f in glob.glob(f"{PARTS}/lumip_mon_*.npz")})
    mon = {}
    for model in models:
        for pair in PAIRS:
            o = pair_monthly(model, pair)
            if o is None or "late" not in o:
                continue
            mon[f"{model}|{pair}"] = o
            L = o["late"]
            s = L["sfcWind"]["eu8"]
            ci = s["ci90_member"] or [np.nan, np.nan]
            g = L.get("G", {}).get("eu8", {}).get("d", np.nan)
            fc = L.get("circulation", {}).get("eu8", {}).get("d_fixed_circ", np.nan)
            lc = "  ".join(f"{v} {L[v]['eu8']:+.2f}" for v in ("cropFrac", "treeFrac", "grassFrac", "lai") if v in L)
            print(f"{model:14s} {pair:7s} n {s['n_A']:2d}/{s['n_B']:2d}  ds {s['d']:+.2f}% [{ci[0]:+.2f}, {ci[1]:+.2f}] boot [{s['boot_p05_p95'][0]:+.2f}, {s['boot_p05_p95'][1]:+.2f}]"
                  f"  same-exp sd {s.get('same_exp', {}).get('sd', np.nan):.2f}  dG {g:+.2f}%  fixed-circ {fc:+.2f}%  | {lc}"
                  f"  | c22 ds {L['sfcWind']['c22']['d']:+.2f}%", flush=True)
            C = L.get("circulation", {})
            print(f"{'':14s} {'':7s} se {s['se']:.2f} ({s['se_source']}, winters-based {s['se_winters']:.2f})  b {C.get('eu8', {}).get('slope_b_mean', np.nan):.2f} r {C.get('eu8', {}).get('corr_mean', np.nan):.2f}"
                  f"  s/G {C.get('eu8', {}).get('d_transfer_sG', np.nan):+.2f}%  dtas {L.get('tas', {}).get('eu8', {}).get('d', np.nan):+.2f} K"
                  f"  | c22: dG {L.get('G', {}).get('c22', {}).get('d', np.nan):+.2f}% fixed {C.get('c22', {}).get('d_fixed_circ', np.nan):+.2f}% se {L['sfcWind']['c22']['se']:.2f} "
                  + " ".join(f"{v} {L[v]['c22']:+.2f}" for v in ("cropFrac", "treeFrac", "grassFrac", "lai") if v in L)
                  + "  ctl crop/tree " + "/".join(f"{L[v]['ctl_B']['eu8']:.1f}" for v in ("cropFrac", "treeFrac") if v in L), flush=True)
            if "early" in o:
                E = o["early"]["sfcWind"]["eu8"]
                print(f"{'':14s} {'early':7s} ds {E['d']:+.2f}% se {E['se_member']:.2f}  " + "  ".join(f"{v} {o['early'][v]['eu8']:+.2f}" for v in ("cropFrac", "treeFrac", "lai") if v in o["early"])
                      + (f"  | change A {o['change_early_late']['A']['eu8']:+.2f}% B {o['change_early_late']['B']['eu8']:+.2f}% DiD {o['change_early_late']['A_minus_B']['eu8']:+.2f}%" if "change_early_late" in o else ""))

    day = daily_part()

    # per-unit land cover: pooled over (model, pair), late window, at fixed circulation where available
    cm = json.load(open(f"{RES}/stability_cmip6.json"))
    unit = {}
    for s in SETS:
        g = "eu8" if s == "eu8" else "cmip22"
        unit[s] = {}
        for v in ("cropFrac", "treeFrac", "lai"):
            pts = []
            for k, o in mon.items():
                L = o["late"]
                if v not in L:
                    continue
                y = L.get("circulation", {}).get(s, {}).get("d_fixed_circ", L["sfcWind"][s]["d"])
                se = L["sfcWind"][s]["se"]
                pts.append((k, L[v][s], y, se, L["sfcWind"][s]["d"]))
            pts = [p for p in pts if np.isfinite(p[3])]
            if len(pts) < 3:
                continue
            x = np.array([p[1] for p in pts])
            y = np.array([p[2] for p in pts])
            se = np.array([p[3] for p in pts])
            b, sb, sb_infl, chi2 = wls_origin(x, y, se)
            ols = np.polyfit(x, y, 1)
            # leave-one-model-out
            mods = sorted({p[0].split("|")[0] for p in pts})
            lomo = {}
            for mm in mods:
                keep = np.array([p[0].split("|")[0] != mm for p in pts])
                if keep.sum() >= 2:
                    lomo[mm] = wls_origin(x[keep], y[keep], se[keep])[0]
            ref = cm["correlations"][g][f"within_G~d_{v}"]["cmip20"]
            unit[s][v] = dict(points={p[0]: dict(dx=p[1], d_fixed_circ=p[2], se=p[3], ds=p[4]) for p in pts},
                              slope_wls_origin=b, se=sb, se_scatter_inflated=sb_infl, reduced_chi2=chi2,
                              ols_slope=float(ols[0]), ols_intercept=float(ols[1]), r=float(np.corrcoef(x, y)[0, 1]),
                              leave_one_model_out=lomo, cmip6_cross_model_slope=ref["slope"], cmip6_cross_model_intercept=ref["intercept"],
                              cmip6_cross_model_r=ref["pearson"], cmip6_n=ref["n"])
            print(f"per unit {s} {v:9s}: LUMIP {b:+.3f} +/- {sb_infl:.3f} (wls through 0, n {len(pts)}, chi2 {chi2:.2f}; OLS {ols[0]:+.3f}, icpt {ols[1]:+.2f}, r {np.corrcoef(x, y)[0, 1]:+.2f})"
                  f"   CMIP6 cross-model {ref['slope']:+.3f} (icpt {ref['intercept']:+.2f}, r {ref['pearson']:+.2f}, n {ref['n']})   LOMO {min(lomo.values()):+.3f}..{max(lomo.values()):+.3f}", flush=True)

    # attribution of the CMIP6 at-fixed-G term
    attr = {}
    for s in SETS:
        g = "eu8" if s == "eu8" else "cmip22"
        ens = {n: m for n, m in cm["models"].items() if not m["wind_only"]}
        wg = np.array([m[g]["within_G"] for m in ens.values()])
        attr[s] = dict(cmip6_within_G_median=float(np.median(wg)), cmip6_within_G_mean=float(wg.mean()), n=len(wg))
        for v in ("cropFrac", "treeFrac", "lai"):
            if v not in unit[s]:
                continue
            u = unit[s][v]
            b, sb = u["slope_wls_origin"], u["se_scatter_inflated"]
            rows = {n: (m[g][f"d_{v}"], m[g]["within_G"]) for n, m in ens.items() if np.isfinite(m[g].get(f"d_{v}", np.nan))}
            dx = np.array([a for a, _ in rows.values()])
            wgs = np.array([c for _, c in rows.values()])
            pred = b * dx
            draws = (b + sb * RNG.standard_normal(4000))[:, None] * dx[None]
            med_draw = np.median(draws, axis=1)
            mean_draw = draws.mean(axis=1)
            share = pred.mean() / wgs.mean()
            attr[s][v] = dict(n_models=len(rows), d_lc_median=float(np.median(dx)), d_lc_mean=float(dx.mean()),
                              within_G_subset_median=float(np.median(wgs)), within_G_subset_mean=float(wgs.mean()),
                              lu_pred_median=float(np.median(pred)), lu_pred_median_p05_p95=[float(q) for q in np.percentile(med_draw, [5, 95])],
                              lu_pred_mean=float(pred.mean()), lu_pred_mean_p05_p95=[float(q) for q in np.percentile(mean_draw, [5, 95])],
                              share_of_subset_mean=float(share),
                              share_p05_p95=[float(q) for q in np.percentile(mean_draw / wgs.mean(), [5, 95])],
                              residual_median=float(np.median(wgs - pred)), residual_mean=float((wgs - pred).mean()),
                              cross_model_pred_mean=float((u["cmip6_cross_model_slope"] * dx).mean()),
                              per_model={n: dict(d_lc=float(a), within_G=float(c), lu_pred=float(b * a), residual=float(c - b * a)) for n, (a, c) in rows.items()})
            a = attr[s][v]
            print(f"attribution {s} {v:9s}: {a['n_models']} models, d_lc median {a['d_lc_median']:+.2f}; within_G subset median {a['within_G_subset_median']:+.2f} mean {a['within_G_subset_mean']:+.2f}; "
                  f"LU part median {a['lu_pred_median']:+.2f} [{a['lu_pred_median_p05_p95'][0]:+.2f}, {a['lu_pred_median_p05_p95'][1]:+.2f}], mean {a['lu_pred_mean']:+.2f} "
                  f"(share of mean {a['share_of_subset_mean']:.0%} [{a['share_p05_p95'][0]:.0%}, {a['share_p05_p95'][1]:.0%}]); residual median {a['residual_median']:+.2f}; "
                  f"cross-model slope would give {a['cross_model_pred_mean']:+.2f}", flush=True)

    # own-model check: LUMIP models that are also in the CMIP6 set
    own = {}
    alias = {"CESM2": "CESM2-WACCM"}
    for k, o in mon.items():
        model, pair = k.split("|")
        cmn = alias.get(model, model)
        if cmn not in cm["models"]:
            continue
        L = o["late"]
        for s in SETS:
            g = "eu8" if s == "eu8" else "cmip22"
            m = cm["models"][cmn][g]
            for v in ("cropFrac", "treeFrac"):
                if v in L and abs(L[v][s]) > 0.5 and np.isfinite(m.get(f"d_{v}", np.nan)):
                    y = L.get("circulation", {}).get(s, {}).get("d_fixed_circ", L["sfcWind"][s]["d"])
                    sl = y / L[v][s]
                    own.setdefault(k, {}).setdefault(s, {})[v] = dict(slope=float(sl), se=float(L["sfcWind"][s]["se"] / abs(L[v][s])),
                                                                     cmip6_name=cmn, cmip6_d_lc=m[f"d_{v}"], cmip6_within_G=m["within_G"],
                                                                     lu_pred=float(sl * m[f"d_{v}"]))
    fits = multi_fits(mon, cm)
    out = dict(meta=dict(pairs=PAIRS, regions_eu8=EU, regions_c22=C22, sign="A minus B: historical minus hist-noLu; ssp370 minus ssp370-ssp126Lu; "
                         "ssp126-ssp370Lu minus ssp126; esm-ssp585 minus esm-ssp585-ssp126Lu", wind="Amon/sfcWind monthly means (winter Sep-Mar)",
                         G="monthly-mean psl, bilinear 2.5 deg, |grad p|/(rho f); a control for circulation differences only"),
               monthly=mon, daily=day, per_unit=unit, attribution=attr, fits=fits, own_model=own)
    for k, o in own.items():
        for s, d in o.items():
            print(f"own model {k:22s} {s}: " + "  ".join(f"{v} slope {x['slope']:+.2f} +/- {x['se']:.2f} %/pp -> {x['cmip6_name']} ssp585 d {x['cmip6_d_lc']:+.2f} pp gives {x['lu_pred']:+.2f}% (its at-fixed-G {x['cmip6_within_G']:+.2f}%)" for v, x in d.items()))
    json.dump(out, open(f"{RES}/lumip_landuse.json", "w"), indent=1, default=float)
    plot(mon, unit, cm)
    print(f"Saved {RES}/lumip_landuse.json")


def plot(mon, unit, cm):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(16, 5))
    cols = {"hist": "tab:blue", "ssp370": "tab:orange", "ssp126": "tab:green", "esm585": "tab:red"}
    ax = axs[0]
    keys = list(mon)
    for i, k in enumerate(keys):
        L = mon[k]["late"]["sfcWind"]["eu8"]
        pair = k.split("|")[1]
        ci = L["ci90_member"]
        ax.errorbar(L["d"], i, xerr=None if ci is None else [[L["d"] - ci[0]], [ci[1] - L["d"]]], fmt="o", color=cols[pair], capsize=2)
        if "same_exp" in L:
            ax.plot(L["same_exp"]["p05_p95"], [i - 0.25] * 2, "-", color="0.6", lw=3, alpha=0.6)
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([k.replace("|", " ") for k in keys], fontsize=7)
    ax.axvline(0, color="0.7", lw=0.5)
    ax.set_xlabel("land-use effect on winter sfcWind, 8 EU regions (%)")
    ax.set_title("A minus B, late window; bars: 90 % member range; grey: same-experiment 5-95 %", fontsize=8)
    for ax, v, lab in ((axs[1], "cropFrac", "d crop fraction (pp)"), (axs[2], "treeFrac", "d tree fraction (pp)")):
        if v not in unit["eu8"]:
            continue
        u = unit["eu8"][v]
        for k, p in u["points"].items():
            ax.errorbar(p["dx"], p["d_fixed_circ"], yerr=1.645 * p["se"], fmt="o", color=cols[k.split("|")[1]], capsize=2)
            ax.annotate(k.split("|")[0].replace("-ESM", "E").replace("1-2-", ""), (p["dx"], p["d_fixed_circ"]), fontsize=6)
        xs = np.array([min(p["dx"] for p in u["points"].values()) - 1, max(p["dx"] for p in u["points"].values()) + 1])
        xm = [m["eu8"].get(f"d_{v}", np.nan) for m in cm["models"].values() if not m["wind_only"]]
        ym = [m["eu8"]["within_G"] for m in cm["models"].values() if not m["wind_only"]]
        ax.plot(xm, ym, "x", color="0.5", ms=5, label="CMIP6 ssp585 at-fixed-G term")
        xx = np.linspace(min(xs[0], np.nanmin(xm)), max(xs[1], np.nanmax(xm)), 10)
        ax.plot(xx, u["slope_wls_origin"] * xx, "-", color="k", lw=1, label=f"LUMIP {u['slope_wls_origin']:+.2f} %/pp")
        ax.plot(xx, u["cmip6_cross_model_slope"] * xx + u["cmip6_cross_model_intercept"], "--", color="0.5", lw=1,
                label=f"CMIP6 cross-model {u['cmip6_cross_model_slope']:+.2f} %/pp")
        ax.axhline(0, color="0.8", lw=0.5)
        ax.axvline(0, color="0.8", lw=0.5)
        ax.set_xlabel(lab)
        ax.set_ylabel("wind change at fixed circulation (%)")
        ax.legend(fontsize=7)
    fig.suptitle("LUMIP land-use effect on winter 10 m wind (blue hist - hist-noLu, orange ssp370 - ssp370-ssp126Lu, green ssp126-ssp370Lu - ssp126, red esm-ssp585 - esm-ssp585-ssp126Lu)", fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{RES}/lumip_landuse.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
