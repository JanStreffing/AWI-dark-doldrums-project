"""
Gradient-wind check of the G weakening (Thomas Jung's hypothesis, 2026-11-19): does the weakening of the geostrophic
speed G under Scandinavian blocking underestimate the weakening of the balanced wind, because broader systems are less
curved and so less supergeostrophic around highs?

Cheap version on the existing daily msl fields (regime grid 2.5 deg, 20-80 N, 80 W-40 E, the fields of the regime and
G analysis), no new data:
  G      |grad p| / (rho f), centred differences on the sphere (as geostrophic.geostrophic, non-periodic)
  kappa  signed curvature of the isobars, (p_y^2 p_xx - 2 p_x p_y p_xy + p_x^2 p_yy) / |grad p|^3, positive cyclonic
         (around a low), negative anticyclonic (around a high); metric terms of the sphere neglected
  V      gradient wind from kappa V^2 + f V = f G, regular branch V = 2 G / (1 + sqrt(1 + 4 kappa G / f)).
         Where the anticyclonic discriminant is negative (1 + 4 kappa G / f < 0) V is set to its limit 2 G and the
         point is counted as 'clipped'.
Per region box (the 8 European regions of the G analysis) and day: mean G, mean V, and both over anticyclonic and over
cyclonic points. Regime labels as in the paper (regime_G_robustness.classify / geostrophic.regime_labels).
Changes: per region 100 (mean_fut / mean_ctl - 1), mean over regions, per K as in regime_G_robustness.
Weak-wind days: days with box-mean G below the control 20th percentile of that region and regime; their mean V/G.
Caveat: this is the balance of the sea-level pressure field, not the wind above the boundary layer; curvature from a
2.5 deg field is smooth by construction; isobar curvature equals streamline curvature only for balanced, steady flow.

Environment: levante, module load python3/2025.01-gcc-13.3.0
Usage: python3 gradient_wind.py
Output: dunkelflaute/results/regimes/gradient_wind.json
"""
import glob
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regimes")
import regime_G_robustness as rr               # noqa: E402  (imports geostrophic, which sets region_events.WIN = 2)
gs, ra, ev = rr.gs, rr.ra, rr.ev

OUT = f"{_ROOT}/dunkelflaute/results/regimes/gradient_wind.json"
EU = gs.EU
RHO, OMEGA, RE = gs.RHO, gs.OMEGA, gs.RE
OURS = ["9km", "31km", "100km", "eerie_ctl"]
NAMES = ["NAO+", "NAO-", "ScBL", "AR"]


def fields(p, lat, lon):
    """p (D, ny, nx) Pa -> G, V, kappa, clipped (D, ny, nx); NaN on the border."""
    phi = np.deg2rad(lat)[None, :, None]
    dy = np.deg2rad(lat[1] - lat[0]) * RE
    dx = np.deg2rad(lon[1] - lon[0]) * RE * np.cos(phi)
    nan = lambda: np.full_like(p, np.nan)                 # noqa: E731
    px, py, pxx, pyy, pxy = nan(), nan(), nan(), nan(), nan()
    px[:, :, 1:-1] = (p[:, :, 2:] - p[:, :, :-2]) / (2 * dx)
    py[:, 1:-1] = (p[:, 2:] - p[:, :-2]) / (2 * dy)
    pxx[:, :, 1:-1] = (p[:, :, 2:] - 2 * p[:, :, 1:-1] + p[:, :, :-2]) / dx ** 2
    pyy[:, 1:-1] = (p[:, 2:] - 2 * p[:, 1:-1] + p[:, :-2]) / dy ** 2
    pxy[:, 1:-1, 1:-1] = (p[:, 2:, 2:] - p[:, 2:, :-2] - p[:, :-2, 2:] + p[:, :-2, :-2]) / (4 * dy * dx[:, 1:-1])
    f = 2 * OMEGA * np.sin(phi)
    gp = np.hypot(px, py)
    G = gp / (RHO * f)
    with np.errstate(divide="ignore", invalid="ignore"):
        kap = (py ** 2 * pxx - 2 * px * py * pxy + px ** 2 * pyy) / gp ** 3
        arg = 1 + 4 * kap * G / f
        clipped = arg < 0
        V = 2 * G / (1 + np.sqrt(np.clip(arg, 0, None)))
    return G, V, kap, clipped


def box_means(G, V, kap, clipped, lat, lon):
    """Per day and region: G, V, G and V over anticyclonic and cyclonic points, clipped share (D, R) each."""
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    out = {k: [] for k in ("G", "V", "G_anti", "V_anti", "G_cyc", "V_cyc", "anti_share", "clipped_share")}
    for r in EU:
        la0, la1, lo0, lo1 = gs.BOX[r]
        sel = (LA >= la0) & (LA <= la1) & (LO >= lo0) & (LO <= lo1)
        if not sel.any():
            sel = (LA >= la0 - 1.25) & (LA <= la1 + 1.25) & (LO >= lo0 - 1.25) & (LO <= lo1 + 1.25)
        g, v, k, c = G[:, sel], V[:, sel], kap[:, sel], clipped[:, sel]
        a, cy = k < 0, k > 0
        with np.errstate(invalid="ignore"):
            out["G"].append(np.nanmean(g, 1)); out["V"].append(np.nanmean(v, 1))
            out["G_anti"].append(np.nanmean(np.where(a, g, np.nan), 1)); out["V_anti"].append(np.nanmean(np.where(a, v, np.nan), 1))
            out["G_cyc"].append(np.nanmean(np.where(cy, g, np.nan), 1)); out["V_cyc"].append(np.nanmean(np.where(cy, v, np.nan), 1))
            out["anti_share"].append(np.mean(a, 1)); out["clipped_share"].append(np.mean(c, 1))
    return {k: np.stack(v, 1) for k, v in out.items()}


def chg(a, b):
    return 100 * (np.nanmean(b, 0) / np.nanmean(a, 0) - 1)          # per region


def per_k(p, dT):
    return 100 * (np.exp(np.log(1 + p / 100) / dT) - 1)


def analyse(Bc, Bf, lc, lf, dT):
    res = {}
    for name, (mc, mf) in [(n, (lc == k, lf == k)) for k, n in enumerate(NAMES)] + [("all", (np.ones(len(lc), bool), np.ones(len(lf), bool)))]:
        row = {}
        for q in ("G", "V", "G_anti", "V_anti", "G_cyc", "V_cyc"):
            c = chg(Bc[q][mc], Bf[q][mf])
            row[q] = dict(change_pct=float(np.nanmean(c)), per_K=float(per_k(np.nanmean(c), dT)))
        rc = np.nanmean(Bc["V"][mc], 0) / np.nanmean(Bc["G"][mc], 0)
        rf = np.nanmean(Bf["V"][mf], 0) / np.nanmean(Bf["G"][mf], 0)
        row["V_over_G"] = dict(ctl=float(np.mean(rc)), fut=float(np.mean(rf)), change_pct=float(np.mean(100 * (rf / rc - 1))))
        row["anti_share"] = [float(np.nanmean(Bc["anti_share"][mc])), float(np.nanmean(Bf["anti_share"][mf]))]
        row["clipped_share"] = [float(np.nanmean(Bc["clipped_share"][mc])), float(np.nanmean(Bf["clipped_share"][mf]))]
        # weak-wind days: box G below the control 20th percentile of the region (this regime)
        thr = np.nanpercentile(Bc["G"][mc], 20, axis=0)
        wc, wf = Bc["G"][mc] < thr, Bf["G"][mf] < thr
        r_c = np.nanmean(np.where(wc, Bc["V"][mc] / Bc["G"][mc], np.nan), 0)
        r_f = np.nanmean(np.where(wf, Bf["V"][mf] / Bf["G"][mf], np.nan), 0)
        row["weak_days_V_over_G"] = dict(ctl=float(np.nanmean(r_c)), fut=float(np.nanmean(r_f)),
                                         change_pct=float(np.nanmean(100 * (r_f / r_c - 1))))
        q10 = lambda X: np.nanpercentile(X, 10, axis=0)        # noqa: E731
        row["q10_change_pct"] = dict(G=float(np.mean(100 * (q10(Bf["G"][mf]) / q10(Bc["G"][mc]) - 1))),
                                     V=float(np.mean(100 * (q10(Bf["V"][mf]) / q10(Bc["V"][mc]) - 1))))
        res[name] = row
    return res


def cmip6_fields(path, pat):
    model = os.path.basename(path)[4:-4]
    Pz, Z = np.load(path), np.load(f"{gs.DC}/daily_cmip6_{model}.npz")
    lat, lon = Pz["lat"], Pz["lon"]
    ii = np.where(lat >= 20.0)[0]
    jj = np.where((lon >= -80.0) & (lon <= 40.0))[0]
    D = {}
    for tag, exp, y0, y1 in (("ctl", "historical", 1951, 1980), ("fut", "ssp585", 2071, 2100)):
        k = Pz[f"{exp}|ymd"] @ np.array([10000, 100, 1])
        w = gs.winter_of(k)
        m = (w >= y0) & (w <= y1)
        D[tag] = (gs.dates_of(k[m]), Pz[f"{exp}|psl"][m][:, ii][:, :, jj].astype(np.float64))
    lc, lf = gs.regime_labels(pat, D["ctl"][0], D["ctl"][1], D["fut"][0], D["fut"][1])
    T = {e: dict(zip(Z[f"{e}|tas_years"].tolist(), Z[f"{e}|tas"].tolist())) for e in ("historical", "ssp585")}
    dT = float(np.mean([T["ssp585"][y] for y in range(2071, 2101) if y in T["ssp585"]])
               - np.mean([T["historical"][y] for y in range(1951, 1981) if y in T["historical"]]))
    B = {t: box_means(*fields(D[t][1], lat[ii], lon[jj]), lat[ii], lon[jj]) for t in D}
    return model, B["ctl"], B["fut"], lc, lf, dT


def main():
    pat = np.load(f"{ra.OUT}/regime_patterns.npz")
    ra.NAMES_ORDER = [str(x) for x in pat["names"]]
    assert ra.NAMES_ORDER == NAMES
    lat, lon = pat["lat"], pat["lon"]
    stored = json.load(open(f"{_ROOT}/dunkelflaute/results/regimes/regime_G_robustness.json"))
    out = dict(method=__doc__.split("Environment")[0], runs={}, cmip6={})
    for name in OURS:
        P = rr.comparison(name, pat)
        mlat = np.load(f"{ra.D}/msl_{ev.COMPARISONS[name][0][0]}.npz")
        la, lo = mlat["lat"], mlat["lon"]
        B = {t: box_means(*fields(P[t]["msl"], la, lo), la, lo) for t in ("ctl", "fut")}
        r = analyse(B["ctl"], B["fut"], P["ctl"]["lab"], P["fut"]["lab"], P["dT"])
        st = stored["runs"][name]["change_per_K"]["ScBL"] if "change_per_K" in stored["runs"][name] else None
        r["dT"] = P["dT"]
        out["runs"][name] = r
        print(f"{name:10s} ScBL per K: G {r['ScBL']['G']['per_K']:+.2f} (stored {st})  V {r['ScBL']['V']['per_K']:+.2f}  "
              f"V/G {r['ScBL']['V_over_G']['ctl']:.3f}->{r['ScBL']['V_over_G']['fut']:.3f}  anti share {r['ScBL']['anti_share'][0]:.2f}  "
              f"clipped {r['ScBL']['clipped_share'][0]:.3f}  weak-day V/G {r['ScBL']['weak_days_V_over_G']['ctl']:.3f}->{r['ScBL']['weak_days_V_over_G']['fut']:.3f}",
              flush=True)
    for path in sorted(glob.glob(f"{gs.DC}/psl_*.npz")):
        model = os.path.basename(path)[4:-4]
        if not os.path.exists(f"{gs.DC}/daily_cmip6_{model}.npz"):
            continue
        model, Bc, Bf, lc, lf, dT = cmip6_fields(path, pat)
        r = analyse(Bc, Bf, lc, lf, dT)
        r["dT"] = dT
        out["cmip6"][model] = r
        print(f"{model:18s} ScBL per K: G {r['ScBL']['G']['per_K']:+.2f}  V {r['ScBL']['V']['per_K']:+.2f}  "
              f"V/G change {r['ScBL']['V_over_G']['change_pct']:+.2f} %  weak-day V/G change {r['ScBL']['weak_days_V_over_G']['change_pct']:+.2f} %",
              flush=True)
    ens = {}
    for reg in NAMES + ["all"]:
        for q in ("G", "V"):
            x = np.array([v[reg][q]["per_K"] for v in out["cmip6"].values()])
            ens[f"{reg}|{q}_per_K"] = dict(median=float(np.median(x)), p10=float(np.percentile(x, 10)), p90=float(np.percentile(x, 90)))
        d = np.array([v[reg]["V"]["per_K"] - v[reg]["G"]["per_K"] for v in out["cmip6"].values()])
        ens[f"{reg}|V_minus_G_per_K"] = dict(median=float(np.median(d)), n_negative=int((d < 0).sum()), n=len(d))
        w = np.array([v[reg]["weak_days_V_over_G"]["change_pct"] for v in out["cmip6"].values()])
        ens[f"{reg}|weak_day_V_over_G_change_pct"] = dict(median=float(np.median(w)), n_negative=int((w < 0).sum()))
    out["cmip6_ensemble"] = ens
    json.dump(out, open(OUT, "w"), indent=1)
    print(json.dumps({k: v for k, v in ens.items() if k.startswith("ScBL") or k.startswith("all")}, indent=1))
    print(OUT)


if __name__ == "__main__":
    main()
