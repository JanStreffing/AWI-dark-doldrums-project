"""
Why does the geostrophic speed G weaken within the Euro-Atlantic regimes (most robustly under Scandinavian blocking),
and why do the models agree?

Definitions are those of regime_G_robustness.py / geostrophic.py / regime_analysis.py, imported rather than copied:
  G       box mean over the 8 European regions of |grad p| / (rho f) from daily msl (2.5 deg)
  labels  stored footing: anomalies from the CONTROL seasonal cycle (both periods), ERA5 EOFs, nearest ERA5 centroid
  change  100 (mean_fut G / mean_ctl G - 1) per region over the regime days, then the mean over the regions
The script first reproduces the per-regime change of regime_G_robustness.json for every run and CMIP6 model.

1  Mean-state / anomaly factorial (2 x 2, day by day, G recomputed from each synthetic field)
   Each period gets its own smooth seasonal cycle C_ctl, C_fut (ra.seasonal_fit: mean + 3 harmonics on Sep-Mar days)
   and its own anomalies a = msl - C (so C + a is exactly the daily field). With D = C_fut - C_ctl (a function of the
   calendar day only):
     G_cc = G(C_ctl + a_ctl)  on control regime days    (= actual control)
     G_fc = G(C_fut + a_ctl)  on control regime days    (background changed, control weather)
     G_cf = G(C_ctl + a_fut)  on future regime days     (control background, future weather)
     G_ff = G(C_fut + a_fut)  on future regime days     (= actual future)
   background = chg(G_cc -> G_fc), anomaly = chg(G_cc -> G_cf), interaction = total - background - anomaly
   (chg as the stored change, mean over regions of per-region %). The regime days are those of the stored footing,
   so the anomaly term also carries any change in which days are labelled with a regime. Sensitivity: future labelled
   with its own seasonal cycle (removes the mean change from the labelling). Winter bootstrap 5-95 % (NB draws).
   % per K: each term scaled by per_k(total) / total, so the terms add up to the stored % per K of the total.
2  Amplitude and scale of the regime anomaly (own-cycle anomalies, so the background change is left out):
   composite peak, its position, area-mean |anomaly|, area above +5 hPa and above half the peak, the centroid of the
   half-peak area, the box means of the composite, G of the composite anomaly and of the composite full field,
   the mean daily G of the anomaly alone G(a) and of the background alone G(C), daily quantiles of the Europe-mean G,
   and the daily projection amplitude of the anomaly on the control composite. Within-regime residual
   e = a - (regime composite of a, same period): mean daily G(e) and rms of e over the box grid points, i.e. whether the
   weakening sits in the stationary regime pattern or in the day-to-day (synoptic) deviations from it.
3  Background gradient: winter-mean (Sep-Mar) msl, sector mean 80 W-40 E (and 10 W-30 E), p(40 N) - p(65 N) in hPa;
   G of the winter-mean field over the boxes; temperature: no gridded temperature exists in the repo, only box means of
   2 m temperature (9 km, 31 km: data/regimes/interannual_ours; 100 km: data/regimes/transfer_stability_parts;
   CMIP6: monthly box tas in data/regimes/interannual_cmip6; none for IFS-FESOM). Proxy for the meridional
   lower-tropospheric temperature gradient: T(SC) - (T(IB) + T(IT)) / 2, winter mean. Cross-model correlations.
4  The factorial for the 20 CMIP6 models (historical 1951-1980 against ssp585 2071-2100, psl on the global 15-80 N band,
   periodic differences, as geostrophic.cmip6()).

Environment: levante, module load python3/2025.01-gcc-13.3.0
Usage: python3 regime_G_attribution.py
Output: dunkelflaute/results/regimes/regime_G_attribution.json, regime_G_attribution/*.png
"""
import glob
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regimes")
import regime_G_robustness as rgr          # noqa: E402  (imports geostrophic, which sets region_events.WIN = 2; no events here)
gs, ra, rg = rgr.gs, rgr.ra, rgr.rg

import matplotlib                          # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt            # noqa: E402
from scipy import stats                    # noqa: E402

OUT = f"{_ROOT}/dunkelflaute/results/regimes"
FIG = f"{OUT}/regime_G_attribution"
RUNS = rgr.RUNS
LAB = {"9km": "9 km", "31km": "31 km", "100km": "100 km", "eerie_ctl": "IFS-FESOM"}
EU = ra.EU
NB = 2000
X_HPA = 5.0
CHUNK = 2000
T2_OURS = {"9km": "interannual_ours", "31km": "interannual_ours", "100km": "transfer_stability_parts"}


# ---------------------------------------------------------------- helpers
def boxG(p, lat, lon, periodic):
    """(D, ny, nx) Pa -> (D, 8) box-mean G, in chunks."""
    return np.concatenate([gs.region_G(gs.geostrophic(p[i:i + CHUNK], lat, lon, periodic), lat, lon, EU)
                           for i in range(0, len(p), CHUNK)])


def clim(dates, coef, shape):
    return (ra.harmonics(dates) @ coef).reshape((len(dates),) + shape)


def chg(Gc, Gf):
    return float(100 * (Gf.mean(0) / Gc.mean(0) - 1).mean())


def per_k(pct, dT):
    return rgr.per_k(pct, dT)


def area_w(lat, lon):
    """Cell areas (1e6 km2) on the regular grid."""
    dl = np.deg2rad(abs(lon[1] - lon[0]))
    dp = np.deg2rad(abs(lat[1] - lat[0]))
    a = (gs.RE / 1e3) ** 2 * dl * dp * np.cos(np.deg2rad(lat))[:, None] * np.ones((1, len(lon)))
    return a / 1e6


def bootstrap_terms(rng, win_c, win_f, lab_c, lab_f, Gcc, Gfc, Gcf, Gff):
    """Winter bootstrap of total / background / anomaly / interaction (B, 5) for the 4 regimes + all days."""
    yc, yf = np.unique(win_c), np.unique(win_f)
    Scc, Nc = rgr.winter_sums(win_c, lab_c, Gcc, yc)
    Sfc, _ = rgr.winter_sums(win_c, lab_c, Gfc, yc)
    Scf, Nf = rgr.winter_sums(win_f, lab_f, Gcf, yf)
    Sff, _ = rgr.winter_sums(win_f, lab_f, Gff, yf)
    wc, wf = rgr.boot_weights(rng, len(yc), NB), rgr.boot_weights(rng, len(yf), NB)
    tot = rgr.change(wc, Scc, Nc, wf, Sff, Nf)
    bg = rgr.change(wc, Scc, Nc, wc, Sfc, Nc)
    an = rgr.change(wc, Scc, Nc, wf, Scf, Nf)
    return dict(total=tot, background=bg, anomaly=an, interaction=tot - bg - an)


def blob(A, lat, lon, aw):
    """Shape metrics of a composite anomaly map A (hPa) on the regime domain."""
    i, j = np.unravel_index(np.argmax(A), A.shape)
    pk = float(A[i, j])
    half = A >= 0.5 * pk
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    w = aw * half * A
    return dict(peak_hPa=pk, peak_lat=float(lat[i]), peak_lon=float(lon[j]), min_hPa=float(A.min()),
                mean_abs_hPa=float((np.abs(A) * aw).sum() / aw.sum()),
                area_above_X_1e6km2=float(aw[A >= X_HPA].sum()), area_above_half_peak_1e6km2=float(aw[half].sum()),
                half_peak_centroid_lat=float((w * LA).sum() / w.sum()), half_peak_centroid_lon=float((w * LO).sum() / w.sum()))


def box_means_map(A, lat, lon):
    return {r: float(np.nanmean(A[box_sel(r, lat, lon)])) for r in EU}


def box_sel(r, lat, lon):
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    la0, la1, lo0, lo1 = gs.BOX[r]
    return (LA >= la0) & (LA <= la1) & (LO >= lo0) & (LO <= lo1)


def gradient_index(M, lat, lon):
    """Winter-mean field M (ny, nx, Pa) -> sector-mean p(40N) - p(65N) hPa for 80W-40E and 10W-30E."""
    i40, i65 = int(np.argmin(np.abs(lat - 40))), int(np.argmin(np.abs(lat - 65)))
    out = {}
    for nm, (a, b) in (("atl_eu", (-80, 40)), ("eu", (-10, 30))):
        jj = (lon >= a) & (lon <= b)
        out[nm] = float((M[i40, jj].mean() - M[i65, jj].mean()) / 100)
    return out


# ---------------------------------------------------------------- core analysis of one comparison
def attribution(c, f, lat, lon, periodic, dom, pat, rng, dT):
    """c, f: dicts with d (dates), msl (D, ny, nx Pa, the field G is computed on), win, lab (stored labels), G.
    dom: (ii, jj) index arrays of the 20-80 N, 80 W-40 E regime domain inside the field."""
    ii, jj = dom
    shp = c["msl"].shape[1:]
    coef_c = ra.seasonal_fit(c["d"], c["msl"])
    coef_f = ra.seasonal_fit(f["d"], f["msl"])
    Cc_c, Cf_c = clim(c["d"], coef_c, shp), clim(c["d"], coef_f, shp)
    Cc_f, Cf_f = clim(f["d"], coef_c, shp), clim(f["d"], coef_f, shp)
    a_c, a_f = c["msl"] - Cc_c, f["msl"] - Cf_f
    G = dict(cc=boxG(c["msl"], lat, lon, periodic), fc=boxG(Cf_c + a_c, lat, lon, periodic),
             cf=boxG(Cc_f + a_f, lat, lon, periodic), ff=boxG(f["msl"], lat, lon, periodic),
             anom_c=boxG(a_c, lat, lon, periodic), anom_f=boxG(a_f, lat, lon, periodic),
             clim_c=boxG(Cc_c, lat, lon, periodic), clim_f=boxG(Cf_f, lat, lon, periodic))
    repro = max(float(np.abs(G["cc"] - c["G"]).max()), float(np.abs(G["ff"] - f["G"]).max()))
    # own-cycle labels for the future (sensitivity)
    lf_own, _ = rgr.classify(pat, f["d"], f["msl"][:, ii][:, :, jj], coef_f.reshape(-1, *shp)[:, ii][:, :, jj].reshape(len(coef_f), -1))
    lab_c, lab_f = c["lab"], f["lab"]
    bs = bootstrap_terms(rng, c["win"], f["win"], lab_c, lab_f, G["cc"], G["fc"], G["cf"], G["ff"])
    aw = area_w(lat[ii], lon[jj])
    res = dict(reproduce_max_abs_G_diff=repro, regimes={})
    for j, nm in enumerate(ra.NAMES_ORDER + ["all"]):
        mc = np.ones(len(lab_c), bool) if nm == "all" else lab_c == j
        mf = np.ones(len(lab_f), bool) if nm == "all" else lab_f == j
        mo = np.ones(len(lab_f), bool) if nm == "all" else lf_own == j
        tot, bg, an = chg(G["cc"][mc], G["ff"][mf]), chg(G["cc"][mc], G["fc"][mc]), chg(G["cc"][mc], G["cf"][mf])
        terms = dict(total=tot, background=bg, anomaly=an, interaction=tot - bg - an)
        # symmetric (other-base) version: background on future weather, anomaly on future background
        terms_alt = dict(background=chg(G["cf"][mf], G["ff"][mf]), anomaly=chg(G["fc"][mc], G["ff"][mf]))
        fac = per_k(tot, dT) / tot if abs(tot) > 1e-3 else 1 / dT
        r = dict(days=[int(mc.sum()), int(mf.sum())], change_pct=terms,
                 change_per_K={k: v * fac for k, v in terms.items()},
                 other_base_pct=terms_alt,
                 boot_p05_p95={k: [float(x) for x in np.percentile(v[:, j], [5, 95])] for k, v in bs.items()},
                 ownclim_labels_pct=dict(total=chg(G["cc"][mc], G["ff"][mo]), anomaly=chg(G["cc"][mc], G["cf"][mo])),
                 G_anomaly_only_change_pct=chg(G["anom_c"][mc], G["anom_f"][mf]),
                 G_clim_only_change_pct=chg(G["clim_c"][mc], G["clim_f"][mf]),
                 G_ctl_mean=float(G["cc"][mc].mean()), G_anom_only_ctl_mean=float(G["anom_c"][mc].mean()),
                 G_clim_only_ctl_mean=float(G["clim_c"][mc].mean()))
        if nm != "all":
            # composites on the regime domain (hPa), own-cycle anomalies
            Ac = a_c[mc][:, ii][:, :, jj].mean(0) / 100
            Af = a_f[mf][:, ii][:, :, jj].mean(0) / 100
            Bc = Cc_c[mc][:, ii][:, :, jj].mean(0) / 100
            Bf = Cf_f[mf][:, ii][:, :, jj].mean(0) / 100
            la, lo = lat[ii], lon[jj]
            gA = [float(gs.region_G(gs.geostrophic(x[None] * 100, la, lo, False), la, lo, EU)[0].mean()) for x in (Ac, Af, Ac + Bc, Af + Bf)]
            # daily Europe-mean G distribution
            ec, ef = G["cc"][mc].mean(1), G["ff"][mf].mean(1)
            q = [5, 10, 25, 50, 75, 90, 95]
            qc, qf = np.percentile(ec, q), np.percentile(ef, q)
            # daily projection amplitude of the anomaly on the control composite (area weights)
            Aw = Ac * aw
            amp_c = np.einsum("dij,ij->d", a_c[mc][:, ii][:, :, jj] / 100, Aw) / (Ac * Aw).sum()
            amp_f = np.einsum("dij,ij->d", a_f[mf][:, ii][:, :, jj] / 100, Aw) / (Ac * Aw).sum()
            r.update(composite=dict(ctl=blob(Ac, la, lo, aw), fut=blob(Af, la, lo, aw),
                                    box_anom_ctl_hPa=box_means_map(Ac, la, lo), box_anom_fut_hPa=box_means_map(Af, la, lo),
                                    G_of_composite_anomaly=[gA[0], gA[1]], G_of_composite_anomaly_change_pct=100 * (gA[1] / gA[0] - 1),
                                    G_of_composite_full=[gA[2], gA[3]], G_of_composite_full_change_pct=100 * (gA[3] / gA[2] - 1),
                                    pattern_corr=float(np.corrcoef((Ac * np.sqrt(aw)).ravel(), (Af * np.sqrt(aw)).ravel())[0, 1])),
                     daily=dict(q=q, europe_G_ctl=qc.tolist(), europe_G_fut=qf.tolist(),
                                q_change_pct=(100 * (qf / qc - 1)).tolist(),
                                frac_below_ctl_p10=[0.1, float(np.mean(ef < qc[1]))], frac_above_ctl_p90=[0.1, float(np.mean(ef > qc[5]))],
                                amp_ctl_q=np.percentile(amp_c, q).tolist(), amp_fut_q=np.percentile(amp_f, q).tolist(),
                                amp_mean=[float(amp_c.mean()), float(amp_f.mean())],
                                amp_sd=[float(amp_c.std()), float(amp_f.std())]))
            # within-regime residual e = a - composite (full field): day-to-day deviation from the regime pattern
            ec_, ef_ = a_c[mc] - a_c[mc].mean(0), a_f[mf] - a_f[mf].mean(0)
            Gec, Gef = boxG(ec_, lat, lon, periodic), boxG(ef_, lat, lon, periodic)
            eu = np.zeros((len(ii), len(jj)), bool)
            for rr in EU:
                eu |= box_sel(rr, lat[ii], lon[jj])
            rms = [float(np.sqrt(np.mean((x[:, ii][:, :, jj][:, eu] / 100) ** 2))) for x in (ec_, ef_)]
            rma = [float(np.sqrt(np.mean((x[:, ii][:, :, jj][:, eu] / 100) ** 2))) for x in (a_c[mc], a_f[mf])]
            r["residual"] = dict(G_resid_ctl_mean=float(Gec.mean()), G_resid_change_pct=chg(Gec, Gef),
                                 rms_resid_EU_hPa=rms, rms_resid_EU_change_pct=100 * (rms[1] / rms[0] - 1),
                                 rms_anom_EU_hPa=rma, rms_anom_EU_change_pct=100 * (rma[1] / rma[0] - 1))
            del ec_, ef_
            r["_maps"] = dict(Ac=Ac, Af=Af, Bc=Bc, Bf=Bf)
        res["regimes"][nm] = r
    # background (winter-mean) state
    Mc = c["msl"][:, ii][:, :, jj].mean(0)
    Mf = f["msl"][:, ii][:, :, jj].mean(0)
    la, lo = lat[ii], lon[jj]
    gc, gf = gradient_index(Mc, la, lo), gradient_index(Mf, la, lo)
    GMc = gs.region_G(gs.geostrophic(Mc[None], la, lo, False), la, lo, EU)[0]
    GMf = gs.region_G(gs.geostrophic(Mf[None], la, lo, False), la, lo, EU)[0]
    res["background"] = dict(
        dp40_65_hPa={k: [gc[k], gf[k]] for k in gc},
        dp40_65_change_pct={k: 100 * (gf[k] / gc[k] - 1) for k in gc},
        dp40_65_change_per_K={k: per_k(100 * (gf[k] / gc[k] - 1), dT) for k in gc},
        G_wintermean_field=[float(GMc.mean()), float(GMf.mean())],
        G_wintermean_field_change_pct=float(100 * (GMf / GMc - 1).mean()),
        G_wintermean_field_change_per_K=per_k(float(100 * (GMf / GMc - 1).mean()), dT),
        zonal_profile_hPa=[(Mc.mean(1) / 100).tolist(), (Mf.mean(1) / 100).tolist()])
    res["_maps"] = dict(dM=(Mf - Mc) / 100, Mc=Mc / 100)
    return res


# ---------------------------------------------------------------- temperature contrast proxy
def t_contrast(series_c, series_f):
    """(T_SC - (T_IB + T_IT)/2) winter means [ctl, fut], K."""
    iSC, iIB, iIT = EU.index("SC"), EU.index("IB"), EU.index("IT")
    out = []
    for T in (series_c, series_f):
        out.append(float(np.mean(T[:, iSC] - 0.5 * (T[:, iIB] + T[:, iIT]))))
    return out


def t2_ours(name):
    (cn, c0, c1), (fn, f0, f1) = rgr.ev.COMPARISONS[name]
    if name not in T2_OURS:
        return None
    S = []
    for run, y0, y1 in ((cn, c0, c1), (fn, f0, f1)):
        if T2_OURS[name] == "interannual_ours":
            z = np.load(f"{ra.D}/interannual_ours/{run}_t2.npz")
            assert list(z["regions"]) == EU
            day, t2 = z["day"].astype("datetime64[D]"), z["t2"]
        else:
            parts = sorted(glob.glob(f"{ra.D}/transfer_stability_parts/{run.replace('tco95_', '')}_*.npz"))
            zs = [np.load(p) for p in parts]
            day = np.concatenate([z["day"] for z in zs]).astype("datetime64[D]")
            t2 = np.concatenate([z["t2"] for z in zs])
        k = gs.ymd_of(day)
        w = gs.winter_of(k)
        sel = np.isin(k // 100 % 100, rg.MONTHS) & (w >= y0) & (w <= y1) & np.isfinite(t2).all(1)
        S.append(t2[sel])
    return t_contrast(*S)


def t2_cmip6(model):
    z = np.load(f"{ra.D}/interannual_cmip6/{model}.npz")
    assert list(z["regions"]) == EU
    S = []
    for exp, y0, y1 in (("historical", 1951, 1980), ("ssp585", 2071, 2100)):
        ym, T = z[f"{exp}|tas_ym"], z[f"{exp}|tas"]
        w = np.where(ym[:, 1] >= 9, ym[:, 0] + 1, ym[:, 0])
        sel = np.isin(ym[:, 1], rg.MONTHS) & (w >= y0) & (w <= y1)
        S.append(T[sel])
    return t_contrast(*S)


# ---------------------------------------------------------------- loaders
def load_ours(name, pat):
    P = rgr.comparison(name, pat)
    lat, lon = np.load(f"{ra.D}/msl_{rgr.ev.COMPARISONS[name][0][0]}.npz")["lat"], np.load(f"{ra.D}/msl_{rgr.ev.COMPARISONS[name][0][0]}.npz")["lon"]
    return P, lat, lon


def load_cmip6(model, pat):
    """As regime_G_robustness.cmip6_model, keeping the full-band psl."""
    c0, c1, f0, f1 = 1951, 1980, 2071, 2100
    Pz, Z = np.load(f"{gs.DC}/psl_{model}.npz"), np.load(f"{gs.DC}/daily_cmip6_{model}.npz")
    lat, lon = Pz["lat"], Pz["lon"]
    ii = np.where(lat >= 20.0)[0]
    jj = np.where((lon >= -80.0) & (lon <= 40.0))[0]
    assert len(ii) == 25 and len(jj) == 49 and np.all(np.diff(jj) == 1)
    P = {}
    for tag, exp, y0, y1 in (("ctl", "historical", c0, c1), ("fut", "ssp585", f0, f1)):
        kp = Pz[f"{exp}|ymd"] @ np.array([10000, 100, 1])
        kz = Z[f"{exp}|ymd"] @ np.array([10000, 100, 1])
        common, i_p, _ = np.intersect1d(kp, kz, return_indices=True)
        w = gs.winter_of(common)
        keep = (w >= y0) & (w <= y1)
        common, i_p = common[keep], i_p[keep]
        psl = Pz[f"{exp}|psl"][i_p].astype(np.float64)
        P[tag] = dict(k=common, msl=psl, win=gs.winter_of(common), d=gs.dates_of(common))
        P[tag]["G"] = boxG(psl, lat, lon, True)
    lc, lf = gs.regime_labels(pat, P["ctl"]["d"], P["ctl"]["msl"][:, ii][:, :, jj], P["fut"]["d"], P["fut"]["msl"][:, ii][:, :, jj])
    P["ctl"]["lab"], P["fut"]["lab"] = lc, lf
    T = {e: dict(zip(Z[f"{e}|tas_years"].tolist(), Z[f"{e}|tas"].tolist())) for e in ("historical", "ssp585")}
    P["dT"] = float(np.mean([T["ssp585"][y] for y in range(f0, f1 + 1) if y in T["ssp585"]])
                    - np.mean([T["historical"][y] for y in range(c0, c1 + 1) if y in T["historical"]]))
    return P, lat, lon, (ii, jj)


# ---------------------------------------------------------------- main
def strip(r):
    maps = {"_bg": r.pop("_maps")}
    for nm, v in r["regimes"].items():
        if "_maps" in v:
            maps[nm] = v.pop("_maps")
    return maps


def main():
    os.makedirs(FIG, exist_ok=True)
    pat = np.load(f"{ra.OUT}/regime_patterns.npz")
    ra.NAMES_ORDER = [str(x) for x in pat["names"]]
    stored = json.load(open(f"{OUT}/regime_G_robustness.json"))
    rng = np.random.default_rng(23)
    res = dict(columns=ra.NAMES_ORDER + ["all"], n_boot=NB, X_hPa=X_HPA, runs={}, cmip6={},
               t_proxy="T(SC) - (T(IB)+T(IT))/2, winter-mean box-mean 2 m temperature (no gridded T in the repo)")
    maps = {}
    for name in RUNS:
        P, lat, lon = load_ours(name, pat)
        dom = (np.arange(len(lat)), np.arange(len(lon)))
        r = attribution(P["ctl"], P["fut"], lat, lon, False, dom, pat, rng, P["dT"])
        r["dT"] = P["dT"]
        r["n_winters"] = [len(np.unique(P["ctl"]["win"])), len(np.unique(P["fut"]["win"]))]
        st = stored["runs"][name]["change_pct"]
        r["reproduce_vs_robustness_json_pct"] = float(max(abs(r["regimes"][c]["change_pct"]["total"] - st[c]) for c in res["columns"]))
        tc = t2_ours(name)
        if tc:
            r["t_contrast_K"] = tc
            r["t_contrast_change_K_per_K"] = (tc[1] - tc[0]) / P["dT"]
        maps[name] = strip(r)
        res["runs"][name] = r
        s = r["regimes"]["ScBL"]
        print(f"{name}: repro G {r['reproduce_max_abs_G_diff']:.1e}, vs json {r['reproduce_vs_robustness_json_pct']:.1e} %; ScBL "
              + "  ".join(f"{k} {v:+.2f}" for k, v in s["change_pct"].items()), flush=True)
        del P
    for model in sorted(stored["cmip6"]):
        P, lat, lon, dom = load_cmip6(model, pat)
        r = attribution(P["ctl"], P["fut"], lat, lon, True, dom, pat, rng, P["dT"])
        r["dT"] = P["dT"]
        st = stored["cmip6"][model]["change_pct"]
        r["reproduce_vs_robustness_json_pct"] = float(max(abs(r["regimes"][c]["change_pct"]["total"] - st[c]) for c in res["columns"]))
        tc = t2_cmip6(model)
        r["t_contrast_K"] = tc
        r["t_contrast_change_K_per_K"] = (tc[1] - tc[0]) / P["dT"]
        maps[model] = strip(r)
        res["cmip6"][model] = r
        s = r["regimes"]["ScBL"]
        print(f"  {model:18s} repro {r['reproduce_vs_robustness_json_pct']:.1e} %; ScBL "
              + "  ".join(f"{k} {v:+.2f}" for k, v in s["change_pct"].items()), flush=True)
        del P
    res["cmip6_ensemble"] = ensemble(res)
    res["correlations"] = correlations(res)
    plots(res, maps)
    json.dump(res, open(f"{OUT}/regime_G_attribution.json", "w"), indent=1)
    report(res)


def ensemble(res):
    E = {}
    M = list(res["cmip6"].values())
    for nm in res["columns"]:
        e = {}
        for unit in ("change_pct", "change_per_K"):
            for k in ("total", "background", "anomaly", "interaction"):
                x = np.array([m["regimes"][nm][unit][k] for m in M])
                e[f"{unit}|{k}"] = dict(median=float(np.median(x)), min=float(x.min()), max=float(x.max()),
                                        n_negative=int((x < 0).sum()))
        for k in ("G_anomaly_only_change_pct", "G_clim_only_change_pct"):
            x = np.array([m["regimes"][nm][k] for m in M])
            e[k] = dict(median=float(np.median(x)), min=float(x.min()), max=float(x.max()), n_negative=int((x < 0).sum()))
        if nm != "all":
            for k in ("peak_hPa", "mean_abs_hPa", "area_above_X_1e6km2", "area_above_half_peak_1e6km2",
                      "half_peak_centroid_lat", "half_peak_centroid_lon"):
                a = np.array([m["regimes"][nm]["composite"]["ctl"][k] for m in M])
                b = np.array([m["regimes"][nm]["composite"]["fut"][k] for m in M])
                d = 100 * (b / a - 1) if not k.startswith("half_peak_centroid") else b - a
                e[f"composite|{k}"] = dict(ctl_median=float(np.median(a)), fut_median=float(np.median(b)),
                                           change_median=float(np.median(d)), change_min=float(d.min()), change_max=float(d.max()),
                                           n_negative=int((d < 0).sum()))
            x = np.array([m["regimes"][nm]["composite"]["G_of_composite_anomaly_change_pct"] for m in M])
            e["G_of_composite_anomaly_change_pct"] = dict(median=float(np.median(x)), min=float(x.min()), max=float(x.max()), n_negative=int((x < 0).sum()))
            x = np.array([m["regimes"][nm]["daily"]["q_change_pct"] for m in M])
            e["daily_q_change_pct_median"] = np.median(x, 0).tolist()
            x = np.array([100 * (m["regimes"][nm]["daily"]["amp_mean"][1] / m["regimes"][nm]["daily"]["amp_mean"][0] - 1) for m in M])
            e["amp_mean_change_pct"] = dict(median=float(np.median(x)), min=float(x.min()), max=float(x.max()), n_negative=int((x < 0).sum()))
            for k in ("G_resid_change_pct", "rms_resid_EU_change_pct", "rms_anom_EU_change_pct"):
                x = np.array([m["regimes"][nm]["residual"][k] for m in M])
                e[f"residual|{k}"] = dict(median=float(np.median(x)), min=float(x.min()), max=float(x.max()), n_negative=int((x < 0).sum()))
        E[nm] = e
    for k in ("atl_eu", "eu"):
        x = np.array([m["background"]["dp40_65_change_per_K"][k] for m in M])
        E[f"background|dp40_65_change_per_K|{k}"] = dict(median=float(np.median(x)), min=float(x.min()), max=float(x.max()), n_negative=int((x < 0).sum()))
    x = np.array([m["background"]["G_wintermean_field_change_per_K"] for m in M])
    E["background|G_wintermean_field_change_per_K"] = dict(median=float(np.median(x)), min=float(x.min()), max=float(x.max()), n_negative=int((x < 0).sum()))
    x = np.array([m["t_contrast_change_K_per_K"] for m in M])
    E["t_contrast_change_K_per_K"] = dict(median=float(np.median(x)), min=float(x.min()), max=float(x.max()))
    return E


def predictors(m):
    return {"dp40_65_atl_eu_per_K": m["background"]["dp40_65_change_per_K"]["atl_eu"],
            "dp40_65_eu_per_K": m["background"]["dp40_65_change_per_K"]["eu"],
            "G_wintermean_field_per_K": m["background"]["G_wintermean_field_change_per_K"],
            "t_contrast_K_per_K": m.get("t_contrast_change_K_per_K", np.nan),
            "ScBL_peak_change_pct": 100 * (m["regimes"]["ScBL"]["composite"]["fut"]["peak_hPa"] / m["regimes"]["ScBL"]["composite"]["ctl"]["peak_hPa"] - 1),
            "ScBL_area_half_peak_change_pct": 100 * (m["regimes"]["ScBL"]["composite"]["fut"]["area_above_half_peak_1e6km2"] / m["regimes"]["ScBL"]["composite"]["ctl"]["area_above_half_peak_1e6km2"] - 1)}


def correlations(res):
    out = {}
    targets = {"ScBL_total_per_K": lambda m: m["regimes"]["ScBL"]["change_per_K"]["total"],
               "ScBL_background_per_K": lambda m: m["regimes"]["ScBL"]["change_per_K"]["background"],
               "ScBL_anomaly_per_K": lambda m: m["regimes"]["ScBL"]["change_per_K"]["anomaly"],
               "all_total_per_K": lambda m: m["regimes"]["all"]["change_per_K"]["total"]}
    for setname, M in (("cmip6", list(res["cmip6"].values())), ("cmip6_plus_runs", list(res["cmip6"].values()) + list(res["runs"].values()))):
        o = {}
        for tn, tf in targets.items():
            y = np.array([tf(m) for m in M])
            for pn in predictors(M[0]):
                x = np.array([predictors(m)[pn] for m in M])
                ok = np.isfinite(x)
                if ok.sum() < 5:
                    continue
                pr = stats.pearsonr(x[ok], y[ok])
                sr = stats.spearmanr(x[ok], y[ok])
                sl = np.polyfit(x[ok], y[ok], 1)[0]
                o[f"{tn}~{pn}"] = dict(n=int(ok.sum()), pearson_r=float(pr[0]), p=float(pr[1]), spearman_r=float(sr[0]), slope=float(sl))
        out[setname] = o
    return out


# ---------------------------------------------------------------- plots
TERMS = ("background", "anomaly", "interaction")
TCOL = {"background": "#1f77b4", "anomaly": "#d62728", "interaction": "#7f7f7f", "total": "k"}


def plots(res, maps):
    COLS = res["columns"]
    # 1: decomposition, % per K, runs + CMIP6
    fig, axs = plt.subplots(1, len(COLS), figsize=(15, 4.4), sharey=True)
    ents = RUNS + ["CMIP6"]
    for ax, nm in zip(axs, COLS):
        for i, e in enumerate(ents):
            for t, term in enumerate(TERMS):
                x = i + (t - 1) * 0.22
                if e == "CMIP6":
                    v = np.array([m["regimes"][nm]["change_per_K"][term] for m in res["cmip6"].values()])
                    ax.bar(x, np.median(v), 0.2, color=TCOL[term], alpha=0.5)
                    ax.plot(np.full(len(v), x) + np.linspace(-0.06, 0.06, len(v)), v, ".", color=TCOL[term], ms=3)
                else:
                    r = res["runs"][e]["regimes"][nm]
                    fac = r["change_per_K"]["total"] / r["change_pct"]["total"] if abs(r["change_pct"]["total"]) > 1e-3 else 1 / res["runs"][e]["dT"]
                    b = np.array(r["boot_p05_p95"][term]) * fac
                    ax.bar(x, r["change_per_K"][term], 0.2, color=TCOL[term])
                    ax.plot([x, x], b, color="k", lw=0.8)
            tot = np.median([m["regimes"][nm]["change_per_K"]["total"] for m in res["cmip6"].values()]) if e == "CMIP6" \
                else res["runs"][e]["regimes"][nm]["change_per_K"]["total"]
            ax.plot([i - 0.35, i + 0.35], [tot, tot], color="k", lw=1.6)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xticks(range(len(ents)), [LAB.get(e, e) for e in ents], rotation=30, fontsize=8)
        ax.set_title(nm if nm != "all" else "all days", fontsize=10)
    axs[0].set_ylabel("change in mean G (% per K)")
    hs = [plt.Rectangle((0, 0), 1, 1, color=TCOL[t]) for t in TERMS] + [plt.Line2D([], [], color="k", lw=1.6)]
    fig.legend(hs, ["background (seasonal cycle)", "anomaly (weather)", "interaction", "total"], loc="lower center", ncol=4, fontsize=8, frameon=False)
    fig.text(0.5, 0.065, "runs: bars with winter-bootstrap 5-95 %; CMIP6: median bar, dots = 20 models", ha="center", fontsize=7.5)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(f"{FIG}/decomposition_per_K.png", dpi=150)
    plt.close(fig)

    # 2: ScBL composite maps
    import cartopy.crs as ccrs
    pat = np.load(f"{ra.OUT}/regime_patterns.npz")
    lat, lon = pat["lat"], pat["lon"]
    rows = RUNS + ["CMIP6 mean"]
    cm = {k: np.mean([maps[m]["ScBL"][k] for m in res["cmip6"]], 0) for k in ("Ac", "Af", "Bc", "Bf")}
    cm_dM = np.mean([maps[m]["_bg"]["dM"] for m in res["cmip6"]], 0)
    fig, axs = plt.subplots(len(rows), 3, figsize=(12, 2.5 * len(rows)), subplot_kw=dict(projection=ccrs.PlateCarree()))
    lv = np.arange(-20, 21, 2.5)
    for i, e in enumerate(rows):
        mp = cm if e.startswith("CMIP6") else maps[e]["ScBL"]
        dM = cm_dM if e.startswith("CMIP6") else maps[e]["_bg"]["dM"]
        dT = np.median([m["dT"] for m in res["cmip6"].values()]) if e.startswith("CMIP6") else res["runs"][e]["dT"]
        panels = ((mp["Ac"], None, "control ScBL anomaly (hPa)"),
                  (mp["Af"] - mp["Ac"], mp["Ac"], "ScBL anomaly change; contours: control +-5, 10"),
                  (dM, None, "winter-mean msl change (hPa)"))
        for k, (fld, cont, tt) in enumerate(panels):
            ax = axs[i, k]
            ax.set_extent([-80, 40, 20, 80], ccrs.PlateCarree())
            if k == 0:
                im0 = ax.contourf(lon, lat, fld, lv, cmap="RdBu_r", extend="both")
            else:
                im1 = ax.contourf(lon, lat, fld, np.arange(-4, 4.1, 0.5), cmap="PuOr_r", extend="both")
            if cont is not None:
                ax.contour(lon, lat, cont, [-10, -5, 5, 10], colors="k", linewidths=0.6)
            ax.coastlines(lw=0.4)
            for r in EU:
                la0, la1, lo0, lo1 = gs.BOX[r]
                ax.plot([lo0, lo1, lo1, lo0, lo0], [la0, la0, la1, la1, la0], color="g", lw=0.6)
            if i == 0:
                ax.set_title(tt, fontsize=9)
            if k == 0:
                ax.text(-0.04, 0.5, f"{LAB.get(e, e)}\n(dT {dT:.1f} K)", transform=ax.transAxes, rotation=90, va="center", ha="right", fontsize=9)
    fig.subplots_adjust(left=0.05, right=0.98, top=0.96, bottom=0.09, wspace=0.05, hspace=0.08)
    c0 = fig.add_axes([0.08, 0.045, 0.26, 0.012])
    cb = fig.colorbar(im0, cax=c0, orientation="horizontal")
    cb.ax.tick_params(labelsize=7)
    cb.set_label("ScBL msl anomaly, own seasonal cycle (hPa)", fontsize=8)
    c1 = fig.add_axes([0.45, 0.045, 0.45, 0.012])
    cb = fig.colorbar(im1, cax=c1, orientation="horizontal")
    cb.ax.tick_params(labelsize=7)
    cb.set_label("future minus control (hPa); green: European boxes", fontsize=8)
    fig.savefig(f"{FIG}/scbl_composites.png", dpi=130)
    plt.close(fig)

    # 3: ScBL diagnostics: daily quantile change, amplitude/area changes
    fig, axs = plt.subplots(1, 3, figsize=(16, 5), gridspec_kw=dict(width_ratios=[1, 1.5, 0.9]))
    ax = axs[0]
    q = res["runs"][RUNS[0]]["regimes"]["ScBL"]["daily"]["q"]
    for m in res["cmip6"].values():
        ax.plot(q, m["regimes"]["ScBL"]["daily"]["q_change_pct"], color="0.8", lw=0.7)
    ax.plot(q, res["cmip6_ensemble"]["ScBL"]["daily_q_change_pct_median"], color="0.3", lw=2, label="CMIP6 median")
    for i, e in enumerate(RUNS):
        ax.plot(q, res["runs"][e]["regimes"]["ScBL"]["daily"]["q_change_pct"], "o-", color=f"C{i}", ms=3, label=LAB[e])
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("quantile of daily Europe-mean G under ScBL (%)")
    ax.set_ylabel("future vs control (%)")
    ax.legend(fontsize=7)
    ax.set_title("daily distribution", fontsize=9)
    # composite metrics change (%), ScBL
    keys = [("peak_hPa", "peak"), ("mean_abs_hPa", "mean |a|"), ("area_above_X_1e6km2", f"area > {X_HPA:.0f} hPa"),
            ("area_above_half_peak_1e6km2", "area > peak/2")]
    ax = axs[1]
    for k, (key, lb) in enumerate(keys):
        v = np.array([100 * (m["regimes"]["ScBL"]["composite"]["fut"][key] / m["regimes"]["ScBL"]["composite"]["ctl"][key] - 1) for m in res["cmip6"].values()])
        ax.boxplot(v, positions=[k], widths=0.5, showfliers=False)
        ax.plot(np.full(len(v), k) + np.linspace(-0.15, 0.15, len(v)), v, ".", color="0.5", ms=3)
        for i, e in enumerate(RUNS):
            c = res["runs"][e]["regimes"]["ScBL"]["composite"]
            ax.plot(k + 0.3, 100 * (c["fut"][key] / c["ctl"][key] - 1), "o", color=f"C{i}", ms=5)
    extras = [(("composite", "G_of_composite_anomaly_change_pct"), "G(composite a)"), ((None, "G_anomaly_only_change_pct"), "mean daily G(a)"),
              (("residual", "G_resid_change_pct"), "mean daily G(e)"), (("residual", "rms_resid_EU_change_pct"), "rms e, boxes")]
    for extra, ((grp, key), lb) in enumerate(extras):
        k = len(keys) + extra
        get = (lambda m: m["regimes"]["ScBL"][key]) if grp is None else (lambda m: m["regimes"]["ScBL"][grp][key])
        v = np.array([get(m) for m in res["cmip6"].values()])
        rv = [get(res["runs"][e]) for e in RUNS]
        ax.boxplot(v, positions=[k], widths=0.5, showfliers=False)
        ax.plot(np.full(len(v), k) + np.linspace(-0.15, 0.15, len(v)), v, ".", color="0.5", ms=3)
        for i, x in enumerate(rv):
            ax.plot(k + 0.3, x, "o", color=f"C{i}", ms=5)
    ax.set_xticks(range(len(keys) + len(extras)), [lb for _, lb in keys] + [lb for _, lb in extras], rotation=35, ha="right", fontsize=8)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_ylabel("future vs control (%)")
    ax.set_title("ScBL anomaly a (own cycle); residual e = a - composite", fontsize=9)
    ax.legend([plt.Line2D([], [], marker="o", ls="", color=f"C{i}") for i in range(len(RUNS))] + [plt.Line2D([], [], marker=".", ls="", color="0.5")],
              [LAB[e] for e in RUNS] + ["CMIP6"], fontsize=7, loc="upper right", ncol=2)
    # centroid displacement
    ax = axs[2]
    for m in res["cmip6"].values():
        c = m["regimes"]["ScBL"]["composite"]
        ax.annotate("", (c["fut"]["half_peak_centroid_lon"], c["fut"]["half_peak_centroid_lat"]),
                    (c["ctl"]["half_peak_centroid_lon"], c["ctl"]["half_peak_centroid_lat"]),
                    arrowprops=dict(arrowstyle="->", color="0.6", lw=0.8))
    for i, e in enumerate(RUNS):
        c = res["runs"][e]["regimes"]["ScBL"]["composite"]
        ax.annotate("", (c["fut"]["half_peak_centroid_lon"], c["fut"]["half_peak_centroid_lat"]),
                    (c["ctl"]["half_peak_centroid_lon"], c["ctl"]["half_peak_centroid_lat"]),
                    arrowprops=dict(arrowstyle="->", color=f"C{i}", lw=1.6))
        ax.plot([], [], color=f"C{i}", label=LAB[e])
    ax.plot([], [], color="0.6", label="CMIP6")
    for r in EU:
        la0, la1, lo0, lo1 = gs.BOX[r]
        ax.plot([lo0, lo1, lo1, lo0, lo0], [la0, la0, la1, la1, la0], color="g", lw=0.6)
        if 8 <= lo0 + 0.3 <= 21 and 57 <= la1 - 0.3 <= 63.5:
            ax.text(lo0 + 0.3, la1 - 0.3, r, fontsize=7, ha="left", va="top", color="g", clip_on=True)
    ax.set_xlim(8, 22)
    ax.set_ylim(56, 63.5)
    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    ax.legend(fontsize=7, loc="lower right")
    ax.set_title("ScBL high: half-peak-area centroid, ctl -> fut", fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{FIG}/scbl_amplitude_scale.png", dpi=150)
    plt.close(fig)

    # 4: cross-model scatter
    prs = [("dp40_65_atl_eu_per_K", "p(40N)-p(65N), 80W-40E (% per K)"), ("G_wintermean_field_per_K", "G of winter-mean field, boxes (% per K)"),
           ("t_contrast_K_per_K", "T(SC) - T(IB,IT) change (K per K)")]
    fig, axs = plt.subplots(2, 3, figsize=(13, 7.6))
    for row, (tn, tl) in enumerate((("total", "ScBL G change, total (% per K)"), ("background", "ScBL G change, background term (% per K)"))):
        for ax, (pn, pl) in zip(axs[row], prs):
            x = [predictors(m)[pn] for m in res["cmip6"].values()]
            y = [m["regimes"]["ScBL"]["change_per_K"][tn] for m in res["cmip6"].values()]
            ax.plot(x, y, "o", color="0.5", ms=4, label="CMIP6")
            for i, e in enumerate(RUNS):
                xv = predictors(res["runs"][e])[pn]
                if np.isfinite(xv):
                    ax.plot(xv, res["runs"][e]["regimes"]["ScBL"]["change_per_K"][tn], "s", color=f"C{i}", ms=7, label=LAB[e])
            cr = res["correlations"]["cmip6"].get(f"ScBL_{tn}_per_K~{pn}")
            cr2 = res["correlations"]["cmip6_plus_runs"].get(f"ScBL_{tn}_per_K~{pn}")
            ax.set_title(f"CMIP6 r={cr['pearson_r']:+.2f} (p={cr['p']:.2f}); +runs r={cr2['pearson_r']:+.2f} (n={cr2['n']})", fontsize=8.5)
            ax.axhline(0, color="k", lw=0.4)
            ax.axvline(0, color="k", lw=0.4)
            ax.set_xlabel(pl, fontsize=8.5)
            if ax is axs[row, 0]:
                ax.set_ylabel(tl, fontsize=8.5)
    axs[0, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(f"{FIG}/cross_model_background.png", dpi=150)
    plt.close(fig)

    # 5: sector-mean winter msl profile change
    fig, ax = plt.subplots(figsize=(6, 4))
    for m in res["cmip6"].values():
        a, b = np.array(m["background"]["zonal_profile_hPa"])
        ax.plot(lat, (b - a) / m["dT"], color="0.8", lw=0.7)
    ax.plot(lat, np.median([(np.subtract(*m["background"]["zonal_profile_hPa"][::-1])) / m["dT"] for m in res["cmip6"].values()], 0),
            color="0.3", lw=2, label="CMIP6 median")
    for i, e in enumerate(RUNS):
        a, b = np.array(res["runs"][e]["background"]["zonal_profile_hPa"])
        ax.plot(lat, (b - a) / res["runs"][e]["dT"], color=f"C{i}", lw=1.5, label=LAB[e])
    ax.axhline(0, color="k", lw=0.5)
    for la in (40, 65):
        ax.axvline(la, color="k", lw=0.4, ls=":")
    ax.set_xlabel("latitude")
    ax.set_ylabel("winter-mean msl change, 80W-40E mean (hPa per K)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(f"{FIG}/sector_msl_profile_change.png", dpi=150)
    plt.close(fig)


def report(res):
    print("\nScBL G change: total = background + anomaly + interaction   (% ; % per K)")
    for name, r in list(res["runs"].items()) + [("CMIP6 median", None)]:
        if r is None:
            e = res["cmip6_ensemble"]["ScBL"]
            print(f"  {name:13s} " + "  ".join(f"{k} {e['change_pct|' + k]['median']:+6.2f} [{e['change_pct|' + k]['min']:+.1f},{e['change_pct|' + k]['max']:+.1f}]"
                                               for k in ("total", "background", "anomaly", "interaction")))
            print(f"  {'':13s} " + "  ".join(f"{k} {e['change_per_K|' + k]['median']:+6.2f} [{e['change_per_K|' + k]['min']:+.2f},{e['change_per_K|' + k]['max']:+.2f}] n<0 {e['change_per_K|' + k]['n_negative']}"
                                               for k in ("total", "background", "anomaly", "interaction")))
            continue
        s = r["regimes"]["ScBL"]
        print(f"  {name:13s} " + "  ".join(f"{k} {s['change_pct'][k]:+6.2f} [{s['boot_p05_p95'][k][0]:+.1f},{s['boot_p05_p95'][k][1]:+.1f}]" for k in s["change_pct"]))
        print(f"  {'':13s} " + "  ".join(f"{k} {s['change_per_K'][k]:+6.2f}" for k in s["change_per_K"])
              + f"   other base bg {s['other_base_pct']['background']:+.2f} an {s['other_base_pct']['anomaly']:+.2f}; own labels tot {s['ownclim_labels_pct']['total']:+.2f} an {s['ownclim_labels_pct']['anomaly']:+.2f}")
    for nm in res["columns"]:
        print(f"\n{nm}: " + "; ".join(f"{n} bg {r['regimes'][nm]['change_per_K']['background']:+.2f} an {r['regimes'][nm]['change_per_K']['anomaly']:+.2f} int {r['regimes'][nm]['change_per_K']['interaction']:+.2f}"
                                    for n, r in res["runs"].items())
              + f"; CMIP6 med bg {res['cmip6_ensemble'][nm]['change_per_K|background']['median']:+.2f} an {res['cmip6_ensemble'][nm]['change_per_K|anomaly']['median']:+.2f}")
    print("\nScBL composite (own-cycle anomaly), ctl -> fut")
    for name, r in res["runs"].items():
        c, d = r["regimes"]["ScBL"]["composite"], r["regimes"]["ScBL"]["daily"]
        print(f"  {name:10s} peak {c['ctl']['peak_hPa']:.1f}->{c['fut']['peak_hPa']:.1f} at ({c['ctl']['peak_lat']:.1f},{c['ctl']['peak_lon']:.1f})->({c['fut']['peak_lat']:.1f},{c['fut']['peak_lon']:.1f})"
              f"  mean|a| {c['ctl']['mean_abs_hPa']:.2f}->{c['fut']['mean_abs_hPa']:.2f}  area>X {c['ctl']['area_above_X_1e6km2']:.1f}->{c['fut']['area_above_X_1e6km2']:.1f}"
              f"  area>half {c['ctl']['area_above_half_peak_1e6km2']:.1f}->{c['fut']['area_above_half_peak_1e6km2']:.1f}"
              f"  centroid ({c['ctl']['half_peak_centroid_lat']:.1f},{c['ctl']['half_peak_centroid_lon']:.1f})->({c['fut']['half_peak_centroid_lat']:.1f},{c['fut']['half_peak_centroid_lon']:.1f})"
              f"  G(compA) {c['G_of_composite_anomaly_change_pct']:+.1f}%  G(compFull) {c['G_of_composite_full_change_pct']:+.1f}%  r {c['pattern_corr']:.2f}"
              f"  G(a) {r['regimes']['ScBL']['G_anomaly_only_change_pct']:+.1f}%  G(C) {r['regimes']['ScBL']['G_clim_only_change_pct']:+.1f}%"
              f"  amp {d['amp_mean'][0]:.2f}->{d['amp_mean'][1]:.2f} sd {d['amp_sd'][0]:.2f}->{d['amp_sd'][1]:.2f}  q% {['%+.1f' % v for v in d['q_change_pct']]}"
              f"  <p10 {d['frac_below_ctl_p10'][1]:.3f} >p90 {d['frac_above_ctl_p90'][1]:.3f}")
        rs = r["regimes"]["ScBL"]["residual"]
        print(f"  {'':10s} residual e: G(e) {rs['G_resid_change_pct']:+.1f}% (ctl {rs['G_resid_ctl_mean']:.2f} m/s of G(a) {r['regimes']['ScBL']['G_anom_only_ctl_mean']:.2f})"
              f"  rms e boxes {rs['rms_resid_EU_hPa'][0]:.2f}->{rs['rms_resid_EU_hPa'][1]:.2f} ({rs['rms_resid_EU_change_pct']:+.1f}%)"
              f"  rms a boxes {rs['rms_anom_EU_hPa'][0]:.2f}->{rs['rms_anom_EU_hPa'][1]:.2f} ({rs['rms_anom_EU_change_pct']:+.1f}%)")
        print(f"  {'':10s} box anom ctl " + " ".join(f"{k} {v:+.1f}" for k, v in c["box_anom_ctl_hPa"].items())
              + " | fut " + " ".join(f"{k} {v:+.1f}" for k, v in c["box_anom_fut_hPa"].items()))
    e = res["cmip6_ensemble"]["ScBL"]
    print("  CMIP6 " + "; ".join(f"{k.split('|')[1]} ctl {v['ctl_median']:.2f} fut {v['fut_median']:.2f} chg med {v['change_median']:+.2f} [{v['change_min']:+.1f},{v['change_max']:+.1f}] n<0 {v['n_negative']}"
                                  for k, v in e.items() if k.startswith("composite|")))
    print(f"  CMIP6 G(composite a) {e['G_of_composite_anomaly_change_pct']}; G(a) {e['G_anomaly_only_change_pct']}; G(C) {e['G_clim_only_change_pct']}; amp {e['amp_mean_change_pct']}")
    print(f"  CMIP6 daily q change median {['%+.1f' % v for v in e['daily_q_change_pct_median']]}")
    print("  CMIP6 " + "; ".join(f"{k}: {v}" for k, v in e.items() if k.startswith("residual|")))
    for nm in ra.NAMES_ORDER:
        print(f"  {nm:5s} residual G(e) %: " + " ".join(f"{n} {r['regimes'][nm]['residual']['G_resid_change_pct']:+.1f}" for n, r in res["runs"].items())
              + f" CMIP6 med {res['cmip6_ensemble'][nm]['residual|G_resid_change_pct']['median']:+.1f}; G(composite a) %: "
              + " ".join(f"{n} {r['regimes'][nm]['composite']['G_of_composite_anomaly_change_pct']:+.1f}" for n, r in res["runs"].items())
              + f" CMIP6 med {res['cmip6_ensemble'][nm]['G_of_composite_anomaly_change_pct']['median']:+.1f}")
    print("\nBackground: p40-p65 (hPa) ctl->fut, % per K (atl_eu, eu); G of winter-mean field % per K; T contrast K, K per K")
    for name, r in res["runs"].items():
        b = r["background"]
        print(f"  {name:10s} atl_eu {b['dp40_65_hPa']['atl_eu'][0]:.2f}->{b['dp40_65_hPa']['atl_eu'][1]:.2f} ({b['dp40_65_change_per_K']['atl_eu']:+.2f}%/K)"
              f"  eu {b['dp40_65_hPa']['eu'][0]:.2f}->{b['dp40_65_hPa']['eu'][1]:.2f} ({b['dp40_65_change_per_K']['eu']:+.2f}%/K)  G(mean) {b['G_wintermean_field_change_per_K']:+.2f}%/K"
              + (f"  T {r['t_contrast_K'][0]:.2f}->{r['t_contrast_K'][1]:.2f} ({r['t_contrast_change_K_per_K']:+.3f} K/K)" if "t_contrast_K" in r else "  T -"))
    E = res["cmip6_ensemble"]
    print("  CMIP6 " + "; ".join(f"{k}: {v}" for k, v in E.items() if k.startswith("background") or k.startswith("t_contrast")))
    print("\nCorrelations")
    for s, o in res["correlations"].items():
        for k, v in o.items():
            print(f"  {s:16s} {k:45s} n {v['n']:2d} r {v['pearson_r']:+.2f} p {v['p']:.3f} rho {v['spearman_r']:+.2f} slope {v['slope']:+.3f}")
    print(f"\nSaved {OUT}/regime_G_attribution.json and {FIG}/*.png")


if __name__ == "__main__":
    main()
