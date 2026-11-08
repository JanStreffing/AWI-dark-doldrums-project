"""
Geostrophic against surface wind: is the increase in low-wind events a weakening of the large-scale
pressure gradients (circulation) or a change in how much of the free-atmosphere flow reaches the
surface (transfer: boundary-layer stability, ageostrophic flow)?

Per region and day:
  G   region mean of the geostrophic speed |grad p| / (rho f) from daily-mean msl on the 2.5 deg grid
      (centred differences, rho = 1.25 kg m-3), grid points inside the region box (DE: DE_BOX)
  s   daily-mean region-mean wind speed (our runs: log-law 100 m speed; CMIP6: sfcWind)
  a_m control-period transfer ratio sum(s) / sum(G) per calendar month
  s*  = a_m G: the surface speed the region would have if the transfer stayed at its control value
Events on the daily footing (48 h = 2 days, Mockert window expansion, Sep-Mar winters), each series
with its own threshold matched to TARGET events per control winter, in speed (not CF) space.
  R_s    future / control events of s       (the actual change, speed space)
  R_G    future / control events of s*      (circulation part)
  R_s / R_G                                 (transfer part)
Field means of ln R over regions; share of the change from circulation = mean ln R_G / mean ln R_s.

Within the Euro-Atlantic regimes (ERA5 EOFs and centroids from regime_analysis.py, anomalies from
the control seasonal cycle): per regime, the change in mean G, in the transfer ratio mean(s)/mean(G)
and in the frequency of low days (daily s, and s*, below the control 10th percentile).

Sources:
  ours   msl from data/regimes (20-80 N, 80 W-40 E; European regions only), speeds from data/regions
  CMIP6  psl from extract_psl.py (global 15-80 N, so every region from 29 N up), sfcWind from
         daily_cmip6_*.npz; historical 1951-1980 against ssp585 2071-2100

Usage: python3 geostrophic.py
Output: dunkelflaute/results/regimes/geostrophic_summary.json
"""
import glob
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("regions", "analysis", "regimes", "cmip6"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/{p}")
import regions as rg                       # noqa: E402
import region_events as ev                 # noqa: E402
import regime_analysis as ra               # noqa: E402
import region_per_K as rpk                 # noqa: E402

ev.WIN = 2
OUT = f"{_ROOT}/dunkelflaute/results/regimes"
DC = f"{_ROOT}/dunkelflaute/data/cmip6"
TARGET = 5.0
MIN_DAYS = 200
RHO, OMEGA, RE = 1.25, 7.292e-5, 6.371e6
EU = ra.EU
BOX = dict(rg.BOXES, DE=rg.DE_BOX)
CMIP_REGIONS = [r for r in rg.REGIONS if BOX[r][0] >= 29.0]
OURS = {"9km": "dart9km", "31km": "tco319", "100km": "tco95", "eerie_ctl": "eerie", "eerie_hist": "eerie"}
LOWQ = 10


def ymd_of(days):
    y = days.astype("datetime64[Y]").astype(int) + 1970
    m = days.astype("datetime64[M]").astype(int) % 12 + 1
    d = (days - days.astype("datetime64[M]")).astype(int) + 1
    return y * 10000 + m * 100 + d


def winter_of(k):
    y, m = k // 10000, k // 100 % 100
    return np.where(m >= 9, y + 1, y)


def dates_of(k):
    """ymd ints -> datetime64[D] for the seasonal harmonics (360-day calendars: day clipped)."""
    y, m, d = k // 10000, k // 100 % 100, k % 100
    ml = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])[m - 1]
    return np.array([np.datetime64(f"{a:04d}-{b:02d}-{min(c, e):02d}") for a, b, c, e in zip(y, m, d, ml)])


def geostrophic(p, lat, lon, periodic):
    """p (D, ny, nx) Pa -> |V_g| (D, ny, nx) m/s, NaN where a centred difference is unavailable."""
    phi = np.deg2rad(lat)[None, :, None]
    dphi = np.deg2rad(lat[1] - lat[0])
    dlam = np.deg2rad(lon[1] - lon[0])
    dpy = np.full_like(p, np.nan)
    dpy[:, 1:-1] = (p[:, 2:] - p[:, :-2]) / (2 * dphi * RE)
    if periodic:
        dpx = (np.roll(p, -1, axis=2) - np.roll(p, 1, axis=2)) / (2 * dlam * RE * np.cos(phi))
    else:
        dpx = np.full_like(p, np.nan)
        dpx[:, :, 1:-1] = (p[:, :, 2:] - p[:, :, :-2]) / (2 * dlam * RE * np.cos(phi))
    f = 2 * OMEGA * np.sin(phi)
    return np.hypot(dpx, dpy) / (RHO * f)


def region_G(vg, lat, lon, regions):
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    out = []
    for r in regions:
        la0, la1, lo0, lo1 = BOX[r]
        sel = (LA >= la0) & (LA <= la1) & (LO >= lo0) & (LO <= lo1)
        if not sel.any():          # box smaller than the grid: points within half a cell
            sel = (LA >= la0 - 1.25) & (LA <= la1 + 1.25) & (LO >= lo0 - 1.25) & (LO <= lo1 + 1.25)
        out.append(np.nanmean(vg[:, sel], axis=1))
    return np.stack(out, axis=1)


def winters(k, y0, y1):
    w = winter_of(k)
    return [w == y for y in range(y0, y1 + 1) if (w == y).sum() >= MIN_DAYS]


def matched_ratio(xc, kc, c0, c1, xf, kf, f0, f1):
    """Events per winter at the control-matched threshold, future/control."""
    mc, mf = winters(kc, c0, c1), winters(kf, f0, f1)
    grid = np.linspace(*np.nanpercentile(xc, [0.1, 40]), 400)
    cc = np.stack([ev.counts_grid(xc[m], grid) for m in mc]).mean(0)
    j = int(np.argmax(cc >= TARGET)) if (cc >= TARGET).any() else int(np.argmax(cc))
    thr = grid[j]
    a = np.array([ev.counts_grid(xc[m], np.array([thr]))[0] for m in mc])
    b = np.array([ev.counts_grid(xf[m], np.array([thr]))[0] for m in mf])
    return float(b.mean() / a.mean()), float(a.mean()), float(b.mean())


def analyse(name, kc, sc, Gc, kf, sf, Gf, c0, c1, f0, f1, regions, dT, labels=None):
    """s*, G (D, R) aligned on ymd ints kc / kf."""
    mon_c, mon_f = kc // 100 % 100, kf // 100 % 100
    per = {}
    for i, r in enumerate(regions):
        ok_c = np.isfinite(sc[:, i]) & np.isfinite(Gc[:, i])
        if ok_c.mean() < 0.99:
            continue
        a = {m: np.nansum(sc[mon_c == m, i]) / np.nansum(Gc[mon_c == m, i]) for m in rg.MONTHS}
        xc = np.array([a[m] for m in mon_c]) * Gc[:, i]
        xf = np.array([a[m] for m in mon_f]) * Gf[:, i]
        Rs, ec, ef = matched_ratio(sc[:, i], kc, c0, c1, sf[:, i], kf, f0, f1)
        RG, gc, gf = matched_ratio(xc, kc, c0, c1, xf, kf, f0, f1)
        sel_c = np.isin(winter_of(kc), range(c0, c1 + 1))
        sel_f = np.isin(winter_of(kf), range(f0, f1 + 1))
        ms_c, ms_f = np.nanmean(sc[sel_c, i]), np.nanmean(sf[sel_f, i])
        mg_c, mg_f = np.nanmean(Gc[sel_c, i]), np.nanmean(Gf[sel_f, i])
        row = dict(ratio_s=Rs, ratio_G=RG, ratio_transfer=Rs / RG, events_ctl=ec, events_fut=ef,
                   mean_s_change_pct=100 * (ms_f / ms_c - 1), mean_G_change_pct=100 * (mg_f / mg_c - 1),
                   transfer_change_pct=100 * ((ms_f / mg_f) / (ms_c / mg_c) - 1),
                   corr_s_G_ctl=float(np.corrcoef(sc[sel_c, i], Gc[sel_c, i])[0, 1]))
        if labels is not None and r in EU:
            lc, lf = labels
            qs = np.nanpercentile(sc[sel_c, i], LOWQ)
            qx = np.nanpercentile(xc[sel_c], LOWQ)
            reg = {}
            for k, nm in enumerate(ra.NAMES_ORDER):
                cc_, ff_ = sel_c & (lc == k), sel_f & (lf == k)
                if cc_.sum() < 30 or ff_.sum() < 30:
                    continue
                reg[nm] = dict(days_ctl=int(cc_.sum()), days_fut=int(ff_.sum()),
                               mean_G_change_pct=100 * (np.nanmean(Gf[ff_, i]) / np.nanmean(Gc[cc_, i]) - 1),
                               transfer_change_pct=100 * ((np.nanmean(sf[ff_, i]) / np.nanmean(Gf[ff_, i])) / (np.nanmean(sc[cc_, i]) / np.nanmean(Gc[cc_, i])) - 1),
                               low_s_frac_ctl=float(np.mean(sc[cc_, i] < qs)), low_s_frac_fut=float(np.mean(sf[ff_, i] < qs)),
                               low_sstar_frac_ctl=float(np.mean(xc[cc_] < qx)), low_sstar_frac_fut=float(np.mean(xf[ff_] < qx)))
            row["regimes"] = reg
        per[r] = row
    V = list(per)
    ls = np.array([np.log(per[r]["ratio_s"]) for r in V])
    lg = np.array([np.log(per[r]["ratio_G"]) for r in V])
    field = dict(n_regions=len(V), mean_ratio_s=float(np.exp(ls.mean())), mean_ratio_G=float(np.exp(lg.mean())),
                 mean_ratio_transfer=float(np.exp((ls - lg).mean())),
                 share_circulation=float(lg.mean() / ls.mean()) if ls.mean() else np.nan,
                 mean_G_change_pct=float(np.mean([per[r]["mean_G_change_pct"] for r in V])),
                 transfer_change_pct=float(np.mean([per[r]["transfer_change_pct"] for r in V])),
                 corr_s_G_ctl=float(np.mean([per[r]["corr_s_G_ctl"] for r in V])))
    if np.isfinite(dT) and dT > 0.3:
        field["pct_per_K_s"] = float(100 * (np.exp(ls.mean() / dT) - 1))
        field["pct_per_K_G"] = float(100 * (np.exp(lg.mean() / dT) - 1))
    eu = [r for r in V if r in EU]
    if eu:
        les = np.array([np.log(per[r]["ratio_s"]) for r in eu])
        leg = np.array([np.log(per[r]["ratio_G"]) for r in eu])
        field["europe"] = dict(n_regions=len(eu), mean_ratio_s=float(np.exp(les.mean())), mean_ratio_G=float(np.exp(leg.mean())),
                               share_circulation=float(leg.mean() / les.mean()) if les.mean() else np.nan)
        if labels is not None:
            agg = {}
            for nm in ra.NAMES_ORDER:
                rows = [per[r]["regimes"][nm] for r in eu if nm in per[r].get("regimes", {})]
                if rows:
                    agg[nm] = {k: float(np.mean([x[k] for x in rows])) for k in rows[0] if k not in ("days_ctl", "days_fut")}
            field["europe"]["regimes"] = agg
    return dict(dT=dT, field=field, regions=per)


def regime_labels(pat, dates_c, msl_c, dates_f, msl_f):
    coef = ra.seasonal_fit(dates_c, msl_c)
    w = pat["wts"].ravel()
    lc = ra.assign(ra.anomalies(dates_c, msl_c, coef) * w, pat["eofs"], pat["cent"])
    lf = ra.assign(ra.anomalies(dates_f, msl_f, coef) * w, pat["eofs"], pat["cent"])
    return lc, lf


def ours(T, pat):
    out = {}
    for name, mod in OURS.items():
        (cn, c0, c1), (fn, f0, f1) = ev.COMPARISONS[name]
        S = {}
        for run in (cn, fn):
            z = np.load(f"{rg.OUT_DIR}/{run}.npz")
            assert list(z["regions"]) == rg.REGIONS
            t = z["time"].astype("datetime64[s]")
            day = (t - np.timedelta64(1, "s")).astype("datetime64[D]")
            u, inv, cnt = np.unique(day, return_inverse=True, return_counts=True)
            s = np.zeros((len(u), len(rg.REGIONS)))
            np.add.at(s, inv, z["spd100"])
            full = cnt == cnt.max()
            s = s[full] / cnt[full, None]
            m = np.load(f"{ra.D}/msl_{run}.npz")
            md = m["dates"].astype("datetime64[D]")
            msl = m["msl"].astype(np.float64)
            k_s, k_m = ymd_of(u[full]), ymd_of(md)
            common, i_s, i_m = np.intersect1d(k_s, k_m, return_indices=True)
            keep = np.isin(common // 100 % 100, rg.MONTHS)
            common, i_s, i_m = common[keep], i_s[keep], i_m[keep]
            G = region_G(geostrophic(msl[i_m], m["lat"], m["lon"], False), m["lat"], m["lon"], EU)
            S[run] = (common, s[i_s][:, [rg.REGIONS.index(r) for r in EU]], G, md[i_m], msl[i_m])
        kc, sc, Gc, dc, mc = S[cn]
        kf, sf, Gf, df_, mf = S[fn]
        sel_c = np.isin(winter_of(kc), range(c0, c1 + 1))
        sel_f = np.isin(winter_of(kf), range(f0, f1 + 1))
        lab_c = np.full(len(kc), -1)
        lab_f = np.full(len(kf), -1)
        a, b = regime_labels(pat, dc[sel_c], mc[sel_c], df_[sel_f], mf[sel_f])
        lab_c[sel_c], lab_f[sel_f] = a, b
        dT = rpk.mean_years(T, fn, f0, f1)[0] - rpk.mean_years(T, cn, c0, c1)[0]
        out[name] = analyse(name, kc, sc, Gc, kf, sf, Gf, c0, c1, f0, f1, EU, dT, labels=(lab_c, lab_f))
        print(f"  {name}: done", flush=True)
    return out


def cmip6(pat):
    out = {}
    c0, c1, f0, f1 = 1951, 1980, 2071, 2100
    for fn in sorted(glob.glob(f"{DC}/psl_*.npz")):
        model = os.path.basename(fn)[4:-4]
        P = np.load(fn)
        Z = np.load(f"{DC}/daily_cmip6_{model}.npz")
        lat, lon = P["lat"], P["lon"]
        ri = [rg.REGIONS.index(r) for r in CMIP_REGIONS]
        ii = np.where(lat >= 20.0)[0]
        jj = np.where((lon >= -80.0) & (lon <= 40.0))[0]
        S = {}
        for exp, (y0, y1) in (("historical", (c0, c1)), ("ssp585", (f0, f1))):
            kp = P[f"{exp}|ymd"] @ np.array([10000, 100, 1])
            kz = Z[f"{exp}|ymd"] @ np.array([10000, 100, 1])
            common, i_p, i_z = np.intersect1d(kp, kz, return_indices=True)
            w = winter_of(common)
            keep = (w >= y0) & (w <= y1)
            common, i_p, i_z = common[keep], i_p[keep], i_z[keep]
            psl = P[f"{exp}|psl"][i_p].astype(np.float64)
            G = region_G(geostrophic(psl, lat, lon, True), lat, lon, CMIP_REGIONS)
            S[exp] = (common, Z[f"{exp}|spd10"][i_z][:, ri], G, psl[:, ii][:, :, jj])
        kc, sc, Gc, mc = S["historical"]
        kf, sf, Gf, mf = S["ssp585"]
        lab = regime_labels(pat, dates_of(kc), mc, dates_of(kf), mf)
        T = {e: dict(zip(Z[f"{e}|tas_years"].tolist(), Z[f"{e}|tas"].tolist())) for e in ("historical", "ssp585")}
        dT = np.mean([T["ssp585"][y] for y in range(f0, f1 + 1) if y in T["ssp585"]]) - np.mean([T["historical"][y] for y in range(c0, c1 + 1) if y in T["historical"]])
        out[model] = analyse(model, kc, sc, Gc, kf, sf, Gf, c0, c1, f0, f1, CMIP_REGIONS, float(dT), labels=lab)
        f = out[model]["field"]
        print(f"  {model:18s} dT {dT:.2f}  R_s {f['mean_ratio_s']:.2f}  R_G {f['mean_ratio_G']:.2f}  share {f['share_circulation']:.2f}", flush=True)
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    pat = np.load(f"{ra.OUT}/regime_patterns.npz")
    ra.NAMES_ORDER = [str(x) for x in pat["names"]]
    res = dict(ours=ours(rpk.global_t(), pat), cmip6=cmip6(pat))
    ens = {}
    # wind-only models (EC-Earth3, the IFS-family check) are reported but kept out of the 20-model ensemble
    inv = json.load(open(f"{DC}/inventory.json"))
    res["cmip6_windonly"] = {m: res["cmip6"].pop(m) for m in list(res["cmip6"]) if inv.get(m, {}).get("wind_only")}
    if res["cmip6"]:
        F = [v["field"] for v in res["cmip6"].values()]
        for key in ("mean_ratio_s", "mean_ratio_G", "mean_ratio_transfer", "share_circulation", "mean_G_change_pct",
                    "transfer_change_pct", "pct_per_K_s", "pct_per_K_G", "corr_s_G_ctl"):
            x = np.array([f.get(key, np.nan) for f in F])
            ens[key] = dict(median=float(np.nanmedian(x)), p10_p90=[float(v) for v in np.nanpercentile(x, [10, 90])])
        x = np.array([f["europe"]["share_circulation"] for f in F])
        ens["europe_share_circulation"] = dict(median=float(np.nanmedian(x)), p10_p90=[float(v) for v in np.nanpercentile(x, [10, 90])])
        for nm in ra.NAMES_ORDER:
            for key in ("mean_G_change_pct", "transfer_change_pct"):
                x = np.array([f["europe"]["regimes"][nm][key] for f in F if nm in f["europe"].get("regimes", {})])
                ens[f"europe_{nm}_{key}"] = float(np.median(x))
            for key in ("low_s", "low_sstar"):
                a = np.array([f["europe"]["regimes"][nm][f"{key}_frac_ctl"] for f in F if nm in f["europe"].get("regimes", {})])
                b = np.array([f["europe"]["regimes"][nm][f"{key}_frac_fut"] for f in F if nm in f["europe"].get("regimes", {})])
                ens[f"europe_{nm}_{key}_frac"] = [float(np.median(a)), float(np.median(b))]
    res["cmip6_ensemble"] = ens
    res["cmip6_regions"] = CMIP_REGIONS
    json.dump(res, open(f"{OUT}/geostrophic_summary.json", "w"), indent=1)

    print("\nOur runs, European regions: R_s (actual), R_G (circulation), share; dG mean, d transfer; corr(s,G)")
    for n, v in res["ours"].items():
        f = v["field"]
        print(f"  {n:11s} dT {v['dT']:.2f}  R_s {f['mean_ratio_s']:.2f}  R_G {f['mean_ratio_G']:.2f}  transfer {f['mean_ratio_transfer']:.2f}  "
              f"share {f['share_circulation']:.2f}  dG {f['mean_G_change_pct']:+.1f}%  dtransfer {f['transfer_change_pct']:+.1f}%  r {f['corr_s_G_ctl']:.2f}")
        for nm, x in f["europe"].get("regimes", {}).items():
            print(f"      {nm:5s} dG {x['mean_G_change_pct']:+5.1f}%  dtransfer {x['transfer_change_pct']:+5.1f}%  "
                  f"low s {x['low_s_frac_ctl']:.3f}->{x['low_s_frac_fut']:.3f}  low s* {x['low_sstar_frac_ctl']:.3f}->{x['low_sstar_frac_fut']:.3f}")
    if ens:
        print(f"\nCMIP6 ({len(res['cmip6'])} models, {len(CMIP_REGIONS)} regions): " + "; ".join(
            f"{k} {v['median']:.2f} [{v['p10_p90'][0]:.2f}, {v['p10_p90'][1]:.2f}]" for k, v in ens.items() if isinstance(v, dict)))
        for nm in ra.NAMES_ORDER:
            print(f"  Europe {nm:5s} dG {ens[f'europe_{nm}_mean_G_change_pct']:+5.1f}%  dtransfer {ens[f'europe_{nm}_transfer_change_pct']:+5.1f}%  "
                  f"low s {ens[f'europe_{nm}_low_s_frac'][0]:.3f}->{ens[f'europe_{nm}_low_s_frac'][1]:.3f}  low s* {ens[f'europe_{nm}_low_sstar_frac'][0]:.3f}->{ens[f'europe_{nm}_low_sstar_frac'][1]:.3f}")
    print(f"Saved {OUT}/geostrophic_summary.json")


if __name__ == "__main__":
    main()
