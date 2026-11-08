"""
Artifact checks for the rising transfer (mean s / mean G) in our runs, geostrophic.py (Section 13b).

transfer change = 100 ((mean s_f / mean G_f) / (mean s_c / mean G_c) - 1) per region, field mean over
regions, exactly as geostrophic.analyse (regions kept only if >= 99% of control days are finite).

  A1  CMIP6 restricted to the 8 European regions (per-region values of geostrophic_summary.json)
  A2  G construction: the box-mean regime msl (extract_msl.py, current) against msl bilinearly
      interpolated to the CMIP6 psl grid (transfer_checks_msl.py), first on the regime sector
      without periodicity (interpolation effect alone), then on the global grid as CMIP6
  A3  all 22 CMIP_REGIONS for our runs with the global bilinear G, against CMIP6 on the same regions
  A4  wind quantity: raw 10 m speed (spd10) against the log-law 100 m speed (spd100)
  A5  native 100 m speed of TCO95 (spd100_native, 100u/100v)

Inputs: data/regions/<run>.npz (spd100, spd10, spd100_native), data/regions/
transfer_checks_dart9km_<run>_spd10.npz (DART 10 m), data/regimes/msl_<run>.npz (current G),
data/regimes/transfer_checks_msl_<run>.npz (CMIP6-style G), results/regimes/geostrophic_summary.json.

Usage: python3 transfer_checks.py      (legacy python, login node, about 1 min)
Output: dunkelflaute/results/regimes/transfer_checks.json
"""
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regimes")
import geostrophic as gs     # noqa: E402  (imports regions, region_events, regime_analysis)

rg, ev, ra = gs.rg, gs.ev, gs.ra
D = f"{_ROOT}/dunkelflaute/data"
OUT = f"{_ROOT}/dunkelflaute/results/regimes/transfer_checks.json"
EU, CR = gs.EU, gs.CMIP_REGIONS
EU_FLAT = ["DE", "FR", "UK", "PL", "SC"]   # EU without the regions with high terrain (IB, IT, BK)
SECTOR = dict(lat=(20.0, 80.0), lon=(-80.0, 40.0))


def daily(t, x):
    """(T,) stamps, (T, R) -> unique days (only complete ones), daily means; a stamp belongs to t - 1 s."""
    day = (t.astype("datetime64[s]") - np.timedelta64(1, "s")).astype("datetime64[D]")
    u, inv, cnt = np.unique(day, return_inverse=True, return_counts=True)
    s = np.zeros((len(u), x.shape[1]))
    np.add.at(s, inv, x)
    full = cnt == cnt.max()
    return gs.ymd_of(u[full]), s[full] / cnt[full, None]


def speeds(run, var):
    if var == "spd10_cmip":          # our run through the CMIP6 daily extractor (extract_daily.py)
        fn = f"{D}/cmip6/daily_{run}.npz"
        if not os.path.exists(fn):
            return None
        z = np.load(fn)
        assert list(z["regions"]) == rg.REGIONS
        return z["run|ymd"] @ np.array([10000, 100, 1]), z["run|spd10"].astype(np.float64)
    if var == "spd10" and run.startswith("dart9km"):
        fn = f"{D}/regions/transfer_checks_dart9km_{run.split('_')[1]}_spd10.npz"
        if not os.path.exists(fn):
            return None
        z = np.load(fn)
        return daily(z["time"], z["spd10"])
    z = np.load(f"{D}/regions/{run}.npz")
    assert list(z["regions"]) == rg.REGIONS
    if var not in z.files:
        return None
    return daily(z["time"], z[var])


def G_of(run, how, regions):
    """-> (ymd, G (D, R)) for 'box' (current), 'bilin_sector', 'bilin_global'."""
    if how == "box":
        m = np.load(f"{ra.D}/msl_{run}.npz")
        lat, lon, msl, per = m["lat"], m["lon"], m["msl"], False
    else:
        fn = f"{ra.D}/transfer_checks_msl_{run}.npz"
        if not os.path.exists(fn):
            return None
        m = np.load(fn)
        lat, lon, msl, per = m["lat"], m["lon"], m["msl"], True
        if how == "bilin_sector":
            ii = np.where((lat >= SECTOR["lat"][0]) & (lat <= SECTOR["lat"][1]))[0]
            jj = np.where((lon >= SECTOR["lon"][0]) & (lon <= SECTOR["lon"][1]))[0]
            lat, lon, msl, per = lat[ii], lon[jj], msl[:, ii][:, :, jj], False
    k = gs.ymd_of(m["dates"].astype("datetime64[D]"))
    keep = np.isin(k // 100 % 100, rg.MONTHS)
    G = gs.region_G(gs.geostrophic(msl[keep].astype(np.float64), lat, lon, per), lat, lon, regions)
    return k[keep], G


def transfer(kc, sc, Gc, kf, sf, Gf, c0, c1, f0, f1, regions):
    sel_c = np.isin(gs.winter_of(kc), range(c0, c1 + 1))
    sel_f = np.isin(gs.winter_of(kf), range(f0, f1 + 1))
    per = {}
    for i, r in enumerate(regions):
        ok_c = np.isfinite(sc[:, i]) & np.isfinite(Gc[:, i])
        if ok_c.mean() < 0.99:
            continue
        ms_c, ms_f = np.nanmean(sc[sel_c, i]), np.nanmean(sf[sel_f, i])
        mg_c, mg_f = np.nanmean(Gc[sel_c, i]), np.nanmean(Gf[sel_f, i])
        per[r] = dict(mean_s_change_pct=100 * (ms_f / ms_c - 1), mean_G_change_pct=100 * (mg_f / mg_c - 1),
                      transfer_change_pct=100 * ((ms_f / mg_f) / (ms_c / mg_c) - 1))
    return per


def field(per, regions=None):
    V = [r for r in per if regions is None or r in regions]
    return dict(n_regions=len(V), **{k: float(np.mean([per[r][k] for r in V]))
                                     for k in ("transfer_change_pct", "mean_G_change_pct", "mean_s_change_pct")})


def ours(summary):
    out = {}
    for name, mod in gs.OURS.items():
        (cn, c0, c1), (fn, f0, f1) = ev.COMPARISONS[name]
        res = {}
        for var in ("spd100", "spd10", "spd100_native", "spd10_cmip"):
            S = {run: speeds(run, var) for run in (cn, fn)}
            if any(v is None for v in S.values()):
                continue
            for how, regions in (("box", EU), ("bilin_sector", EU), ("bilin_global", CR)):
                GG = {run: G_of(run, how, regions) for run in (cn, fn)}
                if any(v is None for v in GG.values()):
                    continue
                A = {}
                for run in (cn, fn):
                    (ks, s), (kg, G) = S[run], GG[run]
                    common, i_s, i_g = np.intersect1d(ks, kg, return_indices=True)
                    keep = np.isin(common // 100 % 100, rg.MONTHS)
                    common, i_s, i_g = common[keep], i_s[keep], i_g[keep]
                    A[run] = (common, s[i_s][:, [rg.REGIONS.index(r) for r in regions]], G[i_g])
                per = transfer(*A[cn], *A[fn], c0, c1, f0, f1, regions)
                key = f"{var}|{how}"
                res[key] = dict(eu=field(per, EU), eu_flat=field(per, EU_FLAT), regions=per,
                                days_ctl=int(np.isin(gs.winter_of(A[cn][0]), range(c0, c1 + 1)).sum()),
                                days_fut=int(np.isin(gs.winter_of(A[fn][0]), range(f0, f1 + 1)).sum()))
                if how == "bilin_global":
                    res[key]["cmip22"] = field(per, CR)
                print(f"  {name:10s} {key:26s} EU dT {res[key]['eu']['transfer_change_pct']:+5.2f}% dG {res[key]['eu']['mean_G_change_pct']:+5.2f}%"
                      + (f"   22 reg ({res[key]['cmip22']['n_regions']}) dT {res[key]['cmip22']['transfer_change_pct']:+5.2f}% "
                         f"dG {res[key]['cmip22']['mean_G_change_pct']:+5.2f}%" if "cmip22" in res[key] else ""), flush=True)
        res["dT"] = float(summary["ours"][name]["dT"])
        out[name] = res
    return out


def ens(vals):
    x = np.array(vals, float)
    return dict(median=float(np.nanmedian(x)), p10_p90=[float(v) for v in np.nanpercentile(x, [10, 90])], n=int(np.isfinite(x).sum()))


def cmip6(summary):
    """Field means over EU and over all CMIP_REGIONS from the per-region values (A1)."""
    out = dict(eu_in_cmip_regions=all(r in summary["cmip6_regions"] for r in EU), models={})
    for grp in ("cmip6", "cmip6_windonly"):
        for m, v in summary[grp].items():
            per = v["regions"]
            out["models"][m] = dict(group=grp, dT=v["dT"], eu=field(per, EU), eu_flat=field(per, EU_FLAT), cmip22=field(per, CR),
                                    stored_field_transfer_change_pct=v["field"]["transfer_change_pct"])
    M = [m for m, v in out["models"].items() if v["group"] == "cmip6"]
    for reg in ("eu", "eu_flat", "cmip22"):
        for k in ("transfer_change_pct", "mean_G_change_pct", "mean_s_change_pct"):
            out[f"ensemble_{reg}_{k}"] = ens([out["models"][m][reg][k] for m in M])
    return out


CHECKS = {   # table row -> (key in ours[...], region set, CMIP6 footing)
    "current (spd100, box G, EU)":           ("spd100|box", "eu", "eu"),
    "A2a bilinear G, sector, EU":            ("spd100|bilin_sector", "eu", "eu"),
    "A2b bilinear G, global grid, EU":       ("spd100|bilin_global", "eu", "eu"),
    "A3 bilinear G, 22 regions":             ("spd100|bilin_global", "cmip22", "cmip22"),
    "A2c box G, EU without IB/IT/BK":         ("spd100|box", "eu_flat", "eu_flat"),
    "A2d bilinear G, EU without IB/IT/BK":    ("spd100|bilin_sector", "eu_flat", "eu_flat"),
    "A4 spd10, box G, EU":                   ("spd10|box", "eu", "eu"),
    "A4 spd10, bilinear G, 22 regions":      ("spd10|bilin_global", "cmip22", "cmip22"),
    "A4b spd10 CMIP6-style daily, bilinear G, 22 regions": ("spd10_cmip|bilin_global", "cmip22", "cmip22"),
    "A5 native 100 m, box G, EU":            ("spd100_native|box", "eu", "eu"),
    "A5 native 100 m, bilinear G, 22 regions": ("spd100_native|bilin_global", "cmip22", "cmip22"),
}


def table(res):
    """Each check x comparison, with the share of the 20 CMIP6 models whose transfer change is at least ours."""
    M = [m for m, v in res["cmip6"]["models"].items() if v["group"] == "cmip6"]
    out = {}
    for row, (key, reg, creg) in CHECKS.items():
        x = np.array([res["cmip6"]["models"][m][creg]["transfer_change_pct"] for m in M])
        xk = np.array([res["cmip6"]["models"][m][creg]["transfer_change_pct"] / res["cmip6"]["models"][m]["dT"] for m in M])
        r = dict(cmip6_transfer=res["cmip6"][f"ensemble_{creg}_transfer_change_pct"],
                 cmip6_G=res["cmip6"][f"ensemble_{creg}_mean_G_change_pct"],
                 cmip6_transfer_per_K=ens(xk),
                 ec_earth3=res["cmip6"]["models"]["EC-Earth3"][creg])
        for name, v in res["ours"].items():
            if key in v:
                f = v[key][reg]
                r[name] = dict(transfer_change_pct=f["transfer_change_pct"], mean_G_change_pct=f["mean_G_change_pct"],
                               mean_s_change_pct=f["mean_s_change_pct"], n_regions=f["n_regions"],
                               cmip6_share_at_or_above=float(np.mean(x >= f["transfer_change_pct"])),
                               transfer_per_K=f["transfer_change_pct"] / v["dT"],
                               cmip6_share_at_or_above_per_K=float(np.mean(xk >= f["transfer_change_pct"] / v["dT"])))
        out[row] = r
    return out


def main():
    summary = json.load(open(f"{gs.OUT}/geostrophic_summary.json"))
    res = dict(cmip6=cmip6(summary), ours=ours(summary))
    res["table"] = table(res)
    res["current_ours"] = {n: dict(transfer_change_pct=v["field"]["transfer_change_pct"],
                                   mean_G_change_pct=v["field"]["mean_G_change_pct"]) for n, v in summary["ours"].items()}
    c = res["cmip6"]
    print("\nCMIP6 (20 models): EU transfer {0[median]:+.2f} [{0[p10_p90][0]:+.2f}, {0[p10_p90][1]:+.2f}], dG {1[median]:+.2f} "
          "[{1[p10_p90][0]:+.2f}, {1[p10_p90][1]:+.2f}]; 22 reg transfer {2[median]:+.2f} [{2[p10_p90][0]:+.2f}, {2[p10_p90][1]:+.2f}], "
          "dG {3[median]:+.2f} [{3[p10_p90][0]:+.2f}, {3[p10_p90][1]:+.2f}]".format(
              c["ensemble_eu_transfer_change_pct"], c["ensemble_eu_mean_G_change_pct"],
              c["ensemble_cmip22_transfer_change_pct"], c["ensemble_cmip22_mean_G_change_pct"]))
    e = c["models"]["EC-Earth3"]
    print(f"EC-Earth3: EU transfer {e['eu']['transfer_change_pct']:+.2f} dG {e['eu']['mean_G_change_pct']:+.2f}; "
          f"22 reg transfer {e['cmip22']['transfer_change_pct']:+.2f} dG {e['cmip22']['mean_G_change_pct']:+.2f}")
    json.dump(res, open(OUT, "w"), indent=1)
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
