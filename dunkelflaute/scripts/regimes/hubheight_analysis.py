"""
Hub-height test, analysis of hubheight_extract.py output. Per CMIP6 model and region (8 European,
22 with G), winter (Sep-Mar) period means, historical 1951-1980 against ssp585 2071-2100, changes in %:

  all cell-days     sfc (day/sfcWind), s10v (|daily-mean uas, vas|), G (geostrophic.py, day/psl),
                    V850 (850 hPa >= HMIN above ground), V100ml (CFday model levels at 100 m, 5 models)
  m925 cell-days    V925 and sfc where 925 hPa is >= HMIN above ground (Eday, 10 models)
  m1000 cell-days   V1000, sfc, s10v, V925 where 1000 hPa is >= HMIN above ground; G weighted by the
                    daily valid fraction of the region (so G is sampled on the same days)
  ratios            s10/G, s10/V1000, V1000/G (m1000), s10/V925, V925/G (m925), s10/V100ml, V100ml/G
  V100pl            power law between 10 m and 1000 hPa on m1000 cell-days, period means:
                    sfc (V1000/sfc)^(ln 10 / ln(h/10)); removes the height change of 1000 hPa
A region enters the m1000 quantities of a model if its mean valid fraction is >= FMIN in both periods.
Field means: mean of the per-region % changes over the regions valid for that model, and over the
regions valid in every model ("common").

Correlations across models with the land-cover changes of stability_cmip6.json; crude event-statistic
scaling with the fit per K = 3.0 - 2.07 x (mean CF change %/K) of methods notes section 14.
TCO95 from transfer_stability.json (B4 chain: native 100 m, 10 m, level 85, G).

Usage: python3 hubheight_analysis.py      (legacy python, login node)
Output: dunkelflaute/results/regimes/hubheight_summary.json
"""
import glob
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("regions", "regimes", "analysis", "cmip6"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/{p}")
import regions as rg                     # noqa: E402
import geostrophic as gs                 # noqa: E402

PARTS = f"{_ROOT}/dunkelflaute/data/regimes/hubheight_parts"
DC = f"{_ROOT}/dunkelflaute/data/cmip6"
RES = f"{_ROOT}/dunkelflaute/results"
OUT = f"{RES}/regimes/hubheight_summary.json"
FMIN = float(os.environ.get("FMIN", 0.05))
EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
C22 = gs.CMIP_REGIONS
SETS = {"eu8": EU, "cmip22": C22}
EXP = ("historical", "ssp585")
# ssp585 day/ua of KACE-1-0-G is corrupt: the levels labelled 1000 and 850 hPa have mean |u| of 14-18 m/s
# (max 108 m/s) against 5.8 / 6.3 m/s in the historical file; its pressure-level quantities are dropped
BAD_PL = {"KACE-1-0-G"}


def wmean(x, f):
    ok = np.isfinite(x) & np.isfinite(f) & (f > 0)
    return float(np.sum(x[ok] * f[ok]) / np.sum(f[ok])) if ok.any() else np.nan


def pct(a, b):
    return 100.0 * (b / a - 1.0)


def model(name):
    P = np.load(f"{PARTS}/hubheight_pl_{name}.npz")
    mlf = f"{PARTS}/hubheight_ml_{name}.npz"
    M = np.load(mlf) if os.path.exists(mlf) else None
    S = np.load(f"{DC}/psl_{name}.npz")
    regs = list(P["regions"])
    per = {e: {} for e in EXP}
    for e in EXP:
        kp = P[f"{e}|ymd"] @ np.array([10000, 100, 1])
        ks = S[f"{e}|ymd"] @ np.array([10000, 100, 1])
        com, ip, is_ = np.intersect1d(kp, ks, return_indices=True)
        G = gs.region_G(gs.geostrophic(S[f"{e}|psl"][is_].astype(np.float64), S["lat"], S["lon"], True), S["lat"], S["lon"], C22)
        if M is not None:
            km = M[f"{e}|ymd"] @ np.array([10000, 100, 1])
            _, ip2, im = np.intersect1d(com, km, return_indices=True)
            assert len(ip2) > 0.98 * len(com), f"{name}: ml days"
        for r in C22:
            i, g = regs.index(r), C22.index(r)
            x = lambda k: P[f"{e}|{k}"][ip, i].astype(np.float64)   # noqa: E731
            Gr = G[:, g]
            f1, f9, f8 = x("frac_m1000"), x("frac_m925"), x("frac850")
            d = dict(sfc=np.nanmean(x("sfc_all")), s10v=np.nanmean(x("s10v_all")), G=np.nanmean(Gr),
                     V850=wmean(x("V850_all"), f8), frac1000=np.nanmean(f1), frac925=np.nanmean(f9),
                     h1000_all=np.nanmean(x("h1000_all")), zgh1000_all=np.nanmean(x("zgh1000_all")), ps=np.nanmean(x("ps_all")),
                     V1000=wmean(x("V1000_m1000"), f1), sfc1=wmean(x("sfc_m1000"), f1), s10v1=wmean(x("s10v_m1000"), f1),
                     V925_1=wmean(x("V925_m1000"), f1), G1=wmean(Gr, f1), h1=wmean(x("h1000_m1000"), f1),
                     zgh1=wmean(x("zgh1000_m1000"), f1), V925=wmean(x("V925_m925"), f9), sfc9=wmean(x("sfc_m925"), f9),
                     G9=wmean(Gr, f9), h925=wmean(x("h925_m925"), f9))
            d["V100pl"] = d["sfc1"] * (d["V1000"] / d["sfc1"]) ** (np.log(10.0) / np.log(d["h1"] / 10.0)) if d["h1"] > 20 else np.nan
            if M is not None:
                v = M[f"{e}|V100ml_all"][im, i].astype(np.float64)
                sf = P[f"{e}|sfc_all"][ip[ip2], i].astype(np.float64)
                d["V100ml"], d["sfc_ml"], d["G_ml"] = np.nanmean(v), np.nanmean(sf), np.nanmean(Gr[ip2])
                d["zlow"], d["z2"] = np.nanmean(M[f"{e}|zlow_all"][im, i]), np.nanmean(M[f"{e}|z2_all"][im, i])
                d["Vlow"] = np.nanmean(M[f"{e}|Vlow_all"][im, i])
            per[e][r] = d
    out = {}
    for r in C22:
        c, f = per["historical"][r], per["ssp585"][r]
        ok1 = c["frac1000"] >= FMIN and f["frac1000"] >= FMIN
        ok9 = c["frac925"] >= 0.5 and f["frac925"] >= 0.5
        o = dict(ctl={k: float(v) for k, v in c.items()}, fut={k: float(v) for k, v in f.items()}, ok1000=bool(ok1), ok925=bool(ok9))
        o["d_sfc"] = pct(c["sfc"], f["sfc"])
        o["d_s10v"] = pct(c["s10v"], f["s10v"])
        o["d_G"] = pct(c["G"], f["G"])
        o["d_s10/G"] = pct(c["sfc"] / c["G"], f["sfc"] / f["G"])
        o["d_V850"] = pct(c["V850"], f["V850"])
        o["d_V850/G"] = pct(c["V850"] / c["G"], f["V850"] / f["G"])
        o["d_s10/V850"] = pct(c["sfc"] / c["V850"], f["sfc"] / f["V850"])
        if ok1:
            o["d_V1000"] = pct(c["V1000"], f["V1000"])
            o["d_sfc_m1000"] = pct(c["sfc1"], f["sfc1"])
            o["d_s10/V1000"] = pct(c["sfc1"] / c["V1000"], f["sfc1"] / f["V1000"])
            o["d_s10v/V1000"] = pct(c["s10v1"] / c["V1000"], f["s10v1"] / f["V1000"])
            o["d_V1000/G"] = pct(c["V1000"] / c["G1"], f["V1000"] / f["G1"])
            o["d_s10/G_m1000"] = pct(c["sfc1"] / c["G1"], f["sfc1"] / f["G1"])
            o["d_G_m1000"] = pct(c["G1"], f["G1"])
            o["d_h1000_m"] = f["h1"] - c["h1"]
            o["d_V100pl"] = pct(c["V100pl"], f["V100pl"])
            o["d_s10/V100pl"] = pct(c["sfc1"] / c["V100pl"], f["sfc1"] / f["V100pl"])
            o["d_V100pl/G"] = pct(c["V100pl"] / c["G1"], f["V100pl"] / f["G1"])
            o["d_V925_m1000"] = pct(c["V925_1"], f["V925_1"])
            o["d_V1000/V925_m1000"] = pct(c["V1000"] / c["V925_1"], f["V1000"] / f["V925_1"])
        if ok9 and np.isfinite(c["V925"]):
            o["d_V925"] = pct(c["V925"], f["V925"])
            o["d_s10/V925"] = pct(c["sfc9"] / c["V925"], f["sfc9"] / f["V925"])
            o["d_V925/G"] = pct(c["V925"] / c["G9"], f["V925"] / f["G9"])
        if "V100ml" in c:
            o["d_V100ml"] = pct(c["V100ml"], f["V100ml"])
            o["d_s10/V100ml"] = pct(c["sfc_ml"] / c["V100ml"], f["sfc_ml"] / f["V100ml"])
            o["d_V100ml/G"] = pct(c["V100ml"] / c["G_ml"], f["V100ml"] / f["G_ml"])
        if name in BAD_PL:
            o = {k: v for k, v in o.items() if not any(t in k for t in ("V1000", "V850", "V925", "V100pl", "_m1000", "h1000_m"))}
        out[r] = {k: (float(v) if not isinstance(v, (dict, bool)) else v) for k, v in o.items()}
    return out, json.loads(str(P["info"])), (json.loads(str(M["info"])) if M is not None else None)


KEYS = ["d_sfc", "d_s10v", "d_G", "d_s10/G", "d_V850", "d_V850/G", "d_s10/V850",
        "d_V1000", "d_sfc_m1000", "d_s10/V1000", "d_s10v/V1000", "d_V1000/G", "d_s10/G_m1000", "d_G_m1000", "d_h1000_m",
        "d_V100pl", "d_s10/V100pl", "d_V100pl/G", "d_V925_m1000", "d_V1000/V925_m1000",
        "d_V925", "d_s10/V925", "d_V925/G", "d_V100ml", "d_s10/V100ml", "d_V100ml/G"]


def field(per_region, regions):
    o = {}
    for k in KEYS:
        v = [per_region[r][k] for r in regions if k in per_region[r] and np.isfinite(per_region[r][k])]
        o[k] = float(np.mean(v)) if len(v) == len(regions) and v else (float(np.mean(v)) if v and k in ("d_sfc", "d_G") else None)
    return o


def field_valid(per_region, regions):
    """Field mean over the regions where the quantity exists for this model (count reported)."""
    o = {}
    for k in KEYS:
        v = [per_region[r][k] for r in regions if k in per_region[r] and np.isfinite(per_region[r][k])]
        o[k] = float(np.mean(v)) if v else None
        o[f"n|{k}"] = len(v)
    return o


def corr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 4:
        return dict(n=int(ok.sum()), r=None)
    r = float(np.corrcoef(a[ok], b[ok])[0, 1])
    sl, ic = np.polyfit(b[ok], a[ok], 1)
    return dict(n=int(ok.sum()), r=r, slope=float(sl), intercept=float(ic))


def ens(vals):
    v = np.array([x for x in vals if x is not None and np.isfinite(x)])
    if len(v) == 0:
        return None
    return dict(n=int(len(v)), median=float(np.median(v)), p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)))


def main():
    inv = json.load(open(f"{DC}/inventory.json"))
    stab = json.load(open(f"{RES}/regimes/stability_cmip6.json"))
    evs = json.load(open(f"{RES}/cmip6/daily_events_summary.json"))["comparisons"]
    names = [n for n in inv if os.path.exists(f"{PARTS}/hubheight_pl_{n}.npz")]
    per, infos = {}, {}
    for n in names:
        per[n], ip, im = model(n)
        infos[n] = dict(pl=ip, ml=im)
        print(f"  loaded {n}", flush=True)
    # common region sets for the m1000 quantities
    common = {s: [r for r in R if all(per[n][r]["ok1000"] for n in names if n not in BAD_PL)] for s, R in SETS.items()}
    res = dict(meta=dict(periods="historical 1951-1980 vs ssp585 2071-2100, winters Sep-Mar", HMIN_m=infos[names[0]]["pl"]["HMIN"], FMIN=FMIN,
                         common_regions_m1000=common, notes="d_x in %; d_h1000_m in m; s10 = day/sfcWind; m1000 = cell-days with 1000 hPa >= HMIN above ground",
                         excluded_pressure_level=sorted(BAD_PL)),
               models={}, ensemble={}, correlations={}, event_scaling={})
    for n in names:
        lc = stab["models"].get(n, {})
        res["models"][n] = dict(wind_only=bool(inv[n].get("wind_only")), info=infos[n], per_region=per[n],
                                field={s: field_valid(per[n], R) for s, R in SETS.items()},
                                field_common={s: field(per[n], common[s]) for s in SETS},
                                landcover={s: {k: lc.get(s, {}).get(k) for k in ("d_cropFrac", "d_treeFrac", "d_lai_pct", "within_G", "total", "dG_pct")} for s in SETS})
    ens_names = [n for n in names if not inv[n].get("wind_only")]
    for s in SETS:
        res["ensemble"][s] = {}
        for foot in ("field", "field_common"):
            res["ensemble"][s][foot] = {k: ens([res["models"][n][foot][s].get(k) for n in ens_names]) for k in KEYS}
        # correlations with land cover and with the at-fixed-G term
        cs = {}
        for k in ("d_s10/G", "d_s10/V1000", "d_s10v/V1000", "d_V1000/G", "d_s10/V100pl", "d_V100pl/G", "d_s10/V925", "d_V925/G",
                  "d_s10/V850", "d_V850/G", "d_s10/V100ml", "d_V100ml/G", "d_sfc", "d_V1000", "d_V100ml", "d_V925", "d_s10/G_m1000"):
            for foot in ("field", "field_common"):
                y = [res["models"][n][foot][s].get(k) for n in ens_names]
                y = [np.nan if v is None else v for v in y]
                for lk in ("d_cropFrac", "d_treeFrac", "d_lai_pct", "within_G"):
                    x = [res["models"][n]["landcover"][s].get(lk) for n in ens_names]
                    x = [np.nan if v is None else v for v in x]
                    cs[f"{foot}|{k}~{lk}"] = corr(y, x)
        res["correlations"][s] = cs
    # crude event-statistic scaling (24-region fit of section 14, applied with eu8 / cmip22 shear changes)
    ev = {}
    for n in names:
        c = evs.get(f"{n}|2071-2100")
        if c is None:
            continue
        ev[n] = dict(dT=c["dT"], per_K=c["wind"]["pct_per_K"], cf_change_pct=c["wind"]["mean_change_pct"],
                     d_sfc24=float(np.nanmean([per[n][r]["d_sfc"] for r in C22])))
    xs = np.array([ev[n]["d_sfc24"] for n in ens_names if n in ev])
    ys = np.array([ev[n]["cf_change_pct"] for n in ens_names if n in ev])
    elast = float((xs * ys).sum() / (xs * xs).sum())
    sc = dict(elasticity_cf_per_speed=elast, note="crude: CF elasticity from a regression through the origin of the 22-region mean CF change on the "
              "sfcWind change across models; shift in per K = 2.07 x elasticity x d(s10/V_hub) / dT, i.e. the mean CF change at hub height is "
              "taken as elasticity x (d sfcWind - d(s10/V_hub)), and the section 14 slope -2.07 converts it to events per K", models={})
    for n in ev:
        row = dict(ev[n])
        for s in SETS:
            for k, foot in (("d_s10/V1000", "field"), ("d_s10/V100pl", "field"), ("d_s10/V100ml", "field"), ("d_s10/V925", "field"),
                            ("d_s10/V1000", "field_common"), ("d_s10/V100pl", "field_common")):
                v = res["models"][n][foot][s].get(k)
                if v is not None:
                    row[f"shift_per_K|{s}|{foot}|{k}"] = 2.07 * elast * v / ev[n]["dT"]
                    row[f"per_K_hub|{s}|{foot}|{k}"] = ev[n]["per_K"] + row[f"shift_per_K|{s}|{foot}|{k}"]
        sc["models"][n] = row
    keys = sorted({k for n in ens_names if n in sc["models"] for k in sc["models"][n] if k.startswith(("shift", "per_K"))})
    sc["ensemble"] = {k: ens([sc["models"][n].get(k, np.nan) for n in ens_names if n in sc["models"]]) for k in keys}
    res["event_scaling"] = sc
    # TCO95
    ts = json.load(open(f"{RES}/regimes/transfer_stability.json"))
    res["tco95"] = {c: dict(B4={k: ts[c]["B4"][k] for k in ("s10/G", "s100n/G", "s100n/s10", "V85/G", "V40/G", "s10/V40", "V40/V85", "s100n/V85")},
                            dG=ts[c]["B1"]["mean_G_change_pct"], ds10=ts[c]["B1"]["mean_s10_change_pct"], level_heights_m=ts[c]["level_heights_m"])
                    for c in ("100km", "hist_fut")}
    json.dump(res, open(OUT, "w"), indent=1)
    # ---------------- printout
    def f_(v, w=6, p=2):
        return f"{v:+{w}.{p}f}" if v is not None and np.isfinite(v) else " " * (w - 2) + "--"
    for s in SETS:
        print(f"\n=== {s}: common m1000 regions {common[s]}")
        cols = ["d_sfc", "d_V1000", "d_V925", "d_V100ml", "d_G", "d_s10/G", "d_s10/V1000", "d_V1000/G", "d_s10/V100pl", "d_s10/V925", "d_V925/G", "d_s10/V100ml", "d_V100ml/G"]
        print(f"{'model':18s} " + " ".join(f"{c[2:]:>11s}" for c in cols) + "   crop  withinG")
        for n in names:
            F = res["models"][n]["field"][s]
            lc = res["models"][n]["landcover"][s]
            print(f"{n:18s} " + " ".join(f"{f_(F.get(c), 8)}({F.get('n|' + c, 0):2d})" for c in cols) + f" {f_(lc.get('d_cropFrac'))} {f_(lc.get('within_G'))}")
        E = res["ensemble"][s]["field"]
        print(f"{'median':18s} " + " ".join(f"{f_(E[c]['median'] if E[c] else None, 8)}({E[c]['n'] if E[c] else 0:2d})" for c in cols))
        print(f"{'p10':18s} " + " ".join(f"{f_(E[c]['p10'] if E[c] else None, 8)}    " for c in cols))
        print(f"{'p90':18s} " + " ".join(f"{f_(E[c]['p90'] if E[c] else None, 8)}    " for c in cols))
        for k, v in res["correlations"][s].items():
            if v.get("r") is not None and k.startswith("field|") and ("crop" in k or "within" in k):
                print(f"  {k:45s} r={v['r']:+.2f} n={v['n']} slope={v['slope']:+.3f}")
    for s in SETS:
        print(f"\n=== {s}, common m1000 regions {common[s]} (all quantities on these regions)")
        cols = ["d_sfc", "d_sfc_m1000", "d_V1000", "d_V100pl", "d_V925_m1000", "d_G", "d_G_m1000", "d_s10/G", "d_s10/G_m1000",
                "d_s10/V1000", "d_s10v/V1000", "d_V1000/G", "d_s10/V100pl", "d_V100pl/G", "d_V1000/V925_m1000", "d_h1000_m"]
        print(f"{'model':18s} " + " ".join(f"{c[2:]:>10s}" for c in cols))
        for n in names:
            F = res["models"][n]["field_common"][s]
            print(f"{n:18s} " + " ".join(f"{f_(F.get(c), 10)}" for c in cols))
        E = res["ensemble"][s]["field_common"]
        for q in ("median", "p10", "p90"):
            print(f"{q:18s} " + " ".join(f"{f_(E[c][q] if E[c] else None, 10)}" for c in cols))
    print("\n=== heights over the 8 European regions, control (m): h1000 all cell-days (negative = below ground) | "
          "valid fraction | h1000 on valid cell-days (zg1000 - orog) | h925 | lowest model levels")
    for n in names:
        pr = res["models"][n]["per_region"]
        row = []
        for r in EU:
            c = pr[r]["ctl"]
            row.append(f"{r} {c['h1000_all']:5.0f}|{c['frac1000']:.2f}|{c['h1']:4.0f}({c['zgh1']:4.0f})|{c['h925']:4.0f}")
        ml = infos[n]["ml"]
        print(f"{n:18s} " + "  ".join(row) + (f"  ml {ml['historical|z_levels_median'][:2]}" if ml else ""))
    print(f"\nCF elasticity {elast:.2f}")
    for k, v in sc["ensemble"].items():
        if v:
            print(f"  {k:55s} median {v['median']:+.2f} [{v['p10']:+.2f}, {v['p90']:+.2f}] n={v['n']}")
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
