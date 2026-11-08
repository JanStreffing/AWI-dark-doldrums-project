"""
Euro-Atlantic weather regimes from daily msl, and how much of the change in Dunkelflaute events
comes from regime frequency versus a change within regimes.

Regimes (reference = ERA5 Sep-Mar 1979/80-2024/25):
  - daily msl on the 2.5 deg grid 20-80 N, 80 W-40 E (extract_msl.py)
  - anomalies from a smooth seasonal cycle (mean + 3 annual harmonics, least squares on Sep-Mar days)
  - weights sqrt(cos lat); EOFs by SVD; first NPC principal components
  - k-means, K clusters (n_init 50); clusters named from the centroid anomaly pattern
Models: anomalies from the comparison's CONTROL seasonal cycle (both periods), projected on the
ERA5 EOFs, assigned to the nearest ERA5 centroid. Sensitivity: each period's own seasonal cycle
(removes the forced mean msl change, keeps changes in variability and sequence).

Events: as region_events.py (6-hourly, 48 h, window expansion, base-rate-matched threshold from
region_events_summary.json), for the European regions. An event's regime = the most frequent regime
over its calendar days.
Decomposition per region (per winter, regime k): events e_k = f_k * lam_k with f_k regime days and
lam_k events per regime day. de = sum_k lam_k^ctl df_k (frequency) + sum_k f_k^fut dlam_k (within).

Usage: python3 regime_analysis.py
Output: dunkelflaute/results/regimes/regime_summary.json, regime_patterns.npz
"""
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/analysis")
import region_events as ev                                                    # noqa: E402
from robustness_mockert_events import mockert_events, trailing_mean          # noqa: E402

D = f"{_ROOT}/dunkelflaute/data/regimes"
OUT = f"{_ROOT}/dunkelflaute/results/regimes"
NPC = 10
K = 4
EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
NAMES = {"NAO+": (65, -20, -1, 40, -20, +1), "NAO-": (65, -20, +1, 40, -20, -1)}


def load_msl(name):
    z = np.load(f"{D}/msl_{name}.npz")
    return z["dates"].astype("datetime64[D]"), z["msl"].astype(np.float64), z["lat"], z["lon"]


def harmonics(dates):
    doy = (dates - dates.astype("datetime64[Y]")).astype(int)
    w = 2 * np.pi * doy / 365.25
    return np.column_stack([np.ones(len(dates))] + [f(k * w) for k in (1, 2, 3) for f in (np.cos, np.sin)])


def seasonal_fit(dates, x):
    A = harmonics(dates)
    coef, *_ = np.linalg.lstsq(A, x.reshape(len(dates), -1), rcond=None)
    return coef


def anomalies(dates, x, coef):
    return x.reshape(len(dates), -1) - harmonics(dates) @ coef


def winter_of_day(dates):
    m = (dates.astype("datetime64[M]").astype(int) % 12) + 1
    y = dates.astype("datetime64[Y]").astype(int) + 1970
    return np.where(m >= 9, y + 1, y)


def name_clusters(cent_maps, lat, lon):
    """Name centroids from where their anomaly peaks: Iceland/Azores dipole for NAO+/-, positive
    centre over Scandinavia for blocking, over the mid-Atlantic for the ridge."""
    LA, LO = np.meshgrid(lat, lon, indexing="ij")

    def at(m, la, lo, r=7.5):
        sel = (np.abs(LA - la) <= r) & (np.abs(LO - lo) <= 2 * r)
        return m[sel].mean()

    scores = {}
    for k, m in enumerate(cent_maps):
        scores[k] = dict(nao=at(m, 40, -25) - at(m, 65, -25), scand=at(m, 62, 15), ridge=at(m, 52, -30))
    names = {}
    left = set(scores)
    for nm, key, sign in (("NAO+", "nao", 1), ("NAO-", "nao", -1), ("ScBL", "scand", 1), ("AR", "ridge", 1)):
        k = max(left, key=lambda j: sign * scores[j][key])
        names[k] = nm
        left.discard(k)
    return [names[k] for k in range(len(cent_maps))], scores


def kmeans(X, k, n_init=50, seed=0, iters=300):
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(n_init):
        c = X[rng.choice(len(X), k, replace=False)]
        for _ in range(iters):
            lab = np.argmin(((X[:, None, :] - c[None]) ** 2).sum(-1), axis=1)
            c2 = np.stack([X[lab == j].mean(0) if (lab == j).any() else c[j] for j in range(k)])
            if np.allclose(c2, c):
                break
            c = c2
        inertia = ((X - c[lab]) ** 2).sum()
        if best is None or inertia < best[0]:
            best = (inertia, c, lab)
    return best[1], best[2]


def assign(anom_w, eofs, cent):
    pcs = anom_w @ eofs.T
    return np.argmin(((pcs[:, None, :] - cent[None]) ** 2).sum(-1), axis=1)


def events_with_days(d, kind, region, y0, y1, thr):
    """Per winter: list of events as arrays of calendar days."""
    x = ev.metric(d, kind)[:, region]
    out = {}
    for y, m in ev.winters(d, y0, y1):
        t = d["time"][m]
        _, evs = mockert_events(trailing_mean(x[m], ev.WIN), ev.WIN, thr)
        out[y] = [np.unique((t[a:b] - np.timedelta64(1, "s")).astype("datetime64[D]")) for a, b in evs]
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    dates, x, lat, lon = load_msl("era5_-")
    win = winter_of_day(dates)
    keep = (win >= 1980) & (win <= 2025)
    dates, x = dates[keep], x[keep]
    coef = seasonal_fit(dates, x)
    wts = np.sqrt(np.cos(np.deg2rad(lat)))[:, None] * np.ones((1, len(lon)))
    a = anomalies(dates, x, coef) * wts.ravel()
    U, s, Vt = np.linalg.svd(a - a.mean(0), full_matrices=False)
    eofs = Vt[:NPC]
    expl = (s ** 2 / (s ** 2).sum())[:NPC]
    pcs = a @ eofs.T
    cent, lab = kmeans(pcs, K)
    cmaps = [(a[lab == k].mean(0) / wts.ravel()).reshape(len(lat), len(lon)) / 100.0 for k in range(K)]   # hPa
    names, scores = name_clusters(cmaps, lat, lon)
    order = [names.index(n) for n in ("NAO+", "NAO-", "ScBL", "AR")]
    cent = cent[order]
    names = [names[k] for k in order]
    cmaps = [cmaps[k] for k in order]
    lab = np.argsort(order)[lab]
    out = dict(npc=NPC, k=K, explained_variance=expl.tolist(), names=names,
               era5=dict(winters=[1980, 2025], freq={n: float((lab == k).mean()) for k, n in enumerate(names)}),
               models={}, events={})
    np.savez(f"{OUT}/regime_patterns.npz", lat=lat, lon=lon, centroid_hPa=np.stack(cmaps), names=np.array(names),
             eofs=eofs, cent=cent, coef=coef, wts=wts)
    print(f"ERA5: {len(dates)} days, EOF 1-{NPC} explain {expl.sum():.0%}; regimes "
          + ", ".join(f"{n} {(lab == k).mean():.1%}" for k, n in enumerate(names)), flush=True)

    S = json.load(open(f"{_ROOT}/dunkelflaute/results/regions/region_events_summary.json"))
    R = S["regions"]
    for name, ((cn, c0, c1), (fn, f0, f1)) in ev.COMPARISONS.items():
        if name == "eerie_drift" or not all(os.path.exists(f"{D}/msl_{n}.npz") for n in (cn, fn)):
            print(f"skip {name}")
            continue
        dc, xc, *_ = load_msl(cn)
        df_, xf, *_ = load_msl(fn)
        wc, wf = winter_of_day(dc), winter_of_day(df_)
        mc, mf = (wc >= c0) & (wc <= c1), (wf >= f0) & (wf <= f1)
        dc, xc, wc = dc[mc], xc[mc], wc[mc]
        df_, xf, wf = df_[mf], xf[mf], wf[mf]
        cc = seasonal_fit(dc, xc)
        regs = {}
        for tag, cf in (("ctlclim", cc), ("ownclim", seasonal_fit(df_, xf))):
            lc = assign(anomalies(dc, xc, cc) * wts.ravel(), eofs, cent)
            lf = assign(anomalies(df_, xf, cf) * wts.ravel(), eofs, cent)
            regs[tag] = (lc, lf)
        lc, lf = regs["ctlclim"]
        freq = {}
        for tag, (a_, b_) in regs.items():
            freq[tag] = {n: dict(control=float((a_ == k).mean()), future=float((b_ == k).mean())) for k, n in enumerate(names)}
        # winter bootstrap of the frequency change (ctl clim)
        rng = np.random.default_rng(3)
        ua, ub = np.unique(wc), np.unique(wf)
        fa = np.array([[np.mean(lc[wc == y] == k) for k in range(K)] for y in ua])
        fb = np.array([[np.mean(lf[wf == y] == k) for k in range(K)] for y in ub])
        bs = np.array([fb[rng.integers(0, len(fb), len(fb))].mean(0) - fa[rng.integers(0, len(fa), len(fa))].mean(0) for _ in range(5000)])
        for k, n in enumerate(names):
            freq["ctlclim"][n]["change_p05_p95"] = [float(v) for v in np.percentile(bs[:, k], [5, 95])]
        out["models"][name] = dict(freq=freq, msl_bias_vs_era5_hPa=float(np.sqrt(np.mean(((harmonics(dc[:1]) @ cc) - (harmonics(dc[:1]) @ coef)) ** 2)) / 100))

        # events by regime, European regions, base-rate-matched thresholds
        regday_c = dict(zip(dc, lc))
        regday_f = dict(zip(df_, lf))
        sc, sf = ev.load(cn), ev.load(fn)
        rows = {}
        for kind in ("comb", "wind"):
            for r in EU:
                if r not in S["valid_regions"][kind]:
                    continue
                i = R.index(r)
                thr = S["comparisons"][name][kind][r]["threshold"]
                part = {}
                for per, src, y0, y1, rd, labs, wins in (("ctl", sc, c0, c1, regday_c, lc, wc), ("fut", sf, f0, f1, regday_f, lf, wf)):
                    E = events_with_days(src, kind, i, y0, y1, thr)
                    nw = len(E)
                    ek = np.zeros(K)
                    for y, evs in E.items():
                        for days in evs:
                            rr = [rd[dd] for dd in days if dd in rd]
                            if rr:
                                ek[np.bincount(rr, minlength=K).argmax()] += 1
                    used = np.isin(wins, list(E))
                    fk = np.array([np.sum(labs[used] == k) for k in range(K)]) / max(nw, 1)     # regime days per winter
                    part[per] = dict(events_per_winter=(ek / nw).tolist(), regime_days_per_winter=fk.tolist(), n_winters=nw)
                ec, ef = np.array(part["ctl"]["events_per_winter"]), np.array(part["fut"]["events_per_winter"])
                fc, ff = np.array(part["ctl"]["regime_days_per_winter"]), np.array(part["fut"]["regime_days_per_winter"])
                lamc = np.divide(ec, fc, out=np.zeros(K), where=fc > 0)
                lamf = np.divide(ef, ff, out=np.zeros(K), where=ff > 0)
                freq_term = float((lamc * (ff - fc)).sum())
                within = float((ff * (lamf - lamc)).sum())
                rows[f"{kind}|{r}"] = dict(part, total_change=float(ef.sum() - ec.sum()), frequency_term=freq_term,
                                           within_regime_term=within, events_per_regime_day_ctl=lamc.tolist(),
                                           events_per_regime_day_fut=lamf.tolist())
        out["events"][name] = rows

    json.dump(out, open(f"{OUT}/regime_summary.json", "w"), indent=1)

    print("\nRegime frequency, control -> future (anomalies from the control climatology; [own climatology]):")
    for name, m in out["models"].items():
        f = m["freq"]
        print(f"  {name:15s} " + "  ".join(f"{n} {f['ctlclim'][n]['control']:.2f}->{f['ctlclim'][n]['future']:.2f} [{f['ownclim'][n]['future']:.2f}]" for n in names))
    print("  ERA5 1980-2025  " + "  ".join(f"{n} {out['era5']['freq'][n]:.2f}" for n in names))
    for kind in ("wind", "comb"):
        print(f"\nEvent change per winter split into regime frequency vs within-regime ({kind}, European regions):")
        for name, rows in out["events"].items():
            tot = [v["total_change"] for k, v in rows.items() if k.startswith(kind)]
            fr = [v["frequency_term"] for k, v in rows.items() if k.startswith(kind)]
            wi = [v["within_regime_term"] for k, v in rows.items() if k.startswith(kind)]
            if tot:
                print(f"  {name:15s} mean over {len(tot)} regions: total {np.mean(tot):+.2f}, frequency {np.mean(fr):+.2f}, within {np.mean(wi):+.2f}   "
                      + "  DE: " + (f"{rows[f'{kind}|DE']['total_change']:+.2f} = {rows[f'{kind}|DE']['frequency_term']:+.2f} + {rows[f'{kind}|DE']['within_regime_term']:+.2f}" if f"{kind}|DE" in rows else "-"))
    print(f"Saved {OUT}/regime_summary.json")


if __name__ == "__main__":
    main()
