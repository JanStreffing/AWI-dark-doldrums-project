"""
Regime decomposition of the change in low-wind events for the CMIP6 models (Jan 2026-11-10, option C): does the
CMIP6 increase also come from within the regimes, or partly from more Scandinavian blocking?

Same decomposition as regime_analysis.py / regime_decomposition_robustness.py for our runs, on the CMIP6 daily footing
of cmip6/daily_events.py:
  events   onshore wind CF (log-law 100 m daily-mean wind), 48 h = 2 daily steps, Mockert window expansion, Sep-Mar
           winters with at least 200 days, thresholds from results/cmip6/daily_events_summary.json (matched so the
           historical 1951-1980 control has 5 events per winter); historical 1951-1980 against ssp585 2071-2100
  regimes  daily psl (data/cmip6/psl_<model>.npz) classified exactly as geostrophic.regime_labels (anomalies from the
           historical seasonal cycle for both periods, ERA5 EOFs, nearest ERA5 centroid). Sensitivity: the future
           classified against its own seasonal cycle.
  event regime  the most frequent regime over the event's days (days without a label are skipped)
  terms    per region e = sum_k f_k lam_k; de = sum lam_k^ctl df_k (frequency) + sum f_k^fut dlam_k (within);
           summed over the eight European regions and divided by dT (events per winter per K), as fig:mechanism c.
Checks: the events per winter reproduce the stored control and future rates of daily_events_summary.json.
Winter bootstrap (independent in control and future, the same winters in every region) per model.

Environment: levante, module load python3/2025.01-gcc-13.3.0
Usage: python3 regime_decomposition_cmip6.py
Output: dunkelflaute/results/regimes/regime_decomposition_cmip6.json
"""
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regimes")
import geostrophic as gs                      # noqa: E402  (sets region_events.WIN = 2, which is the daily footing here)
import regime_analysis as ra                  # noqa: E402
import regions as rg                          # noqa: E402
from robustness_mockert_events import mockert_events, trailing_mean   # noqa: E402

DC = f"{_ROOT}/dunkelflaute/data/cmip6"
OUT = f"{_ROOT}/dunkelflaute/results/regimes/regime_decomposition_cmip6.json"
EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
WIN, MIN_DAYS, K, NB = 2, 200, 4, 2000
C0, C1, F0, F1 = 1951, 1980, 2071, 2100


def winters_of(ymd):
    return np.where(ymd[:, 1] >= 9, ymd[:, 0] + 1, ymd[:, 0])


def key(ymd):
    return ymd @ np.array([10000, 100, 1])


def labels(model, pat):
    Pz = np.load(f"{DC}/psl_{model}.npz")
    lat, lon = Pz["lat"], Pz["lon"]
    ii = np.where(lat >= 20.0)[0]
    jj = np.where((lon >= -80.0) & (lon <= 40.0))[0]
    D = {}
    for tag, exp, y0, y1 in (("ctl", "historical", C0, C1), ("fut", "ssp585", F0, F1)):
        k = key(Pz[f"{exp}|ymd"])
        w = gs.winter_of(k)
        m = (w >= y0) & (w <= y1)
        D[tag] = (k[m], gs.dates_of(k[m]), Pz[f"{exp}|psl"][m][:, ii][:, :, jj].astype(np.float64))
    lc, lf = gs.regime_labels(pat, D["ctl"][1], D["ctl"][2], D["fut"][1], D["fut"][2])
    a = ra.anomalies(D["fut"][1], D["fut"][2], ra.seasonal_fit(D["fut"][1], D["fut"][2])) * pat["wts"].reshape(1, -1)
    lo = ra.assign(a, pat["eofs"], pat["cent"])
    return dict(ctl=dict(zip(D["ctl"][0].tolist(), lc.tolist())),
                fut=dict(zip(D["fut"][0].tolist(), lf.tolist())),
                fut_own=dict(zip(D["fut"][0].tolist(), lo.tolist())))


def per_winter(model, exp, y0, y1, thr, lab):
    """Per winter: events per regime (R, K) and regime days (K,), for the European regions."""
    z = np.load(f"{DC}/daily_cmip6_{model}.npz")
    ymd, wind = z[f"{exp}|ymd"], z[f"{exp}|wind"]
    win = winters_of(ymd)
    ri = [list(z["regions"]).index(r) for r in EU]
    E, F, totals = [], [], []
    for y in range(y0, y1 + 1):
        m = win == y
        if m.sum() < MIN_DAYS:
            continue
        ks = key(ymd[m])
        lk = np.array([lab.get(k, -1) for k in ks.tolist()])
        F.append(np.array([np.sum(lk == j) for j in range(K)], float))
        er, tr = np.zeros((len(EU), K)), np.zeros(len(EU))
        for n, i in enumerate(ri):
            _, evs = mockert_events(trailing_mean(wind[m][:, i], WIN), WIN, thr[EU[n]])
            tr[n] = len(evs)
            for a, b in evs:
                ll = lk[a:b][lk[a:b] >= 0]
                if len(ll):
                    er[n, np.bincount(ll, minlength=K).argmax()] += 1
        E.append(er)
        totals.append(tr)
    return np.array(E), np.array(F), np.array(totals)          # (W, R, K), (W, K), (W, R)


def terms(Ec, Fc, Ef, Ff):
    """Frequency and within-regime terms summed over regions, events per winter."""
    ec, ef = Ec.mean(0), Ef.mean(0)                    # (R, K)
    fc, ff = Fc.mean(0), Ff.mean(0)                    # (K,)
    lc = np.divide(ec, fc, out=np.zeros_like(ec), where=fc > 0)
    lf = np.divide(ef, ff, out=np.zeros_like(ef), where=ff > 0)
    freq_k = (lc * (ff - fc)).sum(0)
    within_k = (ff * (lf - lc)).sum(0)
    return freq_k, within_k, float(ef.sum() - ec.sum())


def main():
    pat = np.load(f"{ra.OUT}/regime_patterns.npz")
    names = [str(x) for x in pat["names"]]
    S = json.load(open(f"{_ROOT}/dunkelflaute/results/cmip6/daily_events_summary.json"))["regions"]
    rng = np.random.default_rng(5)
    out = dict(names=names, regions=EU, n_boot=NB, models={})
    for comp, v in sorted(S.items()):
        model, span = comp.split("|")
        if span != "2071-2100" or not os.path.exists(f"{DC}/psl_{model}.npz"):
            continue
        thr = {r: v["metrics"]["wind"][r]["threshold"] for r in EU}
        lab = labels(model, pat)
        Ec, Fc, Tc = per_winter(model, "historical", C0, C1, thr, lab["ctl"])
        Ef, Ff, Tf = per_winter(model, "ssp585", F0, F1, thr, lab["fut"])
        Eo, Fo, _ = per_winter(model, "ssp585", F0, F1, thr, lab["fut_own"])
        stored = np.array([[v["metrics"]["wind"][r]["control"], v["metrics"]["wind"][r]["future"]] for r in EU])
        repro = float(np.abs(np.stack([Tc.mean(0), Tf.mean(0)], 1) - stored).max())
        dT = v["dT"]
        fk, wk, tot = terms(Ec, Fc, Ef, Ff)
        fo, wo, _ = terms(Ec, Fc, Eo, Fo)
        bs = []
        for _ in range(NB):
            ic, i_f = rng.integers(0, len(Ec), len(Ec)), rng.integers(0, len(Ef), len(Ef))
            a, b, t = terms(Ec[ic], Fc[ic], Ef[i_f], Ff[i_f])
            bs.append([t, a.sum(), b.sum()])
        bs = np.array(bs) / dT
        lost = 1 - (Ec.sum() + Ef.sum()) / (Tc.sum() + Tf.sum())
        out["models"][model] = dict(
            group=v["group"], dT=dT, n_winters=[len(Ec), len(Ef)], reproduce_max_abs_rate=repro,
            events_without_regime_share=float(lost),
            regime_days_per_winter=dict(ctl=Fc.mean(0).tolist(), fut=Ff.mean(0).tolist(), fut_own=Fo.mean(0).tolist()),
            event_share_ctl=(Ec.sum((0, 1)) / Ec.sum()).tolist(),
            events_per_100_regime_days_ctl=(100 * Ec.mean(0).mean(0) / np.maximum(Fc.mean(0), 1e-9)).tolist(),
            total_perK=tot / dT, frequency_perK=float(fk.sum() / dT), within_perK=float(wk.sum() / dT),
            frequency_k_perK=(fk / dT).tolist(), within_k_perK=(wk / dT).tolist(),
            ownclim=dict(frequency_perK=float(fo.sum() / dT), within_perK=float(wo.sum() / dT)),
            boot_p05_p95=dict(total=np.percentile(bs[:, 0], [5, 95]).tolist(),
                              frequency=np.percentile(bs[:, 1], [5, 95]).tolist(),
                              within=np.percentile(bs[:, 2], [5, 95]).tolist()))
        m = out["models"][model]
        print(f"{model:18s} {m['group']:15s} repro {repro:.3f}  total {m['total_perK']:+.3f} = freq {m['frequency_perK']:+.3f}"
              f" [{bs[:, 1].min():+.2f}..] + within {m['within_perK']:+.3f}; own-cycle freq {m['ownclim']['frequency_perK']:+.3f};"
              f" ScBL days {Fc.mean(0)[2]:.1f}->{Ff.mean(0)[2]:.1f}; share ScBL {m['event_share_ctl'][2]:.2f}", flush=True)
    E = {k: v for k, v in out["models"].items() if v["group"] == "cmip6"}
    arr = lambda f: np.array([f(v) for v in E.values()])           # noqa: E731
    tot, fr, wi = arr(lambda v: v["total_perK"]), arr(lambda v: v["frequency_perK"]), arr(lambda v: v["within_perK"])
    fo = arr(lambda v: v["ownclim"]["frequency_perK"])
    frac = fr / tot
    out["ensemble"] = dict(
        n_models=len(E), total_median=float(np.median(tot)), frequency_median=float(np.median(fr)),
        within_median=float(np.median(wi)), frequency_share_of_total_median=float(np.median(frac)),
        frequency_share_p10_p90=np.percentile(frac, [10, 90]).tolist(),
        n_within_positive=int((wi > 0).sum()), n_within_gt_frequency=int((wi > fr).sum()),
        n_frequency_positive=int((fr > 0).sum()),
        n_frequency_boot_above_zero=int(sum(v["boot_p05_p95"]["frequency"][0] > 0 for v in E.values())),
        n_within_boot_above_zero=int(sum(v["boot_p05_p95"]["within"][0] > 0 for v in E.values())),
        ownclim_frequency_median=float(np.median(fo)),
        scbl_event_share_ctl_median=float(np.median(arr(lambda v: v["event_share_ctl"][2]))),
        scbl_event_share_ctl_range=[float(arr(lambda v: v["event_share_ctl"][2]).min()),
                                    float(arr(lambda v: v["event_share_ctl"][2]).max())])
    json.dump(out, open(OUT, "w"), indent=1)
    print(json.dumps(out["ensemble"], indent=1))
    print(OUT)


if __name__ == "__main__":
    main()
