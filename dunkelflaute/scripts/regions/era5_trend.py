"""
Is the forced increase in winter low-wind events already visible in ERA5? (Jan 2026-11-27)

Trend of NH-land low-wind events in ERA5, winters 1980-2026, against global-mean 2 m temperature, on the
daily footing of the CMIP6 analysis and the 5 x 5 deg land boxes of Fig. maps.

  ERA5    dunkelflaute/data/regions/cooccur_era5_wind.npz (day footing: CF of the mean speed of the
          00/06/12/18 UTC instants, log law with ERA5 z0, onshore curve), the 420 grid5 boxes.
  GMST    ERA5 monthly 2t (/pool/data/ERA5/E5/sf/an/1M/167, reduced Gaussian N320, weight
          cos(lat)/points per ring), winter value = mean of July (y-1) to June (y).
  models  AWI-CM3 100 km daily grid5 series (dunkelflaute/data/cmip6_grid5/daily_tco95_*.npz):
          hi1950g (historical, SSP5-8.5 from 2015) over the same winters, with its own global 2t;
          ctl1950g (constant 1950 forcing, 147 winters) as the null.

Estimator, identical for every series: thresholds matched per box to 5 events per winter over the
47-winter window (48 h = 2 days, Mockert window expansion); field mean of events per winter over the
boxes valid everywhere; OLS slope of ln(field mean) on the winter GMST, as 100 (exp(b) - 1) % per K.
Null: every 47-winter window of the control, regressed on the ERA5 GMST series (placebo regressor),
thresholds matched within each window. This keeps the decadal variability of the control.

Usage: python3 era5_trend.py            (legacy python with eccodes; login node, a few minutes)
Output: dunkelflaute/results/regions/era5_trend.json
"""
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/cmip6")
import region_events as ev   # noqa: E402

ev.WIN = 2
TARGET = 5.0
W0, W1 = 1980, 2026
GT_CACHE = f"{_ROOT}/dunkelflaute/data/regions/global_t2m_era5_monthly.json"
OUT = f"{_ROOT}/dunkelflaute/results/regions/era5_trend.json"


def era5_gmst_monthly():
    if os.path.exists(GT_CACHE):
        return {tuple(map(int, k.split("-"))): v for k, v in json.load(open(GT_CACHE)).items()}
    import eccodes
    out = {}
    w = None
    for y in range(W0 - 2, W1 + 1):
        fn = f"/pool/data/ERA5/E5/sf/an/1M/167/E5sf00_1M_{y}_167.grb"
        if not os.path.exists(fn):
            continue
        with open(fn, "rb") as f:
            while True:
                h = eccodes.codes_grib_new_from_file(f)
                if h is None:
                    break
                try:
                    if w is None:
                        pl = np.asarray(eccodes.codes_get_array(h, "pl"))
                        lats = np.asarray(eccodes.codes_get_array(h, "distinctLatitudes"))
                        w = np.repeat(np.cos(np.deg2rad(lats)) / pl, pl)
                    d = eccodes.codes_get(h, "dataDate")
                    v = eccodes.codes_get_values(h)
                    out[(d // 10000, d // 100 % 100)] = float((v * w).sum() / w.sum())
                finally:
                    eccodes.codes_release(h)
    json.dump({f"{a}-{b}": t for (a, b), t in sorted(out.items())}, open(GT_CACHE, "w"), indent=1)
    return out


def winter_T_era5(mon):
    T = {}
    for y in range(W0, W1 + 1):
        ks = [(y - 1, m) for m in range(7, 13)] + [(y, m) for m in range(1, 7)]
        if all(k in mon for k in ks):
            T[y] = float(np.mean([mon[k] for k in ks]))
    return T


def winter_T_model(ann):
    return {y: 0.5 * (ann[y - 1] + ann[y]) for y in ann if y - 1 in ann}


def box_counts(x, winter, y0, y1):
    """x (T,) daily CF of one box -> events per winter (W,) at the threshold matched over y0..y1."""
    ms = [winter == y for y in range(y0, y1 + 1)]
    ok = all(m.sum() >= 200 for m in ms)
    if not ok:
        return None
    grid = np.linspace(*np.nanpercentile(x[np.isin(winter, range(y0, y1 + 1))], [0.1, 40]), 300)
    c = np.stack([ev.counts_grid(x[m], grid) for m in ms])          # (W, K)
    mean = c.mean(axis=0)
    if not (mean >= 0.9 * TARGET).any():
        return None
    k = int(np.argmax(mean >= TARGET)) if (mean >= TARGET).any() else int(np.argmax(mean))
    return c[:, k].astype(float)


def field(wind, winter, cols, y0, y1):
    rows, used = [], []
    for j, c in enumerate(cols):
        r = box_counts(wind[:, c], winter, y0, y1)
        if r is not None:
            rows.append(r)
            used.append(j)
    return np.array(rows), used                                      # (B, W)


def slope(lnE, T):
    T = np.asarray(T, float)
    b = np.polyfit(T - T.mean(), lnE, 1)[0]
    return 100 * (np.exp(b) - 1)


def winters_of(ymd):
    return np.where(ymd[:, 1] >= 9, ymd[:, 0] + 1, ymd[:, 0])


def main():
    e = np.load(f"{_ROOT}/dunkelflaute/data/regions/cooccur_era5_wind.npz", allow_pickle=True)
    ereg = [str(r) for r in e["regions"]]
    day = e["day"].astype("datetime64[D]")
    ymd = np.stack([day.astype("datetime64[Y]").astype(int) + 1970,
                    day.astype("datetime64[M]").astype(int) % 12 + 1], axis=1)
    ewin = winters_of(ymd)
    ew = e["wind_day"]

    D = f"{_ROOT}/dunkelflaute/data/cmip6_grid5"
    mods = {r: np.load(f"{D}/daily_tco95_{r}.npz") for r in ("ctl1950g", "hi1950g")}
    mreg = [str(r) for r in mods["hi1950g"]["regions"]]
    boxes = [r for r in ereg if r.startswith("g") and r in mreg]
    gt = json.load(open(f"{_ROOT}/dunkelflaute/data/regions/global_t2m_tco95.json"))

    # GMST
    T_era = winter_T_era5(era5_gmst_monthly())
    years = [y for y in range(W0, W1 + 1) if y in T_era]
    print(f"ERA5 winters with GMST: {years[0]}-{years[-1]} ({len(years)}); "
          f"GMST change first to last decade {np.mean([T_era[y] for y in years[-10:]]) - np.mean([T_era[y] for y in years[:10]]):+.2f} K")

    # boxes valid in ERA5, hi and every control window would be ideal; use ERA5 and hi over the window
    y0, y1 = years[0], years[-1]
    E_rows, E_used = field(ew, ewin, [ereg.index(b) for b in boxes], y0, y1)
    hi = mods["hi1950g"]
    hwin = winters_of(hi["run|ymd"])
    H_rows, H_used = field(hi["run|wind"], hwin, [mreg.index(b) for b in boxes], y0, y1)
    common = sorted(set(E_used) & set(H_used))
    ei = [E_used.index(j) for j in common]
    hi_i = [H_used.index(j) for j in common]
    lnE_era = np.log(E_rows[ei].mean(axis=0))
    lnE_hi = np.log(H_rows[hi_i].mean(axis=0))
    T_hi = winter_T_model({int(k): v for k, v in gt["tco95_hi1950g"].items()})
    s_era = slope(lnE_era, [T_era[y] for y in years])
    s_hi = slope(lnE_hi, [T_hi[y] for y in years])
    dT_hi = np.polyfit(years, [T_hi[y] for y in years], 1)[0] * (years[-1] - years[0])
    dT_era = np.polyfit(years, [T_era[y] for y in years], 1)[0] * (years[-1] - years[0])
    # events against year as well
    tr_era = 100 * (np.exp(np.polyfit(years, lnE_era, 1)[0] * (years[-1] - years[0])) - 1)
    tr_hi = 100 * (np.exp(np.polyfit(years, lnE_hi, 1)[0] * (years[-1] - years[0])) - 1)
    print(f"boxes: {len(common)} of {len(boxes)}")
    print(f"ERA5  {s_era:+.2f} %/K; linear-trend warming {dT_era:.2f} K; events change over the record {tr_era:+.1f} %")
    print(f"100km hist {s_hi:+.2f} %/K; warming {dT_hi:.2f} K; events change {tr_hi:+.1f} %")

    # null from the control: all 47-winter windows, placebo regressor = ERA5 GMST
    ctl = mods["ctl1950g"]
    cwin = winters_of(ctl["run|ymd"])
    cy = sorted(set(cwin.tolist()))
    cy = [y for y in cy if (cwin == y).sum() >= 200]
    n = len(years)
    Tx = [T_era[y] for y in years]
    cols = [mreg.index(boxes[j]) for j in common]
    null = []
    for a in range(cy[0], cy[-1] - n + 2):
        R, used = field(ctl["run|wind"], cwin, cols, a, a + n - 1)
        if len(used) < 0.9 * len(cols):
            continue
        lnE = np.log(R.mean(axis=0))
        null.append((slope(lnE, Tx), 100 * (np.exp(np.polyfit(range(n), lnE, 1)[0] * (n - 1)) - 1)))
        print(f"  null window {a}-{a + n - 1}: {null[-1][0]:+.2f} %/K", flush=True)
    null = np.array(null)
    p_era = float((null[:, 0] >= s_era).mean())
    p_hi = float((null[:, 0] >= s_hi).mean())
    q = np.percentile(null[:, 0], [5, 50, 95])
    print(f"null %/K 5/50/95: {q[0]:+.2f} {q[1]:+.2f} {q[2]:+.2f} ({len(null)} windows, overlapping)")
    print(f"one-sided p: ERA5 {p_era:.3f}, 100km hist {p_hi:.3f}")
    res = dict(winters=[years[0], years[-1]], n_winters=n, n_boxes=len(common),
               era5=dict(per_K=s_era, warming_K=dT_era, events_change_pct=tr_era, p_one_sided=p_era,
                         lnE=lnE_era.tolist(), T=Tx),
               hi1950g=dict(per_K=s_hi, warming_K=dT_hi, events_change_pct=tr_hi, p_one_sided=p_hi,
                            lnE=lnE_hi.tolist(), T=[T_hi[y] for y in years]),
               null=dict(per_K_5_50_95=q.tolist(), n_windows=len(null), per_K=null[:, 0].tolist(),
                         change_pct=null[:, 1].tolist()))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(res, open(OUT, "w"), indent=1)
    print("Saved", OUT)


if __name__ == "__main__":
    main()
