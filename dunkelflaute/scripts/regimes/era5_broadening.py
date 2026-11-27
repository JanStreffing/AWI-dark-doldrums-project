"""
Do the day-to-day pressure systems over Europe already broaden in ERA5? (Jan 2026-11-27)

The models keep the strength of the daily highs and lows nearly constant while their slopes weaken
(Supplementary Note S1). This tests the same pair of quantities as a trend over the ERA5 winters 1980-2025.

  field     daily msl, 2.5 deg Euro-Atlantic grid (data/regimes/msl_<source>.npz, Sep-Mar days, Pa)
  e         day-to-day deviation = msl - smooth seasonal cycle of the window (ra.seasonal_fit) - that winter's
            mean deviation, so slow background changes (and their trends) are left out
  strength  winter mean of the daily area-weighted rms of e over 35-70 N, 10 W-30 E (hPa)
  slope     winter mean of the daily geostrophic speed of e, box mean over the 8 European regions (m/s)
  scale     strength / slope, a length (relative change = broadening)
Trend: OLS slope of ln(winter value) on the winter GMST, 100 (exp(b) - 1) % per K. GMST: ERA5 July-June
(data/regions/global_t2m_era5_monthly.json, from era5_trend.py); models: mean of the two calendar years.
Expected: the same quantities from the 100 km run, constant-1950 control 1953-2099 against hi1950g 2071-2099.
Null: all 46-winter windows of the 100 km control, regressed on the ERA5 GMST series.

Environment: module load python3/2025.01-gcc-13.3.0
Usage: python3 era5_broadening.py
Output: dunkelflaute/results/regimes/era5_broadening.json
"""
import json
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regimes")
import regime_G_attribution as att      # noqa: E402

ra, gs = att.ra, att.gs
D = f"{_ROOT}/dunkelflaute/data/regimes"
OUT = f"{_ROOT}/dunkelflaute/results/regimes/era5_broadening.json"
W0, W1 = 1980, 2025


def load(name):
    z = np.load(f"{D}/msl_{name}.npz", allow_pickle=True)
    d = z["dates"].astype("datetime64[D]")
    mon = d.astype("datetime64[M]").astype(int) % 12 + 1
    yr = d.astype("datetime64[Y]").astype(int) + 1970
    win = np.where(mon >= 9, yr + 1, yr)
    return dict(d=d, msl=z["msl"].astype(np.float64), win=win, lat=z["lat"], lon=z["lon"])


def per_winter(s, y0, y1):
    """(years, strength hPa, slope m/s) for the winters y0..y1 of source s."""
    m = (s["win"] >= y0) & (s["win"] <= y1)
    d, p, win = s["d"][m], s["msl"][m], s["win"][m]
    coef = ra.seasonal_fit(d, p)
    e = p - att.clim(d, coef, p.shape[1:])
    for y in np.unique(win):
        k = win == y
        e[k] -= e[k].mean(axis=0)
    lat, lon = s["lat"], s["lon"]
    ii = (lat >= 35) & (lat <= 70)
    jj = (lon >= -10) & (lon <= 30)
    w = np.cos(np.deg2rad(lat[ii]))[:, None] * np.ones((1, jj.sum()))
    sub = e[:, ii][:, :, jj]
    rms = np.sqrt((sub ** 2 * w).sum(axis=(1, 2)) / w.sum()) / 100.0
    G = att.boxG(e, lat, lon, False).mean(axis=1)
    ys = np.unique(win)
    keep = [y for y in ys if (win == y).sum() >= 200]
    return (np.array(keep), np.array([rms[win == y].mean() for y in keep]),
            np.array([G[win == y].mean() for y in keep]))


def slope(x, T):
    T = np.asarray(T, float)
    return 100 * (np.exp(np.polyfit(T - T.mean(), np.log(x), 1)[0]) - 1)


def main():
    mon = {tuple(map(int, k.split("-"))): v for k, v in json.load(open(f"{_ROOT}/dunkelflaute/data/regions/global_t2m_era5_monthly.json")).items()}
    T_era = {y: float(np.mean([mon[(y - 1, m)] for m in range(7, 13)] + [mon[(y, m)] for m in range(1, 7)]))
             for y in range(W0, W1 + 1)}
    gt = json.load(open(f"{_ROOT}/dunkelflaute/data/regions/global_t2m_tco95.json"))
    Tm = {r: {int(y): 0.5 * (v[str(int(y) - 1)] + v[y]) for y in v if str(int(y) - 1) in v}
          for r, v in ((k[6:], gt[k]) for k in ("tco95_ctl1950g", "tco95_hi1950g"))}

    res = {}
    era = load("era5_-")
    ys, st, sl = per_winter(era, W0, W1)
    T = [T_era[y] for y in ys]
    res["era5"] = dict(winters=[int(ys[0]), int(ys[-1])], n=len(ys), strength_per_K=slope(st, T),
                       slope_per_K=slope(sl, T), scale_per_K=slope(st / sl, T),
                       strength_hPa=float(st.mean()), slope_ms=float(sl.mean()))
    print("ERA5", {k: round(v, 2) if isinstance(v, float) else v for k, v in res["era5"].items()}, flush=True)

    ctl, hi = load("tco95_ctl1950g"), load("tco95_hi1950g")
    # expected: period difference, control 1953-2099 against hi 2071-2099 (the paper's 100 km comparison)
    yc, stc, slc = per_winter(ctl, 1953, 2099)
    yf, stf, slf = per_winter(hi, 2071, 2099)
    dT = np.mean([Tm["hi1950g"][y] for y in yf]) - np.mean([Tm["ctl1950g"][y] for y in yc])
    pk = lambda a, b: 100 * ((b.mean() / a.mean()) ** (1 / dT) - 1)          # noqa: E731
    res["100km_expected"] = dict(dT=float(dT), strength_per_K=pk(stc, stf), slope_per_K=pk(slc, slf),
                                 scale_per_K=pk(stc / slc, stf / slf))
    print("100 km expected", {k: round(v, 2) for k, v in res["100km_expected"].items()}, flush=True)
    yh, sth, slh = per_winter(hi, W0, W1)
    Th = [Tm["hi1950g"][y] for y in yh]
    res["100km_hist_window"] = dict(strength_per_K=slope(sth, Th), slope_per_K=slope(slh, Th), scale_per_K=slope(sth / slh, Th))
    print("100 km same winters", {k: round(v, 2) for k, v in res["100km_hist_window"].items()}, flush=True)

    # null: windows of the control with the full-control seasonal cycle per window, placebo regressor
    n = len(ys)
    null = []
    for a in range(1953, 2099 - n + 2):
        y_, s_, g_ = per_winter(ctl, a, a + n - 1)
        if len(y_) != n:
            continue
        null.append((slope(s_, T), slope(g_, T), slope(s_ / g_, T)))
    null = np.array(null)
    res["null"] = {k: dict(p5_50_95=np.percentile(null[:, i], [5, 50, 95]).tolist(), sd=float(null[:, i].std()))
                   for i, k in enumerate(("strength_per_K", "slope_per_K", "scale_per_K"))}
    res["null"]["n_windows"] = len(null)
    for i, k in enumerate(("strength_per_K", "slope_per_K", "scale_per_K")):
        res["era5"][f"p_{k}"] = float((null[:, i] <= res["era5"][k]).mean() if k == "slope_per_K"
                                      else (null[:, i] >= res["era5"][k]).mean())
    print("null", json.dumps(res["null"]), flush=True)
    print("ERA5 one-sided p (slope: weaker, strength/scale: larger)",
          {k: round(res["era5"][f"p_{k}"], 3) for k in ("strength_per_K", "slope_per_K", "scale_per_K")})
    json.dump(res, open(OUT, "w"), indent=1)
    print("Saved", OUT)


if __name__ == "__main__":
    main()
