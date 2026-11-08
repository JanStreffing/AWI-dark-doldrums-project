"""
Shared helpers for the solar-only analysis (solar_region_events.py, cmip6/solar_daily_events.py,
solar_grid_per_K.py). Nothing here changes the wind/comb machinery; it only adds:

  day_index(t)          Sep-Mar calendar-day index 0..211 of every step (non-leap reference winter;
                        the day of a step is the day of (stamp - 1 s), as winter_of in region_events.py;
                        29 Feb -> 28 Feb for climatologies, flagged so alignment can drop it)
  climatology(x, di, m) control calendar-day climatology of x (mean over the steps of each calendar day
                        in the masked control winters), smoothed by a 31-day running mean (truncated at
                        the ends of the Sep-Mar season)
  deseason(x, di, clim) x / clim[day], clim floored at CLIM_FLOOR (5 W m-2 as a daily mean)
  flagged(x, thr, win)  the Mockert window-expansion mask (same rule as region_events.counts_grid)
  counts_rows(X, thr)   events per row of a (N, L) matrix of aligned winter series, vectorised
  matched_rows(X, ...)  lowest threshold on a grid at which the mean over rows reaches TARGET
  pairings(n_a, n_b)    deterministic decoupled pairings of wind winters with solar winters

The deseasonalized metric divides the CF by the control climatology of the same calendar day, so a
48 h mean of it is the fraction of the normal irradiance for the time of year; the future is divided
by the CONTROL climatology, so a mean change in solar shows up in it.
"""
import numpy as np

TARGET = 5.0
CLIM_HALF = 15                 # 31-day running mean
CLIM_FLOOR = 0.005             # CF, i.e. 5 W m-2 as a daily mean
NDAY = 212                     # 1 Sep .. 31 Mar of a non-leap winter
_REF = np.arange(np.datetime64("2001-09-01"), np.datetime64("2002-04-01"))
_MD = {(int(str(d)[5:7]), int(str(d)[8:10])): k for k, d in enumerate(_REF)}


def md_index(month, day):
    """(month, day) arrays -> Sep-Mar day index; 29 Feb -> 28 Feb; days outside Sep-Mar -> -1.
    360-day calendars (30 Feb) map to 28 Feb as well; 31-day months are absent there."""
    out = np.full(len(month), -1, int)
    for i, (m, d) in enumerate(zip(month.tolist(), day.tolist())):
        if m == 2 and d > 28:
            d = 28
        out[i] = _MD.get((m, d), -1)
    return out


def day_index(t):
    """datetime64 stamps (end of window) -> (day index, is_feb29)."""
    t1 = (t - np.timedelta64(1, "s")).astype("datetime64[D]")
    s = t1.astype(str)
    m = np.array([int(x[5:7]) for x in s])
    d = np.array([int(x[8:10]) for x in s])
    return md_index(m, d), (m == 2) & (d == 29)


def climatology(x, di, mask):
    """x (T,) or (T, R); di (T,) day index; mask (T,) control steps -> (NDAY,) or (NDAY, R)."""
    x2 = x if x.ndim == 2 else x[:, None]
    s = np.zeros((NDAY, x2.shape[1]))
    n = np.zeros((NDAY, 1))
    ok = mask & (di >= 0)
    np.add.at(s, di[ok], np.nan_to_num(x2[ok]))
    np.add.at(n, di[ok], 1.0)
    raw = s / np.maximum(n, 1)
    raw[n[:, 0] == 0] = np.nan
    sm = np.full_like(raw, np.nan)
    for k in range(NDAY):
        lo, hi = max(0, k - CLIM_HALF), min(NDAY, k + CLIM_HALF + 1)
        sm[k] = np.nanmean(raw[lo:hi], axis=0)
    return sm if x.ndim == 2 else sm[:, 0]


def deseason(x, di, clim):
    c = np.maximum(clim, CLIM_FLOOR)
    return x / (c[di] if x.ndim == 1 else c[di, :])


def trailing(X, win):
    """Trailing mean along the last axis (NaN where the window is incomplete or contains NaN)."""
    X = np.asarray(X, float)
    c = np.concatenate([np.zeros(X.shape[:-1] + (1,)), np.cumsum(X, axis=-1)], axis=-1)
    out = np.full(X.shape, np.nan)
    out[..., win - 1:] = (c[..., win:] - c[..., :-win]) / win
    return out


def flagged(x, thr, win):
    """Window-expansion mask of one series (same rule as region_events.counts_grid)."""
    roll = trailing(x, win)
    below = np.nan_to_num(roll, nan=np.inf) < thr
    n = len(x)
    c = np.concatenate([[0], np.cumsum(below)])
    hi = np.minimum(np.arange(n) + win, n)
    return (c[hi] - c[:n]) > 0


def counts_rows(X, thr, win):
    """X (N, L) aligned winter series -> events per row at threshold thr."""
    return counts_rows_from_roll(trailing(X, win), thr, win)


def matched_rows(X, win, ref=None, ngrid=400):
    """Lowest threshold on a grid between the 0.1st and 40th percentiles (of ref, default X) at which the
    mean events per row reach TARGET; returns (threshold, rate)."""
    v = (X if ref is None else ref)
    v = v[np.isfinite(v)]
    grid = np.linspace(*np.percentile(v, [0.1, 40]), ngrid)
    roll = trailing(X, win)
    best = (grid[-1], -1.0)
    for g in grid:
        c = counts_rows_from_roll(roll, g, win).mean()
        if c > best[1]:
            best = (g, c)
        if c >= TARGET:
            return float(g), float(c)
    return float(best[0]), float(best[1])


def counts_rows_from_roll(roll, thr, win):
    below = (np.nan_to_num(roll, nan=np.inf) < thr).astype(np.int32)
    n = below.shape[1]
    c = np.concatenate([np.zeros((below.shape[0], 1), np.int32), np.cumsum(below, axis=1)], axis=1)
    hi = np.minimum(np.arange(n) + win, n)
    f = (c[:, hi] - c[:, :n]) > 0
    return f[:, 0].astype(int) + (f[:, 1:] & ~f[:, :-1]).sum(axis=1)


def pairings(n_a, n_b, same, ns=16):
    """Pairs (i, j): wind winter i of a set of n_a with solar winter j of a set of n_b, ns solar winters
    per wind winter spread evenly over the set (never i == j when both come from the same set)."""
    k = min(ns, n_b - 1 if same else n_b)
    stride = max(1, n_b // max(k, 1))
    out = []
    for i in range(n_a):
        base = int(round(i * n_b / n_a))
        js, s = [], 1 if same else 0
        while len(js) < k:
            j = (base + s * stride + (s * stride) // n_b) % n_b
            if not (same and j == i) and j not in js:
                js.append(j)
            s += 1
            if s > 4 * n_b:
                break
        out += [(i, j) for j in js]
    return np.array(out)
