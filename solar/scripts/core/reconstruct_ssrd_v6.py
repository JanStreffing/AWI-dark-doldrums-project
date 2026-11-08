"""v6: reduced-feature, 6-hourly variant of the v4/v5 method, for the TCo319 (31 km) and
TCO95 (~100 km) AWI-CM3 runs, which write only tcc and lcc (no tsr, mcc, hcc) and only
every 6 hours (instantaneous samples at 00/06/12/18 UTC, time_instant axis).

Target is unchanged in kind: kt = ssrd / tisr, but over a 6-hour window ending at the
stamp (DT6 = 21600 s). Features are what those runs carry at the window's two ends:
  tcc, lcc            cloud cover at the window END (the stamp itself)
  tcc_prev, lcc_prev  cloud cover at the window START (stamp - 6 h)
  mu                  tisr6 / (DT6 * S0), the window-mean cos-zenith proxy
tisr stays analytic (solar_geometry.py), so it transfers to any grid.

From ERA5's 3-hourly files (clouds instantaneous every 3 h, ssrd as genuine 3-hour sums
in data/era5_3h_mean/), a 6 h window ending at t uses ssrd(t) + ssrd(t - 3 h) and clouds
at t and t - 6 h; the 3 h cloud sample in between is deliberately NOT used, since the
6-hourly model output does not have it.

v7 is the same method with all four cloud layers (tcc, lcc, mcc, hcc at both window ends),
for runs that do write mcc/hcc 6-hourly (the TCo319-DART d runs do; FEATURES_V7).
"""
import os as _os, sys as _sys
_root = _os.path.dirname(_os.path.abspath(__file__))
while not _os.path.isdir(_os.path.join(_root, "data")):
    _root = _os.path.dirname(_root)
_os.chdir(_root)
_sys.path.insert(0, _os.path.join(_root, "scripts", "core"))

import numpy as np

from solar_geometry import toa_irradiance_accumulated, SOLAR_CONSTANT

DT6 = 21600.0
DAY_MU = 0.02                   # same daylight cut as v4/v5, on the 6 h window mean
CLOUDS_V6 = ("tcc", "lcc")
FEATURES_V6 = ["tcc", "lcc", "tcc_prev", "lcc_prev", "mu"]
# v7: the same with all four cloud layers, for runs that write mcc/hcc 6-hourly (TCo319-DART)
CLOUDS_V7 = ("tcc", "lcc", "mcc", "hcc")
FEATURES_V7 = list(CLOUDS_V7) + [f"{c}_prev" for c in CLOUDS_V7] + ["mu"]


def to_6h_fields(times3h, clouds, ssrd3h=None):
    """Collapse 3-hourly ERA5 fields onto 6 h windows ending at 00/06/12/18 UTC.

    times3h (nt,) datetime64; clouds {name: (nt, ...)} instantaneous; ssrd3h (nt, ...) 3-hour
    sums ending at each stamp, or None. Only windows whose stamp, stamp - 3 h and stamp - 6 h
    all exist are kept (drops at most the first window of a month).
    Returns dict(times, <name>, <name>_prev for every cloud field[, ssrd]) on the 6 h stamps."""
    t = np.asarray(times3h, dtype="datetime64[s]")
    pos = {v: i for i, v in enumerate(t)}
    h3, h6 = np.timedelta64(3, "h"), np.timedelta64(6, "h")
    hour = (t - t.astype("datetime64[D]")) // np.timedelta64(1, "h")
    end = [i for i in np.where(hour % 6 == 0)[0] if (t[i] - h3) in pos and (t[i] - h6) in pos]
    mid = [pos[t[i] - h3] for i in end]
    start = [pos[t[i] - h6] for i in end]
    out = dict(times=t[end])
    for name, x in clouds.items():
        out[name], out[f"{name}_prev"] = x[end], x[start]
    if ssrd3h is not None:
        out["ssrd"] = ssrd3h[end].astype("float64") + ssrd3h[mid]
    return out


def to_6h(times3h, tcc, lcc, ssrd3h=None):
    """v6 form of to_6h_fields: tcc and lcc only."""
    return to_6h_fields(times3h, dict(tcc=tcc, lcc=lcc), ssrd3h)


def tisr_6h(lat, lon, times6):
    """Analytic TOA irradiance accumulated over the 6 h window ending at each stamp, J m-2."""
    return toa_irradiance_accumulated(lat, lon, times6, dt_seconds=DT6, n_substeps=24)


def build_features_6h(w6, tisr6, features):
    """Flattened daylight rows. Returns X (n, len(features)) float32 in that order, the day mask
    over the flattened (time, space) array, and kt (clipped [0, 1.5]) if w6 has ssrd."""
    tisr = tisr6.ravel()
    mu = tisr / (DT6 * SOLAR_CONSTANT)
    day = mu > DAY_MU
    cols = {c: (mu if c == "mu" else w6[c].ravel())[day] for c in features}
    X = np.column_stack([cols[c] for c in features]).astype("float32")
    kt = np.clip(w6["ssrd"].ravel()[day] / tisr[day], 0, 1.5) if "ssrd" in w6 else None
    return X, day, kt


def build_features_v6(w6, tisr6):
    return build_features_6h(w6, tisr6, FEATURES_V6)
