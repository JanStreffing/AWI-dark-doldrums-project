"""
Germany combined capacity factor + Mockert et al. (2023) Dunkelflaute events for
TCo1279-DART 1950C and 2080C, Sep-Mar extended winter.

DART-side counterpart of era5/compute_germany_dunkelflaute_2015_2026_mockert.py --
same method, so ERA5 and DART numbers are directly comparable:
  - real Germany land polygon (onshore) + EEZ polygon (offshore), not a bbox.
    Land/sea flag on the DART grid is the ERA5-lsm regrid already used by the
    wind reconstruction (core/domain.py); the polygons then clip it to Germany.
  - solar CF: reconstructed DART ssrd (model_v5 + exact monthly rescale), /10800 s,
    /1000, clip [0,1]; cos(lat)-weighted mean over onshore cells (as ERA5 side).
  - wind CF: 100 m u/v (log-law reconstruction), cubic power curve, hub height
    100 m for both onshore and offshore; mean over onshore / offshore cells.
  - weights: Mockert's own 44% solar / 50% onshore / 6% offshore.
  - events: 48 h (16-sample) trailing running mean of combined CF, threshold
    0.06, every below-threshold sample pulls its whole 16-sample window into the
    event (Mockert Sec. 2.4), per winter block so the mean never bridges Apr-Aug.

Winters are labelled by the year of their Jan-Mar (Sep Y-1 .. Mar Y). A winter is
usable only if BOTH calendar years are outside the spin-up (first 2 years of each
run, proposal Sec. 4.3) and all 7 months exist: 1950C -> 1953..1969 (17),
2080C -> 2083..2092 (10). Feb/Mar 1952 solar inputs are missing from the DART
archive, but 1952 belongs to the spin-up-excluded winters, so nothing usable is lost.

Usage: python compute_germany_dunkelflaute_dart_mockert.py {1950c|2080c}
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
os.chdir(_ROOT)
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts")
from core.capacity_factor import ONSHORE, OFFSHORE, SSRD_ACCUM_SECONDS_DART, solar_cf, wind_cf, wind_speed  # noqa: E402
from core.domain import BBOX_GERMANY, DomainMasks  # noqa: E402

RUNS = {
    "1950c": dict(years=range(1950, 1970), spinup={1950, 1951},
                  wind_dir="wind/data/dart_1950c_100m",
                  solar_dir="solar/data/dart_reconstructed_1950c/{year}"),
    "2080c": dict(years=range(2080, 2093), spinup={2080, 2081},
                  wind_dir="wind/data/dart_2080c_100m",
                  solar_dir="solar/data/dart_reconstructed_2080c/{year}"),
}
SOLAR_FILE = "ssrd_reduced_3h_{year}{month:02d}-{year}{month:02d}.nc"
WIND_FILE = "u100v100_{year}{month:02d}.nc"
LAND_PATH = "boundaries/germany/germany_land.geojson"
EEZ_PATH = "boundaries/germany/germany_eez.geojson"
RESULTS = "dunkelflaute/results"   # one subfolder per period: validation_era5_2015_2026 / 1950c / 2080c

WEIGHTS = dict(solar=0.44, onshore=0.50, offshore=0.06)   # Mockert et al. (2023)
WIN = 16          # 16 * 3 h = 48 h
THRESH = 0.06
DT_H = 3.0
WINTER_MONTHS = [9, 10, 11, 12, 1, 2, 3]


def usable_winters(years, spinup):
    out = []
    for label in range(min(years) + 1, max(years) + 1):
        if {label - 1, label} & spinup:
            continue
        months = [(label - 1, m) for m in (9, 10, 11, 12)] + [(label, m) for m in (1, 2, 3)]
        if all(y in years for y, _ in months):
            out.append((label, months))
    return out


def germany_masks():
    import shapely.geometry as sgeom
    import shapely.prepared
    dm = DomainMasks(BBOX_GERMANY)
    with open(LAND_PATH) as f:
        land_geom = shapely.prepared.prep(sgeom.shape(json.load(f)["features"][0]["geometry"]))
    with open(EEZ_PATH) as f:
        eez_geom = shapely.prepared.prep(sgeom.shape(json.load(f)["features"][0]["geometry"]))
    lon = np.where(dm.lon > 180, dm.lon - 360, dm.lon)
    onshore = np.zeros(dm.lat.shape, bool)
    offshore = np.zeros(dm.lat.shape, bool)
    for i in np.where(dm.bbox)[0]:
        pt = sgeom.Point(lon[i], dm.lat[i])
        if dm.land[i]:
            onshore[i] = land_geom.contains(pt)
        else:
            offshore[i] = eez_geom.contains(pt)
    lo, hi = dm.lo, dm.hi
    return dict(lo=lo, hi=hi, on=onshore[lo:hi], off=offshore[lo:hi],
                w_on=np.cos(np.deg2rad(dm.lat[lo:hi]))[onshore[lo:hi]])


_M = None


def _init():
    global _M
    _M = germany_masks()


def _month(args):
    run, year, month = args
    t0 = time.time()
    cfg = RUNS[run]
    fs = nc.Dataset(f"{cfg['solar_dir'].format(year=year)}/{SOLAR_FILE.format(year=year, month=month)}")
    fw = nc.Dataset(f"{cfg['wind_dir']}/{WIND_FILE.format(year=year, month=month)}")
    tv = fs["time_counter"]
    times = nc.num2date(tv[:], tv.units, calendar=getattr(tv, "calendar", "standard"),
                        only_use_cftime_datetimes=False)
    times = np.array([np.datetime64(f"{t.year:04d}-{t.month:02d}-{t.day:02d}T{t.hour:02d}:{t.minute:02d}:{t.second:02d}")
                      for t in times])
    lo, hi = _M["lo"], _M["hi"]
    ssrd = fs["ssrd"][:, lo:hi][:, _M["on"]]
    u = fw["u100"][:, lo:hi]
    v = fw["v100"][:, lo:hi]
    fs.close(); fw.close()
    if len(times) != u.shape[0]:
        raise RuntimeError(f"{run} {year}-{month:02d}: solar has {len(times)} steps, wind {u.shape[0]}")
    spd = wind_speed(u, v)
    w = _M["w_on"]
    cf_solar = (solar_cf(ssrd, SSRD_ACCUM_SECONDS_DART) * w).sum(axis=1) / w.sum()
    cf_on = wind_cf(spd[:, _M["on"]], **ONSHORE).mean(axis=1)
    cf_off = wind_cf(spd[:, _M["off"]], **OFFSHORE).mean(axis=1)
    return year, month, dict(time=times, solar=cf_solar, on=cf_on, off=cf_off), time.time() - t0


def mockert_events(roll, dt_h=DT_H):
    """Mockert Sec. 2.4 on ONE winter block: expand each below-threshold rolling-mean
    sample to its full 16-sample window, merge -> events. Returns (flagged, events)
    with events = [(start_idx, end_idx_exclusive)]."""
    n = len(roll)
    flagged = np.zeros(n, bool)
    for i in np.where(np.nan_to_num(roll, nan=np.inf) < THRESH)[0]:
        flagged[max(0, i - WIN + 1):i + 1] = True
    events, start = [], None
    for i in range(n):
        if flagged[i] and start is None:
            start = i
        elif not flagged[i] and start is not None:
            events.append((start, i)); start = None
    if start is not None:
        events.append((start, n))
    return flagged, events


def trailing_mean(x):
    m = np.full(len(x), np.nan)
    if len(x) >= WIN:
        c = np.cumsum(np.insert(x, 0, 0.0))
        m[WIN - 1:] = (c[WIN:] - c[:-WIN]) / WIN
    return m


def main():
    run = sys.argv[1].lower()
    cfg = RUNS[run]
    winters = usable_winters(cfg["years"], cfg["spinup"])
    print(f"{run}: {len(winters)} usable winters {winters[0][0]}..{winters[-1][0]}", flush=True)
    ym = sorted({m for _, months in winters for m in months})
    t0 = time.time()
    data = {}
    n_workers = int(os.environ.get("N_WORKERS", 32))
    with ProcessPoolExecutor(max_workers=n_workers, initializer=_init) as ex:
        for k, (y, m, d, dt_s) in enumerate(ex.map(_month, [(run, y, m) for y, m in ym]), 1):
            data[(y, m)] = d
            print(f"  [{k}/{len(ym)}] {y}-{m:02d} {dt_s:.0f}s (elapsed {time.time()-t0:.0f}s)", flush=True)

    cols = {k: [] for k in ("time", "winter", "solar", "on", "off", "combined", "roll", "flagged")}
    ev_rows = []
    for label, months in winters:
        ch = [data[m] for m in months]
        cat = {k: np.concatenate([c[k] for c in ch]) for k in ("time", "solar", "on", "off")}
        comb = WEIGHTS["solar"] * cat["solar"] + WEIGHTS["onshore"] * cat["on"] + WEIGHTS["offshore"] * cat["off"]
        roll = trailing_mean(comb)
        flagged, events = mockert_events(roll)
        for s, e in events:
            ev_rows.append((label, cat["time"][s], cat["time"][e - 1], (e - s) * DT_H,
                            float(comb[s:e].mean()), float(np.nanmin(roll[s:e]))))
        for k, v in (("time", cat["time"]), ("solar", cat["solar"]), ("on", cat["on"]), ("off", cat["off"]),
                     ("combined", comb), ("roll", roll), ("flagged", flagged)):
            cols[k].append(v)
        cols["winter"].append(np.full(len(comb), label))
        print(f"  winter {label}: {len(events)} events, {int(flagged.sum())*DT_H:.0f} h flagged, "
              f"mean CF {comb.mean():.3f}", flush=True)

    os.makedirs(f"{RESULTS}/{run}", exist_ok=True)
    out = f"{RESULTS}/{run}/combined_cf_{run}_sepmar_mockertweights.npz"
    ev = np.array(ev_rows, dtype=[("winter", "i4"), ("start", "datetime64[s]"), ("end", "datetime64[s]"),
                                  ("duration_h", "f8"), ("mean_cf", "f8"), ("min_roll", "f8")])
    np.savez(out, **{k: np.concatenate(v) for k, v in cols.items()}, events=ev,
             winters=np.array([w for w, _ in winters]),
             weights=np.array([WEIGHTS["solar"], WEIGHTS["onshore"], WEIGHTS["offshore"]]))
    print(f"Saved {out}  ({len(ev)} events, {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
