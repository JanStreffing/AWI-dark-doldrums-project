"""
Sub-daily check of the CMIP6 per-K result: region series for the CMIP6 models with 3-hourly wind and
solar (inventory.json 'subdaily': MPI-ESM1-2-HR, MPI-ESM1-2-LR, IPSL-CM6A-LR, MRI-ESM2-0), on three
footings for the SAME winters (historical 1961-1990, ssp585 2071-2100; MRI-ESM2-0 has 3hr output
only from 1960, so the control is 1961-1990 for all four, not 1951-1980 as in daily_events.py):

  6h    our runs' footing: wind CF from the 3hr/uas,vas instants at 00/06/12/18 UTC (log law with
        the ERA5 z0, onshore curve); solar = mean of the two 3hr/rsds windows ending at the instant
  d3h   daily from the same 3-hourly data: CF of the mean speed over the 8 instants in (d 00, d+1 00],
        solar = mean of the 8 windows in that interval
  day   the day/sfcWind and day/rsds files (what daily_events.py uses)
6h vs d3h isolates the footing (same data, same model); d3h vs day isolates daily-mean speed from
instants against the model's own daily mean of every time step.
Regions and land as extract_daily.py (ERA5 land fraction >= 0.5, cos-lat weights).

Usage: python3 extract_3hr.py MODEL      (legacy python)
Output: dunkelflaute/data/cmip6/sub3h_<MODEL>.npz
"""
import glob
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "1")
import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/cmip6")
from core.capacity_factor import ONSHORE, wind_cf   # noqa: E402
import regions as rg                                # noqa: E402
import extract_daily as xd                          # noqa: E402

WINTERS = {"historical": (1961, 1990), "ssp585": (2071, 2100)}
_S = None


def _init(model):
    global _S
    _S = xd.cmip6_setup(model)


def _dates(t, idx=None):
    v = t[:] if idx is None else t[idx]
    ds = nc.num2date(v, t.units, calendar=getattr(t, "calendar", "standard"))
    return np.array([(x.year, x.month, x.day, x.hour, x.minute) for x in ds])


def _winter(ymd):
    return np.where(ymd[:, 1] >= 9, ymd[:, 0] + 1, ymd[:, 0])


def _read(args):
    """One file, one calendar year of it -> (stamps (n,5), values (n, n_lu))."""
    var, fn, year, w0, w1 = args
    S = _S
    with nc.Dataset(fn) as a:
        t = a["time"]
        d = _dates(t)
        w = _winter(d)
        keep = (d[:, 0] == year) & np.isin(d[:, 1], rg.MONTHS) & (w >= w0) & (w <= w1)
        if not keep.any():
            return var, np.zeros((0, 5), int), np.zeros((0, len(S["lu"])))
        idx = np.where(keep)[0]
        sl = slice(int(idx[0]), int(idx[-1]) + 1)
        x = a[var][sl, S["r0"]:S["r1"], :].reshape(sl.stop - sl.start, -1)[:, S["lu"]].astype(np.float64)
        k = keep[sl]
    return var, d[sl][k], x[k]


def collect(model, exp, table, var, w0, w1, n_workers):
    inv = json.load(open(xd.INV))[model]
    jobs = []
    for f in sorted(glob.glob(inv["paths"][exp][f"{table}/{var}"] + "/*.nc")):
        a, b = xd._span(f)
        for y in range(max(a, w0 - 1), min(b, w1) + 1):
            jobs.append((var, f, y, w0, w1))
    got = {}
    with ProcessPoolExecutor(max_workers=n_workers, initializer=_init, initargs=(model,)) as ex:
        for _, st, x in ex.map(_read, jobs):
            for k in range(len(st)):
                got[tuple(st[k])] = x[k]
    keys = sorted(got)
    return np.array(keys), np.stack([got[k] for k in keys])


def to_dt64(st):
    return np.array([np.datetime64(f"{y:04d}-{m:02d}-{d:02d}T{h:02d}:{mi:02d}") for y, m, d, h, mi in st]).astype("datetime64[m]")


def rmean(S, x):
    w = S["w"]
    return np.stack([(x[:, i] * w[i]).sum(axis=1) / w[i].sum() if len(i) else np.full(len(x), np.nan)
                     for i in S["lpos"]], axis=1)


def main():
    model = sys.argv[1]
    nw = int(os.environ.get("N_WORKERS", 32))
    S = xd.cmip6_setup(model)
    out = {}
    for exp, (w0, w1) in WINTERS.items():
        su, u = collect(model, exp, "3hr", "uas", w0, w1, nw)
        sv, v = collect(model, exp, "3hr", "vas", w0, w1, nw)
        sr, r = collect(model, exp, "3hr", "rsds", w0, w1, nw)
        assert np.array_equal(su, sv), f"{model} {exp}: uas and vas stamps differ"
        tw = to_dt64(su)
        tr = to_dt64(sr)                       # window centres
        spd = np.hypot(u, v)
        mon = rg.month_of(tw.astype("datetime64[s]"))
        ll = rg.loglaw_land(spd, S["z0"][mon - 1])
        cf = wind_cf(ll, **ONSHORE)
        sol = np.clip(r / 1000.0, 0, 1)
        rpos = {t: i for i, t in enumerate(tr)}
        h90 = np.timedelta64(90, "m")

        # 6h footing
        hour = (tw - tw.astype("datetime64[D]")) // np.timedelta64(1, "h")
        ends = [i for i in np.where(hour % 6 == 0)[0]
                if (tw[i] - h90) in rpos and (tw[i] - 2 * h90 - h90) in rpos]
        ends = np.array(ends)
        a1 = np.array([rpos[tw[i] - h90] for i in ends])
        a2 = np.array([rpos[tw[i] - 3 * h90] for i in ends])
        out[f"{exp}|time6"] = tw[ends].astype("datetime64[s]")
        out[f"{exp}|wind6"] = rmean(S, cf[ends])
        out[f"{exp}|solar6"] = rmean(S, 0.5 * (sol[a1] + sol[a2]))

        # daily from the 3-hourly instants and windows
        dw = (tw.astype("datetime64[s]") - np.timedelta64(1, "s")).astype("datetime64[D]")
        dr = tr.astype("datetime64[D]")
        days, iw, cw = np.unique(dw, return_inverse=True, return_counts=True)
        sp_d = np.zeros((len(days), spd.shape[1]))
        np.add.at(sp_d, iw, spd)
        rd_pos = {d: i for i, d in enumerate(days)}
        so_d = np.zeros_like(sp_d)
        cr = np.zeros(len(days))
        for k, d in enumerate(dr):
            j = rd_pos.get(d)
            if j is not None:
                so_d[j] += sol[k]
                cr[j] += 1
        ok = (cw == 8) & (cr == 8)
        sp_d, so_d, days = sp_d[ok] / 8, so_d[ok] / 8, days[ok]
        md = rg.month_of((days + np.timedelta64(1, "D")).astype("datetime64[s]"))   # month of day d
        ll_d = rg.loglaw_land(sp_d, S["z0"][md - 1])
        out[f"{exp}|day_d3h"] = days
        out[f"{exp}|wind_d3h"] = rmean(S, wind_cf(ll_d, **ONSHORE))
        out[f"{exp}|solar_d3h"] = rmean(S, so_d)

        # the daily files, same winters
        s1, x1 = collect(model, exp, "day", "sfcWind", w0, w1, nw)
        s2, x2 = collect(model, exp, "day", "rsds", w0, w1, nw)
        common = sorted(set(map(tuple, s1[:, :3])) & set(map(tuple, s2[:, :3])))
        i1 = {tuple(k): i for i, k in enumerate(s1[:, :3])}
        i2 = {tuple(k): i for i, k in enumerate(s2[:, :3])}
        a = np.array([i1[k] for k in common])
        b = np.array([i2[k] for k in common])
        rs = xd.region_series(S, x1[a], x2[b], np.array(common)[:, 1])
        out[f"{exp}|ymd_day"] = np.array(common)
        out[f"{exp}|wind_day"] = rs["wind"]
        out[f"{exp}|solar_day"] = rs["solar"]

        T = xd.global_tas(model, exp, w0, w1)
        out[f"{exp}|tas_years"] = np.array(sorted(T))
        out[f"{exp}|tas"] = np.array([T[y] for y in sorted(T)])
        print(f"  {model} {exp}: {len(ends)} 6-h steps, {len(days)} d3h days ({(~ok).sum()} incomplete dropped), "
              f"{len(common)} daily-file days, tas {len(T)} years", flush=True)
    fn = f"{xd.OUT}/sub3h_{model}.npz"
    np.savez(fn, **out, regions=np.array(rg.REGIONS))
    print(f"Saved {fn}")


if __name__ == "__main__":
    main()
