"""
Daily region series on a CMIP6-like footing: daily-mean 10 m wind speed and daily-mean surface
solar irradiance, Sep-Mar, for the regions of dunkelflaute/scripts/regions/regions.py.

  wind   log law of the DAILY-MEAN 10 m speed to 100 m with the ERA5 fsr land z0 climatology
         (ERA5 land cells within half a model cell), onshore power curve
  solar  daily-mean irradiance / 1000 (solar_cf of the daily total)
Land = ERA5 land fraction >= 0.5 within half a model cell (one definition for every model grid).
Weights cos(lat).

Sources:
  cmip6 MODEL   day/sfcWind, day/rsds (dunkelflaute/data/cmip6/inventory.json); historical 1950-1980
                and ssp585 2040-2100; plus the global annual-mean Amon/tas for those years
  eerie RUN     daily averages mean10ws, msdwswrf (GRIB via the byte index, section 2 stripped);
                hist 1950-2014, ctl 1950-2050, ssp245 2015-2070
  tco95 RUN     the 3-hourly 10u/10v/ssrd aggregated to days: mean speed over the 8 instants in
                (d 00, d+1 00], ssrd summed over the 8 windows ending in that interval; 1952-2099
The EERIE and TCO95 daily versions bridge to the 6-hourly results (same runs, daily footing).

Usage: python3 extract_daily.py {cmip6 MODEL | eerie RUN | tco95 RUN}      (legacy python)
Output: dunkelflaute/data/cmip6/daily_<source>_<name>.npz
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
from core.capacity_factor import ONSHORE, wind_cf   # noqa: E402
import regions as rg                                # noqa: E402

OUT = f"{_ROOT}/dunkelflaute/data/cmip6" + ("_grid5" if rg.REGION_SET == "grid5" else "")
INV = f"{_ROOT}/dunkelflaute/data/cmip6/inventory.json"
RT = "/work/ab0995/a270092/runtime/awicm3-frontiers-xios"
EIDX = f"{rg.BASE_DIR}/eerie_index"
YEARS = {"historical": (1950, 1980), "ssp585": (2040, 2100)}

_S = None


def land_fraction(lat, lon, dlat, dlon):
    with nc.Dataset(rg.ERA5_LSM) as d:
        elat = d["latitude"][:].astype(float)
        elon = d["longitude"][:].astype(float)
        elsm = d["lsm"][:].squeeze().astype(float)
    out = np.empty(len(lat))
    for k, (la, lo) in enumerate(zip(lat, np.asarray(lon) % 360)):
        mi = np.abs(elat - la) <= max(dlat / 2, 0.125)
        mj = np.abs(((elon - lo) + 180) % 360 - 180) <= max(dlon / 2, 0.125)
        out[k] = elsm[np.ix_(mi, mj)].mean()
    return out


def setup_grid(lat, lon, tag):
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    rows = np.where((lat >= rg.LAT_BAND[0] - 1) & (lat <= rg.LAT_BAND[1] + 1))[0]
    r0, r1 = int(rows.min()), int(rows.max()) + 1
    flat = np.repeat(lat[r0:r1], len(lon))
    flon = np.tile(lon, r1 - r0)
    dlat, dlon = float(np.median(np.abs(np.diff(lat)))), float(np.median(np.abs(np.diff(lon))))
    cache = f"{OUT}/static_{tag}.npz"
    if os.path.exists(cache):
        c = np.load(cache, allow_pickle=True)
        land = c["land"]
    else:
        # land fraction only where a region box could need it
        cand = np.zeros(len(flat), bool)
        L = rg.lon180(flon)
        for la0, la1, lo0, lo1 in list(rg.BOXES.values()) + [rg.DE_BOX]:
            cand |= (flat >= la0 - dlat) & (flat <= la1 + dlat) & (L >= lo0 - dlon) & (L <= lo1 + dlon)
        frac = np.zeros(len(flat))
        frac[cand] = land_fraction(flat[cand], flon[cand], dlat, dlon)
        land = frac >= 0.5
    cells = rg.region_cells(flat, flon, land)
    cells.pop("DE_off")
    lu = np.unique(np.concatenate(list(cells.values())))
    if os.path.exists(cache) and "z0" in np.load(cache).files and np.array_equal(np.load(cache)["lu"], lu):
        z0 = np.load(cache)["z0"]
    else:
        z0 = rg.era5_z0_on(flat[lu], flon[lu], dlat, dlon)
        os.makedirs(OUT, exist_ok=True)
        np.savez(cache, land=land, lu=lu, z0=z0)
    return dict(r0=r0, r1=r1, flat=flat, cells=cells, lu=lu, z0=z0, w=np.cos(np.deg2rad(flat[lu])),
                lpos=[np.searchsorted(lu, cells[r]) for r in rg.REGIONS])


def region_series(S, spd10, rsds, month):
    """spd10, rsds (T, n_lu) daily means; month (T,)"""
    ll = rg.loglaw_land(spd10, S["z0"][month - 1])
    w = S["w"]
    def rm(x):
        return np.stack([(x[:, i] * w[i]).sum(axis=1) / w[i].sum() if len(i) else np.full(len(x), np.nan)
                         for i in S["lpos"]], axis=1)
    return dict(wind=rm(wind_cf(ll, **ONSHORE)), solar=rm(np.clip(rsds / 1000.0, 0, 1)), spd10=rm(spd10))


# ---------------- CMIP6
def _cmip6_file(args):
    """One file of one variable -> (ymd of kept days, values on the region land cells)."""
    model, var, fn, y0, y1 = args
    S = _S
    with nc.Dataset(fn) as a:
        t = a["time"]
        dates = nc.num2date(t[:], t.units, calendar=getattr(t, "calendar", "standard"))
        ymd = np.array([(x.year, x.month, x.day) for x in dates])
        keep = np.isin(ymd[:, 1], rg.MONTHS) & (ymd[:, 0] >= y0) & (ymd[:, 0] <= y1)
        if not keep.any():
            return var, np.zeros((0, 3), int), np.zeros((0, len(S["lu"])))
        idx = np.where(keep)[0]
        sl = slice(int(idx[0]), int(idx[-1]) + 1)
        x = a[var][sl, S["r0"]:S["r1"], :].reshape(sl.stop - sl.start, -1)[:, S["lu"]].astype(np.float64)
        k = keep[sl]
    return var, ymd[sl][k], x[k]


def _cmip6_init(model):
    global _S
    _S = cmip6_setup(model)


def cmip6_setup(model):
    inv = json.load(open(INV))[model]
    f = sorted(glob.glob(inv["paths"]["ssp585"]["day/sfcWind"] + "/*.nc"))[0]
    with nc.Dataset(f) as d:
        return setup_grid(d["lat"][:], d["lon"][:], f"cmip6_{model}")


def _span(fn):
    s = os.path.basename(fn).rsplit("_", 1)[-1].replace(".nc", "").split("-")
    return int(s[0][:4]), int(s[1][:4])


def global_tas(model, exp, y0, y1):
    inv = json.load(open(INV))[model]
    out = {}
    for f in sorted(glob.glob(inv["paths"][exp]["Amon/tas"] + "/*.nc")):
        a, b = _span(f)
        if b < y0 or a > y1:
            continue
        with nc.Dataset(f) as d:
            t = d["time"]
            dates = nc.num2date(t[:], t.units, calendar=getattr(t, "calendar", "standard"))
            yrs = np.array([x.year for x in dates])
            w = np.cos(np.deg2rad(d["lat"][:].astype(float)))
            x = d["tas"][:].astype(np.float64)
            g = (x.mean(axis=2) * w[None]).sum(axis=1) / w.sum()
            for y in np.unique(yrs):
                if y0 <= y <= y1 and (yrs == y).sum() == 12:
                    out[int(y)] = float(g[yrs == y].mean())
    return out


def cmip6(model):
    inv = json.load(open(INV))[model]
    S = cmip6_setup(model)
    res = {}
    for exp, (y0, y1) in YEARS.items():
        jobs = []
        wind_only = inv.get("wind_only", False)          # no day/rsds: solar stays NaN
        for var, tv in (("sfcWind", "day/sfcWind"),) + (() if wind_only else (("rsds", "day/rsds"),)):
            for f in sorted(glob.glob(inv["paths"][exp][tv] + "/*.nc")):
                a, b = _span(f)
                if not (b < y0 or a > y1):
                    jobs.append((model, var, f, y0, y1))
        got = {"sfcWind": {}, "rsds": {}}
        with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 16)), initializer=_cmip6_init, initargs=(model,)) as ex:
            for var, ymd, x in ex.map(_cmip6_file, jobs):
                for k in range(len(ymd)):
                    got[var][tuple(ymd[k])] = x[k]
        if wind_only:
            keys, miss = sorted(got["sfcWind"]), 0
        else:
            keys = sorted(set(got["sfcWind"]) & set(got["rsds"]))
            miss = len(set(got["sfcWind"]) ^ set(got["rsds"]))
        ymd = np.array(keys)
        u = np.stack([got["sfcWind"][k] for k in keys])
        r = np.full_like(u, np.nan) if wind_only else np.stack([got["rsds"][k] for k in keys])
        res[exp] = region_series(S, u, r, ymd[:, 1])
        res[exp]["ymd"] = ymd
        res[exp]["tas"] = global_tas(model, exp, y0, y1)
        print(f"  {model} {exp}: {len(keys)} days ({miss} days in only one variable), {keys[0]} .. {keys[-1]}, "
              f"tas {len(res[exp]['tas'])} years", flush=True)
    return S, res


# ---------------- EERIE daily
def _eerie_chunk(args):
    import eccodes
    run, var, rows = args
    ix = np.load(f"{EIDX}/{run}_{var}.npz")
    S = _S
    out_ymd, out_x = [], []
    for k in rows:
        with open(ix["path"][k], "rb") as fh:
            fh.seek(int(ix["offset"][k]))
            h = eccodes.codes_new_from_message(rg.grib2_strip_local(fh.read(int(ix["size"][k]))))
        try:
            d = eccodes.codes_get(h, "dataDate")
            if (d // 100 % 100) not in rg.MONTHS:
                continue
            assert eccodes.codes_get(h, "numberOfDataPoints") == 1440 * 721
            v = eccodes.codes_get_values(h)
            out_ymd.append((d // 10000, d // 100 % 100, d % 100))
            out_x.append(v[S["glob_idx"]])
        finally:
            eccodes.codes_release(h)
    return var, out_ymd, out_x


def _eerie_init():
    global _S
    _S = eerie_setup()


def eerie_setup():
    lat = np.linspace(90, -90, 721)
    lon = np.arange(1440) * 0.25
    S = setup_grid(lat, lon, "eerie_era5grid")
    S["glob_idx"] = S["r0"] * 1440 + S["lu"]          # flat index into the global 721 x 1440 field
    return S


def eerie(run):
    S = eerie_setup()
    data = {}
    for var in ("mean10ws", "msdwswrf"):
        n = len(np.load(f"{EIDX}/{run}_{var}.npz")["path"])
        jobs = [(run, var, list(range(a, min(n, a + 365)))) for a in range(0, n, 365)]
        got = {}
        with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 64)), initializer=_eerie_init) as ex:
            for v, ymd, xs in ex.map(_eerie_chunk, jobs):
                for a, b in zip(ymd, xs):
                    got[a] = b
        data[var] = got
    keys = sorted(set(data["mean10ws"]) & set(data["msdwswrf"]))
    ymd = np.array(keys)
    u = np.stack([data["mean10ws"][k] for k in keys])
    r = np.stack([data["msdwswrf"][k] for k in keys])
    res = region_series(S, u, r, ymd[:, 1])
    res["ymd"] = ymd
    print(f"  eerie {run}: {len(keys)} days, {keys[0]} .. {keys[-1]}", flush=True)
    return S, {"run": res}


# ---------------- TCO95 daily from 3-hourly
def _tco95_year(args):
    run, year = args
    S = _S
    f = {}
    for v in ("10u", "10v", "ssrd"):
        with nc.Dataset(f"{RT}/{run}/outdata/oifs/atm_remapped_3h_{v}_3h_{year}-{year}.nc") as d:
            if v == "10u":
                tv = d["time_instant"]
                times = np.array([np.datetime64(str(x)[:19]) for x in nc.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)])
            f[v] = d[v][:, S["r0"]:S["r1"], :].reshape(len(times), -1)[:, S["lu"]].astype(np.float64)
    day = (times - np.timedelta64(1, "s")).astype("datetime64[D]")
    spd = np.hypot(f["10u"], f["10v"])
    ud, inv, cnt = np.unique(day, return_inverse=True, return_counts=True)
    ok = cnt == 8
    sp = np.zeros((len(ud), spd.shape[1]))
    ss = np.zeros_like(sp)
    np.add.at(sp, inv, spd)
    np.add.at(ss, inv, f["ssrd"])
    sp, ss, ud = sp[ok] / 8, ss[ok] / 86400.0, ud[ok]
    ymd = np.array([(int(str(x)[:4]), int(str(x)[5:7]), int(str(x)[8:10])) for x in ud])
    k = np.isin(ymd[:, 1], rg.MONTHS)
    res = region_series(S, sp[k], ss[k], ymd[k][:, 1])
    res["ymd"] = ymd[k]
    return res


def _tco95_init(run):
    global _S
    _S = tco95_setup(run)


def tco95_setup(run):
    with nc.Dataset(f"{RT}/{run}/outdata/oifs/atm_remapped_3h_10u_3h_1960-1960.nc") as d:
        return setup_grid(d["lat"][:], d["lon"][:], "tco95_grid")


def tco95(run):
    S = tco95_setup(run)
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 24)), initializer=_tco95_init, initargs=(run,)) as ex:
        parts = list(ex.map(_tco95_year, [(run, y) for y in range(1952, 2100)]))
    res = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    print(f"  tco95 {run}: {len(res['ymd'])} days", flush=True)
    return S, {"run": res}


def main():
    src, name = sys.argv[1], sys.argv[2]
    S, res = dict(cmip6=cmip6, eerie=eerie, tco95=tco95)[src](name)
    os.makedirs(OUT, exist_ok=True)
    out = {}
    for exp, r in res.items():
        for k, v in r.items():
            if k == "tas":
                out[f"{exp}|tas_years"] = np.array(sorted(v))
                out[f"{exp}|tas"] = np.array([v[y] for y in sorted(v)])
            else:
                out[f"{exp}|{k}"] = v
    fn = f"{OUT}/daily_{src}_{name}.npz"
    np.savez(fn, **out, regions=np.array(rg.REGIONS), n_cells=np.array([len(S["cells"][r]) for r in rg.REGIONS]))
    print(f"Saved {fn}; cells " + ", ".join(f"{r} {len(S['cells'][r])}" for r in rg.REGIONS))


if __name__ == "__main__":
    main()
