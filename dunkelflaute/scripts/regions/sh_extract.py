"""
Southern-Hemisphere winter (April-September) daily series for the SH test of the geostrophic
mechanism: region-mean daily wind for the regions of sh_regions.py, and daily-mean sea-level
pressure on a 2.5 deg SH grid (65-15 S, global ring), on the daily footing of the CMIP6 analysis.

Wind, per region (land cells, weights cos(lat); DART unweighted, as extract_dart.py):
  spd10   daily-mean 10 m speed
  spd100  log law of the daily-mean 10 m speed to 100 m with the ERA5 fsr land z0 climatology
  wind    onshore power-curve CF of spd100 (as extract_daily.py)
Land: ERA5 land fraction >= 0.5 within half a model cell (extract_daily.land_fraction); DART:
ERA5 lsm by nearest neighbour (tco1279_era5_regrid_index.npz), z0 from the same ERA5 cell.
Sea-level pressure:
  CMIP6   day/psl, bilinear, periodic in longitude (as extract_psl.py)
  ours    box mean of all source points within +-1.25 deg of each 2.5 deg point (as extract_msl.py)

Sources and daily means:
  cmip6 MODEL    day/sfcWind, day/psl; historical 1951-1980, ssp585 2071-2100; global tas copied
                 from daily_cmip6_<MODEL>.npz
  tco95 RUN      ctl1950g / hi1950g 1952-2099: speed = mean of the 8 three-hourly 10u/10v instants
                 in (d 00, d+1 00]; msl from the daily file atm_remapped_1d_msl
  tco319 RUN     ctl1950d 1950-1969 / ssp585d 2080-2099: 6-hourly 10u/10v and msl, the 4 instants
                 in (d 00, d+1 00]
  eerie RUN      ctl / hist / ssp245: daily averages mean10ws and mmsl (GRIB via the byte index,
                 section 2 stripped)
  dart9km RUN    1950c 1950-1969 / 2080c 2080-2092: native 3-hourly 10u/10v on the reduced
                 Gaussian grid (atm_reduced_3h, all months exist; the 9 km wind and solar
                 RECONSTRUCTIONS cover Sep-Mar only, the raw 10 m wind does not); msl from the
                 3-hourly remapped file, mean of the 06/12/18/24 UTC instants

Usage: python3 sh_extract.py {cmip6 MODEL | tco95 RUN | tco319 RUN | eerie RUN | dart9km RUN}  (legacy python)
Output: dunkelflaute/data/sh/sh_<source>_<name>.npz with {exp}|ymd (D, 3), {exp}|spd10, spd100,
        wind (D, R), {exp}|psl (D, 21, 144) float32 Pa; regions, n_cells, plat, plon
"""
import glob
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "1")
import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("", "/regions", "/cmip6"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts{p}")
from core.capacity_factor import ONSHORE, wind_cf   # noqa: E402
import regions as rg                                # noqa: E402
import sh_regions as sh                             # noqa: E402
from extract_daily import land_fraction            # noqa: E402

INV = f"{_ROOT}/dunkelflaute/data/cmip6/inventory.json"
RT = "/work/ab0995/a270092/runtime/awicm3-frontiers-xios"
H = "/work/ab0995/ICCP_AWI_hackthon_2025"
EIDX = f"{rg.BASE_DIR}/eerie_index"
DART_INDEX = f"{_ROOT}/wind/data/results/tco1279_era5_regrid_index.npz"
CMIP_YEARS = {"historical": (1951, 1980), "ssp585": (2071, 2100)}
TCO319 = {"ctl1950d": ("TCo319-DART-ctl1950d-gibbs-charn", range(1950, 1970)),
          "ssp585d": ("TCo319-DART-ssp585d-gibbs-charn", range(2080, 2100))}
DART = {"1950c": ("TCo1279-DART-1950C", range(1950, 1970)), "2080c": ("TCo1279-DART-2080C", range(2080, 2093))}
REG = list(sh.sh_boxes())
NB = len(sh.PLAT) * len(sh.PLON)
_S = None


def _dates(tv):
    return np.array([np.datetime64(str(x)[:19]) for x in nc.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)])


def _ymd(days):
    s = days.astype(str)
    return np.array([(int(x[:4]), int(x[5:7]), int(x[8:10])) for x in s]).reshape(-1, 3)


# ---------------- grids
def setup_regular(lat, lon, tag):
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    rows = np.where((lat >= sh.LAT_BAND[0] - 1) & (lat <= sh.LAT_BAND[1] + 1))[0]
    r0, r1 = int(rows.min()), int(rows.max()) + 1
    flat = np.repeat(lat[r0:r1], len(lon))
    flon = np.tile(lon, r1 - r0)
    dlat, dlon = float(np.median(np.abs(np.diff(lat)))), float(np.median(np.abs(np.diff(lon))))
    cache = f"{sh.OUT_DIR}/static_{tag}.npz"
    if os.path.exists(cache):
        c = np.load(cache)
        land, lu, z0 = c["land"], c["lu"], c["z0"]
        cells = sh.sh_cells(flat, flon, land)
        assert np.array_equal(lu, np.unique(np.concatenate(list(cells.values()))))
    else:
        cand = np.zeros(len(flat), bool)
        L = rg.lon180(flon)
        for la0, la1, lo0, lo1 in sh.sh_boxes().values():
            cand |= (flat >= la0 - dlat) & (flat <= la1 + dlat) & (L >= lo0 - dlon) & (L <= lo1 + dlon)
        frac = np.zeros(len(flat))
        frac[cand] = land_fraction(flat[cand], flon[cand], dlat, dlon)
        land = frac >= 0.5
        cells = sh.sh_cells(flat, flon, land)
        lu = np.unique(np.concatenate(list(cells.values())))
        z0 = sh.era5_z0_sh(flat[lu], flon[lu], dlat, dlon)
        os.makedirs(sh.OUT_DIR, exist_ok=True)
        np.savez(cache, land=land, lu=lu, z0=z0)
    # psl band for box means (our runs)
    prow = np.where((lat >= sh.PLAT[0] - 1.25) & (lat <= sh.PLAT[-1] + 1.25))[0]
    p0, p1 = int(prow.min()), int(prow.max()) + 1
    bid = sh.psl_boxes(np.repeat(lat[p0:p1], len(lon)), np.tile(lon, p1 - p0))
    psel = np.where(bid >= 0)[0]
    return dict(r0=r0, r1=r1, flat=flat, cells=cells, lu=lu, z0=z0, w=np.cos(np.deg2rad(flat[lu])),
                lpos=[np.searchsorted(lu, cells[r]) for r in REG], p0=p0, p1=p1, psel=psel, pbid=bid[psel],
                pcnt=np.bincount(bid[psel], minlength=NB).astype(float))


def region_series(S, spd10, month):
    """spd10 (T, n_lu) daily means; month (T,)"""
    spd10 = np.ma.filled(np.ma.asarray(spd10, dtype=np.float64), np.nan)
    ll = rg.loglaw_land(spd10, S["z0"][month - 1])
    w = S["w"]

    def rm(x):
        return np.stack([(x[:, i] * w[i]).sum(axis=1) / w[i].sum() if len(i) else np.full(len(x), np.nan)
                         for i in S["lpos"]], axis=1)
    return dict(spd10=rm(spd10), spd100=rm(ll), wind=rm(wind_cf(ll, **ONSHORE)))


def box_mean(S, field2d):
    """field (ny_band, nx) on the source band rows p0:p1 -> (21, 144) float32."""
    x = field2d.ravel()[S["psel"]]
    s = np.bincount(S["pbid"], weights=x, minlength=NB)
    assert (S["pcnt"] > 0).all(), "every 2.5 deg box needs source points"
    return (s / S["pcnt"]).reshape(len(sh.PLAT), len(sh.PLON)).astype(np.float32)


# ---------------- CMIP6
def _span(fn):
    s = os.path.basename(fn).rsplit("_", 1)[-1].replace(".nc", "").split("-")
    return int(s[0][:4]), int(s[1][:4])


def cmip6_setup(model):
    inv = json.load(open(INV))[model]
    f = sorted(glob.glob(inv["paths"]["ssp585"]["day/sfcWind"] + "/*.nc"))[0]
    with nc.Dataset(f) as d:
        return setup_regular(d["lat"][:], d["lon"][:], f"cmip6_{model}")


def _cmip6_init(model):
    global _S
    _S = cmip6_setup(model)


def _bilinear(x, lat, lon):
    """x (D, ny, nx) on (lat, lon) -> (D, 21, 144) on the SH psl grid, periodic in longitude."""
    from scipy.interpolate import RegularGridInterpolator
    if lat[0] > lat[-1]:
        lat, x = lat[::-1], x[:, ::-1]
    lon = lon % 360
    o = np.argsort(lon)
    lon, x = lon[o], x[:, :, o]
    lonp = np.concatenate([lon[-1:] - 360, lon, lon[:1] + 360])
    xp = np.concatenate([x[:, :, -1:], x, x[:, :, :1]], axis=2)
    LA, LO = np.meshgrid(sh.PLAT, sh.PLON % 360, indexing="ij")
    pts = np.column_stack([LA.ravel(), LO.ravel()])
    out = np.empty((len(x), len(sh.PLAT), len(sh.PLON)), np.float32)
    for k in range(len(x)):
        out[k] = RegularGridInterpolator((lat, lonp), xp[k], bounds_error=True)(pts).reshape(len(sh.PLAT), len(sh.PLON))
    return out


def _cmip6_job(args):
    var, fn, year = args
    S = _S
    with nc.Dataset(fn) as a:
        t = a["time"]
        ds = nc.num2date(t[:], t.units, calendar=getattr(t, "calendar", "standard"))
        ymd = np.array([(x.year, x.month, x.day) for x in ds])
        keep = (ymd[:, 0] == year) & np.isin(ymd[:, 1], sh.MONTHS)
        if not keep.any():
            return var, np.zeros((0, 3), int), None
        idx = np.where(keep)[0]
        sl = slice(int(idx[0]), int(idx[-1]) + 1)
        k = keep[sl]
        if var == "sfcWind":
            x = a[var][sl, S["r0"]:S["r1"], :].reshape(sl.stop - sl.start, -1)[:, S["lu"]].astype(np.float64)[k]
            return var, ymd[sl][k], x
        lat = a["lat"][:].astype(float)
        rows = np.where((lat >= sh.PLAT[0] - 5) & (lat <= sh.PLAT[-1] + 5))[0]
        r0, r1 = int(rows.min()), int(rows.max()) + 1
        x = a["psl"][sl, r0:r1, :].astype(np.float64)[k]
        lon = a["lon"][:].astype(float)
    return var, ymd[sl][k], _bilinear(x, lat[r0:r1], lon)


def cmip6(model):
    inv = json.load(open(INV))[model]
    S = cmip6_setup(model)
    old = np.load(f"{_ROOT}/dunkelflaute/data/cmip6/daily_cmip6_{model}.npz")
    res = {}
    for exp, (y0, y1) in CMIP_YEARS.items():
        wdir = inv["paths"][exp]["day/sfcWind"]
        pbase = wdir.split("/day/")[0]
        pdir = sorted(x for x in glob.glob(f"{pbase}/day/psl/*/*") if glob.glob(f"{x}/*.nc"))[-1]
        jobs = []
        for var, d in (("sfcWind", wdir), ("psl", pdir)):
            for f in sorted(glob.glob(f"{d}/*.nc")):
                a, b = _span(f)
                jobs += [(var, f, y) for y in range(max(a, y0), min(b, y1) + 1)]
        got = {"sfcWind": {}, "psl": {}}
        with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 32)), initializer=_cmip6_init, initargs=(model,)) as ex:
            for var, ymd, x in ex.map(_cmip6_job, jobs):
                for k in range(len(ymd)):
                    got[var][tuple(ymd[k])] = x[k]
        keys = sorted(set(got["sfcWind"]) & set(got["psl"]))
        miss = len(set(got["sfcWind"]) ^ set(got["psl"]))
        ymd = np.array(keys)
        r = region_series(S, np.stack([got["sfcWind"][k] for k in keys]), ymd[:, 1])
        r["psl"] = np.stack([got["psl"][k] for k in keys])
        r["ymd"] = ymd
        r["tas_years"], r["tas"] = old[f"{exp}|tas_years"], old[f"{exp}|tas"]
        res[exp] = r
        print(f"  {model} {exp}: {len(keys)} days ({miss} in only one variable), {keys[0]} .. {keys[-1]}", flush=True)
    return S, res


# ---------------- TCO95
def tco95_setup(run):
    with nc.Dataset(f"{RT}/{run}/outdata/oifs/atm_remapped_3h_10u_3h_1960-1960.nc") as d:
        return setup_regular(d["lat"][:], d["lon"][:], "tco95_grid")


def _tco95_init(run):
    global _S
    _S = tco95_setup(run)


def _daily_speed(times, u, v, n_per_day):
    day = (times - np.timedelta64(1, "s")).astype("datetime64[D]")
    spd = np.hypot(u, v)
    ud, inv, cnt = np.unique(day, return_inverse=True, return_counts=True)
    sp = np.zeros((len(ud), spd.shape[1]))
    np.add.at(sp, inv, spd)
    ok = cnt == n_per_day
    return ud[ok], sp[ok] / n_per_day


def _month_slice(times, months=sh.MONTHS):
    m = rg.month_of(times)
    idx = np.where(np.isin(m, months))[0]
    return slice(int(idx[0]), int(idx[-1]) + 1) if len(idx) else None


def _tco95_year(args):
    run, year = args
    S = _S
    with nc.Dataset(f"{RT}/{run}/outdata/oifs/atm_remapped_3h_10u_3h_{year}-{year}.nc") as du, \
            nc.Dataset(f"{RT}/{run}/outdata/oifs/atm_remapped_3h_10v_3h_{year}-{year}.nc") as dv:
        t = _dates(du["time_instant"])
        sl = _month_slice(t)
        n = sl.stop - sl.start
        u = du["10u"][sl, S["r0"]:S["r1"], :].reshape(n, -1)[:, S["lu"]].astype(np.float64)
        v = dv["10v"][sl, S["r0"]:S["r1"], :].reshape(n, -1)[:, S["lu"]].astype(np.float64)
    ud, sp = _daily_speed(t[sl], u, v, 8)
    with nc.Dataset(f"{RT}/{run}/outdata/oifs/atm_remapped_1d_msl_1d_{year}-{year}.nc") as d:
        tv = d["time_instant"] if "time_instant" in d.variables else d["time_counter"]
        pd_ = (_dates(tv) - np.timedelta64(1, "s")).astype("datetime64[D]")
        pk = np.where(np.isin(pd_, ud))[0]
        psl = {pd_[k]: box_mean(S, d["msl"][k, S["p0"]:S["p1"], :].astype(np.float64)) for k in pk}
    keep = np.array([x in psl for x in ud])
    ud, sp = ud[keep], sp[keep]
    ymd = _ymd(ud)
    res = region_series(S, sp, ymd[:, 1])
    res["psl"] = np.stack([psl[x] for x in ud])
    res["ymd"] = ymd
    return res


def tco95(run):
    S = tco95_setup(run)
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 32)), initializer=_tco95_init, initargs=(run,)) as ex:
        parts = list(ex.map(_tco95_year, [(run, y) for y in range(1952, 2100)]))
    res = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    print(f"  tco95 {run}: {len(res['ymd'])} days", flush=True)
    return S, {"run": res}


# ---------------- TCo319
def _tco319_path(run, var, year):
    return f"{H}/{TCO319[run][0]}/outdata/oifs/atm_remapped_6h_{var}_6h_{year}-{year}.nc"


def tco319_setup(run):
    with nc.Dataset(_tco319_path(run, "10u", TCO319[run][1][0])) as d:
        return setup_regular(d["lat"][:], d["lon"][:], "tco319_grid")


def _tco319_init(run):
    global _S
    _S = tco319_setup(run)


def _tco319_year(args):
    run, year = args
    S = _S
    with nc.Dataset(_tco319_path(run, "10u", year)) as du, nc.Dataset(_tco319_path(run, "10v", year)) as dv, \
            nc.Dataset(_tco319_path(run, "msl", year)) as dm:
        tv = du["time_instant"] if "time_instant" in du.variables else du["time_counter"]
        t = _dates(tv)
        sl = _month_slice(t)
        n = sl.stop - sl.start
        u = du["10u"][sl, S["r0"]:S["r1"], :].reshape(n, -1)[:, S["lu"]].astype(np.float64)
        v = dv["10v"][sl, S["r0"]:S["r1"], :].reshape(n, -1)[:, S["lu"]].astype(np.float64)
        ud, sp = _daily_speed(t[sl], u, v, 4)
        tm = _dates(dm["time_instant"] if "time_instant" in dm.variables else dm["time_counter"])
        assert np.array_equal(tm, t)
        day = (t - np.timedelta64(1, "s")).astype("datetime64[D]")
        psl = []
        for dd in ud:
            k = np.where(day == dd)[0]
            acc = sum(dm["msl"][i, S["p0"]:S["p1"], :].astype(np.float64) for i in k) / len(k)
            psl.append(box_mean(S, acc))
    ymd = _ymd(ud)
    res = region_series(S, sp, ymd[:, 1])
    res["psl"] = np.stack(psl)
    res["ymd"] = ymd
    return res


def tco319(run):
    S = tco319_setup(run)
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 20)), initializer=_tco319_init, initargs=(run,)) as ex:
        parts = list(ex.map(_tco319_year, [(run, y) for y in TCO319[run][1]]))
    res = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    print(f"  tco319 {run}: {len(res['ymd'])} days", flush=True)
    return S, {"run": res}


# ---------------- EERIE daily
def eerie_setup():
    lat = np.linspace(90, -90, 721)
    lon = np.arange(1440) * 0.25
    S = setup_regular(lat, lon, "eerie_era5grid")
    S["glob_idx"] = S["r0"] * 1440 + S["lu"]
    return S


def _eerie_init():
    global _S
    _S = eerie_setup()


def _eerie_chunk(args):
    import eccodes
    run, var, rows = args
    ix = np.load(f"{EIDX}/{run}_{var}.npz")
    S = _S
    out = {}
    for k in rows:
        with open(ix["path"][k], "rb") as fh:
            fh.seek(int(ix["offset"][k]))
            h = eccodes.codes_new_from_message(rg.grib2_strip_local(fh.read(int(ix["size"][k]))))
        try:
            d = eccodes.codes_get(h, "dataDate")
            if (d // 100 % 100) not in sh.MONTHS:
                continue
            assert eccodes.codes_get(h, "numberOfDataPoints") == 1440 * 721
            v = eccodes.codes_get_values(h)
            key = (d // 10000, d // 100 % 100, d % 100)
            if var == "mean10ws":
                out[key] = v[S["glob_idx"]]
            else:
                out[key] = box_mean(S, v.reshape(721, 1440)[S["p0"]:S["p1"]])
        finally:
            eccodes.codes_release(h)
    return var, out


def eerie(run):
    S = eerie_setup()
    data = {}
    for var in ("mean10ws", "mmsl"):
        n = len(np.load(f"{EIDX}/{run}_{var}.npz")["path"])
        jobs = [(run, var, list(range(a, min(n, a + 365)))) for a in range(0, n, 365)]
        got = {}
        with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 64)), initializer=_eerie_init) as ex:
            for v, out in ex.map(_eerie_chunk, jobs):
                got.update(out)
        data[var] = got
    keys = sorted(set(data["mean10ws"]) & set(data["mmsl"]))
    ymd = np.array(keys)
    res = region_series(S, np.stack([data["mean10ws"][k] for k in keys]), ymd[:, 1])
    res["psl"] = np.stack([data["mmsl"][k] for k in keys])
    res["ymd"] = ymd
    print(f"  eerie {run}: {len(keys)} days, {keys[0]} .. {keys[-1]}", flush=True)
    return S, {"run": res}


# ---------------- DART 9 km
def dart_setup(run):
    d = np.load(DART_INDEX)
    lat, lon = d["tco_lat"].astype(float), d["tco_lon"].astype(float)
    cells = sh.sh_cells(lat, lon, d["land"])
    allc = np.unique(np.concatenate(list(cells.values())))
    lo, hi = int(allc.min()), int(allc.max()) + 1
    cache = f"{sh.OUT_DIR}/static_dart9km.npz"
    if os.path.exists(cache) and np.array_equal(np.load(cache)["allc"], allc):
        z0 = np.load(cache)["z0"]
    else:
        clim = sh.era5_z0_clim_global()
        z0 = clim[:, d["i"][allc], d["j"][allc]]
        os.makedirs(sh.OUT_DIR, exist_ok=True)
        np.savez(cache, allc=allc, z0=z0)
    f = sorted(glob.glob(f"{H}/{DART[run][0]}/outdata/oifs/atm_remapped_3h_msl_3h_*.nc"))[0]
    with nc.Dataset(f) as m:
        S = setup_psl_only(m["lat"][:], m["lon"][:])
    S.update(lo=lo, hi=hi, union=allc - lo, lu=allc, z0=z0, w=np.ones(len(allc)), cells=cells,
             lpos=[np.searchsorted(allc, cells[r]) for r in REG], lat_check=lat[lo:hi][allc - lo])
    return S


def setup_psl_only(lat, lon):
    lat, lon = np.asarray(lat, float), np.asarray(lon, float)
    prow = np.where((lat >= sh.PLAT[0] - 1.25) & (lat <= sh.PLAT[-1] + 1.25))[0]
    p0, p1 = int(prow.min()), int(prow.max()) + 1
    bid = sh.psl_boxes(np.repeat(lat[p0:p1], len(lon)), np.tile(lon, p1 - p0))
    psel = np.where(bid >= 0)[0]
    return dict(p0=p0, p1=p1, psel=psel, pbid=bid[psel], pcnt=np.bincount(bid[psel], minlength=NB).astype(float))


def _dart_init(run):
    global _S
    _S = dart_setup(run)


def _dart_month(args):
    run, year, month = args
    t0 = time.time()
    S = _S
    base = f"{H}/{DART[run][0]}/outdata/oifs"
    fu, fv = (f"{base}/atm_reduced_3h_{v}_3h_{year}{month:02d}-{year}{month:02d}.nc" for v in ("10u", "10v"))
    fm = f"{base}/atm_remapped_3h_msl_3h_{year}{month:02d}-{year}{month:02d}.nc"
    if not all(os.path.exists(x) for x in (fu, fv, fm)):
        return None, f"{year}-{month:02d} missing"
    with nc.Dataset(fu) as du, nc.Dataset(fv) as dv:
        t = _dates(du["time_counter"])
        assert np.array_equal(t, _dates(dv["time_counter"]))
        if month == 4:
            assert np.allclose(du["lat"][S["lo"]:S["hi"]][S["union"]], S["lat_check"], atol=1e-4)
        spd = np.empty((len(t), len(S["union"])))
        for a in range(0, len(t), 31):
            b = min(len(t), a + 31)
            u = du["10u"][a:b, S["lo"]:S["hi"]][:, S["union"]].astype(np.float64)
            v = dv["10v"][a:b, S["lo"]:S["hi"]][:, S["union"]].astype(np.float64)
            spd[a:b] = np.hypot(u, v)
    ud, sp = _daily_speed(t, spd, np.zeros_like(spd), 8)
    with nc.Dataset(fm) as dm:
        tm = _dates(dm["time_counter"])
        hour = (tm - tm.astype("datetime64[D]")) // np.timedelta64(1, "h")
        day = (tm - np.timedelta64(1, "s")).astype("datetime64[D]")
        psl = {}
        for dd in ud:
            k = np.where((day == dd) & (hour % 6 == 0))[0]
            if len(k) != 4:
                continue
            acc = sum(dm["msl"][i, S["p0"]:S["p1"], :].astype(np.float64) for i in k) / 4
            psl[dd] = box_mean(S, acc)
    keep = np.array([x in psl for x in ud])
    ud, sp = ud[keep], sp[keep]
    ymd = _ymd(ud)
    res = region_series(S, sp, ymd[:, 1])
    res["psl"] = np.stack([psl[x] for x in ud])
    res["ymd"] = ymd
    return res, f"{year}-{month:02d} {len(ud)} days {time.time() - t0:.0f}s"


def dart9km(run):
    S = dart_setup(run)
    jobs = [(run, y, m) for y in DART[run][1] for m in sh.MONTHS]
    parts = []
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 24)), initializer=_dart_init, initargs=(run,)) as ex:
        for r, msg in ex.map(_dart_month, jobs):
            print(f"  {msg}", flush=True)
            if r is not None:
                parts.append(r)
    res = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    print(f"  dart9km {run}: {len(res['ymd'])} days", flush=True)
    return S, {"run": res}


def main():
    src, name = sys.argv[1], sys.argv[2]
    t0 = time.time()
    S, res = dict(cmip6=cmip6, tco95=tco95, tco319=tco319, eerie=eerie, dart9km=dart9km)[src](name)
    out = {f"{exp}|{k}": v for exp, r in res.items() for k, v in r.items()}
    fn = f"{sh.OUT_DIR}/sh_{src}_{name}.npz"
    np.savez(fn, **out, regions=np.array(REG), n_cells=np.array([len(S["cells"][r]) for r in REG]),
             plat=sh.PLAT, plon=sh.PLON)
    print(f"Saved {fn} in {time.time() - t0:.0f}s; cells " + ", ".join(f"{r} {len(S['cells'][r])}" for r in REG[:7])
          + f"; grid cells min {min(len(S['cells'][r]) for r in REG[7:])}")


if __name__ == "__main__":
    main()
