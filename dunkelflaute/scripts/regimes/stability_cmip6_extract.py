"""
Question C of the transfer study (methods notes 13d): why is the CMIP6 transfer change at fixed G
negative? Monthly lower-tropospheric stability, surface fluxes, snow and vegetation per CMIP6 model
and region, historical 1951-1980 against ssp585 2071-2100 winters (Sep-Mar, winter = year of its
January, as geostrophic.winter_of).

Land cells: the region land cells of extract_daily.py (ERA5 land fraction >= 0.5 within half a model
cell, cos-lat weights). Where a variable's grid is the day/sfcWind grid, the cached land mask
dunkelflaute/data/cmip6/static_cmip6_<model>.npz is reused read-only, so the cells are those of the
transfer series; other grids get the same land-fraction recipe computed here.

Per cell and month (then region means):
  tas, ts, ps, hfss (positive upward), hfls, sfcWind (Amon), ta850 (Amon/ta at 85000 Pa)
  th850_ths  theta(850) - theta(tas, ps), theta = T (1e5 / p)^0.2857  (K)
  th850_tts  theta(850) - theta(ts, ps)
  dT850s     ta850 - tas
  snc (LImon, %), snw (LImon, kg m-2), lai, treeFrac, cropFrac, grassFrac (Lmon, if present)
  sftlf      model land fraction of the region cells (fx, any member/experiment)
ta850 quantities use only "stab" cells: ta850 valid in every month of both periods and control
winter-mean ps >= PS_MIN (default 90000 Pa, so 850 hPa lies at least ~50 hPa above the ground).
Everything else uses all land cells, and also the stab cells (suffix _st) for a like-for-like ratio.

Usage: python3 stability_cmip6_extract.py [MODEL ...]     (legacy python; all inventory models by default)
Output: dunkelflaute/data/regimes/stability_cmip6_parts/cmip6_<model>.npz
        arrays <exp>|<var> of shape (7 months Sep..Mar, R regions), climatological means over winters
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
for p in ("", "/regions", "/cmip6", "/regimes"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts{p}")
import regions as rg            # noqa: E402
import extract_daily as xd      # noqa: E402

INV = f"{_ROOT}/dunkelflaute/data/cmip6/inventory.json"
DC = f"{_ROOT}/dunkelflaute/data/cmip6"
OUTD = f"{_ROOT}/dunkelflaute/data/regimes/stability_cmip6_parts"
PERIODS = {"historical": (1951, 1980), "ssp585": (2071, 2100)}
MON = [9, 10, 11, 12, 1, 2, 3]
KAP = 0.2857
PS_MIN = float(os.environ.get("PS_MIN", 90000.0))
AMON = ("tas", "ts", "ps", "hfss", "hfls", "sfcWind")
LAND = (("LImon", "snc"), ("LImon", "snw"), ("Lmon", "lai"), ("Lmon", "treeFrac"), ("Lmon", "cropFrac"),
        ("Lmon", "grassFrac"))


def vdir(member, table, var, grid=None):
    """Latest version directory of member/table/var, preferring the given grid label."""
    ds = sorted(glob.glob(f"{member}/{table}/{var}/*/v*"))
    if not ds:
        return None
    if grid and any(d.split("/")[-2] == grid for d in ds):
        ds = [d for d in ds if d.split("/")[-2] == grid]
    return sorted(ds, key=lambda d: d.split("/")[-1])[-1]


def span(fn):
    s = os.path.basename(fn).rsplit("_", 1)[-1].replace(".nc", "").split("-")
    return int(s[0][:4]), int(s[1][:4])


class Grid:
    """Region land cells on one model grid (flat index into the lat band rows r0:r1)."""
    def __init__(self, model, lat, lon):
        lat, lon = np.asarray(lat, float), np.asarray(lon, float)
        rows = np.where((lat >= rg.LAT_BAND[0] - 1) & (lat <= rg.LAT_BAND[1] + 1))[0]
        self.r0, self.r1 = int(rows.min()), int(rows.max()) + 1
        flat = np.repeat(lat[self.r0:self.r1], len(lon))
        flon = np.tile(lon, self.r1 - self.r0)
        st = f"{DC}/static_cmip6_{model}.npz"
        land = None
        if os.path.exists(st):
            f = sorted(glob.glob(json.load(open(INV))[model]["paths"]["ssp585"]["day/sfcWind"] + "/*.nc"))[0]
            with nc.Dataset(f) as d:
                same = np.allclose(d["lat"][:], lat) and np.allclose(d["lon"][:], lon) if (len(d["lat"]) == len(lat) and len(d["lon"]) == len(lon)) else False
            if same:
                land = np.load(st)["land"]
                self.source = "static_cmip6 cache (day/sfcWind grid)"
        if land is None:
            dlat, dlon = float(np.median(np.abs(np.diff(lat)))), float(np.median(np.abs(np.diff(lon))))
            cand = np.zeros(len(flat), bool)
            L = rg.lon180(flon)
            for la0, la1, lo0, lo1 in list(rg.BOXES.values()) + [rg.DE_BOX]:
                cand |= (flat >= la0 - dlat) & (flat <= la1 + dlat) & (L >= lo0 - dlon) & (L <= lo1 + dlon)
            frac = np.zeros(len(flat))
            frac[cand] = xd.land_fraction(flat[cand], flon[cand], dlat, dlon)
            land = frac >= 0.5
            self.source = "ERA5 land fraction computed here"
        cells = rg.region_cells(flat, flon, land)
        cells.pop("DE_off")
        self.lu = np.unique(np.concatenate([cells[r] for r in rg.REGIONS]))
        self.pos = [np.searchsorted(self.lu, cells[r]) for r in rg.REGIONS]
        self.w = np.cos(np.deg2rad(flat[self.lu]))
        self.shape = (len(lat), len(lon))

    def rmean(self, x, mask=None):
        """x (..., n_lu) -> (..., R) weighted region means, NaN-aware, optional cell mask (n_lu,)."""
        out = []
        for p in self.pos:
            w = self.w[p] * (1.0 if mask is None else mask[p])
            v = x[..., p]
            ww = np.where(np.isfinite(v), w, 0.0)
            s = ww.sum(-1)
            out.append(np.where(s > 0, np.nansum(np.where(np.isfinite(v), v, 0.0) * ww, -1) / np.where(s > 0, s, 1), np.nan))
        return np.stack(out, -1)


def read_months(d_path, var, grid_cache, model, plev=None):
    """-> (G, dict exp -> (W, 7, n_lu) cell values per winter and month, NaN-filled) or None."""
    out = {}
    G = None
    for exp, (y0, y1) in PERIODS.items():
        dd = d_path[exp]
        if dd is None:
            return None, None
        fs = [f for f in sorted(glob.glob(dd + "/*.nc")) if not (span(f)[1] < y0 - 1 or span(f)[0] > y1)]
        if not fs:
            return None, None
        W = y1 - y0 + 1
        arr = None
        for f in fs:
            with nc.Dataset(f) as d:
                d.set_auto_mask(True)
                lat, lon = d["lat"][:], d["lon"][:]
                key = (len(lat), len(lon), float(lat[0]), float(lon[0]))
                if key not in grid_cache:
                    grid_cache[key] = Grid(model, lat, lon)
                G = grid_cache[key]
                t = d["time"]
                dates = nc.num2date(t[:], t.units, calendar=getattr(t, "calendar", "standard"))
                ym = np.array([(x.year, x.month) for x in dates])
                win = np.where(ym[:, 1] >= 9, ym[:, 0] + 1, ym[:, 0])
                keep = np.isin(ym[:, 1], MON) & (win >= y0) & (win <= y1)
                if not keep.any():
                    continue
                idx = np.where(keep)[0]
                sl = slice(int(idx[0]), int(idx[-1]) + 1)
                v = d[var]
                if plev is not None:
                    pl = d["plev"][:].astype(float)
                    if np.nanmax(pl) < 2000:            # hPa
                        pl = pl * 100
                    k = int(np.argmin(np.abs(pl - plev)))
                    assert abs(pl[k] - plev) < 1, f"{model} {var}: no {plev} Pa level ({pl})"
                    x = v[sl, k, G.r0:G.r1, :]
                else:
                    x = v[sl, G.r0:G.r1, :]
                x = np.ma.filled(np.ma.masked_invalid(x).astype(np.float64), np.nan)
                x = x.reshape(x.shape[0], -1)[:, G.lu]
                x[np.abs(x) > 1e15] = np.nan
                if arr is None:
                    arr = np.full((W, 7, len(G.lu)), np.nan)
                kk = keep[sl]
                for xi, (yy, mm), ww in zip(x[kk], ym[sl][kk], win[sl][kk]):
                    arr[ww - y0, MON.index(mm)] = xi
        if arr is None:
            return None, None
        out[exp] = arr
    return G, out


def model(name):
    t0 = time.time()
    inv = json.load(open(INV))[name]
    mem = {e: inv["paths"][e]["Amon/tas"].split("/Amon/")[0] for e in PERIODS}
    grid = inv["paths"]["historical"]["Amon/tas"].split("/")[-2]
    cache = {}
    raw, info = {}, {"grid_label": grid, "dirs": {}}
    for table, var, plev in [("Amon", v, None) for v in AMON] + [("Amon", "ta", 85000.0)] + [(t, v, None) for t, v in LAND]:
        dp = {e: vdir(mem[e], table, var, grid) for e in PERIODS}
        info["dirs"][f"{table}/{var}"] = dp
        G, a = read_months(dp, var, cache, name, plev)
        if a is None:
            continue
        raw[("ta850" if var == "ta" else var)] = (G, a)
    # sftlf: any member / experiment of this model
    sft = None
    for e in ("historical", "ssp585", "piControl"):
        base = mem["historical"].split("/historical/")[0]
        for dd in sorted(glob.glob(f"{base.replace('/CMIP/', '/*/')}/{e}/*/fx/sftlf/*/v*")):
            f = sorted(glob.glob(dd + "/*.nc"))
            if f:
                sft = f[0]
                break
        if sft:
            break
    Gt = raw["tas"][0]
    out = {"regions": np.array(rg.REGIONS), "months": np.array(MON)}
    # stab mask on the tas grid
    ok = "ta850" in raw and raw["ta850"][0].shape == Gt.shape and "ps" in raw and raw["ps"][0].shape == Gt.shape
    if ok:
        ta, ps = raw["ta850"][1], raw["ps"][1]
        valid = np.all(np.isfinite(ta["historical"]), axis=(0, 1)) & np.all(np.isfinite(ta["ssp585"]), axis=(0, 1))
        psc = np.nanmean(ps["historical"], axis=(0, 1))
        stab = valid & (psc >= PS_MIN)
        out["n_cells_land"] = np.array([len(p) for p in Gt.pos])
        out["n_cells_stab"] = np.array([int(stab[p].sum()) for p in Gt.pos])
    else:
        stab = None
        info["note_stab"] = "ta850/ps missing or on another grid"
    for exp in PERIODS:
        for v, (G, a) in raw.items():
            x = a[exp]
            out[f"{exp}|{v}"] = np.nanmean(G.rmean(x), axis=0)                    # (7, R)
            out[f"{exp}|{v}|winters"] = np.nanmean(G.rmean(x), axis=1)            # (W, R) winter means (equal month weights)
            if stab is not None and G.shape == Gt.shape and v != "ta850":
                out[f"{exp}|{v}_st"] = np.nanmean(G.rmean(x, stab), axis=0)
        if stab is not None:
            ta, tas, ts, ps = (raw[k][1][exp] for k in ("ta850", "tas", "ts", "ps"))
            th850 = ta * (1e5 / 85000.0) ** KAP
            for k, xx in (("th850_ths", th850 - tas * (1e5 / ps) ** KAP), ("th850_tts", th850 - ts * (1e5 / ps) ** KAP),
                          ("dT850s", ta - tas), ("ta850", ta)):
                out[f"{exp}|{k}"] = np.nanmean(Gt.rmean(xx, stab), axis=0)
                out[f"{exp}|{k}|winters"] = np.nanmean(Gt.rmean(xx, stab), axis=1)
    if sft:
        with nc.Dataset(sft) as d:
            lat, lon = d["lat"][:], d["lon"][:]
            key = (len(lat), len(lon), float(lat[0]), float(lon[0]))
            G = cache.get(key) or Grid(name, lat, lon)
            x = np.ma.filled(d["sftlf"][G.r0:G.r1, :].astype(np.float64), np.nan).reshape(-1)[G.lu]
        out["sftlf"] = G.rmean(x)
        out["sftlf_min"] = np.array([np.nanmin(x[p]) if len(p) else np.nan for p in G.pos])
        info["sftlf_file"] = sft
    info["grid_source"] = Gt.source
    info["vars"] = sorted(raw)
    out["info"] = np.array(json.dumps(info))
    os.makedirs(OUTD, exist_ok=True)
    np.savez(f"{OUTD}/cmip6_{name}.npz", **out)
    return name, time.time() - t0, sorted(raw), (int(out["n_cells_stab"].sum()), int(out["n_cells_land"].sum())) if ok else None


def main():
    names = sys.argv[1:] or list(json.load(open(INV)))
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 6))) as ex:
        for n, dt, vs, nc_ in ex.map(model, names):
            print(f"  {n:18s} {dt:5.0f}s  vars {','.join(vs)}  stab/land cells {nc_}  [{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
