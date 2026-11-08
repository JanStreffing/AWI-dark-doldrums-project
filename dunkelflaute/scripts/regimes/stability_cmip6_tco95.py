"""
The monthly stability diagnostics of stability_cmip6_extract.py for the TCO95 (~100 km) reruns
ctl1950g / hi1950g, so TCO95 can be placed on the CMIP6 cross-model relations (question C, methods
notes 13d). Land cells and cos-lat weights of extract_tco95.py (model lsm >= 0.5), all regions.

Per cell and month, then region means:
  2t, skt, msl, blh, fsr, sd, rsn, slhf, sshf, T850 (pl t at 85000 Pa), z850 (pl z)
  hfss       = -sshf / 10800 s: monthly sshf is the mean of the 3 h accumulations (J m-2), IFS sign
               positive downward; hfss is the CMIP sign (positive upward), W m-2
  snc        IFS 43r3 snow cover, 100 min(1, (1000 sd / rsn) / 0.1 m) in %, from monthly means
  ps_est     msl exp(-g zs / (Rd Tm)), Tm = 2t + 0.0065 zs / 2, zs = ERA5 orography conservatively
             remapped to the model output grid (no surface pressure or orography in the output)
  th850_ths  theta(T850, 850 hPa) - theta(2t, ps_est), theta = T (1e5 / p)^0.2857
  th850_tts  the same with skt
  inv850_z   T850 - 2t + g (z850 / g - zs - 2 m) / cp (dry static energy form, needs no ps)
  dT850s     T850 - 2t
ta850 quantities on "stab" cells only: winter-mean ps_est of ctl1950g 1953-1962 >= PS_MIN (90000 Pa).

ERA5 orography file (made once with cdo, kept in the parts directory):
  cdo -f nc4 -t ecmwf copy /pool/data/ERA5/E5/sf/an/IV/129/E5sf00_IV_INVARIANT_129.grb raw.nc
  cdo -f nc4 remapcon,<any TCO95 1m 2t file> -setgridtype,regular raw.nc era5_z_on_tco95.nc

Usage: python3 stability_cmip6_tco95.py {ctl1950g|hi1950g}      (N_WORKERS, default 8)
Output: dunkelflaute/data/regimes/stability_cmip6_parts/tco95_<run>.npz, arrays (Y, 12, R)
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "1")
import netCDF4 as nc
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("", "/regions", "/regimes"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts{p}")
import regions as rg            # noqa: E402
import extract_tco95 as xt      # noqa: E402

OUTD = f"{_ROOT}/dunkelflaute/data/regimes/stability_cmip6_parts"
ZS = f"{OUTD}/era5_z_on_tco95.nc"
G0, CP, RD, KAP = 9.80665, 1004.7, 287.05, 0.2857
PS_MIN = float(os.environ.get("PS_MIN", 90000.0))
SURF = ("2t", "skt", "msl", "blh", "fsr", "sd", "rsn", "slhf", "sshf")
YEARS = {"ctl1950g": range(1952, 2100), "hi1950g": range(1950, 2100)}
_S = None


def pl_path(run, var, year):
    return f"{xt.RT}/{run}/outdata/oifs/atm_remapped_1m_pl_{var}_1m_pl_{year}-{year}.nc"


def setup(run):
    S = xt.setup(run)
    lu = S["land_u"]
    with nc.Dataset(ZS) as d:
        zs = d["Z"][0][S["r0"]:S["r1"], :].astype(np.float64).ravel()[lu] / G0
    S["zs"] = zs
    S["lu"] = lu
    S["wl"] = S["w"][lu]
    # stab mask from ctl1950g winters 1953-1962 (one mask for both runs)
    ps = []
    for y in range(1952, 1963):
        x = {}
        for v in ("msl", "2t"):
            with nc.Dataset(xt.path("ctl1950g", "1m", v, y)) as d:
                x[v] = d[v][:, S["r0"]:S["r1"], :].reshape(12, -1)[:, lu].astype(np.float64)
        mo = np.arange(1, 13)
        sel = (mo >= 9) if y == 1952 else ((mo <= 3) if y == 1962 else np.isin(mo, rg.MONTHS))
        ps.append(ps_est(x["msl"], x["2t"], zs)[sel])
    S["stab"] = np.concatenate(ps).mean(0) >= PS_MIN
    return S


def ps_est(msl, t2, zs):
    return msl * np.exp(-G0 * zs / (RD * (t2 + 0.0065 * zs / 2)))


def rmean(x, mask=None):
    S = _S
    w = S["wl"] if mask is None else S["wl"] * mask
    return np.stack([(x[..., p] * w[p]).sum(-1) / w[p].sum() if w[p].sum() > 0 else np.full(x.shape[:-1], np.nan)
                     for p in S["lpos"]], -1)


def _init(run):
    global _S
    _S = setup(run)


def _year(args):
    run, year = args
    S = _S
    f = {}
    for v in SURF:
        with nc.Dataset(xt.path(run, "1m", v, year)) as d:
            f[v] = d[v][:, S["r0"]:S["r1"], :].reshape(12, -1)[:, S["lu"]].astype(np.float64)
    for v in ("t", "z"):
        with nc.Dataset(pl_path(run, v, year)) as d:
            k = int(np.argmin(np.abs(d["pressure_levels"][:] - 85000.0)))
            f[v + "850"] = d[v][:, k, S["r0"]:S["r1"], :].reshape(12, -1)[:, S["lu"]].astype(np.float64)
    zs, st = S["zs"], S["stab"].astype(float)
    ps = ps_est(f["msl"], f["2t"], zs)
    th850 = f["t850"] * (1e5 / 85000.0) ** KAP
    out = {}
    for v in SURF:
        out[v] = rmean(f[v])
        out[v + "_st"] = rmean(f[v], st)
    out["hfss"] = -rmean(f["sshf"]) / 10800.0
    out["hfls"] = -rmean(f["slhf"]) / 10800.0
    out["hfss_st"] = -rmean(f["sshf"], st) / 10800.0
    out["snc"] = rmean(100 * np.minimum(1.0, (1000.0 * f["sd"] / np.maximum(f["rsn"], 50.0)) / 0.1))
    out["ps_est"] = rmean(ps)
    out["ta850"] = rmean(f["t850"], st)
    out["tas"], out["tas_st"] = out.pop("2t"), out.pop("2t_st")
    out["th850_ths"] = rmean(th850 - f["2t"] * (1e5 / ps) ** KAP, st)
    out["th850_tts"] = rmean(th850 - f["skt"] * (1e5 / ps) ** KAP, st)
    out["inv850_z"] = rmean(f["t850"] - f["2t"] + G0 * (f["z850"] / G0 - zs - 2.0) / CP, st)
    out["dT850s"] = rmean(f["t850"] - f["2t"], st)
    out["z850_above_zs"] = rmean(f["z850"] / G0 - zs, st)
    return year, out


def main():
    run = sys.argv[1]
    t0 = time.time()
    S = setup(run)
    years = list(YEARS[run])
    res = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 8)), initializer=_init, initargs=(run,)) as ex:
        for y, o in ex.map(_year, [(run, y) for y in years]):
            res[y] = o
    keys = list(res[years[0]])
    out = {k: np.stack([res[y][k] for y in years]) for k in keys}
    out["years"] = np.array(years)
    out["regions"] = np.array(S["land_regions"])
    out["n_cells_land"] = np.array([len(p) for p in S["lpos"]])
    out["n_cells_stab"] = np.array([int(S["stab"][p].sum()) for p in S["lpos"]])
    os.makedirs(OUTD, exist_ok=True)
    np.savez(f"{OUTD}/tco95_{run}.npz", **out)
    print(f"{run}: {len(years)} years, {time.time() - t0:.0f}s; stab/land cells " +
          ", ".join(f"{r} {a}/{b}" for r, a, b in zip(S["land_regions"], out["n_cells_stab"], out["n_cells_land"])))


if __name__ == "__main__":
    main()
