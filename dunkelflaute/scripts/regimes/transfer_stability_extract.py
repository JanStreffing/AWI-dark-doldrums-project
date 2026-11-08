"""
Per region and day diagnostics of the surface-wind transfer in the TCO95 (~100 km) reruns
ctl1950g / hi1950g, Sep-Mar, for transfer_stability.py.

For each (run, calendar year) the Sep-Mar days of that year (Jan-Mar and Sep-Dec), region means
over the land cells of the 8 European regions (regime_analysis.EU; cells, cos-lat weights and the
ERA5 z0 log law exactly as extract_tco95.py: model lsm >= 0.5 from the 1952 monthly lsm), daily
means of the eight 3 h stamps 03..24 UTC (a stamp belongs to the day of stamp - 1 s):

  s10, s100n, s100ll      10 m speed, native 100 m speed, log-law 100 m speed (ERA5 z0)
  V[k]                    model-level speed, levels 85..91 (k = 0..6; 91 lowest)
  T[k]                    model-level temperature
  h[k]                    (z_k - z_91) / g, height above the lowest level (m), k = 0..5
  dth[k]                  theta_k - theta_91 ~ T_k - T_91 + g (z_k - z_91) / cp (dry, no q)
  sh2[k]                  |V_k - V_91|^2 (vector difference), k = 0..5
  ri_st[k]                fraction of cell-steps with bulk Ri(91 -> k) > 0.25
  ri_med[k]               region-day median of the per cell-step bulk Ri(91 -> k)
  t2, skt, blh            2 m temperature, skin temperature, boundary-layer height
  sshf                    3 h sshf (J m-2, accumulated over the 3 h ending at the stamp) / 10800 s,
                          W m-2, IFS sign: positive downward (into the surface) = stable
  sshf_pos                fraction of cell-steps with sshf > 0
  mslbox                  daily-mean msl from the 1d file on the 2.5 deg box grid of extract_msl.py
Monthly (region-month land means, all 12 months): sd, rsn, asn, tsn, fsr, blh, sshf, skt, 2t, stl1,
plus snow cover fractions min(1, sd / 0.015 m) and min(1, (1000 sd / rsn) / 0.1 m).

Usage: python3 transfer_stability_extract.py RUN Y0-Y1[,Y0-Y1...]   (N_WORKERS, default 12)
Output: dunkelflaute/data/regimes/transfer_stability_parts/<run>_<year>.npz
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
import extract_msl as xm        # noqa: E402

EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
OUTD = f"{_ROOT}/dunkelflaute/data/regimes/transfer_stability_parts"
G, CP = 9.80665, 1004.7
NL = 7                           # levels 85..91
SURF = ("10u", "10v", "100u", "100v", "2t", "skt", "blh", "sshf")
ML = ("u", "v", "t", "z")
MON = ("sd", "rsn", "asn", "tsn", "fsr", "blh", "sshf", "skt", "2t", "stl1")
_S = None


def setup(run):
    S = xt.setup(run)
    eu_cells = [S["cells"][r] for r in EU]
    union = np.unique(np.concatenate(eu_cells))
    pos = [np.searchsorted(union, c) for c in eu_cells]
    z0pos = np.searchsorted(S["land_u"], union)
    assert np.array_equal(S["land_u"][z0pos], union)
    w = S["w"][union]
    return dict(r0=S["r0"], r1=S["r1"], union=union, pos=pos, w=w, z0=S["z0"][:, z0pos],
                n_cells=[len(c) for c in eu_cells])


def _init(run):
    global _S
    _S = setup(run)


def rmean(x):
    """x (T, N) over union cells -> (T, R) weighted region means."""
    w = _S["w"]
    return np.stack([(x[:, p] * w[p]).sum(axis=1) / w[p].sum() for p in _S["pos"]], axis=1)


def rmedian(x):
    return np.stack([np.median(x[:, p], axis=1) for p in _S["pos"]], axis=1)


def _times(d):
    tv = d["time_instant"] if "time_instant" in d.variables else d["time_centered"]
    return np.array([np.datetime64(str(x)[:19]) for x in nc.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)])


def _read(d, v, i0, i1, ml):
    S = _S
    if ml:
        x = d[v][i0:i1, :, S["r0"]:S["r1"], :]
        return x.reshape(x.shape[0], x.shape[1], -1)[:, :, S["union"]].astype(np.float64)
    x = d[v][i0:i1, S["r0"]:S["r1"], :]
    return x.reshape(x.shape[0], -1)[:, S["union"]].astype(np.float64)


def _block(run, year, i0, i1, months):
    """Daily region diagnostics for the steps i0:i1 (whole days)."""
    S = _S
    f = {}
    for v in SURF:
        with nc.Dataset(xt.path(run, "3h", v, year)) as d:
            d.set_auto_mask(False)
            f[v] = _read(d, v, i0, i1, False)
    for v in ML:
        with nc.Dataset(f"{xt.RT}/{run}/outdata/oifs/atm_remapped_3h_ml_{v}_3h_ml_{year}-{year}.nc") as d:
            d.set_auto_mask(False)
            f[v] = _read(d, v, i0, i1, True)
    nt = i1 - i0
    out = {}
    s10 = np.hypot(f["10u"], f["10v"])
    z0 = S["z0"][months - 1]                            # (T, N)
    out["s10"] = rmean(s10)
    out["s100n"] = rmean(np.hypot(f["100u"], f["100v"]))
    out["s100ll"] = rmean(rg.loglaw_land(s10, z0))
    for v in ("2t", "skt", "blh"):
        out[v.replace("2t", "t2")] = rmean(f[v])
    sh = f["sshf"] / 10800.0
    out["sshf"] = rmean(sh)
    out["sshf_pos"] = rmean((sh > 0).astype(float))
    V = np.hypot(f["u"], f["v"])                        # (T, 7, N)
    out["V"] = np.stack([rmean(V[:, k]) for k in range(NL)], axis=1)
    out["T"] = np.stack([rmean(f["t"][:, k]) for k in range(NL)], axis=1)
    dz = (f["z"][:, :NL - 1] - f["z"][:, NL - 1:]) / G  # (T, 6, N)
    dth = f["t"][:, :NL - 1] - f["t"][:, NL - 1:] + G * dz / CP
    sh2 = (f["u"][:, :NL - 1] - f["u"][:, NL - 1:]) ** 2 + (f["v"][:, :NL - 1] - f["v"][:, NL - 1:]) ** 2
    ri = G / f["t"][:, NL - 1:] * dth * dz / np.maximum(sh2, 1e-2)
    out["h"] = np.stack([rmean(dz[:, k]) for k in range(NL - 1)], axis=1)
    out["dth"] = np.stack([rmean(dth[:, k]) for k in range(NL - 1)], axis=1)
    out["sh2"] = np.stack([rmean(sh2[:, k]) for k in range(NL - 1)], axis=1)
    out["ri_st"] = np.stack([rmean((ri[:, k] > 0.25).astype(float)) for k in range(NL - 1)], axis=1)
    out["ri_med"] = np.stack([rmedian(ri[:, k]) for k in range(NL - 1)], axis=1)
    nd = nt // 8
    return {k: v.reshape(nd, 8, *v.shape[1:]).mean(axis=1) for k, v in out.items()}


def _year(args):
    run, year = args
    t0 = time.time()
    fn = f"{OUTD}/{run}_{year}.npz"
    if os.path.exists(fn):
        return year, 0.0, "exists"
    with nc.Dataset(xt.path(run, "3h", "10u", year)) as d:
        times = _times(d)
    day = (times - np.timedelta64(1, "s")).astype("datetime64[D]")
    mon = (day.astype("datetime64[M]").astype(int) % 12) + 1
    assert len(times) % 8 == 0 and (day[::8] == day[7::8]).all(), f"{run} {year}: stamps not 8 per day"
    res = []
    # month blocks keep memory per worker at a few 100 MB
    for m in sorted(rg.MONTHS):
        idx = np.where(mon == m)[0]
        if len(idx) == 0:
            continue
        i0, i1 = int(idx[0]), int(idx[-1]) + 1
        assert i1 - i0 == len(idx) and len(idx) % 8 == 0
        b = _block(run, year, i0, i1, mon[i0:i1])
        b["day"] = day[i0:i1:8]
        res.append(b)
    out = {k: np.concatenate([b[k] for b in res]) for k in res[0]}
    # daily msl, 2.5 deg boxes (as extract_msl.py: box mean of the 1d daily-average file)
    fm = f"{xt.RT}/{run}/outdata/oifs/atm_remapped_1d_msl_1d_{year}-{year}.nc"
    gs = xm._reg_setup(fm)
    with nc.Dataset(fm) as d:
        tm = _times(d)
        dm = tm.astype("datetime64[D]")
        ii = np.where(np.isin(dm, out["day"]))[0]
        x = d["msl"][ii, gs["r0"]:gs["r1"], :].reshape(len(ii), -1)[:, gs["sel"]].astype(np.float64)
    out["mslbox"] = np.stack([xm.box_mean(x[k], gs["bid"], gs["cnt"]) for k in range(len(ii))])
    out["mslday"] = dm[ii]
    # monthly land means
    for v in MON:
        with nc.Dataset(xt.path(run, "1m", v, year)) as d:
            x = d[v][:, _S["r0"]:_S["r1"], :]
            x = x.reshape(x.shape[0], -1)[:, _S["union"]].astype(np.float64)
        out[f"m_{v}"] = rmean(x)
        if v == "sd":
            sd = x
        if v == "rsn":
            rsn = x
    out["m_scf_swe"] = rmean(np.minimum(1.0, sd / 0.015))
    out["m_scf_depth"] = rmean(np.minimum(1.0, (1000.0 * sd / np.maximum(rsn, 50.0)) / 0.1))
    os.makedirs(OUTD, exist_ok=True)
    np.savez(fn + ".tmp.npz", **out)
    os.replace(fn + ".tmp.npz", fn)
    return year, time.time() - t0, "ok"


def main():
    run = sys.argv[1]
    years = []
    for part in sys.argv[2].split(","):
        a, b = part.split("-")
        years += list(range(int(a), int(b) + 1))
    S = setup(run)
    print(f"{run}: {len(years)} years, EU cells " + ", ".join(f"{r} {n}" for r, n in zip(EU, S["n_cells"])), flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 12)), initializer=_init, initargs=(run,)) as ex:
        for y, dt, st in ex.map(_year, [(run, y) for y in years]):
            print(f"  {run} {y} {st} {dt:.0f}s [{time.time() - t0:.0f}s]", flush=True)
    print(f"{run} done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
