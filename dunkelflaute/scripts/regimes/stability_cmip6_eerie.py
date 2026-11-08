"""
The monthly stability diagnostics of stability_cmip6_extract.py for EERIE IFS-FESOM2-SR (IFS 48r1,
0.25 deg output), so EERIE can be placed on the CMIP6 cross-model relations (question C, methods
notes 13d). Comparisons as region_events.COMPARISONS: eerie_ctl (ctl winters 1972-2050 against
ssp245 2041-2070) and eerie_hist (hist 1971-2000 against ssp245 2041-2070).

Data: monthly-mean GRIB messages located through the kerchunk parquet references (as eerie_index.py,
path / offset / size only) and decoded with eccodes after dropping the local section 2
(regions.grib2_strip_local). Every decoded message is checked for its date (and level).
  2D_monthly_avg   mean2t, mskt, msp, msshf (W m-2, ECMWF sign positive downward), msd_m (m w.e.)
  3D_monthly_avg   mt at 850 hPa
Cells: the region land cells of extract_daily.eerie_setup (ERA5 land fraction >= 0.5, cos-lat).
Derived as in stability_cmip6_extract.py: th850_ths = theta(mt850) - theta(mean2t, msp), th850_tts
with mskt, dT850s = mt850 - mean2t, hfss = -msshf (CMIP sign, positive upward). ta850 quantities on
cells whose control winter-mean msp >= 90000 Pa. Snow cover snc = 100 min(1, sd / 0.025 m w.e.)
(IFS 0.1 m depth threshold at an assumed 250 kg m-3 snow density: no monthly snow density stored).

Usage: python3 stability_cmip6_eerie.py       (legacy python: eccodes, pandas; N_WORKERS, default 32)
Output: dunkelflaute/data/regimes/stability_cmip6_parts/eerie_summary.json (stability_cmip6_analysis.py
        reads it) and eerie_<run>.npz (winter-month region means)
"""
import glob
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import pandas as pd

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("", "/regions", "/cmip6", "/regimes"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts{p}")
import regions as rg            # noqa: E402
import extract_daily as xd      # noqa: E402

K = "/work/bm1344/k202193/Kerchunk/IFS-FESOM_sr"
RUNS = {"hist": ("hist-1950/v20240304", 1950), "ctl": ("control-FESOM", 1950), "ssp245": ("ssp245_FESOM", 2015)}
V2 = ("mean2t", "mskt", "msp", "msshf", "msd_m")
NLEV = 19
COMPS = {"eerie_ctl": (("ctl", 1972, 2050), ("ssp245", 2041, 2070)),
         "eerie_hist": (("hist", 1971, 2000), ("ssp245", 2041, 2070))}
MON = [9, 10, 11, 12, 1, 2, 3]
OUTD = f"{_ROOT}/dunkelflaute/data/regimes/stability_cmip6_parts"
KAP, PS_MIN = 0.2857, 90000.0
EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
_S = None


def refs(run, store, var):
    fs = sorted(glob.glob(f"{K}/{RUNS[run][0]}/atmos/gr025/{store}.parq/{var}/refs.*.parq"), key=lambda x: int(x.split(".")[-2]))
    d = pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True)
    return d["path"].fillna("").values.astype(str), d["offset"].fillna(-1).values.astype(np.int64), d["size"].fillna(-1).values.astype(np.int64)


def months_needed():
    need = {}
    for (a, c0, c1), (b, f0, f1) in COMPS.values():
        for run, y0, y1 in ((a, c0, c1), (b, f0, f1)):
            for w in range(y0, y1 + 1):
                for m in MON:
                    need.setdefault(run, set()).add((w - 1 if m >= 9 else w, m))
    return {r: sorted(v) for r, v in need.items()}


def _init():
    global _S
    _S = xd.eerie_setup()


def _decode(job):
    import eccodes
    run, var, ym, path, off, size, lev = job
    with open(path, "rb") as fh:
        fh.seek(int(off))
        h = eccodes.codes_new_from_message(rg.grib2_strip_local(fh.read(int(size))))
    try:
        d = eccodes.codes_get(h, "dataDate")
        lv = eccodes.codes_get(h, "level") if lev else None
        v = eccodes.codes_get_values(h)
    finally:
        eccodes.codes_release(h)
    assert (d // 10000, d // 100 % 100) == tuple(ym), f"{run} {var} {ym}: message date {d}"
    if lev:
        assert lv == 850, f"{run} {var} {ym}: level {lv}"
    return run, var, tuple(ym), v[_S["glob_idx"]]


def main():
    t0 = time.time()
    S = xd.eerie_setup()
    need = months_needed()
    jobs = []
    for run, yms in need.items():
        y_start = RUNS[run][1]
        for var in V2:
            p, o, s = refs(run, "2D_monthly_avg", var)
            for y, m in yms:
                k = (y - y_start) * 12 + (m - 1)
                jobs.append((run, var, (y, m), p[k], o[k], s[k], False))
        p, o, s = refs(run, "3D_monthly_avg", "mt")
        # level index of 850 hPa from the first time step
        _init()
        lk = None
        for li in range(NLEV):
            import eccodes
            with open(p[li], "rb") as fh:
                fh.seek(int(o[li]))
                h = eccodes.codes_new_from_message(rg.grib2_strip_local(fh.read(int(s[li]))))
            lv = eccodes.codes_get(h, "level")
            eccodes.codes_release(h)
            if lv == 850:
                lk = li
                break
        assert lk is not None
        for y, m in yms:
            k = ((y - y_start) * 12 + (m - 1)) * NLEV + lk
            jobs.append((run, "mt850", (y, m), p[k], o[k], s[k], True))
    print(f"{len(jobs)} messages, setup {time.time() - t0:.0f}s", flush=True)
    got = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 32)), initializer=_init) as ex:
        for i, (run, var, ym, v) in enumerate(ex.map(_decode, jobs, chunksize=8)):
            got[(run, var, ym)] = v
            if i % 1000 == 0:
                print(f"  {i} [{time.time() - t0:.0f}s]", flush=True)
    w = S["w"]
    lpos = S["lpos"]
    R = rg.REGIONS

    def rmean(x, mask=None):
        ww = w if mask is None else w * mask
        return np.array([(x[p] * ww[p]).sum() / ww[p].sum() if ww[p].sum() > 0 else np.nan for p in lpos])

    # stab mask: control winters of eerie_ctl (ctl) and eerie_hist (hist), mean msp >= PS_MIN, one per comparison
    summ = {}
    for comp, ((cn, c0, c1), (fn, f0, f1)) in COMPS.items():
        def yms(y0, y1):
            return [(w_ - 1 if m >= 9 else w_, m) for w_ in range(y0, y1 + 1) for m in MON]
        psc = np.mean([got[(cn, "msp", ym)] for ym in yms(c0, c1)], axis=0)
        stab = (psc >= PS_MIN).astype(float)
        per = {}
        for tag, run, y0, y1 in (("ctl", cn, c0, c1), ("fut", fn, f0, f1)):
            acc = {}
            for ym in yms(y0, y1):
                g = {v: got[(run, v, ym)] for v in V2 + ("mt850",)}
                th850 = g["mt850"] * (1e5 / 85000.0) ** KAP
                q = dict(tas=rmean(g["mean2t"]), tas_st=rmean(g["mean2t"], stab), ts=rmean(g["mskt"]), ps=rmean(g["msp"]),
                         hfss=-rmean(g["msshf"]), hfss_st=-rmean(g["msshf"], stab), sd=rmean(g["msd_m"]),
                         snc=rmean(100 * np.minimum(1.0, g["msd_m"] / 0.025)), ta850=rmean(g["mt850"], stab),
                         th850_ths=rmean(th850 - g["mean2t"] * (1e5 / g["msp"]) ** KAP, stab),
                         th850_tts=rmean(th850 - g["mskt"] * (1e5 / g["msp"]) ** KAP, stab),
                         dT850s=rmean(g["mt850"] - g["mean2t"], stab))
                for k, v in q.items():
                    acc.setdefault(k, []).append(v)
            per[tag] = {k: np.mean(v, axis=0) for k, v in acc.items()}
        np.savez(f"{OUTD}/eerie_{comp}.npz", regions=np.array(R), n_cells_stab=np.array([int(stab[p].sum()) for p in lpos]),
                 n_cells_land=np.array([len(p) for p in lpos]), **{f"{t}|{k}": v for t, d in per.items() for k, v in d.items()})
        out = {}
        C22 = [r for r in R if dict(rg.BOXES, DE=rg.DE_BOX)[r][0] >= 29.0]
        for nm, regs in (("eu8", EU), ("cmip22", C22)):
            ii = [R.index(r) for r in regs]
            o = {}
            for k in per["ctl"]:
                h, f = per["ctl"][k], per["fut"][k]
                o[f"d_{k}"] = float(np.nanmean((f - h)[ii]))
                o[f"ctl_{k}"] = float(np.nanmean(h[ii]))
                if k in ("ps", "sd"):
                    o[f"d_{k}_pct"] = float(np.nanmean((100 * (f / h - 1))[ii]))
            o["ratio_dtas_dta850"] = o["d_tas_st"] / o["d_ta850"]
            out[nm] = o
        summ[comp] = out
        e = out["eu8"]
        print(f"{comp}: eu8 d th850-ths {e['d_th850_ths']:+.2f} K, dT850s {e['d_dT850s']:+.2f}, hfss {e['d_hfss']:+.2f}, "
              f"snc {e['d_snc']:+.1f}, tas {e['d_tas']:+.2f}, ta850 {e['d_ta850']:+.2f}", flush=True)
    json.dump(summ, open(f"{OUTD}/eerie_summary.json", "w"), indent=1)
    print(f"done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
