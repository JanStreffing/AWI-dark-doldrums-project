"""
Global, area-weighted annual-mean 2 m temperature per year for every run used in the multi-region
test, so event ratios can be normalised per K of global warming.

  dart9km    TCo1279-DART 1950C/2080C, 3-hourly reduced-Gaussian 2t. Every 9th step (9 and 8 are
             coprime, so all eight 3-hourly slots are sampled equally); cell weight cos(lat)/n_ring.
  tco319     ctl1950d: 6-hourly 2t (1950-1969); ssp585d: monthly 2t (2015-2099).
  tco95      ctl1950g / hi1950g: monthly 2t, 1950-2099.
  eerie      hist / ctl / ssp245: monthly mean2t (0.25 deg), decoded from GRIB via the byte index.
Annual means from complete years only (12 months).

Usage: python3 global_t2m.py MODEL          (dart9km | tco319 | tco95 | eerie; legacy python)
Output: dunkelflaute/data/regions/global_t2m_<MODEL>.json  {run: {year: K}}
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
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
import regions as rg   # noqa: E402

H = "/work/ab0995/ICCP_AWI_hackthon_2025"
RT = "/work/ab0995/a270092/runtime/awicm3-frontiers-xios"
STRIDE = 9


def _wreg(lat, nlon):
    w = np.cos(np.deg2rad(lat))
    return w[:, None] * np.ones((1, nlon)) / w.sum() / nlon


# ---------------- DART 9 km
_W = None


def _dart_init():
    global _W
    d = np.load(f"{_ROOT}/wind/data/results/tco1279_era5_regrid_index.npz")
    lat = d["tco_lat"].astype(np.float64)
    u, inv, cnt = np.unique(lat, return_inverse=True, return_counts=True)
    w = np.cos(np.deg2rad(lat)) / cnt[inv]
    _W = w / w.sum()


def _dart_month(args):
    run, y, m = args
    f = f"{H}/TCo1279-DART-{run}/outdata/oifs/atm_reduced_3h_2t_3h_{y}{m:02d}-{y}{m:02d}.nc"
    with nc.Dataset(f) as d:
        v = d["2t"]
        n = v.shape[0]
        idx = np.arange((y * 12 + m) % STRIDE, n, STRIDE)
        s = 0.0
        for i in idx:
            s += float(np.dot(v[i].astype(np.float64), _W))
    return run, y, m, s / len(idx)


def dart():
    out = {}
    jobs = [(r, y, m) for r, ys in (("1950C", range(1950, 1970)), ("2080C", range(2080, 2093)))
            for y in ys for m in range(1, 13)]
    mon = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", 32)), initializer=_dart_init) as ex:
        for r, y, m, t in ex.map(_dart_month, jobs):
            mon[(r, y, m)] = t
            print(f"  {r} {y}-{m:02d} {t - 273.15:.2f} C", flush=True)
    for r in ("1950C", "2080C"):
        key = f"dart9km_{r.lower()}"
        out[key] = {}
        for y in sorted({k[1] for k in mon if k[0] == r}):
            ms = [mon[(r, y, m)] for m in range(1, 13) if (r, y, m) in mon]
            if len(ms) == 12:
                out[key][y] = float(np.mean(ms))
    return out


# ---------------- TCo319 31 km
def tco319():
    out = {"tco319_ctl1950d": {}, "tco319_ssp585d": {}}
    for y in range(1950, 1970):
        with nc.Dataset(f"{H}/TCo319-DART-ctl1950d-gibbs-charn/outdata/oifs/atm_remapped_6h_2t_6h_{y}-{y}.nc") as d:
            w = _wreg(d["lat"][:].astype(float), len(d["lon"]))
            v = d["2t"]
            s = [float((v[a:a + 100].astype(np.float64) * w[None]).sum(axis=(1, 2)).sum()) for a in range(0, v.shape[0], 100)]
            out["tco319_ctl1950d"][y] = sum(s) / v.shape[0]
    for f in sorted(glob.glob(f"{H}/TCo319-DART-ssp585d-gibbs-charn/outdata/oifs/atm_remapped_1m_2t_1m_*.nc")):
        y = int(f.split("_")[-1][:4])
        with nc.Dataset(f) as d:
            w = _wreg(d["lat"][:].astype(float), len(d["lon"]))
            x = d["2t"][:].astype(np.float64)
            if x.shape[0] == 12:
                out["tco319_ssp585d"][y] = float((x * w[None]).sum(axis=(1, 2)).mean())
    return out


# ---------------- TCO95 100 km
def tco95():
    out = {}
    for run in ("ctl1950g", "hi1950g"):
        out[f"tco95_{run}"] = {}
        for y in range(1950, 2100):
            f = f"{RT}/{run}/outdata/oifs/atm_remapped_1m_2t_1m_{y}-{y}.nc"
            if not os.path.exists(f):
                continue
            with nc.Dataset(f) as d:
                w = _wreg(d["lat"][:].astype(float), len(d["lon"]))
                x = d["2t"][:].astype(np.float64)
                if x.shape[0] == 12:
                    out[f"tco95_{run}"][y] = float((x * w[None]).sum(axis=(1, 2)).mean())
    return out


# ---------------- EERIE
def eerie():
    import eccodes
    lat = np.linspace(90, -90, 721)
    w = (np.cos(np.deg2rad(lat))[:, None] * np.ones((1, 1440))).ravel()
    w /= w.sum()
    out = {}
    for run in ("hist", "ctl", "ssp245"):
        ix = np.load(f"{rg.BASE_DIR}/eerie_index/{run}_mean2t.npz")
        mon = {}
        for k in range(len(ix["path"])):
            with open(ix["path"][k], "rb") as f:
                f.seek(int(ix["offset"][k]))
                h = eccodes.codes_new_from_message(rg.grib2_strip_local(f.read(int(ix["size"][k]))))
            try:
                d = eccodes.codes_get(h, "dataDate")
                assert eccodes.codes_get(h, "numberOfDataPoints") == 1440 * 721
                mon[(d // 10000, d // 100 % 100)] = float(np.dot(eccodes.codes_get_values(h), w))
            finally:
                eccodes.codes_release(h)
        out[f"eerie_{run}"] = {y: float(np.mean([mon[(y, m)] for m in range(1, 13)]))
                               for y in sorted({k[0] for k in mon}) if all((y, m) in mon for m in range(1, 13))}
        print(f"  eerie {run}: {len(mon)} months, years {min(out[f'eerie_{run}'])}-{max(out[f'eerie_{run}'])}", flush=True)
    return out


def main():
    model = sys.argv[1]
    out = dict(dart9km=dart, tco319=tco319, tco95=tco95, eerie=eerie)[model]()
    fn = f"{rg.BASE_DIR}/global_t2m_{model}.json"
    json.dump({k: {str(y): t for y, t in v.items()} for k, v in out.items()}, open(fn, "w"), indent=1)
    for k, v in out.items():
        ys = sorted(v)
        print(f"{k}: {ys[0]}-{ys[-1]}, first {v[ys[0]] - 273.15:.2f} C, last {v[ys[-1]] - 273.15:.2f} C")
    print(f"Saved {fn}")


if __name__ == "__main__":
    main()
