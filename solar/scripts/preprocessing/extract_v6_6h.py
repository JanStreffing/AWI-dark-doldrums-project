"""Training rows for v6 (6-hourly, tcc/lcc only) or v7 (6-hourly, all four cloud layers),
Sep-Mar 2015-2026, stride 8 -- same grid thinning, same lat/elevation columns and same true
3-hour-sum ssrd source as the v4 cache, but collapsed onto 6 h windows by
reconstruct_ssrd_v6.to_6h_fields().

Usage: python3 extract_v6_6h.py [--variant v6|v7] [N_MONTHS_LIMIT]
Cache: data/training_cache_<variant>_6h_stride8/<stamp>.npz
       (X float32 in FEATURES_V6 / FEATURES_V7 order, y=kt over 6 h, lat, elev int16)."""

import os as _os, sys as _sys
_root = _os.path.dirname(_os.path.abspath(__file__))
while not _os.path.isdir(_os.path.join(_root, "data")):
    _root = _os.path.dirname(_root)
_os.chdir(_root)
_sys.path.insert(0, _os.path.join(_root, "scripts", "core"))
_sys.path.insert(0, _os.path.join(_root, "scripts", "preprocessing"))

import os, sys
import numpy as np
import xarray as xr
from joblib import Parallel, delayed
from reconstruct_ssrd_v4 import CLOUDS_DIR, MEAN3H_DIR
from reconstruct_ssrd_v6 import CLOUDS_V6, CLOUDS_V7, FEATURES_V6, FEATURES_V7, to_6h_fields, tisr_6h, build_features_6h
from extract_all_v4 import available_stamps, orography

STRIDE = 8
MONTHS = (9, 10, 11, 12, 1, 2, 3)
VARIANTS = dict(v6=(CLOUDS_V6, FEATURES_V6), v7=(CLOUDS_V7, FEATURES_V7))
VARIANT = sys.argv[sys.argv.index("--variant") + 1] if "--variant" in sys.argv else "v6"
CLOUDS, FEATURES = VARIANTS[VARIANT]
CACHE = f"data/training_cache_{VARIANT}_6h_stride8"


def load_3h(stamp, stride):
    sub = dict(latitude=slice(None, None, stride), longitude=slice(None, None, stride))
    clouds = xr.open_dataset(f"{CLOUDS_DIR}/era5_clouds_3h_{stamp}.nc").isel(**sub)
    ssrd = xr.open_dataset(f"{MEAN3H_DIR}/era5_ssrd_3hmean_{stamp}.nc").isel(**sub)
    ct = clouds.valid_time.values.astype("datetime64[s]")
    st = ssrd.valid_time.values.astype("datetime64[s]")
    common = np.intersect1d(ct, st)
    clouds = clouds.isel(valid_time=np.isin(ct, common))
    ssrd = ssrd.isel(time=np.isin(st, common))
    out = dict(times=common, lat=clouds.latitude.values, lon=clouds.longitude.values,
               clouds={c: clouds[c].values.astype("float32") for c in CLOUDS},
               ssrd=ssrd.ssrd.values.astype("float32"))
    clouds.close(); ssrd.close()
    return out


def extract(stamp):
    path = f"{CACHE}/{stamp}.npz"
    if os.path.exists(path):
        return stamp, None
    m = load_3h(stamp, STRIDE)
    w6 = to_6h_fields(m["times"], m["clouds"], m["ssrd"])
    tisr6 = tisr_6h(m["lat"], m["lon"], w6["times"])
    X, day, kt = build_features_6h(w6, tisr6, FEATURES)
    nlat, nlon = len(m["lat"]), len(m["lon"])
    cell = np.nonzero(day)[0] % (nlat * nlon)
    ilat, ilon = cell // nlon, cell % nlon
    elev = orography()[ilat * STRIDE, ilon * STRIDE]
    np.savez(path + ".tmp.npz", X=X, y=kt.astype("float32"), lat=m["lat"][ilat].astype("float32"),
             elev=np.round(elev).astype("int16"))
    os.replace(path + ".tmp.npz", path)
    return stamp, len(kt)


if __name__ == "__main__":
    os.makedirs(CACHE, exist_ok=True)
    stamps = [s for s in available_stamps() if int(s[4:6]) in MONTHS]
    args = [a for a in sys.argv[1:] if a not in ("--variant", VARIANT)]
    limit = int(args[0]) if args else 0
    if limit:
        stamps = stamps[:limit]
    print(f"{VARIANT}: {len(stamps)} Sep-Mar months, features {FEATURES}", flush=True)
    for s, n in Parallel(n_jobs=int(os.environ.get("NJOBS", "16")), backend="loky")(delayed(extract)(s) for s in stamps):
        if n: print(f"  {s}: {n:,} rows", flush=True)
