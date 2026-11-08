"""
Byte index of the EERIE IFS-FESOM2-SR 6-hourly GRIB messages on Levante, from the kerchunk parquet
references (only path, offset and size are used; decoding is done with eccodes in
extract_eerie.py because the parquet stores do not open with zarr 3).

One row per time step, in the store's order, missing chunks kept as empty paths so that row k is
time step k of the store. The time of row k is checked against the decoded message in
extract_eerie.py.

Usage: python3 eerie_index.py [VAR ...]         (module python3/2025.01: needs pandas + pyarrow)
Output: dunkelflaute/data/regions/eerie_index/<run>_<var>.npz
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

K = "/work/bm1344/k202193/Kerchunk/IFS-FESOM_sr"
RUNS = {"hist": "hist-1950/v20240304", "ctl": "control-FESOM", "ssp245": "ssp245_FESOM"}
VARS = {"10u": "2D_6hourly_instant", "10v": "2D_6hourly_instant", "2t": "2D_6hourly_instant",
        "ssrd": "2D_6hourly_accumulated", "tisr": "2D_6hourly_accumulated",
        "mmsl": "2D_daily_avg", "mean2t": "2D_monthly_avg",
        "mean10ws": "2D_daily_avg", "msdwswrf": "2D_daily_avg"}
OUT = "/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/data/regions/eerie_index"


def main():
    os.makedirs(OUT, exist_ok=True)
    only = set(sys.argv[1:])
    for run, rp in RUNS.items():
        for v, store in VARS.items():
            if only and v not in only:
                continue
            fs = sorted(glob.glob(f"{K}/{rp}/atmos/gr025/{store}.parq/{v}/refs.*.parq"),
                        key=lambda x: int(x.split(".")[-2]))
            d = pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True)
            path = d["path"].fillna("").values.astype(str)
            off = d["offset"].fillna(-1).values.astype(np.int64)
            size = d["size"].fillna(-1).values.astype(np.int64)
            last = int(np.where(path != "")[0][-1]) + 1     # refs pages are padded to 100000 rows
            path, off, size = path[:last], off[:last], size[:last]
            n_missing = int((path == "").sum())
            np.savez(f"{OUT}/{run}_{v}.npz", path=path, offset=off, size=size)
            print(f"{run} {v}: {len(path)} rows, {n_missing} missing")


if __name__ == "__main__":
    main()
