"""PRODUCTION reconstruction: v5-area model (trained on genuine 3-hour-sum ERA5
ssrd/tsr, era5_3h_mean -- NOT the earlier 1h-accum-sampled-every-3h data v3 used)
+ the exact constrained monthly rescale (rescale.py) against DART's own real
monthly ssrd, for one DART year's Jan/Feb/Mar/Sep/Oct/Nov/Dec (the "Sep-Mar
extended winter" months within that calendar year).
Writes full 3-hourly ssrd netCDF files, one per month, at native resolution.

The exact-sum constraint (memory rule: rescaling is part of the model) is verified
and logged for every month.

Usage: python3 dart_reconstruct_year.py <run_dir> <year> [--cell-limit N]
Output: data/dart_reconstructed_<1950c|2080c>/<year>/ssrd_reduced_3h_<stamp>-<stamp>.nc
"""
import os as _os, sys as _sys
_root = _os.path.dirname(_os.path.abspath(__file__))
while not _os.path.isdir(_os.path.join(_root, "data")):
    _root = _os.path.dirname(_root)
_os.chdir(_root)
_sys.path.insert(0, _os.path.join(_root, "scripts", "core"))

import argparse
import os
import time

import joblib
import netCDF4 as nc
import numpy as np
import xarray as xr

from reconstruct_ssrd_dart import (DAY_MU, DT_DART, MONTHLY_ANCHOR_DIVISOR, TEMPLATE_DIR,
                                   is_leap, log, nearest_regular_index)
from rescale import constrained_rescale

MODEL_PATH = "models/model_v5/model_kt_v5.joblib"
FEATURES = ["T", "tcc", "hcc", "mcc", "lcc", "mu"]
MONTHS = [1, 2, 3, 9, 10, 11, 12]     # Sep-Mar extended winter within THIS calendar year
KMAX = 1.1
TB = 8
# Same exclusion the model was TRAINED under (train_v5_model.py: LAT_MIN, ELEV_MAX =
# -60.0, 3000) -- predictions outside this mask are extrapolation, never validated.
# Not relevant to this project's own domains (Germany, South Korea both sit well
# inside it), but writing a real-looking number there anyway is misleading for
# anyone who doesn't already know to discount those cells. Written as NaN instead.
LAT_MIN, ELEV_MAX = -60.0, 3000.0
REGRID_INDEX = "/work/ab0246/a270092/software/dark-doldrums-project/wind/data/results/tco1279_era5_regrid_index.npz"
ERA5_OROGRAPHY = "data/static/era5_geopotential_surface.nc"


def dart_elevation_mask(lat, sl):
    """True where this DART cell is INSIDE the trained region (lat>=LAT_MIN and
    elevation<=ELEV_MAX). Elevation isn't a DART field -- brought in via the
    wind team's precomputed ERA5->TCo1279 nearest-neighbour regrid index (verified
    bit-for-bit identical cell ordering to this script's own lat/lon, same grid).
    sl is the same cell slice main() applies everywhere (identity unless
    --cell-limit is set, e.g. for a quick test run on a cell subset)."""
    idx = np.load(REGRID_INDEX, allow_pickle=True)
    i, j, tco_lat = idx["i"][sl], idx["j"][sl], idx["tco_lat"][sl]
    assert len(tco_lat) == len(lat) and np.max(np.abs(tco_lat - lat)) == 0, \
        "regrid index cell ordering does not match this run's lat array"
    z = xr.open_dataset(ERA5_OROGRAPHY)["z"].squeeze().values.astype("float32") / 9.80665
    elev = z[i, j]
    return (lat >= LAT_MIN) & (elev <= ELEV_MAX)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir"); ap.add_argument("year", type=int)
    ap.add_argument("--cell-limit", type=int, default=0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    y = args.year
    d = f"{args.run_dir}/outdata/oifs"
    run_tag = "2080c" if "2080" in os.path.basename(args.run_dir.rstrip("/")).lower() else "1950c"
    out_dir = args.out or f"data/dart_reconstructed_{run_tag}/{y}"
    os.makedirs(out_dir, exist_ok=True)
    N = args.cell_limit or None
    sl = slice(0, N)

    def month_ok(mth):
        need = [f"{d}/atm_reduced_3h_{v}_3h_{y}{mth:02d}-{y}{mth:02d}.nc" for v in ("tcc", "hcc", "mcc", "lcc", "tsr")]
        need.append(f"{d}/atm_remapped_1m_ssrd_1m_{y}{mth:02d}-{y}{mth:02d}.nc")
        miss = [f for f in need if not os.path.exists(f)]
        if miss:
            log(f"  SKIPPING {y}{mth:02d}: missing source file(s) {[os.path.basename(m) for m in miss]}")
            return False
        return True

    months_avail = [m for m in MONTHS if month_ok(m)]
    if not months_avail:
        log(f"year {y}: no months available (all 5 missing source data), nothing to do")
        return

    model = joblib.load(MODEL_PATH)["model"]
    first = xr.open_dataset(f"{d}/atm_reduced_3h_tcc_3h_{y}{months_avail[0]:02d}-{y}{months_avail[0]:02d}.nc")
    lat = first.lat.values[sl].astype("float32"); lon = first.lon.values[sl].astype("float32"); first.close()
    ncell = lat.size
    log(f"year {y}: {ncell:,} cells, months {months_avail}" +
        (f" (skipped {sorted(set(MONTHS) - set(months_avail))}, missing source data)" if len(months_avail) < len(MONTHS) else ""))

    valid_mask = dart_elevation_mask(lat, sl)
    log(f"year {y}: {100*(1-valid_mask.mean()):.2f}% of cells outside the trained region "
        f"(lat<{LAT_MIN} or elevation>{ELEV_MAX:.0f}m) -- written as NaN, not extrapolated")

    # ---- tisr template, read ONCE for the whole year (all 5 months share it) ----
    year_start = np.datetime64(f"{y}-01-01T03:00:00", "s")
    tpl = f"{TEMPLATE_DIR}/tisr_template_dart_{'366' if is_leap(y) else '365'}day.nc"
    ds = xr.open_dataset(tpl)
    nt_year = ds.sizes["time_counter"]
    idx_by_month = {}
    for mth in months_avail:
        tt = xr.open_dataset(f"{d}/atm_reduced_3h_tcc_3h_{y}{mth:02d}-{y}{mth:02d}.nc").time_counter.values
        ix = np.round((tt.astype("datetime64[s]") - year_start) / np.timedelta64(3, "h")).astype(int)
        assert ix.min() >= 0 and ix.max() < nt_year, f"{mth}: idx out of range"
        idx_by_month[mth] = ix
    all_idx = np.unique(np.concatenate(list(idx_by_month.values())))
    tisr_year = np.empty((len(all_idx), ncell), "float32")
    t0 = time.time()
    for s_ in range(0, ncell, 100_000):
        e_ = min(s_ + 100_000, ncell)
        tisr_year[:, s_:e_] = ds.tisr.isel(cell=slice(s_, e_)).values[all_idx]
    ds.close()
    log(f"template read in {(time.time()-t0)/60:.1f} min")

    worst_rel_err = 0.0
    for mth in months_avail:
        tm = time.time()
        stamp = f"{y}{mth:02d}"; tag = f"{stamp}-{stamp}"
        dsv = {v: xr.open_dataset(f"{d}/atm_reduced_3h_{v}_3h_{tag}.nc") for v in ("tcc", "hcc", "mcc", "lcc", "tsr")}
        monthly = xr.open_dataset(f"{d}/atm_remapped_1m_ssrd_1m_{tag}.nc")
        times = dsv["tcc"].time_counter.values
        ntime = len(times)
        rows_t = np.searchsorted(all_idx, idx_by_month[mth])
        tisr_m = tisr_year[rows_t]

        ssrd_pred = np.zeros((ntime, ncell), "float64")
        for t0_ in range(0, ntime, TB):
            t1_ = min(t0_ + TB, ntime)
            f = {v: dsv[v][v].isel(time_counter=slice(t0_, t1_), cell=sl).values.astype("float32") for v in dsv}
            tisr = tisr_m[t0_:t1_]
            mu = tisr / (DT_DART * 1361.0)
            day = mu > DAY_MU
            if day.any():
                X = np.column_stack([np.clip(f["tsr"][day] / tisr[day], 0, 1.5), f["tcc"][day], f["hcc"][day],
                                     f["mcc"][day], f["lcc"][day], mu[day]])
                p = np.clip(model.predict(X), 0, 1.5) * tisr[day]
                full = np.zeros(day.shape, "float64"); full[day] = p
                ssrd_pred[t0_:t1_] = full

        # ---- DART's own monthly ssrd -> native cells (nearest cell on the remapped grid) ----
        mlat, mlon = monthly.lat.values, monthly.lon.values
        ilat, ilon = nearest_regular_index(lat, lon, mlat, mlon)
        actual_wm2 = (monthly.ssrd.squeeze().values / MONTHLY_ANCHOR_DIVISOR)[ilat, ilon].astype("float64")
        target = actual_wm2 * ntime * DT_DART     # J/m2 monthly sum DART implies

        # ---- exact constrained rescale (rescale.py) -- monthly sum is a HARD constraint ----
        tisr_full = tisr_m.astype("float64")
        q, diag = constrained_rescale(ssrd_pred, tisr_full, target, kmax=KMAX)
        log(f"  {stamp} rescale check: max_rel_sum_err={diag['max_rel_sum_err']:.2e} "
            f"ceiling_raised={diag['n_ceiling_raised']} twilight={diag['n_twilight']} "
            f"max_kt={diag['max_kt_normal']:.3f} night_max={diag['night_max_normal']:.1e}")
        worst_rel_err = max(worst_rel_err, diag["max_rel_sum_err"])

        # Mask applied to the OUTPUT only, after the rescale's own exact-sum check has
        # already run over the full array -- doesn't touch the rescale numerics for
        # valid cells, just blanks the never-validated ones before anything is written.
        q_out = np.where(valid_mask[None, :], q, np.nan)

        out_path = f"{out_dir}/ssrd_reduced_3h_{tag}.nc"
        wds = nc.Dataset(out_path, "w", format="NETCDF4")
        wds.createDimension("time_counter", ntime); wds.createDimension("cell", ncell)
        secs = ((times.astype("datetime64[s]") - np.datetime64(f"{y}-01-01T00:00:00", "s")) /
                np.timedelta64(1, "s")).astype("int64")
        tv = wds.createVariable("time_counter", "i8", ("time_counter",)); tv[:] = secs
        tv.units = f"seconds since {y}-01-01 00:00:00"; tv.calendar = "proleptic_gregorian"
        wds.createVariable("lat", "f4", ("cell",))[:] = lat
        wds.createVariable("lon", "f4", ("cell",))[:] = lon
        v = wds.createVariable("ssrd", "f4", ("time_counter", "cell"), zlib=True, complevel=4,
                               chunksizes=(ntime, min(100_000, ncell)), fill_value=np.nan)
        v[:] = q_out.astype("float32")
        v.units = "J m-2"
        v.long_name = "Reconstructed surface solar radiation downwards, 3-hour accumulation"
        v.comment = ("v5-area-weighted ML model (kt=ssrd/tisr from cloud fractions, tsr, solar geometry) "
                     "+ exact constrained monthly rescale to DART's own monthly ssrd (data/dart_analytical_tisr "
                     "templates, rescale.py). Monthly sum matches DART's real value exactly for valid cells "
                     f"(see log). kt capped at {KMAX}; night is exactly zero. "
                     "NaN south of 60S or above 3000m elevation -- outside the region the model was "
                     "trained and validated on; see solar/README.md Key Finding 5.")
        wds.createVariable("valid_training_region", "i1", ("cell",))[:] = valid_mask.astype("i1")
        wds.model = MODEL_PATH; wds.features = ", ".join(FEATURES)
        wds.close()
        for ds_ in list(dsv.values()) + [monthly]:
            ds_.close()
        log(f"  {stamp}: wrote {out_path}  ({(time.time()-tm)/60:.1f} min)")

    log(f"year {y} DONE. worst relative monthly-sum error across all months: {worst_rel_err:.2e}")


if __name__ == "__main__":
    main()
