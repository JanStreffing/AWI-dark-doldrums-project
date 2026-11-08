"""
Log-law reconstruction of 100 m wind for one month of a TCo1279-DART run.
Supports both 1950C (1950-1969) and 2080C (2080-2092) -- same grid, same method,
only the source directory / output directory / land-z0 cycling base year differ.

Streams through the month in small time-chunks (default 40 steps ~ 5 days) rather
than loading the whole month at once -- the first version of this script loaded
everything into memory and was OOM-killed by the login node's per-user memory
cgroup at ~50 GB resident. Chunked, peak memory is roughly chunk_size * n_cells *
4 bytes * (a handful of arrays) ~ a few GB, safe on any node.

z0 sources:
  land:   ERA5 monthly fsr climatology (2016-2025), year-cycled --
          DART year Y uses ERA5 year 2016 + (Y - CYCLE_BASE) % 10, same calendar month.
          CYCLE_BASE is the run's first year (1950 for 1950C, 2080 for 2080C).
          (precomputed once in tco1279_land_z0_2016_2025.npz, reused for both runs)
  ocean:  Charnock relation, solved iteratively from THIS month's own TCo1279 u10/v10.
          z0 = alpha * u*^2/g + 0.11*nu/u*,  u* = kappa*speed10/ln(10/z0)

Land/sea split and the ERA5->TCo1279 nearest-neighbour regrid index are both
precomputed once in tco1279_era5_regrid_index.npz (TCo1279's grid is identical
across DART runs at this resolution -- verified 1950C vs 2080C bit-for-bit).

Output matches the original DART files' own convention: dims (time_counter, cell),
time_counter carrying the real timestamps, lat/lon as (cell,) coords, float32.

Usage: python reconstruct_month.py YEAR MONTH [--run 1950c|2080c] [--chunk N] [--pilot-steps N]
  --run NAME       : which DART run (default 1950c)
  --chunk N        : timesteps per chunk (default 40, ~5 days at 3-hourly)
  --pilot-steps N  : only process the first N timesteps (quick correctness check),
                     still writes a (small) output file so the format can be checked
"""

import os
import sys
import time
import numpy as np
import netCDF4 as nc

YEAR = int(sys.argv[1])
MONTH = int(sys.argv[2])
RUN = sys.argv[sys.argv.index("--run") + 1] if "--run" in sys.argv else "1950c"
CHUNK = int(sys.argv[sys.argv.index("--chunk") + 1]) if "--chunk" in sys.argv else 40
PILOT_N = int(sys.argv[sys.argv.index("--pilot-steps") + 1]) if "--pilot-steps" in sys.argv else None

RUNS = {
    "1950c": dict(dart_name="TCo1279-DART-1950C", cycle_base=1950, out_dir="data/dart_1950c_100m"),
    "2080c": dict(dart_name="TCo1279-DART-2080C", cycle_base=2080, out_dir="data/dart_2080c_100m"),
}
cfg = RUNS[RUN]
CYCLE_BASE = cfg["cycle_base"]

KAPPA, G, NU, ALPHA_CH = 0.4, 9.81, 1.5e-5, 0.018
DART_DIR = f"/work/ab0995/ICCP_AWI_hackthon_2025/{cfg['dart_name']}/outdata/oifs"
OUT_DIR = cfg["out_dir"]
os.makedirs(OUT_DIR, exist_ok=True)

t0 = time.time()

idx = np.load("data/results/tco1279_era5_regrid_index.npz")
land = idx["land"]  # (n_cells,) bool
n_cells = land.size
lz = np.load("data/results/tco1279_land_z0_2016_2025.npz")
land_z0_table = lz["land_z0"]  # (10 era5-years, 12 months, n_land_cells)

era5_year = 2016 + (YEAR - CYCLE_BASE) % 10
land_z0_month = land_z0_table[era5_year - 2016, MONTH - 1].astype(np.float32)  # (n_land_cells,)
print(f"DART {YEAR}-{MONTH:02d} -> land z0 from ERA5 {era5_year}-{MONTH:02d}", flush=True)

tag = f"{YEAR}{MONTH:02d}-{YEAR}{MONTH:02d}"
u10_in = nc.Dataset(f"{DART_DIR}/atm_reduced_3h_10u_3h_{tag}.nc")
v10_in = nc.Dataset(f"{DART_DIR}/atm_reduced_3h_10v_3h_{tag}.nc")
nt_full = u10_in.dimensions["time_counter"].size
nt = min(PILOT_N, nt_full) if PILOT_N else nt_full
time_counter_vals = u10_in.variables["time_counter"][:nt]
assert u10_in.dimensions["cell"].size == n_cells
print(f"month has {nt_full} timesteps, processing {nt}  ({time.time()-t0:.1f}s)", flush=True)

outfile = f"{OUT_DIR}/u100v100_{YEAR}{MONTH:02d}.nc"
out = nc.Dataset(outfile, "w")
out.createDimension("time_counter", nt)
out.createDimension("cell", n_cells)
tvar = out.createVariable("time_counter", "f8", ("time_counter",))
tvar.units = u10_in.variables["time_counter"].units if hasattr(u10_in.variables["time_counter"], "units") else ""
tvar[:] = time_counter_vals
latvar = out.createVariable("lat", "f4", ("cell",))
lonvar = out.createVariable("lon", "f4", ("cell",))
latvar[:] = idx["tco_lat"]
lonvar[:] = idx["tco_lon"]
u100var = out.createVariable("u100", "f4", ("time_counter", "cell"), zlib=True, complevel=1,
                              chunksizes=(min(CHUNK, nt), n_cells))
v100var = out.createVariable("v100", "f4", ("time_counter", "cell"), zlib=True, complevel=1,
                              chunksizes=(min(CHUNK, nt), n_cells))
out.title = f"log_law_100m_{YEAR}{MONTH:02d}"
out.description = f"Log-law-reconstructed 100 m wind for {cfg['dart_name']}"
out.Conventions = "CF-1.6"
out.method = "log law: v100 = v10 * ln(100/z0) / ln(10/z0)"
out.land_z0_source = (f"ERA5 {era5_year}-{MONTH:02d} monthly fsr climatology, regridded onto TCo1279 "
                       f"(DART year {YEAR} -> ERA5 year {era5_year})")
out.ocean_z0_source = "Charnock relation, solved iteratively from this month's own DART u10/v10"
out.land_sea_mask_source = "ERA5 lsm regridded onto TCo1279 via nearest-neighbour lookup"

sum_ratio_land = sum_ratio_ocean = 0.0
sum_spd10_land = sum_spd100_land = sum_spd10_ocean = sum_spd100_ocean = 0.0
n_land_pts = n_ocean_pts = 0
z0_ocean_min, z0_ocean_max = np.inf, -np.inf

for c0 in range(0, nt, CHUNK):
    c1 = min(c0 + CHUNK, nt)
    u10 = u10_in.variables["10u"][c0:c1].astype(np.float32)
    v10 = v10_in.variables["10v"][c0:c1].astype(np.float32)
    nchunk = c1 - c0
    land2d = np.broadcast_to(land, (nchunk, n_cells))

    spd10 = np.hypot(u10, v10)
    spd_safe = np.maximum(spd10, 0.5)

    z0 = np.where(land2d, 1.0, 2.0e-4).astype(np.float32)
    for _ in range(8):
        ustar = KAPPA * spd_safe / np.log(10.0 / z0)
        z0_new = ALPHA_CH * ustar ** 2 / G + 0.11 * NU / ustar
        z0 = np.where(land2d, z0, z0_new)
    del ustar, z0_new

    z0[:, land] = land_z0_month[None, :]
    ratio = np.log(100.0 / z0) / np.log(10.0 / z0)
    u100 = (u10 * ratio).astype(np.float32)
    v100 = (v10 * ratio).astype(np.float32)

    u100var[c0:c1] = u100
    v100var[c0:c1] = v100

    sum_ratio_land += ratio[land2d].sum(); sum_ratio_ocean += ratio[~land2d].sum()
    sum_spd10_land += spd10[land2d].sum(); sum_spd10_ocean += spd10[~land2d].sum()
    spd100 = np.hypot(u100, v100)
    sum_spd100_land += spd100[land2d].sum(); sum_spd100_ocean += spd100[~land2d].sum()
    n_land_pts += land2d.sum(); n_ocean_pts += (~land2d).sum()
    z0_ocean_min = min(z0_ocean_min, z0[~land2d].min()); z0_ocean_max = max(z0_ocean_max, z0[~land2d].max())

    del u10, v10, spd10, spd_safe, z0, ratio, u100, v100, spd100, land2d
    print(f"  chunk {c0}:{c1}  ({time.time()-t0:.1f}s)", flush=True)

out.close()
u10_in.close(); v10_in.close()

print(f"\nsanity: mean ratio land={sum_ratio_land/n_land_pts:.3f}  ocean={sum_ratio_ocean/n_ocean_pts:.3f}")
print(f"mean spd10 land={sum_spd10_land/n_land_pts:.2f}  mean spd100 land={sum_spd100_land/n_land_pts:.2f}")
print(f"mean spd10 ocean={sum_spd10_ocean/n_ocean_pts:.2f}  mean spd100 ocean={sum_spd100_ocean/n_ocean_pts:.2f}")
print(f"z0 ocean: min {z0_ocean_min:.2e} max {z0_ocean_max:.2e}")
print(f"\nsaved {outfile}  total elapsed {time.time()-t0:.1f}s")
