# TCO95 rerun for the Dunkelflaute study: design spec

Draft 2026-10-28. Status: for review by Jan before anyone sets the run up.

## Purpose

The 9 km (TCo1279-DART) and 31 km (TCo319-DART) slices lack 100 m wind and 3-hourly surface solar, so the Dunkelflaute analysis reconstructs both (log law from 10 m wind, ML model from clouds). Nothing in the future climate can check those reconstructions. A cheap TCO95L91-CORE2 rerun that writes the missing fields directly gives:

1. **Perfect-model test of the method.** Run the reconstructions on TCO95's own 10 m wind and clouds and score them against TCO95's own 100 m wind and `ssrd`, separately in the 1950 control and in the 2080s. This tests whether the reconstruction error changes with climate, which is the main caveat on the 9 km 2080C result.
2. **Event change without any reconstruction.** Mockert events from native 100 m wind and native `ssrd`, compared with events from the reconstructed fields, for 1950C and 2080C.
3. **Natural-variability baseline.** A long 1950 control answers how often two 10-winter windows differ by as much as the 9 km 1950C to 2080C change (+58%) by chance alone.
4. **Time evolution.** A continuous transient shows whether the change emerges gradually or late.
5. **Hub-height sensitivity.** Winds at 100, 120 and 140 m, to see how much taller turbines change event frequency and its climate response. The main study stays at a fixed 100 m.

## Model and configuration

- AWI-CM3 v3.1 (OpenIFS 43r3 + FESOM2 on the CORE2 mesh), TCO95L91, the same configuration as the existing TCO95L91-CORE2 runs `ctl1950d`, `hi1950d` and `ssp585d` (the 20-year subsets in `/work/ab0995/ICCP_AWI_hackthon_2025/`).
- The physics, tuning, forcing and initial state must match the existing runs exactly. **Only the output changes.**
  - Superseded 2026-11-02: the TCO319 tuning of the d runs did not work at TCO95 (strong drift). Production runs `ctl1950g` (E1) and `hi1950g` (E2) in `/work/ab0995/a270092/runtime/awicm3-frontiers-xios/` use release-style parameters and start from the spun-up `own1950g` state. See PLAN.md, "Notes for the paper".
- Output is configured in esm_tools' OpenIFS XIOS templates for cycle 43r3 (`esm_tools/namelists/oifs/43r3/xios/`), not in the generated `file_def.xml` in `work/`, which is regenerated every run. See https://awi-esm3.readthedocs.io/en/latest/how_to/xios_output.html (applies to AWI-CM3 v3.1 and later).
- Field availability, checked against `field_def_ifs.xml` in the AWI-CM3 v3.1 / 43r3 tree: `10u`, `10v`, `100u`, `100v`, `ssrd`, `tsr`, `tisr`, `tcc`, `lcc`, `mcc`, `hcc`, `2t`, `msl`, `sshf`, `blh`, `skt`, `fsr` and `sr` are all defined. Model-level (`3D_ml`) and pressure-level (`3D_pl`) groups exist. **There is no height-level group and no 120 m or 140 m diagnostic**, so those heights come from model-level winds (see below).

## Experiments

| ID | Forcing | Years | Starts from | Minimum |
|---|---|---|---|---|
| E1 `ctl1950` | constant 1950 | 100+ | same initial state as `ctl1950d` | 1950-1969, matching the 9 km 1950C |
| E2 `hist-ssp585` | historical 1950-2014, then SSP5-8.5 2015-2100 | 151 | same initial state as `hi1950d` | 2080-2099, matching the 9 km 2080C (2080-2092) |

The first two years of each are spin-up for the analysis, as for the 9 km slices. The heavy output below should be on for the whole of both runs, so the variability baseline and the time evolution come for free.

## Output

All atmosphere fields on the same regular grid as the existing TCO95 output (`atm_remapped_*`, 400 x 192), so old and new files compare directly.

**3-hourly, 2D, global:**

| Field | Kind | Used for |
|---|---|---|
| `10u`, `10v` | instantaneous | input to the wind reconstruction |
| `100u`, `100v` | instantaneous | truth for the wind reconstruction |
| `ssrd` | accumulated over 3 h | truth for the solar reconstruction |
| `tsr`, `tisr` | accumulated over 3 h | inputs to v5 solar; `tisr` checks the analytic geometry |
| `tcc`, `lcc`, `mcc`, `hcc` | instantaneous | inputs to v5 (all four) and v6 (`tcc`, `lcc`) solar |
| `2t` | instantaneous | temperature during events, demand |
| `msl` | instantaneous | weather regimes, composites |
| `sshf`, `blh`, `skt` | `sshf` accumulated, others instantaneous | stability, to explain the log-law error |

Accumulated fields must be accumulated over the 3 h output interval and say so in `cell_methods`. The pipeline converts J/m2 to W/m2 by dividing by the accumulation window, and nothing in the units attribute distinguishes 3 h from other windows.

**3-hourly, model levels, Europe only:**

- `u`, `v` and `z` (geopotential) on the lowest N model levels, with N chosen so that the highest level is above about 250 m everywhere. Take the level heights from the run's own vertical coordinate (`NAMVV1` in `fort.4` or the init file), not from a generic L91 table.
- If `z` cannot be written on model levels, write `t` and `q` on the same levels plus `sp`, and compute heights hypsometrically offline.
- Restrict this output to a European zoom domain (about 30-75 N, 20 W-40 E). Global model-level output at 3 h is not needed and dominates the storage.

**Daily:** `z` at 500 hPa (weather regimes).

**Monthly:** `ssrd` (the anchor the solar rescale uses at 9 and 31 km), `fsr` and `sr` (the roughness the model itself uses, to compare with the ERA5 `fsr` that the log law borrows).

**Once (static):** `lsm`, surface geopotential, `sdor`.

**Rough storage, uncompressed:**
- 2D 3-hourly: 0.9 GB per field per year (from the existing 6-hourly files, 0.45 GB), about 14 GB/yr for the 16 fields above.
- Model levels on the Europe zoom: small by comparison.
- Both runs together: roughly 250 years, about 3.6 TB uncompressed, likely 2-3x less with netCDF compression.

## How 120 m and 140 m are derived

- Interpolate `u` and `v` linearly in ln(height above ground) between the two model levels that bracket the target height.
- Validate the procedure first at 100 m against the native `100u`/`100v`. If the interpolated 100 m wind does not reproduce the native one, the 120/140 m numbers are not trustworthy.
- Turbine power curves stay generic. Only the hub height changes, so the sensitivity isolates the wind profile.

## Acceptance checks after the first leg

Following the AWI-ESM3 practice of checking output, not the queue:

- Every new file exists in `outdata/` with a non-zero `time_counter`, and the date file has only advanced because output was written.
- 3-hourly stamps fall at 00, 03, ..., 21 UTC, and accumulated fields carry a 3 h `cell_methods` interval.
- The monthly mean of 3-hourly `ssrd` / 10800 equals monthly `ssrd` divided by its own accumulation window.
- Native `100u` speeds exceed `10u` speeds almost everywhere over land and sea.
- If the binary, namelists and restarts are identical to the existing runs, fields the old runs also wrote (6-hourly `10u`, `tcc`, `msl`) should reproduce the existing files exactly at matching times, since output settings should not feed back on the model. Any difference means the configuration is not the same.

## Open items for Jan

- [ ] Which existing TCO95 runs to branch from: the recent one by a270221, the earlier one by a270210, or the subset in the hackathon archive? Where are their run directories, restarts and runscripts?
- [ ] Physics parity across resolutions: the TCo319 runs are named `gibbs-charn`, while the TCO95 and TCo1279 names carry no such suffix. If this is a Charnock change, it alters sea-surface roughness and therefore offshore 10 m to 100 m shear, which is exactly what the perfect-model test is about. Compare the Charnock setting in each run's `fort.4` before the rerun, and match the TCo1279 run if they differ.
- [ ] Length of E1 (100 years or more) and whether E2 runs the full 1950-2100 transient, within the budget.
- [ ] Machine, account and who runs it.
