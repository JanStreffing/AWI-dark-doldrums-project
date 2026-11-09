# Methods notes for the Dunkelflaute paper

A running log, started 2026-11-05. Every method decision, parameter, convention check and validation number goes here when it is made, with the script and output that produced it. Results and next steps live in [`PLAN.md`](../PLAN.md); this file is the Methods section in waiting. Items still to decide or verify are marked **open**.

## 1. Model data

### 1.1 TCo1279-DART (~9 km), the main dataset

- AWI-CM3 (OpenIFS 43r3 + FESOM2 on the DART mesh, 4-25 km ocean), time slices from Moon et al. (2025), run on Aleph (IBS, Busan). Copy on Levante: `/work/ab0995/ICCP_AWI_hackthon_2025/TCo1279-DART-{1950C,2080C}/outdata/oifs`.
- 1950C: constant 1950 forcing, 1950-1969. 2080C: SSP5-8.5, 2080-2092.
- Native reduced Gaussian grid, 6,599,680 cells, 3-hourly: `10u`, `10v`, `2t`, `cp`, `lsp`, `tcc`, `lcc`, `mcc`, `hcc`, `tsr`, `w850`. Monthly on a 5136 x 2560 grid: `ssrd` and other radiation fields. `msl` 3-hourly on the same remapped grid.
- Not written: 100 m wind and sub-monthly `ssrd`. Both are reconstructed (Sections 5-7).
- `sr` (param 173) in the TCo1279 init files is a flat 1e-4 m placeholder, so it cannot serve as a roughness length.
- Accumulated fields are 3-hour sums ending at the stamp. Monthly `ssrd` is the mean of 3-hour accumulations (`cell_methods = "time: mean (interval: 3 h)"`), so it divides by 10800 s.
- Spin-up: the first 2 years of each slice are excluded (proposal Sec. 4.3). A winter is usable only if both calendar years are outside the spin-up and all 7 months exist: 1950C winters 1953-1969 (17), 2080C winters 2083-2092 (10). `tcc`/`tsr` for Feb/Mar 1952 are missing from the archive, which only touches a spin-up winter.

### 1.2 TCo319-DART (~31 km), d runs

- Same model family and Aleph production configuration as 9 km, runs named `*-gibbs-charn`. Levante copy: `/work/ab0995/ICCP_AWI_hackthon_2025/TCo319-DART-{ctl1950d,ssp585d}-gibbs-charn/outdata/oifs`.
- `ctl1950d`: constant 1950, 6-hourly output 1950-1969. `ssp585d`: continuous SSP5-8.5 transient 2015-2099, with 6-hourly output only for 2080-2099 (monthly `10u`/`10v`/`2t`/`z500` for 2015-2099).
- Regular remapped grid, 1312 x 640. 6-hourly instantaneous fields at 00/06/12/18 UTC (`time_instant`): `10u`, `10v`, `2t`, `tcc`, `lcc`, `mcc`, `hcc`, `msl`, `w850`. Every year complete (1460/1464 records, checked 2026-11-05).
- Monthly `ssrd` is the mean of 6-hour accumulations (`cell_methods = "time: mean (interval: 6 h)"`): /21600 s gives a global mean of 184.7 W/m2 (1955), /10800 would give 369.
- No 6-hourly `tsr`, so solar uses v7 (Section 7.3), not v5.
- No land-sea mask or roughness in the output. Land-sea mask from `/work/ab0246/a270092/input/oifs-43r3/TCO319L137/ICMGGa0iwINIT_DART` (DART-matched coastline, 1030 of 421,120 cells differ from the plain init, land fraction 0.286 vs 0.284), mapped to the output grid by nearest neighbour. `sr` in that init file is again a flat 1e-4 m placeholder (one distinct value globally), so land roughness is borrowed from ERA5 as at 9 km.
- **Open**: spin-up exclusion for `ctl1950d`. The `ssp585d` window needs none (continuous transient).

### 1.3 TCO95 (~100 km), Levante reruns for the perfect-model test

- Spec: [`docs/tco95_rerun_spec.md`](tco95_rerun_spec.md). Runs: `/work/ab0995/a270092/runtime/awicm3-frontiers-xios/{ctl1950g,hi1950g}`, both 1950-2099 (E2 lacks 2100).
- Configuration: AWI-CM3 frontiers-xios (OpenIFS frontiers with the yearly aerosol-scaling branch `awicm-3-aerscaleyearly-xios`, FESOM2 `aleph_prod` 0374b8c5), CORE2 ocean mesh, TCO95L91. Release-style AWI-CM3 v3.1 parameters (`namelists/fesom_own_gm`), no TCO319 tuning. Both runs branch from the 20-year spin-up `own1950g` at 1974-12-31, relabelled to 1950-01-01. E1 holds 1950 forcing (`NCMIPFIXYR 1950`), E2 is historical then SSP5-8.5.
- **TCO95 is a separate configuration, not just a third resolution.** The TCO319 tuning of the d runs (RPRCON 1e-3, ENTRORG 1.25e-3, GGAUSSB -0.45) did not work at TCO95: the tuned control `ctl1950df` drifted, 2t from 14.36 to 12.99 C over 1950-68. Any 100 vs 31/9 km difference therefore mixes resolution with tuning.
- Output used here (regular 400 x 192 grid), 3-hourly: `10u`, `10v`, native `100u`/`100v`, `ssrd`, `tsr`, `tisr`, `tcc`, `lcc`, `mcc`, `hcc`, `2t`, `msl`, `sshf`, `blh`, `skt`; model levels 85-91 (`u`, `v`, `z`, `t`, `q`, `sp`). Daily `z` at 500 hPa. Monthly `ssrd`, `fsr`, `sr`, `lsm`.
- Conventions, checked against the output 2026-11-05: instantaneous fields stamped at the instant, stamps 03:00 ... 00:00. Accumulated fields cover the 3 h ending at the stamp: model `tisr` against analytic `tisr` for that window, r = 1.0000 and mean ratio 0.9997 (centred window r = 0.92, following window 0.72). A stamp belongs to the month of (stamp - 1 s), as on DART.
- Model-level heights over the open North Sea (January 1960 mean of `z`/g): 32, 55, 88, 131, 184, 247, 321 m for levels 91 to 85, so 120 and 140 m are bracketed. **Open**: subtract surface geopotential before interpolating; the lowest level at 32 m rather than about 10 m suggests nonzero spectral orography over the sea.

## 2. Reanalysis and observations

- **ERA5**, 0.25 deg regular grid. Wind training data: 3-hourly `u10`/`v10`/`u100`/`v100`, full months 2015-2026 and Oct-Feb winters 1979-2014 (`wind/data/era5/training_data/`). Solar: 3-hourly clouds (instantaneous) and genuine 3-hour sums of `ssrd`/`tsr` built from hourly data (`solar/data/era5_3h_mean/`), plus monthly `ssrd` (mean of daily accumulations, /86400). `fsr` monthly 2015-2025, `lsm`.
- Accumulation divisors, J/m2 to W/m2: ERA5 monthly 86400; ERA5 hourly 3600; ERA5 3-hour sums 10800; DART monthly and 3-hourly 10800; TCo319 monthly 21600; TCO95 3-hourly and monthly 10800. The units attribute is `J m**-2` in all cases, so the window has to be read from `cell_methods` or verified.
- **SMARD** (Bundesnetzagentur), Germany, hourly "Realisierte Erzeugung" for solar, wind onshore and wind offshore, 2015-2026. Resampled to 3 h with `label="right", closed="right"`, matching the backward-looking ERA5 accumulation. An earlier forward-looking bin had depressed the solar correlation from about 0.97 to 0.75.
- SMARD installed capacity: Jan-1 anchors per year, interpolated linearly in time, not held flat. Offshore capacity tripled within 2015 (993 to 3283 MW), and a flat denominator gave CF up to 297%. CFs clipped at 1.02. 2015-17 are flagged as offshore capacity-step years; excluding them raises offshore r against SMARD from 0.70 to 0.86 (reconstruction) and 0.87 (ERA5).
- SMARD is used as an independent sanity check only, never as a calibration target. Curtailment lowers SMARD CFs in high-output hours, but does not touch Dunkelflaute, which is a shortage regime.

## 3. Domain and spatial aggregation

- Germany onshore: land cells whose centre lies inside the Natural Earth land polygon (`boundaries/germany/germany_land.geojson`). Offshore: sea cells inside the Marine Regions EEZ polygon (`germany_eez.geojson`). Land means `lsm >= 0.5`.
- On DART the land/sea flag is ERA5 `lsm` regridded by nearest neighbour onto TCo1279 (`wind/data/results/tco1279_era5_regrid_index.npz`). On TCo319 it is the DART init `lsm` (Section 1.2), and on TCO95 the model's own monthly `lsm`.
- Cell counts: ERA5 729 land, 123 EEZ; DART 4005 land, 430 EEZ; TCO95 55 land, 10 EEZ.
- Weighting: cos(lat) on regular grids (ERA5, TCo319, TCO95) for wind and solar. On DART, solar is cos(lat)-weighted and wind unweighted. DART cells are near-equal-area, so the cos(lat) weight there is a small inconsistency (Germany spans cos(47 N)/cos(55 N) = 1.19).

## 4. Capacity factors and combination

- Wind: generic cubic power curve, applied per sample before any averaging. Zero below cut-in, (v^3 - v_in^3)/(v_rated^3 - v_in^3) up to rated, 1 to cut-out, 0 above. Onshore (3, 12, 25) m/s, offshore (3, 13, 25) m/s. Hub height 100 m on land and at sea. A fixed hub height is a deliberate choice: fleet realism (hub heights, park age, per country) is out of scope. Hub-height sensitivity at 120/140 m only from TCO95 model levels.
- Solar: CF = clip(G / 1000, 0, 1) with G the window-mean irradiance in W/m2 (STC 1000 W/m2). No temperature or orientation correction.
- Combined CF = w_s CF_solar + w_on CF_on + w_off CF_off, with Mockert et al. (2023) weights, Germany's 2018 mix (IRENA): 0.44 / 0.50 / 0.06. Sensitivity: Germany 2024 (BNetzA) 0.577 / 0.369 / 0.054.

## 5. Wind reconstruction: 10 m to 100 m

- Log law, per component: v100 = v10 ln(100/z0) / ln(10/z0), no stability correction, no displacement height.
- Land z0 from ERA5 `fsr`. On DART: monthly fields 2016-2025, year-cycled onto the run (DART year Y uses ERA5 year 2016 + (Y - first run year) mod 10, same month), regridded by nearest neighbour. On TCO95 and TCo319: monthly climatology 2016-2025, geometric mean of the ERA5 land cells inside each model cell.
- Ocean z0 from the Charnock relation, solved by 8 fixed-point iterations from that step's own 10 m wind: z0 = alpha u*^2/g + 0.11 nu/u*, u* = kappa v10 / ln(10/z0), with kappa 0.4, g 9.81, nu 1.5e-5, alpha 0.018. First guess 2e-4 m, 10 m speed floored at 0.5 m/s inside the iteration.
- Code: `wind/scripts/methods/dart_loglaw/reconstruct_month.py` (DART), `dunkelflaute/scripts/tco95/perfect_model_series.py` (TCO95).
- Against ERA5 native 100 m wind: RMSE 1.19 m/s onshore and 0.60 m/s offshore, r 0.97 and 0.997. Below 3 m/s, MAPE 28-31% onshore and 15% offshore; the log law beat an ERA5-fitted regression on MAPE and was chosen because it carries no grid-specific coefficients.
- Known bias: on ERA5 2016-26 the log law gives onshore CF 0.252 against 0.207 from native 100 m wind (+22%), and 46% too high in calm hours. The reconstruction counts 17% fewer events than native ERA5 (5.64 vs 6.82 per winter).
- Not addressed (limitations): fetch (a 9 km cell is at or below the about 10 km fetch a 100 m wind needs to equilibrate), displacement height over forest, stability.

## 6. Wind calibration: quantile mapping v1

- Script `wind/scripts/methods/calibration/calibrate_quantile_mapping.py`, samples from `extract_germany_era5_samples.py`, table `wind/models/quantile_mapping_v1/qm_table.npz`, evaluation `dunkelflaute/results/wind_calibration/qm_evaluation.json`.
- Target: ERA5 native 100 m wind speed, Germany cells.
- Bins use only quantities every model grid has: month (Sep-Mar, 7) x 3-hourly UTC slot (8, by stamp hour) x surface class (4: offshore; onshore `fsr` < 0.1, 0.1-0.5, > 0.5 m, fixed edges). 105 quantiles (0.001-0.999, finer in both tails). Linear interpolation between quantiles; beyond the outermost quantile the end offset is held constant; output floored at 0.
- Validation:
  - 2-fold cross-validation by winter on Sep-Mar 2016-26: onshore CF 0.252 to 0.207 (native 0.207); calm-hour bias 1.46x to 1.10x; events 5.64 to 6.18 per winter (native 7.00); per-winter event correlation 0.84 to 0.92. All three use the reconstructed solar CF (`cf_solar_recon`), so only the wind differs; with native solar as well, native ERA5 gives 6.82 (Section 5, fig_validation). Checked 2026-11-08.
  - Independent Oct-Feb 1980-2014 winters, never fitted: onshore CF 0.281 to 0.236 (native 0.233).
- Applied to DART (`dunkelflaute/scripts/dart/apply_wind_calibration_dart.py`): events 1950C 3.47 to 3.59, 2080C 5.50 to 6.00, ratio 1.58 to 1.67, threshold curve unchanged in shape.
- On TCO95 the onshore class uses the annual-mean ERA5 climatology z0 of each cell, not the month's value. **Open**: align with the DART application (month value), minor.

## 7. Solar reconstruction: monthly to sub-daily ssrd

### 7.1 Analytic top-of-atmosphere irradiance

- `solar/scripts/core/solar_geometry.py`: Spencer (1971) / Iqbal (1983) for eccentricity, declination and equation of time; S0 = 1361 W/m2; cos z clipped at 0.
- Accumulated over the window ending at the stamp by the midpoint rule: 12 sub-steps on regular grids for 3 h, 8 on the DART template, 24 for 6 h windows.
- Model `tisr` reproduces it (TCO95: r = 1.0000).

### 7.2 v5 (3-hourly, used on DART)

- Target kt = ssrd / tisr, clipped to [0, 1.5]. Features: T = tsr/tisr (clipped [0, 1.5]), tcc, hcc, mcc, lcc, mu = tisr / (dt S0). Daylight rows only (mu > 0.02).
- `HistGradientBoostingRegressor` (sklearn 1.1.1), 250 iterations, random_state 0, sample weights cos(lat).
- Training data: ERA5 Sep-Mar 2015-2026 (80 months), every 8th grid point, cells south of 60 S and above 3000 m excluded, 32 M rows subsampled. Three year-block CV folds (2015-18, 2019-22, 2023-26) plus a final model on all years.
- Exact constrained monthly rescale (`solar/scripts/core/rescale.py`): the prediction is scaled per cell so each calendar month's sum equals the true monthly total, with kt capped at 1.1 by iterative water-filling. Twilight cells (analytic tisr zero all month but a nonzero monthly total) get a weight peaking at local noon.
- Validation, Germany, ERA5: r = 0.998 against native `ssrd`; against SMARD r = 0.931, equal to ERA5's own 0.931. Before rescale R2 = 0.994; held-out folds RMSE 8.5 W/m2, flat across years.

### 7.3 v6 and v7 (6-hourly, for TCo319 and TCO95)

- Same target, model, mask, weights, training data and folds as v5, on 6 h windows ending at 00/06/12/18 UTC. ERA5 6 h `ssrd` = ssrd(t) + ssrd(t - 3 h). Clouds are taken only at the window's two ends, as 6-hourly output provides; the 3 h sample in between is deliberately unused. mu over the 6 h window.
  - v6 features: tcc, lcc at both ends, mu.
  - v7 features: tcc, lcc, mcc, hcc at both ends, mu (the TCo319 d runs write all four).
- Code `solar/scripts/core/reconstruct_ssrd_v6.py`, `preprocessing/extract_v6_6h.py --variant v6|v7`, `training/train_v6_model.py --variant v6|v7`, `jobs/train_v{6,7}_model.sbatch`.
- Pilot check: mean 6 h kt for January 2015 equals v4's 3 h value (0.539).
- v6, Germany, ERA5 held-out folds, 6 h windows (`solar/scripts/evaluation/evaluate_v6_germany.py`): r = 0.982, RMSE 19.8 W/m2, bias +0.2 W/m2 (v5 aggregated to 6 h: 0.996, 9.1). With native ERA5 wind, events per winter at 0.06: native solar 6.73, v5 6.91, v6 6.73; within 0.5 events at every threshold 0.04-0.10.
- Without the monthly rescale v6 is +17 W/m2 too bright over German winters and undercounts events by 15%, so the rescale is essential.
- v7, same evaluation: r = 0.989, RMSE 15.6 W/m2, bias +0.2 W/m2; before rescale +14.3 W/m2. Events with native ERA5 wind: 6.91 per winter at 0.06 (native solar 6.73), within 0.2 at every threshold 0.04-0.10. Training job 28748419, 30 min.

## 8. Event detection

- Mockert et al. (2023), Sec. 2.4: trailing running mean of the combined CF over 48 h (16 steps at 3 h, 8 at 6 h), threshold 0.06. Every time step contributing to a below-threshold mean is flagged (window expansion), and contiguous flagged steps form one event, so the minimum duration is 48 h.
- Extended winter Sep-Mar, labelled by the year of its Jan-Mar. The running mean is computed per winter block and never bridges Apr-Aug. A winter needs at least 1600 3-hourly steps (complete). Code `dunkelflaute/scripts/dart/compute_germany_dunkelflaute_dart_mockert.py`, `era5/compute_germany_dunkelflaute_2015_2026_mockert.py`, shared functions in `scripts/analysis/robustness_mockert_events.py`.
- Sensitivity: thresholds 0.04-0.10, windows 24/48/72 h, 2024 weights.
- Li-style comparison: capacity-weighted wind CF < 0.2 and solar CF < 0.2, instantaneous, runs longer than 24 h. At night this reduces to the wind condition.
- Mockert's own reference: ERA5 1979-2018, July-June years, about 4.3 events per year, median 3.2 d, longest 8 d, almost all in SON/DJF, CFs from atlite bias-corrected against observed generation. Our ERA5 rate (5.64, native 6.82) is above that; the uncorrected power curves are the likely cause.

## 9. Statistics

- Event counts per winter compared with Welch t, Mann-Whitney, an exact conditional Poisson test, and a year-block bootstrap of the ratio (20,000 resamples of whole winters). 9 km result: 3.47 to 5.50, ratio 1.58, Poisson p = 0.015, bootstrap 90% interval 1.19-2.12 (`dunkelflaute/results/robustness/robustness_summary.json`).
- These within-run tests assume independent winters and understate decadal variability (Section 10). The block null from the TCO95 control is the more honest test: for Germany at 9 km it gives p = 0.13-0.20 (combined, base-rate matched, Section 10c), so ~~significant at p = 0.015~~ the Germany-only result is not significant on its own.
- 6-hourly sampling (`dunkelflaute/scripts/analysis/sampling_3h_vs_6h.py`): wind sampled at 00/06/12/18 UTC, solar summed over 6 h, 48 h = 8 steps. Events per winter change by 0% (ERA5 native and recon, 2080C), +2% (SMARD), -8% (1950C) at 0.06, per-winter r 0.93-0.99.

## 10. Perfect-model test and natural variability (TCO95)

- Code `dunkelflaute/scripts/tco95/perfect_model_series.py` (Germany series, job `dunkelflaute/jobs/perfect_model_tco95.sbatch`) and `perfect_model_events.py` (events and statistics); summary `dunkelflaute/results/tco95/perfect_model_summary.json`.
- For every year 1952-2099 of both runs, Germany CF from the model's native 100 m wind and `ssrd`, and from the 9 km recipes applied to its 10 m wind and clouds: log law with borrowed ERA5 z0, the same plus calibration, log law with the model's own `fsr`, v5 solar with exact monthly rescale to the model's own monthly total.
- Periods: control = `ctl1950g` winters 1953-2099 (147 winters) and its 1953-1969 subset; future = `hi1950g` winters 2081-2099 and 2083-2092.
- Results:
  - **Reconstruction error vs climate.** Onshore log-law CF / native 1.239 (control) and 1.264 (2081-99); calibrated 1.035 and 1.055; model-`fsr` log law 1.272 and 1.298. Offshore log law 0.889 and 0.885 (calibrated 0.929 and 0.925). v5 solar after rescale 1.000 in both, r = 0.996; before rescale 0.974 and 1.019, so the rescale also absorbs a climate-dependent solar error.
  - **Event ratio preserved.** At 0.06 the reconstruction undercounts (control: native 7.08, calibrated 6.48, uncalibrated 6.00 per winter), but the future/control ratio matches: native 0.99, calibrated 1.03 (2081-99 vs 1953-2099).
  - ~~**No change at TCO95.**~~ **No change at TCO95 for Germany** (Mockert mix). Future/control about 1.0 at all thresholds 0.05-0.10. Over 24 NH regions TCO95 does show the increase (Section 10c), so this is a regional, not a model-wide, null. Decadal means along E2 (6.3-9.2) barely exceed the control's range (4.9-8.9) (corrected 2026-11-08: earlier 'stay inside').
  - **Regional perfect-model test (2026-11-06, paper repo `scripts/perfectmodel/reduce_perfectmodel_regions.py`, from existing series in seconds).** Over the 24 regions, where TCO95 does change, the uncalibrated log law on the 6-hourly footing with matched thresholds reproduces the native ratio: geometric mean 1.30 against 1.30 (wind), 1.19 against 1.21 (combined, 18 regions). Single regions scatter: r of log ratios 0.64 (wind) and 0.87 (combined), rms difference 12 % and 8 %. Native ratio range across thresholds 0.05-0.10 for 2081-2099 is 0.99-1.10 (0.95 only for the 2083-2092 span).
  - Note: the perfect-model test used v5 solar (TCO95 writes 3-hourly `tsr`); the 31 km v7 recipe is untested in perfect-model mode.
  - **Natural variability.** All placements of a 10-winter block against a disjoint 17-winter block in the 147-winter control (14,762 pairs): at 0.06 (7 events/winter) the ratio reaches 1.58 in 0.3%. At thresholds giving a 9 km-like base rate it reaches 1.58 in 2.3-5.3% (0.040: 2.85 events/winter, 2.9%; 0.045: 3.78, 4.8%; 0.050: 4.67, 5.3%; 0.055: 5.85, 2.3%). The 5-95% range at 0.045 is 0.60-1.58, wider than Poisson noise (about 0.66-1.34), so decadal variability matters.

## 10b. 31 km pipeline and the three-resolution comparison

- Germany series: `dunkelflaute/scripts/tco319/compute_germany_tco319.py`. Box 37 x 45 cells; Germany 593 land and 88 EEZ cells (land/sea from the DART init `lsm`, Section 1.2). ERA5 z0 per onshore cell, median 0.77 m in January (ERA5 itself and TCO95: 0.75).
- Wind: log law and calibration as in Sections 5-6, applied to the 6-hourly instantaneous 10 m wind. UTC slots are then 0, 2, 4, 6 of the calibration's 8. The onshore class follows the month's z0, as in the DART application.
- Solar: v7 on 6 h windows ending at the stamp. Clouds at the stamp and 6 h earlier come from consecutive 6-hourly records; the window ending 06 UTC on 1 January needs the previous year's file and is not built (night over Germany). Exact monthly rescale with target = monthly-mean `ssrd` x 4 x days in month.
- Events: `dunkelflaute/scripts/tco319/tco319_events.py`. All three resolutions on a 6-hourly footing (48 h = 8 steps). 9 km and 100 km 3 h series are degraded as in Section 9: wind sampled at 00/06/12/18 UTC, solar the mean of the two 3 h values. 9 km uses v5 solar with uncalibrated and calibrated wind (`results/wind_calibration/dart_calibrated_*.npz`, same time axis as the stored series, checked). Ratios with an exact Poisson test and a winter bootstrap; output `dunkelflaute/results/tco319/tco319_events_summary.json`.
- Job: `dunkelflaute/jobs/tco319_pipeline.sbatch` (also re-scores v5/v6/v7 in `evaluate_v6_germany.py`). Run 2026-11-05, job 28749088, 3 min.
- Independence check: per-winter event counts at 0.06, 9 km vs 31 km in the same years, r = -0.35 (1953-69, n = 17) and 0.04 (2083-92, n = 10); winter-mean onshore CF r = -0.12 and +0.47. The two resolutions are independent samples at winter scale.
- Results are in PLAN.md (TCO95 and 31 km sections).

### 1.4 EERIE IFS-FESOM2-SR (~9 km), second high-resolution model

- IFS TCo1279 atmosphere with FESOM2 NG5 ocean (EERIE project, AWI/ECMWF), member r1i1p1f1, version v20240304. Runs: `hist-1950` 1970-2014, `eerie-control-1950` 1971-2050 (constant 1950 forcing), `highres-future-ssp245` 2015-2070.
- Used: the 0.25 deg regridded 6-hourly output, which is ERA5's own grid (1440 x 721, 90 N to 90 S, 0-359.75 E): instantaneous `10u`, `10v`, `2t`, and `ssrd` accumulated over the 6 h ending at the validity time (step 0-6; global mean /21600 s = 181 W/m2).
- Access: raw GRIB on `/work/bm1344/a270228/fdb_{hist1950,ctl1950,...}/REGULARLL_grids/`, located through the kerchunk parquet references under `/work/bm1344/k202193/Kerchunk/IFS-FESOM_sr/<run>/atmos/gr025/2D_6hourly_{instant,accumulated}.parq` (only path, offset and size are used: `dunkelflaute/scripts/regions/eerie_index.py`). The parquet stores do not open with zarr 3, so messages are decoded with eccodes (legacy python env), about 0.13 s per field.
- The reference pages are padded with empty rows to 100,000; trailing empties are dropped. Row k is time step k, except in the control's instantaneous store, which misses three 6 h steps (rows 6815, 20210, 24589). Times of every row are resolved by bisection on decoded validity times and spot-checked every 997 rows; every decoded message is checked against them.
- Land/sea: ERA5 `lsm`. Land z0: the ERA5 fsr climatology of the same cell, so no regridding of z0, and the QM calibration (fitted on this grid) applies unchanged.
- Daily (`2D_daily_avg`: `mmsl`) and monthly (`2D_monthly_avg`: `mean2t`) stores start in 1950 for hist and control (65 and 101 years), unlike the 6-hourly ones. Their GRIB2 messages carry a 52-byte local-use section 2 where eccodes (2.43 in the Python bindings, 2.46 on the command line) expects 57 bytes for that local definition, and every later key then fails. Section 2 is optional in GRIB2, so it is dropped before decoding (`regions.grib2_strip_local`); the decoded fields are then normal regular_ll 1440 x 721 (`avg_msl`, paramId 235151, 24 h average from 00 UTC of dataDate).
- Forcing differs from the AWI-CM3 runs: SSP2-4.5 to 2070, against SSP5-8.5 to 2099. A ratio from EERIE is therefore for less warming.
- Supplement restructured 2026-11-09 (Jan): Notes ordered by first citation, number chains moved into Supplementary Tables S1-S7 (transfer checks, cross-model correlates, LUMIP effects, height levels, reconstruction checks, IFS-FESOM baseline, footings), region coordinate table replaced by Supplementary Fig. S1. Table values were read from the result files, not the old prose. Corrections found on the way: CMIP6 V1000/G is +1.95 % (shown as +1.9, formerly '2.0'; 18 of 19 positive confirmed); the 2 m temperature change correlates with the fixed-G term at r = -0.34 (Europe) and -0.40 (22 regions), now listed; 100 km native 100 m wind -2.1 % against -1.3 % at 10 m; 'no CMIP6 model archives a 100 m wind' is wrong (EC-Earth3 ua100m, 6hrPlev, ESGF).
- Paper presentation, 2026-11-09 (Jan): the model is called "IFS-FESOM" (as in Ghosh et al. 2026; "SR" is not defined in any source we found), labelled "IFS-FESOM 9 km" in the figures, and enters the main text and figures only as SSP2-4.5 2041-2070 against the control 1972-2050, like every AWI-CM3 comparison. The historical baseline (1971-2000) moved to Supplement app:baseline. Analysis outputs still carry both (`eerie_ctl`, `eerie_hist`); the paper's `scripts/common/dfroles.py` RUN_ORDER drops `eerie_hist`.
  - Paper numbers that changed: wind fit for our runs +5.37 -> +5.34 %/K, mix +3.85 -> +3.83 %/K; wind pattern r with CMIP6 0.31-0.63 -> 0.41-0.63; warming 2.0-2.3 -> 2.3 K; mean G over Europe -1.2..-4.0 -> -1.7..-4.0 %; transfer +1.3..+3.4 -> +2.0..+3.4 % (bilinear G), +1.6..+4.1 -> +2.3..+4.1 % (box G); at fixed G +0.97..+2.02 -> +1.38..+2.02 % (Europe), +0.28..+0.80 -> +0.53..+0.80 % (22 regions); ScBL events per 100 regime days 3.4-3.9 -> 3.4-3.8; CMIP6 models at or above our lowest transfer change 3 -> 2 (Europe), 8 -> 2 (22 regions).
  - Corrected at the same time: the winter-mean wind CF change "-0.3 to -4.2 %" came from 31 km on the 9 km spans; on the four main comparisons it is -0.6 to -4.2 %. "In every variant at or above the CMIP6 90th percentile" was already too broad for the control comparison: in the CMIP6-style daily extraction on 22 regions EERIE ctl gives +1.23 % against a p90 of 1.31 % (3 of 20 models at or above); per K it is above (+0.53 against +0.34 %/K).
  - Historical baseline against control: dT 2.02 / 2.30 K; wind +5.88 / +5.75, solar +3.43 / +4.49, mix +4.37 / +3.17 %/K; wind map area mean +4.74 / +4.79 %/K, FDR boxes 99 / 211, pattern r 0.31 / 0.41; within-regime term positive in 5 / 7 of 8; transfer Europe (bilinear, 10 m) +1.30 / +1.95 %, 22 regions (A3) +0.75 / +1.30 % (8 / 2 CMIP6 models at or above); per K over Europe above every CMIP6 model in both.

## 10c. Multi-region test (all four models)

- Purpose: whether the 9 km increase is a resolution effect (then 9 km should exceed 31 km in most regions) or sampling (then 9 km minus 31 km scatters around zero). Outcome: neither framing held. All resolutions show an increase in most regions; 9 km minus 31 km is small and not significant.
- Regions: `dunkelflaute/scripts/regions/regions.py`. 24 Northern-Hemisphere land regions: DE (land polygon) plus 23 lat/lon boxes (FR, IB, UK, PL, SC, IT, BK, UA, RU, TR, KZ, US_NP, US_SP, US_MW, US_NE, US_CA, CA_PR, CN_NE, CN_N, CN_E, KR, JP, IN_N), land cells (`lsm >= 0.5`) only. Plus DE_off (sea cells in the German EEZ) to rebuild the Mockert Germany mix. Southern Hemisphere left out: the 9 km wind and solar reconstructions exist for Sep-Mar only.
- Cell counts (DE, DE_off): DART 4006, 587; TCo319 593, 88; TCO95 55, 10; EERIE/ERA5 grid 729, 121.
- **DART Germany offshore in the original pipeline is clipped.** `compute_germany_dunkelflaute_dart_mockert.py` takes candidate cells from `BBOX_GERMANY` (47-55 N, 6-15 E), while the EEZ polygon spans 52.9-55.9 N, 3.35-14.75 E. It keeps 430 of 587 EEZ cells and drops the EEZ north of 55 N and west of 6 E (and a sliver of land north of 55 N and east of 15 E: 4005 vs 4006). ERA5, TCo319 and TCO95 Germany series use the full polygons. Offshore carries 6% weight, so the effect on events is small, but the 9 km Germany offshore CF is not like-for-like. The region extraction uses the full polygons on every grid.
- Extractors (`dunkelflaute/scripts/regions/extract_{dart,tco319,tco95,eerie}.py`, job `dunkelflaute/jobs/regions_extract.sbatch`), per region cos(lat)-weighted means (DART unweighted):
  - DART 9 km: the existing global reconstructions (log-law 100 m wind, v5 solar with rescale), 3-hourly. Check against the Germany pipeline for January 2085: solar 0.0335 vs 0.0333 (the pipeline cos(lat)-weights DART solar), onshore wind 0.4647 vs 0.4647.
  - TCo319 31 km: the Germany recipe (Section 10b) on every region's land cells; ERA5 z0 from the ERA5 land cells within half a model cell. Check for 2085: DE solar, uncalibrated and calibrated onshore wind and calibrated offshore wind identical to `compute_germany_tco319.py` to 1e-10.
  - TCO95 100 km: uncalibrated log law (same recipe as the others), native 100 m wind, native 3 h `ssrd`.
  - EERIE: uncalibrated log law and QM, native 6 h `ssrd`, `2t`.
- Wind in the region comparison is the uncalibrated log law in every model: the QM table was fitted over Germany only, and at 9 km the Germany ratio is the same either way (1.73 vs 1.76).
- Events (`region_events.py`): 6-hourly footing, 48 h = 8 steps, Mockert window expansion. Metrics: comb = 0.44 solar + 0.56 wind, and wind alone. The vectorised counter matches `mockert_events` exactly (0 mismatches over 114 winter x threshold cases).
- Thresholds are matched to base rate: per comparison, region and metric, the lowest threshold at which that comparison's control has 5 events per winter. This removes mean-CF differences between models and regions; regions with little winter wind (TR, US_CA, CN_E, IN_N) then measure solar-driven lows. Fixed 0.06 is reported as well.
- Validity: a region enters a field test for a metric only if every comparison's control and both null controls reach at least 90% of the target rate. For the combined metric this leaves out TR, US_CA, CN_N, CN_E, KR and IN_N (18 regions remain): in sunny, low-wind regions low combined CF comes in long spells that merge, so 5 separate events per winter are never reached at any threshold. Wind alone: all 24 regions valid.
- Run 2026-11-05: extraction jobs 28757727-28757899 (DART 13 and 9 min, TCo319 6 min, TCO95 3 min, EERIE 9-15 min per run); `region_events.py` 1 min. Results in PLAN.md (Phase 2b).
- Tests: Poisson exact per region; per-region block null (all disjoint 17/10-winter placements) from the TCO95 constant-1950 control (147 winters) and the EERIE control (79 winters). Field tests: mean over regions of log ratio, and for 9 km against 31 km the mean of log(r9/r31), with the null drawn from block placements using the same winters in every region, so the spatial correlation of the real samples is kept.

## 10d. Normalisation per K of global warming

- Global annual-mean 2t per year (`dunkelflaute/scripts/regions/global_t2m.py`, outputs `dunkelflaute/data/regions/global_t2m_<model>.json`), complete years only, area-weighted:
  - DART 9 km: 3-hourly reduced-Gaussian `2t`, every 9th step (9 and 8 coprime, so all eight 3-hourly slots are sampled equally), weight cos(lat)/points-per-ring.
  - TCo319: `ctl1950d` 6-hourly `2t`; `ssp585d` monthly `2t` (2015-2099). TCO95: monthly `2t`. EERIE: monthly `mean2t` (GRIB, section 2 stripped, Section 1.4).
- dT per comparison = mean of the future years minus mean of the control years (calendar years matching the winter spans). Sensitivity per region = ln(ratio)/dT, quoted as 100 (exp(ln(ratio)/dT) - 1) %/K. Field mean over the valid regions; 5-95% from 5000 winter resamples, the same winter indices in every region (`region_per_K.py`). This interval is within-run and omits decadal variability beyond the sampled winters.
- The DART 1950C control cools by 0.76 K over 1950-69 (13.25 to 12.49 C, 5-year means), a drift that dT includes because it uses the compared years.
- Mean vs tail: the future series is rescaled to the control's Sep-Mar mean (x * mean_ctl / mean_fut, per region) and events recounted at the same threshold; the share of ln(ratio) that disappears is the share explained by the mean change.

## 13. Weather regimes

- msl: daily means on a 2.5 deg Euro-Atlantic grid, 20-80 N, 80 W-40 E (25 x 49), by box mean of all source points within +-1.25 deg of each box centre (`dunkelflaute/scripts/regimes/extract_msl.py`, outputs `dunkelflaute/data/regimes/msl_<source>_<run>.npz`). Source points per box: ERA5 N320 9-81, TCO95 4-9.
  - ERA5: `/pool/data/ERA5/E5/sf/an/1D/151`, daily average msl, reduced Gaussian N320, Sep-Mar 1979-2025.
  - DART 9 km: 3-hourly msl on the remapped 5136 x 2560 grid (uncompressed, 13 GB per month). TCo319: 6-hourly. Day = mean of the 06, 12, 18 and 24 UTC instants, i.e. the window (d 00, d+1 00]; a first version used 00-18 UTC and lost the 1st of every month on DART, because a monthly file starts at 03 UTC.
  - TCO95: daily file `atm_remapped_1d_msl`. EERIE: daily average `mmsl` (avg over the 24 h from 00 UTC of dataDate), 1950-2050 (control), 1950-2014 (hist).
- Regimes (`regime_analysis.py`): ERA5 Sep-Mar 1979/80-2024/25 (9,976 days). Anomalies from a mean + 3 annual harmonics (least squares on Sep-Mar days). Weights sqrt(cos lat). EOFs by SVD; 10 PCs explain 82% (EOF1-4: 20, 15, 14, 8%). k-means, k = 4, 50 initialisations (own implementation, numpy). Clusters named from the centroid pattern (Iceland-Azores dipole, Scandinavian and mid-Atlantic highs): NAO+ (low -16 hPa at 60 N 15 W), NAO- (high +14 hPa at 65 N 35 W), ScBL (high +14 hPa at 60 N 10 E), AR (high +13 hPa at 52 N 25 W). ERA5 frequencies 27/24/26/23%.
- Models: anomalies from the seasonal cycle of the comparison's control period, applied to both periods, projected on the ERA5 EOFs (ERA5 weights), assigned to the nearest ERA5 centroid in 10-PC space. Sensitivity: each period's own seasonal cycle (removes the forced mean msl change).
- Events by regime: European regions DE, FR, IB, UK, PL, SC, IT, BK; same events and base-rate-matched thresholds as Section 10c. An event's regime is the most frequent regime over its calendar days (a day = the 6-hourly steps whose (stamp - 1 s) falls on it). Decomposition per region with f_k regime days per winter and lam_k events per regime day: de = sum lam_k^ctl df_k (frequency term) + sum f_k^fut dlam_k (within-regime term); the two terms add up to the total exactly.
- Per-K form for the paper (fig:mechanism, 2026-11-07, `paper/scripts/mechanism/reduce_mechanism.py`): the same decomposition summed over the 8 regions and divided by the comparison's dT (9 km 5.50, 31 km 5.68, 100 km 4.90, EERIE ctl 2.30, hist 2.02 K), in events per winter per K. Total = frequency + within: 9 km +0.246 = -0.036 + 0.282; 31 km +0.233 = -0.024 + 0.257; 100 km +0.220 = -0.014 + 0.234; EERIE ctl +0.135 = -0.028 + 0.163; hist +0.144 = -0.027 + 0.172. Within-regime term positive in 8, 8, 8, 7, 5 of 8 regions. Mean G change per regime (fig:mechanism d) in %/K: CMIP6 median ScBL -1.13 (20 of 20 negative), AR -1.14; EC-Earth3 below all 20 CMIP6 models in ScBL (-2.89) and AR (-2.19), and over Europe as a whole -9.2 % or -1.85 %/K (next HadGEM3-GC31-MM -1.30 %/K). The earlier "G falls by 11.5 and 7.6 %" (9 and 31 km, ScBL) is -2.21 and -1.38 %/K.
- Name of the blocking regime, checked 2026-11-07 against the ERA5 centroid (`paper/scripts/ed_centroids/`, fig:ed_centroids): the high peaks at +13.9 hPa at 60 N 10 E (southern Norway, Skagerrak); +12.3 hPa over northern Scandinavia, +6.3 hPa over Germany; low at 65 N 40 W. Single blocking days peak between 50-70 N and 27.5 W-40 E; south of 57.5 N on 28 % of days; anomaly larger over Germany than over northern Scandinavia on 34 %; Germany positive on 79 %. Reassigning ERA5 days reproduces the frequencies (27.0, 23.9, 26.5, 22.6 %). This is the fourth regime of the classical four-regime set, "often referred to as Scandinavian blocking" (Cassou, ECMWF seminar 2010; ECMWF Newsletter 165 "European/Scandinavian blocking"); the 7-regime set of Grams et al. (2017) splits it into European and Scandinavian blocking, so our cluster takes in days of both. Name kept (Jan: no new names).

## 14. CMIP6 on a daily footing, and the daily bridges

- Inventory (`dunkelflaute/scripts/cmip6/inventory.py`, output `dunkelflaute/data/cmip6/inventory.json`): models under `/pool/data/CMIP6/data/ScenarioMIP/*/*/ssp585` with a historical counterpart; the first member (r1i1p1* preferred) that has `day/sfcWind`, `day/rsds` and `Amon/tas` in both experiments; latest version directory. 20 models: AWI-ESM-1-REcoM, FGOALS-g3, CanESM5, CMCC-CM2-SR5, ACCESS-CM2, MPI-ESM1-2-HR, INM-CM4-8, INM-CM5-0, IPSL-CM6A-LR, HadGEM3-GC31-LL (r1i1p1f3), HadGEM3-GC31-MM (r1i1p1f3), UKESM1-0-LL (r1i1p1f2), MPI-ESM1-2-LR, MRI-ESM2-0, CESM2-WACCM (r2i1p1f1), NorESM2-LM, NorESM2-MM, KACE-1-0-G, GFDL-CM4, GFDL-ESM4. Four also have 3-hourly `uas`, `vas`, `rsds` in both: MPI-ESM1-2-HR, MPI-ESM1-2-LR, IPSL-CM6A-LR, MRI-ESM2-0. All grids are 1-D regular lat/lon; calendars include noleap, 365_day and 360_day, so dates are carried as (year, month, day) integers.
- Extraction (`extract_daily.py`): historical 1950-1980, ssp585 2040-2100, Sep-Mar days. Land = ERA5 land fraction >= 0.5 within half a model cell; land z0 = ERA5 fsr climatology as in Section 10c. Wind CF = onshore curve of the log-law 100 m speed of the daily-mean 10 m speed; solar CF = daily-mean rsds / 1000. Region means cos(lat)-weighted; a region with no land cells gives NaN (CanESM5 has 1 Korea cell). Warming from `Amon/tas`, global cos(lat)-weighted annual means, complete years.
- Bridges on the same footing: EERIE daily averages `mean10ws` and `msdwswrf` (GRIB, section 2 stripped); TCO95 3-hourly `10u`/`10v`/`ssrd` aggregated to days (mean speed over the 8 instants in (d 00, d+1 00], ssrd summed over the 8 windows ending in that interval).
- Events (`daily_events.py`): 48 h = 2 daily steps, window expansion, winters need >= 200 days (360-day calendars have 210). Comparisons: CMIP6 historical winters 1951-1980 against ssp585 2071-2100 and 2041-2070; EERIE and TCO95 as in Sections 10c-10d, plus TCO95 hi1950g 1953-1980 against 2071-2099 in CMIP6 style. Threshold matched to 5 events/winter in each comparison's control. Field means over the regions valid in every comparison (all 24 on the daily footing, both metrics). Null for the field mean: TCO95 daily constant-1950 control, all disjoint 30/30-winter placements, the same winters in every region. Ensemble statistics: CV of ln(ratio) against CV of ln(ratio)/dT, correlation and regression of ln(ratio) on dT across the 20 models.
- Run 2026-11-05: extraction jobs 28760189-28760196 (four CMIP6 batches, EERIE, TCO95), CanESM5 run on the login node as a test; analysis about 10 min. A first summary of Germany across models counted 21 entries because the TCO95 CMIP-style label also ends in "2071-2100"; filter on the `group` field instead.
- What sets the inter-model spread (2026-11-05, from `daily_events_summary.json`, wind only, 24-region field means, 2071-2100). This compares per K against the change in winter mean wind CF per K.
  - Across the 20 CMIP6 models the correlation is -0.90. A linear fit gives per K = 3.0 - 2.07 x (mean CF change, %/K). The median mean change is -1.9 %/K.
  - The intercept, +3.0 %/K, is the increase a model would show with no mean decline.
  - The bridges lose much less mean wind: -0.3 to -0.75 %/K (EERIE, TCO95). They sit above the CMIP6 line, at 5.4-6.0 %/K against a predicted 3.7-4.6. So in AWI-CM3 and EERIE the increase comes mostly from a broader low tail, while in CMIP6 it comes mostly from the mean decline.
  - Resolution pairs in CMIP6: HadGEM3-GC31 MM/LL +4.7/+5.1, MPI-ESM1-2 HR/LR +6.9/+7.2, NorESM2 MM/LM +9.9/+10.4 %/K. In each pair the finer grid is slightly lower, and the differences are small against the spread between models.
- IFS-family check, 2026-11-05: EC-Earth3 r1i1p1f1, whose atmosphere is IFS cy36r4.
  - Setup:
    - It has no `day/rsds`, so it runs wind only. `inventory.py --add-wind-only EC-Earth3` gives the entry `wind_only: true`, and `extract_daily.py` leaves solar NaN.
    - `daily_events.py` puts it in group `cmip6_windonly`. It is left out of the common-region sets and the ensemble statistics, so the 20-model numbers are unchanged.
    - `geostrophic.py` reports it under `cmip6_windonly`.
    - Jobs 28764548 (daily wind) and 28764549 (psl).
  - Result, wind only:
    - Field-mean ratio 1.81, +12.6 %/K at dT 4.99 K, 24 of 24 regions up, mean CF -16.9 %, 45 % of the increase from the mean.
    - It is the highest of all 21 CMIP6 models and lies on or above the line from the inter-model fit (predicted +10.0).
    - Mid-century: ratio 1.71 at dT 3.05, +19.3 %/K. Most of the change has arrived by 2041-2070.
  - Geostrophic, 22 regions: R_s 1.81, R_G 1.49, share from circulation 0.67, mean G -7.9 % (beyond the CMIP6 10th percentile of -6.9 %), transfer -0.2 %. Within the European regimes, G changes by -7.6 % (NAO+), -7.0 % (NAO-), -13.6 % (ScBL) and -10.5 % (AR). The transfer changes by -3.0 to +1.0 %.
  - So the small mean decline and the rising transfer of AWI-CM3 and EERIE (+1.6 to +4.1 %) do not come from an IFS atmosphere as such. The older IFS cycle in EC-Earth3 behaves like the rest of CMIP6 (flat transfer) and has the strongest gradient weakening. What remains in common to our runs: the newer IFS cycles (OpenIFS 43r3 in AWI-CM3, IFS 48r1 in EERIE), the FESOM2 ocean, and the AWI coupled setup.

### 14b. Sub-daily check with the 3-hourly CMIP6 models

- Models: MPI-ESM1-2-HR, MPI-ESM1-2-LR, IPSL-CM6A-LR, MRI-ESM2-0 (`3hr/uas`, `3hr/vas`, `3hr/rsds` in both experiments). Winters 1961-1990 against 2071-2100, because MRI-ESM2-0 3hr output starts in 1960. That differs from the 1951-1980 control of Section 14, so the daily footing is recomputed for the same winters.
- `uas`/`vas` are instants: 03-24 UTC (MPI, IPSL) or 00-21 UTC (MRI). `rsds` is the 3 h mean stamped at the window centre (01:30, 04:30, ...).
- Three footings from one extraction (`dunkelflaute/scripts/cmip6/extract_3hr.py`, output `dunkelflaute/data/cmip6/sub3h_<model>.npz`), same regions, land, z0 and CF curves as Section 14:
  - 6h: wind CF at the 00/06/12/18 UTC instants; solar = mean of the two 3 h windows ending at the instant (as `to_6h` for our runs); 48 h = 8 steps.
  - d3h: CF of the mean speed over the 8 instants in (d 00, d+1 00]; solar = mean of the 8 windows centred on day d. Days with fewer than 8 instants are dropped (one per winter, 31 March, whose last instant falls in April).
  - day: the `day/` files, as Section 14.
- Events and thresholds (`subdaily_events.py`, output `dunkelflaute/results/cmip6/subdaily_summary.json`): matched per footing, region and metric. Field means over the regions valid on every footing in every model: 24 for wind, 21 for comb.
- Result (wind, %/K on 6h / d3h / day):

  | model | 6h | d3h | day |
  |---|---|---|---|
  | IPSL-CM6A-LR | +7.5 | +7.9 | +7.6 |
  | MPI-ESM1-2-HR | +7.9 | +7.7 | +7.7 |
  | MPI-ESM1-2-LR | +7.1 | +7.0 | +7.2 |
  | MRI-ESM2-0 | +4.9 | +5.5 | +5.6 |

  - The regional patterns of ln(ratio) correlate at r = 0.91-0.99 between 6h and day.
  - Comb behaves the same way: within 0.6 %/K. MRI-ESM2-0 decreases on every footing (-2.2 to -2.8 %/K).
  - So the daily footing does not bias the CMIP6 per-K values, consistent with the EERIE and TCO95 bridges of Section 14.
- Run 2026-11-05: jobs 28762356-28762359, about 10 min each.

## 13b. Geostrophic against surface wind

- Question: is the increase in low-wind events a weakening of the large-scale pressure gradients, or a change in how much of the free-atmosphere flow reaches the surface (stability, ageostrophic flow)?
- G: daily geostrophic speed |grad p| / (rho f) from daily-mean msl on the 2.5 deg grid, rho = 1.25 kg m-3, centred differences. Region G = mean over the grid points inside the region box (DE: the DE_BOX; a box smaller than the grid takes points within half a cell).
  - Our runs: the regime msl of Section 13 (Euro-Atlantic only, so the 8 European regions).
  - CMIP6: `day/psl` bilinearly interpolated (periodic in longitude) to a global 2.5 deg grid, 15-80 N (`dunkelflaute/scripts/cmip6/extract_psl.py`, output `dunkelflaute/data/cmip6/psl_<model>.npz`). The grid's 20-80 N, 80 W-40 E part is the regime grid. Bilinear rather than box mean, because most CMIP6 grids are about as coarse as 2.5 deg. The 22 regions with a southern edge at or north of 29 N are kept; CN_E and IN_N are dropped because f is small there.
- s: daily mean of the region-mean speed. Our runs use the log-law 100 m speed; CMIP6 uses `sfcWind`. The log-law factor is fixed per cell and month, so it cancels in relative changes of the transfer.
- Transfer a_m = sum(s) / sum(G) over control days of calendar month m. Pseudo-surface speed s* = a_m G in both periods: the surface speed if the transfer kept its control value.
- Events in speed space on the daily footing (48 h, window expansion), each series with its own control-matched threshold (5 events/winter).
  - R_s = actual ratio; R_G = ratio for s* (circulation part); R_s / R_G = transfer part.
  - Field means of ln R; share from circulation = mean ln R_G / mean ln R_s.
- Within regimes (ERA5 EOFs and centroids, anomalies from the control seasonal cycle; CMIP6 days assigned the same way): per regime and region, the change in mean G, in the transfer mean(s)/mean(G), and in the share of days below the control 10th percentile of s and of s*.
- Script `dunkelflaute/scripts/regimes/geostrophic.py`, output `dunkelflaute/results/regimes/geostrophic_summary.json`; about 1 min on the login node.
- Correlation of daily s with G in the control: 0.88-0.91 (our runs), 0.84 median (CMIP6).
- Result, field means:

  | comparison | R_s | R_G | share from circulation | mean G | transfer |
  |---|---|---|---|---|---|
  | 9 km | 1.26 | 1.35 | 1.32 | -4.6% | +4.1% |
  | 31 km | 1.20 | 1.21 | 1.06 | -2.2% | +3.2% |
  | 100 km | 1.25 | 1.13 | 0.54 | -3.7% | +2.6% |
  | EERIE ctl | 1.07 | 1.08 | 1.10 | -2.2% | +2.3% |
  | EERIE hist | 1.04 | 1.07 | 1.55 | -1.4% | +1.6% |
  | CMIP6, 20 models, median [10-90%] | 1.34 [1.17, 1.54] | 1.29 [1.16, 1.48] | 0.91 [0.60, 1.22] | -3.6% [-6.9, -2.3] | 0.0% [-2.5, +1.3] |

  - CMIP6 per K: +7.2 %/K for s, +6.0 %/K for s*.
  - The increase is mainly circulation: weaker pressure gradients. The transfer does not explain it. In our runs the transfer rises by 2-4%, which partly offsets the weakening of G and is consistent with their smaller mean wind decline than CMIP6 (Section 14). In CMIP6 the transfer is flat.
- Within regimes:
  - G weakens in every regime in CMIP6: Europe medians -2.3% (NAO+ and NAO-), -4.5% (ScBL), -5.4% (AR).
  - In our runs the weakening is largest under Scandinavian blocking: -11.5% at 9 km, -7.6% at 31 km, -4.7% for EERIE ctl.
  - Low-s* days increase within every regime, matching the within-regime result of Section 13.
  - Correction (2026-11-06, figure check): G weakens in every regime in the CMIP6 medians and under ScBL and AR in every run, but not in all runs under NAO+/-: 31 km NAO- +3.2 %, 9 km NAO+ -0.1 %, and 3 of 20 CMIP6 models strengthen G under NAO+/-.
  - On the 8 European regions the CMIP6 transfer change is -0.9 % [-4.0, +1.6], not the 0.0 % of the 22 regions. On the CMIP6-style G (bilinear msl) our runs give +1.3 to +3.4 % on 8 regions; the box-mean G above gives +1.6 to +4.1 %.
- Limits: daily-mean msl at 2.5 deg resolves synoptic gradients only. The transfer term lumps stability, friction, ageostrophic flow and sub-2.5 deg pressure structure. Our runs' msl uses a box mean, CMIP6's a bilinear interpolation.

## 13c. Is the rising transfer an artifact?

- Why: in AWI-CM3 and EERIE the transfer of 13b rises (+1.6 to +4.1 %), while in CMIP6 it is flat. Five checks test whether that comes from how it is measured.
- Run 2026-11-06:
  - Scripts:
    - `dunkelflaute/scripts/regimes/transfer_checks_msl.py`: msl the CMIP6 way, bilinear, periodic in longitude, global 2.5 deg, 15-80 N. Output `dunkelflaute/data/regimes/transfer_checks_msl_<run>.npz`.
    - `dunkelflaute/scripts/regimes/transfer_checks_dart10m.py`: 9 km raw 10 m speed. Output `dunkelflaute/data/regions/transfer_checks_dart9km_<run>_spd10.npz`.
    - `dunkelflaute/scripts/regimes/transfer_checks.py`: the analysis.
  - Job 28773133 (`dunkelflaute/jobs/transfer_checks.sbatch`), 21 min, 26 GB, 16 workers. The analysis takes about 1.5 min on the login node.
  - Output `dunkelflaute/results/regimes/transfer_checks.json`, which holds every variant and the per-region values. It reproduces `geostrophic_summary.json` exactly.
- Same regions. All 8 European regions are in `CMIP_REGIONS`.
  - On these 8 regions alone, the CMIP6 transfer change is -0.9 % [10-90 %: -4.0, +1.6], and EC-Earth3 is -1.5 %.
  - The 9 km and 31 km rises (+4.1 %, +3.2 %) exceed every CMIP6 model. The 100 km rise (+2.6 %) equals the highest (CMCC-CM2-SR5, +2.5 %).
- G construction:
  - Box mean against bilinear interpolation changes the European transfer change by 0.06 percentage points at 100 km, and by 0.3-0.6 at 9 km, 31 km and EERIE.
  - Almost all of that is IB and BK. At 31 km the IB change in G is -4.4 % with the box mean and -1.7 % with bilinear interpolation, probably from the msl reduction over high terrain responding to warming.
  - Without IB, IT and BK the two constructions agree within 0.2 percentage points. The rise is then +1.0 to +2.8 %, and no CMIP6 model reaches it, except against EERIE hist (25 % of models).
- 22 regions with bilinear G, transfer change:

  | | 9 km | 31 km | 100 km | EERIE ctl | EERIE hist | CMIP6 median [10-90 %] | EC-Earth3 |
  |---|---|---|---|---|---|---|---|
  | transfer | +2.7 % | +2.3 % | +2.4 % | +1.3 % | +0.8 % | +0.0 % [-2.5, +1.3] | -0.2 % |
  | mean G | -4.6 % | -3.1 % | -3.9 % | -1.6 % | -1.2 % | -3.6 % [-6.9, -2.3] | -7.9 % |

  - The transfer rises in 18-20 of the 22 regions in our runs, against a CMIP6 median of 10. CN_N is the only region where it clearly falls in all our runs (-2.7 to -4.2 %).
  - Per K of global warming: +0.4 to +0.6 %/K in our runs, against 0.0 %/K [-0.5, +0.3] in CMIP6. Per K puts EERIE, with its smaller warming, on the same footing as the others.
- Wind quantity:
  - `spd100` is the uncalibrated log law with z0 fixed per cell and calendar month. QM feeds only the capacity factor, never `spd100`.
  - Replacing `spd100` by the raw 10 m speed changes the transfer change by at most 0.1 percentage points in every run.
  - The 9 km ratio of `spd100` to the 10 m speed is the same in 1950C and 2080C (1.88 in DE). TCO95 through the CMIP6 daily extractor (daily-mean 10 m speed, ERA5 land fraction) gives +2.4 % on the 22 regions, as the 3-hourly series does.
- Native 100 m wind (TCO95 `100u`/`100v`): +1.7 % on the European and on the 22 regions, against +2.4-2.6 % with the log law.
  - So about a third of the rise at 100 m is near-surface shear, which the neutral log law cannot carry. The rest is in the native 100 m wind. CMIP6 has only 10 m wind, so this part is not compared.
- Result: the rising transfer is not an artifact of region choice, G construction or wind quantity. In every variant it lies at or above the CMIP6 90th percentile, except EERIE hist on the 22 regions.
- Corrections to earlier notes:
  - The comment in `extract_psl.py` that box mean and bilinear interpolation agree closely holds for the msl field (daily rms 0.3-0.5 hPa), not for the change in G over IB, IT and BK at high resolution.
  - The transfer of 13b refers to the 10 m wind (the log-law factor cancels). At 100 m it is about a third smaller.
  - EC-Earth3 is an outlier in G, not in transfer.

## 13d. What raises the transfer in TCO95

- Data: the TCO95 rerun (ctl1950g, hi1950g, 1950-2099).
  - 3-hourly `10u`/`10v`, `100u`/`100v`, `2t`, `skt`, `blh`, `sshf`, and model levels 85-91.
  - Monthly `sd`, `fsr`.
  - Levels 85-90 lie 289, 215, 152, 99, 56 and 24 m above level 91.
- Region-day diagnostics: daily means of the eight 3 h stamps over the land cells and cos-lat weights of `extract_tco95.py`, for the 8 European regions, Sep-Mar. `extract_tco95.py` uses each run's own lsm >= 0.5; ERA5 supplies only z0. G is as in 13b.
- Scripts:
  - `dunkelflaute/scripts/regimes/transfer_stability_extract.py` (job script `dunkelflaute/jobs/transfer_stability_extract.sbatch`).
  - `transfer_stability.py` and `transfer_stability_cmip6.py`.
- Outputs:
  - `dunkelflaute/results/regimes/transfer_stability.json` and `transfer_stability_cmip6.json`.
  - Figures `transfer_stability_{classes,regionmonth}_{100km,hist_fut}.png`.
  - Intermediate files in `dunkelflaute/data/regimes/transfer_stability_parts/`.
- Run 2026-11-06: jobs 28773134 (ctl1950g, 15 min) and 28773136 (hi1950g, 10 min), 16 workers each.
- Data notes:
  - The 3 h `sshf` is accumulated over the 3 h ending at the stamp, in J m-2, so divide by 10800 s. It is positive downward. The monthly `sshf` is the mean of the 3 h accumulations.
  - `sd` is in m water equivalent.
  - Over land, OpenIFS 43r3 derives the 10 m wind from the 40 m wind with z0 = 0.03 m and Monin-Obukhov stability (`surfexcdriver_ctl_mod.F90`, `sppcfl_mod.F90`). So the 10 m wind is 9 % above the lowest-level wind.
- Baseline: rebuilt from the raw files, the 100km comparison (ctl1950g 1953-2099 against hi1950g 2081-2099) gives +2.56 %, identical to 13b.
  - Raw 10 m speed: +2.45 %. Native 100 m: +1.68 %.
  - hi1950g 1951-1980 against 2071-2099: +2.28 % (G -3.36 %).
  - Control drift (ctl1950g 1953-1980 against 2071-2099): 0.00 %.
- Method: classes are control quintiles per region of `sshf`, `blh`, 2t - skt, the inversion theta(level 85, about 300 m) - theta(2 m), and G. Joint classes are 5 x 5 pairs.
  - Between-class part: future class frequencies with the control mean transfer per class.
  - Within-class part: the rest.
  - The order of the split changes the parts by less than 0.2 points. The 90 % ranges come from a winter bootstrap.
- Result, transfer change of the log-law 100 m series (`transfer_stability.py`; +2.56 %, against +2.45 % for raw 10 m), mean over 8 regions:
  - Weaker G alone (G quintiles): between-class +1.0 % [0.7, 1.4]. s/G falls from 0.90 in the weakest G quintile to 0.62 in the strongest, so fewer strong-gradient days raise the ratio.
  - G x inversion: between +2.65 % [2.33, 3.00], within -0.09 % [-0.50, 0.32].
  - G x blh: +2.71 % and -0.15 %.
  - So weaker gradients and a less stable boundary layer account for the whole rise: about +1.0 % from G and +1.6 % from stability.
  - The hist_fut sample agrees: G x inversion gives +2.39 / -0.10.
- The destabilization is not a response to the wind. At fixed G and fixed 10 m speed (25 joint classes), the inversion still weakens by 0.30 K (raw 0.31 K), blh deepens by 26 m (+5 %), and 2t - skt falls by 0.13 K (raw 0.15 K).
- The bulk Richardson number from region-day means is circular: its shear denominator scales with the wind. It is not used.
- Height: of the +2.45 % in s10/G, +1.68 % is already in V(level 85)/G, +0.41 % between level 85 and 40 m, and +0.35 % in the 10 m diagnostic.
  - So the whole boundary layer couples more strongly to the free flow; the change is not only in the surface layer.
  - V85/G is fully explained by G x blh classes (between +1.71, within -0.03).
  - Native 100 m / 10 m falls from 1.656 to 1.643 (-0.75 %), less near-surface shear.
- Snow, across 56 region-months:
  - The transfer change correlates with the snow-cover change at r = -0.42, and -0.57 with the G change removed. On G and snow cover together, R^2 = 0.45.
  - Region-months with control snow cover above 0.1 rise by +3.4 %, the others by +1.4 %. The largest are BK January +7.3 % and DE March +5.9 %.
  - Snow is an amplifier, not the only cause. IB (+2.0 %) and FR (+2.6 %) are nearly snow-free and still destabilize (inversion -0.27 and -0.24 K at fixed G and s10).
- Roughness: `fsr` rises by 1.3 % (0.652 to 0.660 m) where snow is lost. That is the wrong sign to raise the wind, and the land 10 m diagnostic uses a fixed z0 anyway.
- CMIP6 with the same G split (8 European regions, 20 models):
  - The between-G part is +0.86 % (median; 10-90 % +0.51 to +2.14), about the same as TCO95's +1.0 %.
  - The at-fixed-G part is -2.0 % (-5.0 to +0.85). EC-Earth3: -4.3 %. TCO95: +1.5 %.
  - So the model difference is entirely in the transfer at fixed pressure gradient. Why the CMIP6 term is negative is the subject of 13e.
- Historical period (hi1950g 1951-1980 against 1991-2020): transfer +0.07 % (bootstrap +/- 0.5 %), G -1.2 %, inversion -0.05 K. No transfer signal is detectable in the period with observations, so the observed trend cannot constrain it. An observational test has to use the interannual sensitivity instead.

## 13e. Why the CMIP6 transfer falls at fixed G

- Question: the at-fixed-G term (13d) is negative in CMIP6 and positive in our runs. Is that stability, or something else?
- Monthly winter (Sep-Mar) diagnostics over the region land cells (the land mask of `extract_daily.py`), historical 1951-1980 against ssp585 2071-2100:
  - stability theta(850 hPa) - theta(tas, ps), theta = T (1e5/p)^0.2857, on cells whose control winter-mean ps is at least 900 hPa (removes most of TR, US_CA and parts of CN_N and IN_N);
  - ta850 - tas, `hfss` (positive upward), `snc`, `lai`, `treeFrac`, `cropFrac`.
  - Use the theta version: it is about 0.2 K less negative than ta850 - tas, because theta(850) warms 1.0475 times T850.
- Scripts in `dunkelflaute/scripts/regimes/`: `stability_cmip6_extract.py`, `stability_cmip6_tco95.py`, `stability_cmip6_eerie.py`, `stability_cmip6_split22.py`, `stability_cmip6_split_ours.py`, `stability_cmip6_analysis.py`.
- Outputs in `dunkelflaute/results/regimes/`: `stability_cmip6.json`, `stability_cmip6_split22.json`, `stability_cmip6_split_ours.json`, `stability_cmip6_scatter_{eu8,cmip22}.png`. Intermediate files in `dunkelflaute/data/regimes/stability_cmip6_parts/`.
- Run 2026-11-06: CMIP6 extraction 3 min on the login node (6 workers), TCO95 1 min, EERIE job 28775770 (33 s).
- Data notes:
  - TCO95 stores no surface pressure or orography. ps is estimated from msl, 2t and ERA5 orography (`era5_z_on_tco95.nc`, conservative remap with cdo). A version built from z850 that needs no ps agrees within 0.01 K.
  - TCO95 monthly `sshf` is the mean of 3 h accumulations, positive downward: hfss = -sshf/10800 s. EERIE `avg_ishf` is positive downward and is negated.
  - Snow cover is not the same quantity everywhere. TCO95: the IFS depth formula on monthly mean `sd` and `rsn`. EERIE has no snow density, so 250 kg m-3 is assumed (cover = min(1, sd/0.025 m)). CMIP6 `snc` is the monthly mean fraction.
  - HadGEM3 and UKESM `ta` carry a `missing_value` that does not cast to the data type; values above 1e15 are masked.
  - Coverage: `snc` for 14 of 20 models plus EC-Earth3; land-cover fields for 10-12. KACE-1-0-G has no `snc`, `lai` or `sftlf`; EC-Earth3 has no Lmon vegetation on /pool. No model outputs roughness length.
  - The 9 km and 31 km runs have no monthly temperature at 850 hPa (TCo319 stores only `pl_z_850`), so they have no stability values.
  - Our runs have zero land-cover change: IFS uses fixed land cover and a monthly LAI climatology. Checked for OpenIFS 43r3.
  - Checked for EERIE IFS 48r1 (2026-11-07):
    - cvh, cvl, tvh, tvl, slt, cl, dl, sdor and sdfor are bit-identical between the 1950 initial file `ifs_inputs/tco1279l137/i2hj/LWDA/1950010100/gfc/ICMGGi2hjINIT` and the 2017 model states `Cont1950_ini_Storyline/ICMGGi2hj+000384` and `Hist_SSP245_ini_StoryLine/ICMGGi2hj+000384` (max abs difference 0 at all 6,599,680 points).
    - The yearly climate files `ifsINIT/ICMCL_1279_*` copy one template year. On 15 January, LAI and albedo are identical between 1950, 1972, 1983, 2014, 2033 and 2050.
    - `arpifs/climate/updclie.F90` updates only albedo, LAI, SST and sea ice during a run, never the cover fields.
    - The vegetation fields are not in the EERIE output (kerchunk stores and FDB indices). The climate files on disk stop at 2052, so 2052-2070 rests on the code and the input structure.
    - Scratch: scratchpad/eerie_lc/initcmp.py and initcmp1.txt.
- At-fixed-G term on the CMIP6 footing (10 m speed, G from msl interpolated bilinearly on the global grid):
  - CMIP6, 22 regions (new): -1.38 % [10-90 %: -3.76, +0.11]. The 8-region rerun reproduces `transfer_stability_cmip6.json` exactly (-2.00 % [-5.00, +0.85]).
  - Ours, 8 European regions: +1.67 % (9 km), +2.02 % (31 km), +1.42 % (100 km), +1.38 % (EERIE ctl), +0.97 % (EERIE hist).
  - Ours, 22 regions: +0.53, +0.80, +0.79, +0.69, +0.28 %. At or above the highest CMIP6 model (+0.54 %).
- Stability does not explain the spread:
  - CMIP6 destabilises too. dtheta on 22 regions -0.42 K [-0.81, -0.10], negative in 19 of 20 models; on the 8 European regions -0.02 K [-0.45, +0.34].
  - TCO95 (-0.15 K on 8, -0.55 K on 22) and EERIE (-0.11, -0.24 K) lie inside the CMIP6 range. So do their changes in hfss (+1.1 to +2.4 W m-2), snow cover, and dtas/dta850 (1.07 to 1.16 in all three).
  - Across 20 models the at-fixed-G term against dtheta gives r = +0.18 (8 regions) and +0.24 (22). The sign is the wrong one for stability control. hfss r = +0.15 / -0.10, snow cover r = +0.16 / +0.14 (n = 14), warming ratio r = -0.14 / -0.10.
  - Pooled over model and region (440 points, region means removed): dtheta r = +0.03, hfss r = -0.04.
  - EC-Earth3 destabilises most (-0.78 K, snow cover -40 %) and has the third most negative term (-4.3 %).
- Land-cover change is the best predictor:
  - Crop fraction change: r = +0.89 (8 regions), +0.81 (22), n = 10; partial r = +0.83 given dtas and dG.
  - Tree fraction change: r = -0.67 (8), -0.38 (22), n = 11. LAI change: r = -0.33 (8), -0.68 (22), n = 12.
  - Within a model, regions with more LAI or tree gain have a more negative term (pooled, model means removed: r = -0.52 for LAI, -0.53 for tree fraction; crop +0.24).
  - INM-CM4-8 and INM-CM5-0 have essentially no land-cover change and terms of -0.1 to -0.3 %.
  - In one model family: AWI-ESM-1-REcoM (crop -2.1 %) has +0.4 %; MPI-ESM1-2-LR and -HR (crop -4.5 and -7.2 %) have -2.4 and -5.2 %.
  - The size is plausible: turning about 7 % of the area from crop (z0 about 0.1 m) into forest (about 1 m) lowers the 10 m to G ratio by roughly 4 %.
  - Sea ice and coast play no clear role: model land fraction in the region cells is 85-99 % (r = +0.20 / +0.39).
- Our runs against the CMIP6 relation at zero land-cover change:
  - Predicted: +0.1 % (8 regions) and -0.7 % (22) from the crop fit; -1.1 % from tree fraction; -1.3 / -0.6 % from LAI.
  - Our runs lie 0.9 to 1.9 points above the crop fit on 8 regions and 1.0 to 1.5 on 22, about 1.4 to 2 residual standard deviations. Against the stability fits they are 3.1 to 3.6 points above at the same dtheta.
  - Snow with dtas and dG controlled: partial r = -0.28 / -0.50 (n = 14). That is the TCO95 sign (more snow loss, higher term), but weak.
- Interpretation: the model difference has two parts.
  - (Revised by the direct LUMIP test in 13g: the land-use part is about 1.4 points, not 2, because the cross-model crop slope overstates the per-pp effect about two-fold.)
  - About 2 points: the CMIP6 at-fixed-G term is pulled down by surface roughening from prescribed SSP5-8.5 land use and vegetation change (cropland loss, more trees, higher LAI). IFS keeps land cover fixed. This is a forcing difference, not a difference in climate physics, and it acts on the 10 m wind that the CMIP6 pipeline extrapolates to 100 m; at hub height the roughness effect would be smaller.
  - About 1 to 1.5 points: at zero land-cover change, CMIP6 transfer still does not respond to the destabilisation it also has, while IFS does (13d, about +1.5 %).
- Limits:
  - Land-cover output exists for only 10-12 models. Two of the three most negative, KACE-1-0-G (-8.8 %) and EC-Earth3 (-4.3 %; MPI-ESM1-2-HR is -5.2 %), cannot be tested. No roughness output, so the roughness step is inferred, not measured.
  - Crop and tree fraction changes are collinear; the regression cannot separate them.
  - Monthly 850 hPa stability is coarser than the daily 300 m inversion of 13d.

## 13f. Interannual sensitivity and stations: no observational constraint

- Question: in warm against cold winters, does the transfer at fixed G rise, as in our runs under warming, or fall, as in the CMIP6 forced change? And does each model's interannual sensitivity predict its forced change? (13d showed no detectable trend in the observed period, so the trend cannot be used.)
- Winter transfer: ln(sum_k mean s_k / sum_k mean G_k) over the 5 G-quintile classes k of each region, so the G mix is held at climatology.
- Sensitivity: the coefficient of winter regional land 2 m temperature in a regression of the winter transfer on that temperature and ln(winter-mean G). All series linearly detrended. In % per K, averaged over the 8 European regions; 90 % ranges from a 1000-draw winter bootstrap.
- Forced comparison: the within-G-quintile transfer change of `transfer_stability_cmip6.py`, divided by the regional winter land warming over the same windows. The CMIP6 values reproduce `transfer_stability_cmip6.json` exactly.
- Scripts in `dunkelflaute/scripts/regimes/`: `interannual_hadisd_download.py` and `interannual_hadisd_process.py` (run with `python3 -I`), `interannual_extract_cmip6.py`, `interannual_extract_ours.py`, `interannual_extract_era5.py`, `interannual_extract_era5_wind.py`, `interannual_transfer.py`. Job scripts `dunkelflaute/jobs/interannual_extract_ours.sbatch`, `interannual_era5_wind.sbatch`.
- Outputs: `dunkelflaute/results/regimes/interannual_transfer.json`, `interannual_transfer_scatter.png`. Intermediate data in `dunkelflaute/data/regimes/interannual_cmip6/`, `interannual_ours/`, `interannual_era5_{t,wind}.npz`; stations in `dunkelflaute/data/obs/` (`hadisd_raw/` 3.8 GB, `hadisd_selection.json`, `hadisd_daily.npz`).
- Run 2026-11-06: job 28775425 (CMIP6 `day/sfcWind`, `day/psl`, `Amon/tas` for historical winters 1951-2014, 21 models, 8 min), 28775530 (31 km and 9 km 2t, 2 min), 28775965 (ERA5 hourly 10 m wind, 1.5 min). Analysis 2.5 min on the login node. HadISD download about 3 min.
- Windows: CMIP6 interannual on historical 1951-2014 (64 winters; 1951-1980 gives similar values). 100 km interannual on ctl1950g (147 winters); forced on hi1950g 1951-1980 against 2071-2099 (the 13d pairing gives +0.28 %/K instead of +0.31). 31 km and 9 km have only 17 winters each.
- Models:
  - Ours, interannual: +0.57 %/K [+0.41, +0.73] (100 km), +0.65 (31 km), +0.88 (9 km), +0.65 (EERIE ctl). Forced at fixed G: +0.28 to +0.54 % per K of regional warming.
  - CMIP6, interannual: median +0.27 %/K [10-90 %: -0.06, +0.69]; forced -0.40 %/K [-0.78, +0.16]. EC-Earth3: +0.48 interannual, -0.56 forced.
  - Almost every model has a positive interannual sensitivity, whatever the sign of its forced change.
- Cross-model relation: interannual against forced, r = -0.27 over 21 CMIP6 models (permutation p = 0.25; Spearman -0.15), r = +0.02 over all 25, pooled models x regions r = +0.07. Variants (daily residual model, direction controls, quadratic detrend, first differences) give r = -0.35 to -0.11. **The interannual sensitivity is not an emergent constraint on the forced transfer change.**
- Observations: HadISD v3.4.3.2025f (Dunn et al. 2012, 2014, 2016; Dunn 2019; Smith et al. 2011), winters 1980-2025 (46 winters).
  - 423 wind stations (DE 31, FR 101, IB 39, UK 92, PL 35, SC 47, IT 48, BK 54) and 430 temperature stations. Criteria: at least 80 % valid days in at least 30 winters, elevation below 800 m.
  - Daily mean of 00, 06, 12, 18 UTC, all four required. Each station's speed is divided by its calendar-month mean before the region mean. G from ERA5 daily msl.
  - Observed sensitivity: -0.14 %/K [-0.94, +0.60] with station temperature, -0.15 with ERA5 t2m, -0.22 per-station detrended. Quadratic detrend +0.04, first differences +0.20 [-0.43, +0.71], with direction controls -0.85 [-1.95, +0.15].
  - The regions disagree: UK -1.5, IT -0.9, DE -0.8, SC -0.8 against BK +1.6, PL +0.8, FR +0.5.
- ERA5 10 m wind through the same methods:
  - Region means: +0.76 %/K [+0.46, +1.03], like our runs. Sampled at the stations on the station days: +0.48 [+0.17, +0.76], so station sampling alone lowers the value by about 0.3.
  - Stations minus ERA5 at the same points, paired: -0.62 %/K [-1.30, -0.02].
- ERA5 inversion (t850 - t2m) as regressor: -0.75 % per K of inversion [-2.9, +1.5], uninformative (inversion and T correlate at only r = -0.25).
- Interpretation:
  - Observed warm winters lack the rise in transfer at fixed G that ERA5 and our runs show. On the like-for-like test (stations against ERA5 at the same points) the difference is significant at about the 90 % level. The observations sit with the lowest CMIP6 models (CESM2-WACCM, FGOALS-g3, CMCC-CM2-SR5, ACCESS-CM2).
  - This points to a process discrepancy in IFS-type boundary layers on interannual time scales. Because interannual and forced sensitivities are unrelated across models, it does not decide between our runs and CMIP6 for the forced change.
- Limits:
  - 46 observed winters, strong regional disagreement.
  - Station transfer anomalies correlate only r = 0.37 with ERA5 at the same points. Exposure changes and instrument steps are not removed by a linear detrend, although quadratic detrending and first differences do not close the gap.
  - Winter temperature is about 56 % explained by circulation, so the T regressor mixes advection and stability.
  - The 9 km and 31 km samples are short, and the 1950c control drifts.

## 13g. Direct test of the land-use attribution with LUMIP

- Question: do experiments that differ only in land use reproduce the cross-model land-use attribution of the CMIP6 at-fixed-G term (13e)?
- Pairs, land-use effect = A minus B: historical minus hist-noLu; ssp370 minus ssp370-ssp126Lu; ssp126-ssp370Lu minus ssp126; esm-ssp585 minus esm-ssp585-ssp126Lu.
  - The two ssp370/ssp126 pairs give SSP3 minus SSP1 land use; the esm585 pair gives SSP5 minus SSP1 land use under emission-driven SSP5-8.5.
- Windows: late 1985-2014 and early 1951-1980 for the historical pair; late 2071-2100 and early 2016-2035 for the scenario pairs.
- Winter = Sep-Mar, equal month weights, labelled by its January. Region land cells as in `stability_cmip6_extract.Grid`.
- Wind: Amon/sfcWind (monthly means), because daily sfcWind is missing in most LUMIP experiments.
- Control members are restricted to the i/p/f variant of the LUMIP members. CanESM5 LUMIP runs are p2, so only p2 controls are used. The UKESM f3 historical member is excluded.
- Monthly G: |grad p|/(rho f) of the monthly-mean psl, bilinear to the global 2.5 deg grid, region boxes as in `geostrophic.py`. Used only as a circulation control.
- At-fixed-circulation effect: 100 (ln s_A/s_B - b ln G_A/G_B), with b the per-region regression slope of monthly ln s on ln G anomalies. b is 0.16-0.22. The s/G version (b = 1) is kept as a variant.
- Uncertainty: member spread of the single-member window mean, pooled over both experiments, se = sigma sqrt(1/n_A + 1/n_B). Where both experiments have one member, the interannual spread divided by sqrt(W). Also a winter bootstrap and differences between disjoint member groups of the same experiment.
- Scripts in `dunkelflaute/scripts/regimes/`: `lumip_inventory.py`, `lumip_extract_monthly.py`, `lumip_extract_daily.py`, `lumip_analysis.py`. Job script `dunkelflaute/jobs/lumip_extract_daily.sbatch`.
- Outputs: `dunkelflaute/results/regimes/lumip_inventory.json`, `lumip_landuse.json`, `lumip_landuse.png`; console log `dunkelflaute/logs/lumip_analysis_login.log`. Intermediate files in `dunkelflaute/data/regimes/lumip_parts/` (94 MB).
- Run 2026-11-06: inventory 37 s and monthly extraction 104 s (16 workers) on the login node; daily extraction job 28777216 (1 min 42 s, 22 GB read); analysis 17 s.
- Data notes:
  - No LUMIP experiment on /pool has a roughness length variable.
  - CanESM5 and GFDL-ESM4 hist-noLu have no sfcWind, uas or vas. CMCC-ESM2, MIROC-ES2L and NorESM2-LM hist-noLu hold only Lmon/nbp.
  - ssp370-ssp126Lu has no wind for CNRM-ESM2-1, IPSL-CM6A-LR, MIROC-ES2L or GFDL-ESM4.
  - BCC-CSM2-MR outputs lai but no crop or tree fraction. IPSL hist-noLu daily sfcWind ends in 2009.
- Usable pairs: historical for BCC, CESM2, CNRM, IPSL, MPI-ESM1-2-LR and UKESM; ssp370 and ssp126 for BCC, CanESM5, MPI-ESM1-2-LR and UKESM; esm585 for BCC and MPI-ESM1-2-LR.
- Land cover in Europe, SSP3 minus SSP1: crop -6.3 to -7.6 pp in CanESM5, MPI-ESM1-2-LR and UKESM. Tree fraction rises by 5.4 pp (CanESM5) and 2.3-2.6 pp (MPI-ESM1-2-LR); in UKESM it falls by 0.8-1.1 pp and grass and pasture take up the crop. This is the direction of the CMIP6 ssp585 change (median crop -5 pp, tree +4 pp).
- Result, 8 European regions, late window, at fixed circulation:
  - ssp370 pair: -2.40 % (BCC), -2.52 % (CanESM5), -2.82 % (MPI-ESM1-2-LR), +0.47 % (UKESM).
  - ssp126 pair: -3.27 % (BCC), -2.23 % (CanESM5), -3.26 % (MPI-ESM1-2-LR), +1.00 % (UKESM).
  - Standard errors 0.5-0.8 %. Member groups of the same experiment differ with a standard deviation of 0.4-1.2 %.
  - esm585 pair: -0.68 % (BCC) and +0.52 % (MPI-ESM1-2-LR), not significant; the European land-cover difference is small (MPI: crop -2.1, tree -1.2, pasture +6.2 pp).
  - Historical pair: late +0.54 % mean (-1.95 to +2.72) at crop +2.3 to +4.7 pp and mixed tree change; early +2.35 % (+0.62 to +3.98) at crop +6.1 to +7.9 pp. So European land-cover change from 1951-1980 to 1985-2014 lowers the 10 m wind by 1.8 % on average (-4.3 to +0.2).
  - Monthly G differs by less than 2 % except CNRM historical (+3.1 %) and BCC ssp126 (-2.6 %). Removing the circulation part changes the effects by 0.5 points or less.
- 22 regions: the historical pair gives +3.1 to +10.7 %, from deforestation outside Europe (crop +12 to +14 pp, tree -4 to -13 pp). The scenario pairs give -0.7 to +0.7 %.
- Per unit land cover, weighted fit through the origin over the 12 late-window points of 6 models:
  - crop +0.25 +/- 0.11 %/pp (8 regions), against the cross-model slope of +0.60 in 13e;
  - tree -0.33 +/- 0.20 %/pp, against the cross-model -0.325;
  - with the historical early windows: crop +0.31 +/- 0.07, tree -0.45 +/- 0.15;
  - 22 regions: crop +0.34 against +0.46, tree -0.67 against -0.29.
  - The fits scatter more than internal variability allows (reduced chi2 10-14): the effect depends on what the crop becomes (forest in CanESM5 and MPI, grass in UKESM). Ranges come from a bootstrap over LUMIP models with scatter-inflated errors.
  - LAI cannot be used: it differs by less than 0.2 between LUMIP pairs, while the CMIP6 ssp585 LAI rise is driven by CO2 and climate.
- Land-use part of the CMIP6 at-fixed-G term, 8 regions, from each model's ssp585 crop and tree change:
  - crop and tree together: -1.46 % [5-95 % bootstrap -2.11, -0.41], 60 % [17, 86] of those 10 models' mean (-2.44 %);
  - all variants (crop only, tree only, with historical early windows, raw wind, s/G): -1.0 to -1.9 %.
  - So about -1.3 of the 20-model median of -2.0 % is land use, if the models with land-cover output are representative.
  - At zero land-cover change the CMIP6 term is about -1.0 % (-0.5 to -1.45). That puts our runs 2.4-3.0 points above CMIP6 (2.0-3.0 including IFS-FESOM2-SR hist), not 0.9-1.9 as in 13e.
  - 22 regions: land-use part -1.06 % [-1.78, -0.94], 59 % of -1.78 %.
  - Own-model check with the scenario pairs: about -1.7 of the -2.36 % of MPI-ESM1-2-LR and about -2.0 to -2.4 of the -4.97 % of CanESM5 are land use. The UKESM historical pair accounts for about all of its -3.07 %.
- Daily footing (day/sfcWind in both experiments, log law with the fixed ERA5 z0, onshore curve, 48 h events, threshold set to 5 per winter in experiment B):
  - mean CF, historical minus hist-noLu: +6.7 % (CESM2), +5.1 % (UKESM), +3.9 % (BCC), -4.1 % (IPSL, 1985-2009);
  - event ratio: 0.82 [0.79, 0.87] (CESM2), 0.94 [0.90, 0.98] (UKESM), 0.92 [0.84, 1.01] (BCC), 1.035 [0.99, 1.09] (IPSL);
  - UKESM ssp370 and ssp126 pairs: event ratios 0.985 and 1.01, not significant.
  - Events change by about -4.5 % per +1 % of wind; mean-removed ratios are 0.91-1.02, so the event change is mostly the mean shift. At that rate the -1.4 % land-use part means roughly 6-7 % more events.
  - The fixed z0 in the 10 m to 100 m step does not see the land-use roughening, so these CF and event effects at hub height are overstated.
- Interpretation: land-cover change lowers the CMIP6 10 m winter wind at fixed circulation, as inferred in 13e, but by about half the cross-model crop estimate. About 1.4 of the 3.4-4.0 point model difference on 8 regions is land use; the rest is not. (13h then shows that most of the remaining 10 m difference is also confined to the surface layer: at about 100 m the CMIP6 and IFS transfer changes at fixed G are close.)
- Limits:
  - Six models carry the crop and tree fits. Every scenario pair rests on one LUMIP member, except UKESM ssp370 with two.
  - Monthly-mean G is only a circulation control.
  - CMIP6 tree change includes dynamic vegetation, so the attribution is to land-cover change of any cause, none of which IFS has.
  - CESM2 stands in for CESM2-WACCM in the own-model check.

## 13h. Hub-height test of the CMIP6 wind decline

- Question: is the CMIP6 decline of the 10 m wind over Europe smaller at hub height, as a near-surface roughness effect of land-cover change (13e, 13g) would make it?
- Scripts: `dunkelflaute/scripts/regimes/hubheight_extract.py`, `hubheight_analysis.py`; job script `dunkelflaute/jobs/hubheight_extract.sbatch`.
- Output: `dunkelflaute/results/regimes/hubheight_summary.json`. Intermediate files in `dunkelflaute/data/regimes/hubheight_parts/` (307 MB, can be deleted; the summary holds everything).
- Run 2026-11-06: job 28777371 (6 array tasks, 6-31 min each, about 1.6 node-h), job 28777752 (6 min, rerun of ACCESS-CM2 and HadGEM3-GC31-MM after adding the remap); analysis 1.7 min on the login node.
- Periods, winters and land cells as in Section 14: historical 1951-1980 against ssp585 2071-2100, Sep-Mar, region land cells of `extract_daily.py`.
- Inputs:
  - No CMIP6 model on /pool has a 100 m wind variable.
  - Daily 1000 hPa wind from `day/ua,va` (plev8) for 19 models; GFDL-ESM4 uses `Eday`.
  - Daily 925 hPa wind from `Eday/ua,va` (plev19) for 10 models. MRI-ESM2-0 and GFDL-CM4 have historical Eday only for 1961-2000; HadGEM3-GC31-MM has no ssp585 Eday.
  - Daily model-level wind from `CFday` for 4 models (AWI-ESM-1-REcoM, MPI-ESM1-2-HR, MPI-ESM1-2-LR, NorESM2-LM), interpolated linearly in ln z to 100 m above ground. Lowest levels at about 31 and 141 m (MPI, AWI) and 61 and 194 m (NorESM2-LM). IPSL historical `CFday` covers only the 1890s, 1910s and 1990s and is not used.
  - Not used: `6hrLev` (too large), `6hrPlevPt` (starts at 925 hPa), Amon plev19 (monthly mean vector understates mean speed).
- Height of 1000 hPa above ground: hypsometric, with daily ps = day/psl x (Amon ps / Amon psl). Agrees with daily zg1000 minus orog within 1-5 m.
  - Over all control winter cell-days 1000 hPa is below ground in most of Europe: DE -110 to -255 m, FR -130 to -265 m, UK -70 to +30 m.
  - It is at least 50 m above ground on 0.17-0.47 of UK cell-days, about 0.1-0.3 in DE, FR, PL, SC, and at most 0.13 in IB, IT, BK. On those cell-days it lies 60-200 m above ground, typically 100-140 m.
  - 925 hPa lies 330-665 m above ground, well above hub height.
- Masks: V1000 quantities use only cell-days with 1000 hPa at least 50 m above ground, and only regions with at least 5 % such cell-days in both periods. Regions valid in every model: DE, FR, UK, PL, SC (8-region set), plus UA, RU, KZ, CN_N (22-region set). sfcWind, G, V925 and model-level V100 use all regions.
  - MPI, AWI, INM, NorESM, CanESM5, CMCC and EC-Earth3 extrapolate winds below ground, so the mask is needed. HadGEM3, UKESM and ACCESS-CM2 additionally mask 9-26 % of the cell-days counted as valid; FGOALS and CESM2-WACCM about 12 %. Those days are dropped for both s10 and V1000.
- Data notes:
  - HadGEM3, UKESM and ACCESS-CM2 u and v are on staggered points and are remapped bilinearly onto the sfcWind cells.
  - KACE-1-0-G ssp585 `day/ua` is corrupt (mean |u| 14-18 m/s, max 108 m/s at the levels labelled 1000 and 850 hPa); its pressure-level values are excluded, its sfcWind and G kept.
  - V100pl: power law between sfcWind and V1000 to 100 m on period means. It removes the +1.9 m median rise of 1000 hPa under warming.
  - |V| comes from daily-mean components while sfcWind is a daily-mean speed. The vector version |uas, vas| gives s10/V1000 -2.1 % instead of -2.5 % on the common regions.
  - Check: s10/G on all days reproduces `transfer_stability_cmip6.json` for each model except FGOALS-g3 (-1.08 against -1.25 %), which moves the median from -0.88 to -0.79 %.
- Result, 8 regions, CMIP6 medians [10-90 %]:
  - sfcWind -5.1 % [-8.8, -1.5]; V1000 -2.2 % [-4.6, -0.5]; V925 -1.4 % [-2.7, +0.6] (10 models); model-level V100 -2.7 % (4 models); G -3.8 % [-6.9, -1.9].
  - s10/V1000 -2.4 % [-4.7, 0.0]; V1000/G +2.0 % [+0.1, +3.5], positive in 18 of 19 models.
  - On the five common regions alone: sfcWind -4.1 %, V100pl -2.1 %.
  - 22 regions: sfcWind -4.4 %, V1000 -1.8 %, s10/V1000 -2.4 %, V1000/G +2.0 %.
- Model-level check: model-level and power-law 100 m agree within 0.2 points in MPI-ESM1-2-LR and NorESM2-LM, about 1 point in MPI-ESM1-2-HR and AWI. In these four, 100 m keeps 11-64 % of the 10 m decline (AWI -0.2 against -1.8 %, MPI-ESM1-2-HR -4.4 against -8.4 %).
- EC-Earth3: sfcWind -10.6 %, V1000 -9.8 %, s10/V1000 -3.2 %, V1000/G +0.9 %. Its decline is in the circulation, not the near-surface layer.
- Land cover, across models (8 regions):
  - s10/V1000 against crop fraction change r = +0.83 (n = 10; 22 regions +0.70), against tree fraction change r = -0.86 (n = 11). s10/V925 against crop r = +0.91 (n = 6).
  - V1000/G against crop r = +0.03 (22 regions -0.21); V925/G against crop r = -0.30.
  - The 13e at-fixed-G term correlates with s10/V1000 at r = +0.74 and with V1000/G at r = +0.30 (n = 19).
  - At zero crop change the fit leaves s10/V1000 at -0.4 % (8 regions) and -1.5 % (22).
- TCO95, 8 regions, native 100 m: s10/V100 +0.76 %, V100/G +1.68 %, V(level 85)/G +1.68 %. Its 100 m wind declines slightly more than its 10 m wind, the opposite of CMIP6.
- Interpretation:
  - The land-cover signal is confined to the layer below about 100 m. Above it the CMIP6 change of V/G resembles the IFS runs (+2.0 % against +1.7 %).
  - So the CMIP6 transfer drop at fixed G (13d, 13e) is a near-surface shear increase, and the pipeline's log-law extrapolation with fixed z0 carries it up to hub height. That roughly doubles the CMIP6 hub-height decline over Europe, by about 2 percentage points.
  - This revises 13g's reading that the physics gap to our runs is 2.4-3.0 points. Measured at about 100 m, the CMIP6 transfer at fixed G (+2.0 %) and ours (+1.7 %) are close; the remaining 10 m difference is the surface layer.
- Jan 2026-11-07: the crude per-K event consequence below is dropped from the paper; only the direct 10 m vs ~100 m wind result is used.
- Crude event consequence (labelled crude; Section 14 daily footing, which uses 10 m wind for our runs too):
  - CF elasticity 2.02 to wind speed (regression through the origin of the 24-region CF change on the 22-region sfcWind change across models). Shift in per K = 2.07 x 2.02 x d(s10/V_hub) / dT, using the Section 14 slope.
  - With the 22-region shear the CMIP6 median per K falls from +7.5 to +5.2 (V1000), +5.6 (V100pl) or +4.9 %/K (model levels, 4 models); with the 8-region shear to +4.6 to +5.0 %/K.
  - EC-Earth3: +12.6 to about +9.9 %/K.
  - TCO95 moves the other way, from +5.4 to about +6.1 %/K.
  - So at hub height the CMIP6 and AWI-CM3 per-K increases roughly meet.
- Limits:
  - Valid 1000 hPa cell-days are higher-pressure, calmer days (DE sfcWind 4.0 against 4.6 m/s). Their sfcWind change tracks the all-day change (-4.6 against -5.1 %), but G falls more on them (-3.3 against -2.3 %).
  - Model-level V100/G on all days is only +0.2 % [-0.4, +2.0] in the four models, so the agreement with IFS above 100 m is less firm than the near-surface result.
  - Land cover exists for 10-12 models and model levels for 4.
  - The event consequence goes through two linear fits and is an estimate of size, not a recomputation of events at hub height.

## 13i. Southern-Hemisphere test of the geostrophic mechanism

- Run 2026-11-07 (Jan: "Yes we shall do this and check").
- Question: if NH winter low-wind events increase because G weakens (13b), SH winter (April-September) events should change with the sign of the G change, also where G strengthens.
- Regions (`dunkelflaute/scripts/regions/sh_regions.py`), land cells (ERA5 land fraction >= 0.5) in lat/lon boxes, equatorward edge at or poleward of 28 S so that f is not small: AU_SE southern and eastern Australia (39-28 S, 138-154 E), AU_SW south-western Australia (35.5-28 S, 114-124 E), ZA southern Africa (35-28 S, 17-32 E), AR_C central Chile and Argentina (38-30 S, 72-57 W), PA Patagonia (52-38 S, 73-63 W), BR_S southern Brazil and Uruguay (35-28 S, 58-48 W), NZ New Zealand (47.5-34 S, 166-179 E). Plus 47 half-open 5 x 5 deg land boxes over 55-25 S with ERA5 land fraction >= 0.25.
- SH winter = April-September of one year, labelled by that year; a winter needs >= 170 days.
- Extraction (`sh_extract.py`, job `dunkelflaute/jobs/sh_extract.sbatch`, jobs 28811049-28811057, 2-6 min each; output `dunkelflaute/data/sh/sh_<source>_<name>.npz`, about 4.4 GB, gitignored): daily-mean 10 m speed, its log-law 100 m speed (ERA5 fsr z0 climatology) and onshore CF per region; daily-mean msl on a 2.5 deg grid, 65-15 S, global ring.
  - CMIP6: the 20 models of Section 14 plus EC-Earth3, `day/sfcWind` and `day/psl` (bilinear, periodic), historical 1951-1980 and ssp585 2071-2100. FGOALS-g3 lacks 2100 in one variable (29 future winters).
  - 100 km: 3-hourly 10u/10v (mean speed of the 8 instants in (d 00, d+1 00]) and the daily msl file, 1952-2099, both runs.
  - 31 km: 6-hourly 10u/10v and msl, 4 instants per day.
  - IFS-FESOM2-SR: daily `mean10ws` and `mmsl`, all three runs.
  - 9 km: native 3-hourly 10u/10v on the reduced grid (`atm_reduced_3h`, all months exist in both runs; only the reconstructions are Sep-Mar), ERA5 land and z0 by nearest neighbour, unweighted means; msl from the 3-hourly remapped file, 06/12/18/24 UTC.
  - Our runs' msl is a box mean, CMIP6's bilinear, as in 13b.
- Analysis (`dunkelflaute/scripts/regimes/sh_geostrophic.py`, about 2.5 min on the login node, output `dunkelflaute/results/sh/sh_geostrophic_summary.json`): G = |grad p| / (rho |f|) as in 13b; s = sfcWind (CMIP6) or log-law 100 m speed (our runs); R_s, R_G (s* = a_m G) and wind-CF events R_cf with control-matched thresholds; G per K = 100 ((G_fut / G_ctl)^(1/dT) - 1). Comparison spans as Section 10c with April-September winters. Winter bootstrap for our runs (5000, same winters in all regions); null from the 147 SH winters of the 100 km control. Matched threshold reaches 90 % of the target rate in all but 2 of 54 regions x 21 models x 3 series. Daily s against G in the control: r = 0.77-0.88 per SH region (NH 0.84).
- Cross-hemisphere test (CMIP6, wind-CF events per K against G per K, NH 22 regions from 13b and Section 14):
  - Per model and region: NH r = -0.42, slope -3.3, intercept +4.3 %/K (n = 440); SH r = -0.50, slope -3.8, intercept +1.7 %/K (n = 140). Same slope, lower SH intercept.
  - Region medians (22 NH + 7 SH): r = -0.66, p = 1e-4; the event change has the opposite sign to the G change in 29 of 29 regions.
  - Within each model across its 29 regions: r < 0 in 20 of 20 (median -0.51).
  - NH fit applied to the SH: r = 0.50, bias -2.3 %/K, sign right in 68 % of model-region pairs.
  - Where G strengthens (38 SH model-region pairs) events fall in 68 % (median -2.6 %/K); where it weakens (102) they rise in 82 % (median +6.4 %/K).
  - Speed-space events give the same picture: NH r = -0.42, SH r = -0.43, slopes -3.3 and -3.4, region medians opposite in sign 29 of 29.
- SH CMIP6 medians per region, G / events (%/K): AU_SE -1.49 / +7.3, AU_SW -1.31 / +2.9, ZA +0.24 / -1.7, AR_C -0.20 / +1.6, PA -0.59 / +9.3, BR_S +0.31 / -2.8, NZ -0.29 / +7.4. G strengthens in 13 of 20 models in ZA and BR_S.
- SH 5 deg grid (CMIP6 medians, 47 boxes): r = -0.67. G strengthens at 30-25 S and 55-50 S and weakens at 45-35 S; events fall in 12 of the 20 boxes where G strengthens. Southern tip of South America: G up, events up.
- SH field (7 regions, CMIP6 median): events +3.4 %/K, s* +3.2 %/K, share from circulation 0.71 (slope of ln R_s on ln R_G 0.73, r 0.65), mean G -2.0 % [-3.9, -1.0], transfer +0.1 % [-1.5, +3.0]. Grid mean events +1.0 %/K.
- Our runs, 7 SH regions, CF events per K [5-95 % bootstrap], mean G, transfer: 9 km +3.4 [1.8, 5.2], -4.4 %, +3.2 %; 31 km +4.9 [3.5, 6.4], -5.3 %, +3.3 %; 100 km +2.2 [0.8, 3.5], -4.9 %, +3.1 %; IFS-FESOM2-SR ctl +1.2 [-0.9, 3.3], -1.4 %, +1.2 %; hist +1.0 [-1.8, 3.8], -2.1 %, +1.3 %. The rising transfer of our runs (13b, 13c) appears in the SH as well. G strengthens in at most 1 of 7 named regions in our runs (on the grid in 11-15 of 47 boxes; events and G opposite in 31-42 boxes). Pooled over runs and regions: r = -0.46, slope -3.4, opposite sign in 80 %.
- Null, 100 km control, 7-region field ratio 5-95 %: 30/30 winters 0.95-1.05, 17/10 winters 0.90-1.08.
- Verdict: the sign of the G change predicts the sign of the event change in both hemispheres, with the same sensitivity (about -3.5 % events per % G). Mean G explains about a quarter of the model-region variance, and the NH has an extra +2.5 %/K at unchanged mean G that the SH does not.
- Against Qu et al. (2025): no overall SH winter decrease over land here; events decrease only where G strengthens. Checked 2026-11-07 on the article page: Qu et al. use hourly wind over the whole year and find wind-drought frequency -5 to -15 % between 30 S and 30 N (SSP5-8.5), +10-20 % at 30-60 N. Their decrease is tropical and subtropical, so it does not conflict with our SH regions poleward of 28 S.
- Limits: daily-mean msl at 2.5 deg; msl reduction over the Andes and the South African plateau (lowest s-G correlation, 0.77, in AR_C and ZA); some coarse CMIP6 grids have 2-6 cells in NZ and AU_SW; the null for IFS-FESOM2-SR comes from the 100 km control.
- Paper figure: `scripts/southern/plot_southern.py` (fig:southern) in the paper repo.

## 15. Gridded per-K maps and the aggregation scale

- Grid: `REGION_SET=grid5` in `dunkelflaute/scripts/regions/regions.py` replaces the 24 named regions with 5 x 5 deg boxes over 20-70 N (half-open edges), kept where the ERA5 land fraction is at least 0.25. That gives 420 boxes, plus DE. The box list is cached in `dunkelflaute/data/regions_grid5/grid5_boxes.json`. Everything else in the region pipeline is unchanged: land mask, z0, log law to 100 m, QM, CF curves, solar v7 for AWI-CM3, native `ssrd` for EERIE.
  - Our runs: the same extraction scripts on the 6-hourly footing, output `dunkelflaute/data/regions_grid5/`. CMIP6: `extract_daily.py` with the same switch on the daily footing, output `dunkelflaute/data/cmip6_grid5/`.
- Per box (`dunkelflaute/scripts/regions/grid_per_K.py`):
  - Threshold matched to 5 events per control winter, searched on 200 levels between the 0.1st and 40th percentiles of the control series. A box is valid if its control reaches 0.9 x 5.
  - Ratio of events per winter, per K = 100 (exp(ln ratio / dT) - 1), and the mean-removed ratio (future rescaled to the control mean, same threshold).
  - Comparisons on the 6-hourly footing (WIN 8): 9 km, 31 km, 100 km, EERIE ctl -> SSP2-4.5 and hist -> SSP2-4.5. Daily footing (WIN 2): 20 CMIP6 models (1951-1980 vs 2071-2100, and 2041-2070), plus the EERIE and TCO95 daily bridges.
- Significance:
  - Our runs: a per-box chance range from the TCO95 constant-1950 control (147 winters). Take all placements of two disjoint blocks with the comparison's own winter counts, and use the 5-95 % range of ln ratio. For the 100 km comparison (147 + 19 winters) too few placements remain, so the control block is cut to W - n_b - 40 winters. A shorter block widens the range, so the test stays conservative.
  - CMIP6: the share of models with ln ratio > 0. Agreement means at least 80 % of models.
- Pattern correlation: cos-lat weighted, over the boxes valid in both fields, between each field and the CMIP6 median. As a reference for how similar two independent models are, each CMIP6 model is also correlated with the median of the other 19.
- Aggregation scale:
  - The box series are combined, cos-lat weighted, into 10, 20 and 40 deg tiles (at least 2, 4 and 8 land boxes; a time step needs at least half the tile weight present). They are also combined into the 24 named regions and Europe (35-70 N, 10 W-30 E), each approximated by the boxes whose centres fall inside it.
  - Events are recounted on each aggregate with its own matched threshold. The result is per K of the area-weighted mean ln ratio over the aggregates.
- Run 2026-11-05:
  - Grid extractions: TCo319 jobs 28762539 and 28762540 (4 workers, because 10 ran out of memory). EERIE jobs 28762760-28762762 (16 workers, because 64 ran out of memory).
  - Analysis: job 28763196, about 1 h, output `dunkelflaute/results/maps/grid_per_K.npz` and `grid_per_K_summary.json`.
  - Figures: `dunkelflaute/scripts/plotting/plot_per_K_maps.py` (module python3/2025.01, cartopy), output `per_K_maps_wind.png`, `per_K_maps_tail.png` and `per_K_vs_scale.png` in `dunkelflaute/results/maps/`.
- Results, wind only, area mean over the 5 deg boxes:

  | comparison | per K (%/K) | mean removed | boxes up | above null p95 | below null p05 |
  |---|---|---|---|---|---|
  | 9 km | +3.5 | +3.1 | 81 % | 35 % | 4 % |
  | 31 km | +2.8 | +3.2 | 79 % | 40 % | 7 % |
  | 100 km | +2.2 | +2.6 | 67 % | 47 % | 15 % |
  | EERIE ctl | +4.8 | +5.1 | 76 % | 46 % | 7 % |
  | EERIE hist | +4.7 | +4.4 | 72 % | 29 % | 5 % |
  | CMIP6 median map | +3.7 | +2.4 | 80 % | at least 80 % of models agree on an increase in 50 % of boxes, on a decrease in 4 % | |

  - Pattern: increases over the US Great Plains and Midwest, Europe, Kazakhstan and North China. Decreases over the Sahara and Arabia, where every model agrees. CMIP6 zonal mean: peak +7.1 %/K at 35-40 N, -0.6 at 20-25 N, +1.0 at 65-70 N.
  - Pattern correlation with the CMIP6 median: 9 km 0.44, 31 km 0.54, 100 km 0.63, EERIE 0.31-0.41. Single CMIP6 model against the rest: median 0.60 (0.27-0.70). So our runs agree with CMIP6 about as well as one CMIP6 model agrees with another; the 9 km and EERIE runs sit at the lower end.
  - Solar + wind: CMIP6 median +2.7 %/K, but -0.3 %/K with the mean removed. At the 5 deg scale the combined increase is almost all a mean shift.
  - Mid-century CMIP6 (2041-2070): +4.6 %/K, about the same as late century, so roughly linear again.
- Aggregation scale, wind only, area-mean per K (mean removed in brackets):

  | | 5 deg | 10 deg | 20 deg | 40 deg | named |
  |---|---|---|---|---|---|
  | CMIP6 median | +4.1 [+2.3] | +5.5 [+2.6] | +7.0 [+2.7] | +8.9 [+3.9] | +8.3 [+4.2] |
  | 9 km | +3.5 [+3.1] | +3.3 [+3.0] | +3.4 [+3.0] | +5.6 [+3.9] | +6.1 [+4.6] |
  | 31 km | +2.8 [+3.2] | +2.4 [+3.2] | +2.5 [+3.4] | +3.0 [+3.4] | +4.3 [+4.0] |
  | 100 km | +2.2 [+2.6] | +1.7 [+2.4] | +1.5 [+2.4] | +2.5 [+2.4] | +4.9 [+4.1] |
  | EERIE ctl | +4.8 [+5.1] | +4.8 [+5.3] | +4.3 [+4.9] | +3.6 [+4.4] | +5.6 [+5.3] |

  - In CMIP6 the change grows with the scale of aggregation, from +4 %/K for a 5 deg box to +9 %/K for a 40 deg tile. The mean-removed part hardly grows (+2.3 to +3.9), so the growth is the mean wind decline. A working explanation, not yet tested: aggregation lowers the variability of the series relative to its mean, so the same relative decline in the mean moves more of the distribution below a fixed threshold.
  - In our runs the change is nearly flat from 5 to 20 deg. They lose less mean wind than CMIP6 (Section 14), consistent with the rising transfer of Section 13b.
  - The named-region values are higher for everyone, partly because of footprint rather than scale. Those regions sit on mid-latitude wind areas and leave out the desert boxes. The CMIP6 named value (+8.3 %/K, box-approximated) matches the +7.5 %/K field mean of Section 14. That accounts for the 5 deg boxes showing about half the Section 14 value.
  - CMIP6 named regions, wind only, median %/K: Europe +9.6, DE +8.3, US_MW +10.6, IN_N +12.2, CN_E +0.6, CA_PR +2.9.
  - Solar + wind: the CMIP6 median grows from +2.9 to +7.2 %/K at 40 deg; the mean-removed part grows from -0.5 to +2.3. In our runs, solar + wind at 40 deg is noisy (14 tiles, fewer valid solar boxes) and drops to +0.5 to +0.6 %/K for 9 km, 31 km and 100 km.
- Notes:
  - The d3h footing of 14b drops 31 March in every winter.
  - The 40 deg tiles cover only 14 aggregates and mix land types. They do not form a clean scale series on their own, so the scale result rests on 5 to 20 deg together with the named regions.

## 15a. Field significance of the per-K maps

- Script `dunkelflaute/scripts/regions/fieldsig_grid.py` (job 28790343, 4 min); output `results/maps/fieldsig_grid.npz`, `fieldsig_grid_summary.json`, `fieldsig_null_counts.npz`.
- Null: TCO95 constant-1950 control (147 winters), 6-hourly, thresholds matched to 5 events per control winter in every box as in `grid_per_K.py`. For each comparison all placements of two disjoint blocks of consecutive winters with the comparison's own winter counts; control block cut to W - n_b - 40 winters where too few placements remain (100 km 88, EERIE ctl 77). Same placement in every box, so spatial correlation is kept (`block_pairs`, `block_logratios` from `region_events.py`). 1722 to 14 762 placements.
- Per box: two-sided p = min(1, 2 min(p_up, p_down)), p_up = (1 + #null >= observed)/(1 + P). FDR after Benjamini and Hochberg (1995), as recommended by Wilks (2016), alpha_FDR = 0.10 (resulting p threshold 0.020-0.059). Gaussian p from the null mean and SD passes 5-15 % fewer boxes.
- Check: null 5-95 % per box equals `null_p05`/`null_p95` in `grid_per_K.npz` exactly.
- Results, wind, FDR-significant boxes of 420 (increase/decrease): 9 km 120 (108/12), 31 km 162 (135/27), 100 km 252 (193/59), EERIE ctl 211 (187/24), hist 99 (90/9). Uncorrected p < 0.05 in 117-248 boxes against 21 by chance.
- Area mean (cos-lat mean of ln ratio): 1.09-1.20 against a null 95th percentile of 1.01-1.02, above every placement in all five runs (p <= 0.0007). Solar + wind: 79-175 significant boxes, every area mean above the whole null.
- The significance test uses the area mean of ln ratio, 0.1-0.2 %/K below the per-box %/K means in the map titles (9 km +3.4 against +3.5).
- Caveat: the null comes from the 100 km model; if the 9 km run or EERIE has more decadal variability, the test is too lenient for them.
- CMIP6 per-model ranges (`fieldsig_cmip6_bootstrap.py`, `results/maps/fieldsig_cmip6_bootstrap.json`): 5000 winter bootstraps, same winters in every region, stored thresholds; reproduces `daily_events_summary.json` exactly; per-model 5-95 % about +-1 %/K.
- Through-origin fits of ln ratio on dT (paper fig:perK): wind CMIP6 +7.2 %/K, our five main comparisons +5.4 %/K; solar + wind +4.0 and +3.9 %/K.

## 15b. Europe-wide co-occurrence and the variance explanation of the scale growth

- Question: do winter low-wind spells hit more of Europe at once under warming, and why does the per-K change grow with aggregation scale in CMIP6 (Section 15)?
- ERA5 input (`dunkelflaute/scripts/regions/cooccur_extract_era5.py`, job 28790362, 3 min): hourly 10u/10v (N320) at 00/06/12/18 UTC, Sep-Mar, winters 1980-2026 (47). Land and z0 as in the model extractions, log law to 100 m, onshore curve. Footings: 6-hourly instantaneous and daily mean of the four instants. Regions: the 8 European regions and the 420 grid5 boxes. Output `dunkelflaute/data/regions/cooccur_era5_wind.npz`.
- Co-occurrence (`cooccur_events.py`, job 28790567, 3 min):
  - Regions DE, FR, IB, UK, PL, SC, IT, BK, wind only. Thresholds exactly as `region_events.py` (6-hourly, our runs) and `daily_events.py` (daily, CMIP6, EC-Earth3, bridges); checked identical.
  - N(t) = number of regions in an event at step t; statistic = share of winter time with N >= k.
  - Future variants: all (control thresholds), mean removed, dependence only (thresholds matched anew in the future, so regional rates are fixed).
  - Null: TCO95 control, 400 random disjoint block placements per comparison; 5-95 % of ln ratio.
- Results, at least 4 of 8 regions at once:
  - Control: ERA5 2.1 % of winter time (6-hourly), 1.5 % (daily); our runs 2.1-2.6 %; CMIP6 median 1.9 % (1.2-3.0). Joint lows 10-14 times more frequent than for independent regions. Pairwise Jaccard pattern against ERA5: r = 0.93-0.98 (ours), 0.80-0.94 (CMIP6). ERA5 halves: 1.45 and 1.53 % (daily).
  - All: +7.7 (9 km), +15.8 (31 km), +12.1 (100 km), +12.1 / +10.7 %/K (EERIE ctl / hist); CMIP6 median +24.9 %/K (+10.4 to +32.0), 19 of 20 up, 18 beyond the null. Days per winter: 9 km 5.6 to 8.4, 31 km 4.4 to 10.2, 100 km 5.3 to 9.3; CMIP6 4.1 to 10.5.
  - Dependence only: -1.6 to +6.2 %/K (ours), CMIP6 median +4.6 (-0.6 to +9.0), 17 of 20 up, 4 beyond the null. At k = 3, four of five of ours and 8 of 20 CMIP6 beyond the null, at +2 to +6 %/K.
  - Share of ln change from co-occurrence itself: CMIP6 0.14 (k = 3), 0.23 (k = 4); ours 0.2-0.5.
  - So Europe-wide spells grow 2-4 times faster than single-region events, mostly because every region has more events; coherence itself rises only slightly.
- Generation check (`cooccur_generation_check.py`): DE (SMARD), PL and SE (energy-charts), 6-hourly, winters 2016-2026, against ERA5 DE, PL, SC: Jaccard DE-PL 0.18 vs 0.21, DE-SE 0.12 vs 0.11, PL-SE 0.10 vs 0.09; at least 2 of 3 at once 4.4 vs 4.2 % of winter time.
- Variance explanation (`cooccur_scale_cv.py`, job 28790631, 3.5 min):
  - NH, CMIP6 median: per K +3.9, +5.5, +7.0, +8.9 (5 to 40 deg, reproduces Section 15); mean removed +2.2 to +3.9; mean-only +1.8 to +4.9. Mean change about -5 % at every scale. CV 0.92, 0.76, 0.59, 0.40; elasticity d ln events / d ln threshold 1.1, 1.7, 2.5, 3.3. Mean-only ln ratio = elasticity x mean decline at r = 0.99 (20 models x 4 scales). ERA5 aggregates the same way (CV 0.90 to 0.42 daily).
  - Europe family (5, 10, 20 deg aggregates in 35-70 N, 10 W-30 E, and Europe itself), CMIP6 median: +5.2, +6.5, +7.7, +9.7 %/K; mean removed +3.3 to +4.7; mean-only +2.0 to +6.2; mean change -7 to -8 %; CV 0.78 to 0.33 (ERA5 daily 0.90 to 0.40).
  - Our runs: +2.3 to +4.0 %/K at 5 deg, +3.6 to +6.0 for Europe; mean change -2 to +3 %, mean-only part within +-1.7 %/K; their aggregates are just as elastic.
  - So the CMIP6 scale growth is the mean decline acting on a less variable aggregate.
- Caveats: Europe as one series is a single aggregate; CMIP6 is less variable than ERA5 at the Europe scale (CV 0.33 vs 0.40), partly the footing; k >= 6 is too rare to use.
- Paper figure: `scripts/europe/plot_europe.py` (fig:europe) in the paper repo.

## 15c. Solar-only events and the attribution of the combined change

- Run 2026-11-07, prompted by Jan's question whether the paper is mostly wind because solar showed no change (it is not: solar had never been analysed on its own).
- Code (analysis repo): `dunkelflaute/scripts/regions/solar_common.py` (helpers), `solar_region_events.py` (24 regions, 6-hourly, our runs, about 10 min on the login node), `dunkelflaute/scripts/cmip6/solar_daily_events.py` (daily, CMIP6 and bridges, job 28811385, 9 min), `dunkelflaute/scripts/regions/solar_grid_per_K.py` (5 deg grid, job 28811157, 13 min). Outputs `results/regions/solar_region_events_summary.json`, `solar_region_events_counts.npz`, `solar_grid_per_K.npz`, `solar_grid_per_K_summary.json`, `results/cmip6/solar_daily_events_summary.json`.
- Series: the solar CF the combined metric already uses (9 km v5 with exact monthly rescale, 31 km v7, 100 km and IFS-FESOM2-SR native ssrd, CMIP6 daily-mean rsds / 1000). Events, winters, matched threshold (5 per control winter), validity (90 %), bootstrap and block nulls as in Sections 10c, 10d, 14, 15, 15a. For solar the threshold levels span the 0.1st to 40th percentiles of the control's 48 h means, because half the raw 6-hourly values are night zeros.
- Two definitions. Absolute: the solar CF itself. Deseasonalized: CF divided by the control calendar-day climatology (31-day running mean, floored at 5 W m-2; the future divided by the control climatology). The 48 h mean on the 6-hourly footing covers whole days, so the diurnal cycle does not enter. Grid boxes whose climatology drops below the floor (polar night) are left out of the deseasonalized maps.
- Seasonality: absolute events put 92-94 % of their time in Nov-Jan (73-78 % in Dec-Jan; uniform 43 %), so they measure the depth of the solstice minimum; Scandinavia never reaches 5 separate events per winter. Deseasonalized events put 48-59 % of their time in Nov-Jan.
- Mean solar CF, Sep-Mar, 24-region mean, our runs: 9 km -0.80, 31 km -0.68, 100 km -0.55, IFS-FESOM2-SR -0.75 (ctl) and -0.57 (hist) %/K. Monthly at 9 km: Sep +1.4, Dec -9.5, Jan -10.3, Mar -6.0 %. Largest losses at high latitude and in the continental interiors (SC, RU, KZ, CA_PR, CN_N, CN_NE).
- Solar events, our runs, absolute (23 regions): 9 km +3.3, 31 km +2.2, 100 km +4.5, IFS-FESOM2-SR +4.5 (ctl) and +3.4 (hist) %/K; up in 18-21 regions; field mean above the TCO95 block null in every comparison (p <= 0.0006); share from the mean 18-44 %. Deseasonalized (24 regions): +7.4, +6.9, +8.0, +11.9, +9.4 %/K, up in 20-22, share from the mean 34-65 %.
- CMIP6 daily, 22 common regions (no SC, RU), absolute: median -0.2 %/K (10-90 % -5.2 to +3.9), 9 of 20 up, 9 beyond the null up and 9 down, origin fit -0.5. Deseasonalized: +2.4 (-7.5 to +9.7), 12 of 20 up. Mean solar CF over the 24 regions: median +0.44 %/K (14 of 20 brighten); Europe +6 to +8 %, CN_E +4 %, high-latitude interiors -3 to -4 %; by month +3 to +5 % in Sep-Nov, -1.5 to -2.9 % in Dec-Feb. Across models solar per K correlates with the mean solar change per K at r = -0.66 (absolute) and -0.81 (deseasonalized).
- Daily bridges agree with the 6-hourly values: TCO95 +5.0 against +4.5 %/K, IFS-FESOM2-SR +4.7 and +3.5 against +4.5 and +3.4 %/K (absolute).
- Grid (5 deg), deseasonalized: area mean per K 9 km +6.6, 31 km +6.0, 100 km +9.5, IFS-FESOM2-SR +14.6 and +12.8; FDR-significant boxes 254, 271, 292, 287, 238 of 331-340 (20-51 decreases); area mean above every null placement. CMIP6 median map +5.3 %/K, per-model area mean median +3.7 (18 of 20 up); at least 80 % of models agree on an increase in 51 % of boxes, on a decrease in 19 %; zonal median -3 %/K at 25-35 N, +6 to +22 %/K at 40-65 N. Pattern correlation with the CMIP6 median 0.76-0.90 (single CMIP6 model against the rest: median 0.84; wind maps 0.31-0.63); per-box per K against per-box mean solar change r = -0.78 to -0.86. Mean solar CF over the grid: our runs -1.1 to -1.35 %/K, CMIP6 median map -0.54 %/K, 19 of 20 models darken.
- Grid, absolute: our runs +1.2 to +5.4 %/K (9 and 31 km +1.3 and +1.2, mean removed -0.4), 90 % of event time in Nov-Jan, no valid boxes north of about 55 N; CMIP6 per-model median +0.6, 13 of 20 up.
- Attribution of the combined change. Wind and solar winters aligned by calendar slot (29 Feb dropped except in 360-day calendars) and recombined across winters, 16 solar winters per wind winter, with the threshold matched in the recombined control: control wind with control solar (base), future wind with control solar, control wind with future solar, both future. Our runs (18 comb regions): wind +1.6 to +4.1, solar +0.1 to +2.3, interaction -0.4 to +0.1 %/K; the coupled change on the same series reproduces the real one (9 km +4.3 against +4.2). CMIP6 (24 regions): median wind +5.6, solar -0.4 %/K; MRI-ESM2-0 and ACCESS-CM2 have solar -2.3 and -2.5 %/K, which turn their comb change negative. Variant with the solar part replaced by its control climatology in both periods: ours +0.4 to +4.2, CMIP6 median +5.2 %/K.
- Why the mix is weaker than wind: during control comb events the wind supplies 48-54 % of the 48 h combined CF in our runs (43 % CMIP6), and the wind sits at the 5th percentile of its own 48 h distribution (7th in CMIP6), where its future quantile falls by 12-17 % (CMIP6 21 %). So a given wind change moves the combined CF about half as much. On the same regions, wind through the mix gives +1.6 to +4.1 against +3.4 to +5.2 %/K for wind alone (CMIP6 median +5.6 against +7.5). Solar adds to the mix in our runs and slightly offsets it in CMIP6.
- Paper: fig:perK (a wind, b solar absolute, c wind + solar, d per K against mean wind CF change, the former fig:meantail a); fig:maps_solar (Extended Data, deseasonalized). Through-origin fits CMIP6 / ours: wind +7.2 / +5.4, solar -0.5 / +3.3, mix +4.0 / +3.9 %/K. Panel d: r = -0.90, fit 3.03 - 2.07x, our runs 1.2-2.2 %/K above the line.
- Aerosol forcing (checked 2026-11-07, answering the earlier open question). Both models prescribe time-varying anthropogenic aerosol, direct effect only.
  - AWI-CM3 / OpenIFS 43r3, OIFS branch `awicm-3-aerscaleyearly-xios` (T. Semmler 2022, tip 1a7cb8d), source `/work/ab0246/a270092/model_codes/awicm3-frontiers-xios/oifs-43r3/src/ifs/`:
    - Baseline is the CAMS climatology `aerosol_cams_climatology_43R3a.nc` (`phys_ec/su_aer_climatology.F90:51`; `NAERMACC=1`, `LNEWAER=.T.` by default in `phys_radi/suecrad.F90:639,651`).
    - Black carbon, organic matter and sulfate (`RMACBC1/2`, `RMACOR1/2`, `RMACSU1`) are multiplied by a 3x3 deg factor map for each year 1850-2085 from `aerosol_scale_1850_2085_r2005.nc` (`su_aer_climatology.F90:177`, `utility/updtim.F90:759-763`). The factor is MACv2 anthropogenic AOD (Kinne gt_ant) over its 2005 value. Dust and sea salt are not scaled. After 2085 the 2085 map is held.
    - `NCMIPFIXYR>0` freezes the year, 0 reloads each 1 January (`updtim.F90:745-746`). The same file is used whatever `SSPNAME` says, and its 2025-2085 scenario is undocumented.
    - Factor over Germany (47-55 N, 6-15 E): 1950 1.42, 1980 1.80, 2020 0.71, 2085+ 0.49. Eastern China 0.07 to 0.39, global mean 0.44 to 0.84.
    - No aerosol-cloud effect. Droplet number comes from 10 m wind (`surf/module/surfrad_ctl_mod.F90:541-572`, `phys_radi/liquid_effective_radius.F90:144-163`); `LAERCCN` false, `NAERCLD` 0.
    - 100 km: fort.4 has `NCMIPFIXYR = 1950` in ctl1950g and 0 in hi1950g; the log prints "UPDTIM AEROSOL SCALING APPLIED" after each reload. 9 and 31 km: no namelists on Levante. A clear-sky fingerprint (clear-sky surface solar / (1 - albedo), future minus control, %) matches the 100 km pattern: C. Europe +0.1 / +0.5 / +0.6, E. US +2.9 / +3.1 / +3.4, E. China -6.7 / -6.9 / -6.3, N. Pacific -3.4 / -3.7 / -3.0 at 9 / 31 / 100 km. Medium-high confidence the scaling was on.
  - IFS-FESOM2-SR (`/work/bm1344/a270228/ifs-bundle-DE_CY48R1/source/raps/bin/model:3725-3813`, Ghosh2026 Sect. 3): hist uses time-varying `aerosols_CONFESS_1950_2015.nc` (`NCMIPFIXYR=-1`, before 1970 repeats 1970); ctl the same file frozen at 1950, i.e. 1970 aerosol; SSP2-4.5 `aerosol_cams-extended-macv2sp-SSP245_..._CONFESS_1951-2100_4D.nc`. C. Europe sulfate 4.3x its 2015 value in 1950-1970, 0.65x in 2050. `NAERCLD` 0, direct effect only.
  - DJF clear-sky surface solar, central Europe: +2.4 % at 100 km, -0.9 % at 31 km, -2.8 % at 9 km, as water vapour rises 35-67 %. So the aerosol decrease (about 70 % in mass) does not brighten our European winter. The missing indirect effect and differences in cloud response are the likely reasons for the gap to CMIP6's +6 to +8 % (about 70 % confidence; a clean test needs clear-sky downward solar, not in the output).
  - Script: scratchpad/clearsky_fingerprint.py.

## 11. Software

- Solar model steps: `/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin/python` (sklearn 1.1.1); the v1-v7 pickles do not load in sklearn 1.6.1. Its `shapely` needs `LD_LIBRARY_PATH=<that env>/lib`.
- Other steps: module `python3/2025.01-gcc-13.3.0` (sklearn 1.6.1).
- GRIB inspection: `grib_ls` / `grib_get_data` from `/work/ab0246/HPC_libraries/...` (eccodes).

## 12. Limitations to state

- Meteorological potential only: no demand, storage or grid.
- Fixed hub height and generic power curves; uncorrected CFs put the ERA5 event rate about 30% above Mockert's.
- Turbine icing and snow on panels make real cold events slightly deeper than the index, and both become rarer under warming, the opposite direction to the projected increase. Rare in practice: icing mostly on low-mountain hilltops, and winter solar CF is only about 0.03-0.08. To quantify: share of event hours with 2t below 0 C.
- Log law without stability, fetch or displacement height. The reconstruction has a cross-resolution transfer that is only partly tested (z0 coarsening for wind; nothing for solar).
- Single model family; SSP5-8.5 only; short 9 km windows (17 + 10 winters).
