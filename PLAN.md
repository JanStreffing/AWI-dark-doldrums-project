# Plan: projected change in winter wind-and-solar droughts across model resolutions

~~Plan: resolution dependence of projected Dunkelflaute change~~ (retitled 2026-11-05: the multi-region test shows the increase at all three resolutions, Phase 2b)

Written 2026-10-28, after taking over the intern project (Jing/Amanda Wang, Sep 8 to Oct 5 2026) and copying all of its data into this folder.

## Current direction (decided 2026-11-05)

### What we found (details in Phase 2b, numbers in `dunkelflaute/results/regions/` and `dunkelflaute/results/regimes/`)

1. **The increase is not a 9 km feature.** Over 24 Northern-Hemisphere land regions, low-wind events increase under SSP5-8.5 at 9, 31 and 100 km (wind-only mean ratio 1.37 / 1.31 / 1.25 on the 9 km spans; null 0.96-1.04). 9 km minus 31 km is small and not significant (ratio of ratios 1.05, p = 0.11-0.29).
2. **Germany is a weak-signal region at 31 and 100 km.** One region with 10-20 winters can miss a common signal. This likely explains why published Germany studies disagree (Strnad et al. flat, Heinrich & Geyer slightly up).
3. **Per K of global warming every run agrees.** About +5-6 % low-wind events per K (+3-4 %/K for the 0.44 solar + 0.56 wind mix), the same at 9, 31 and 100 km, in EERIE IFS-FESOM2-SR, and under SSP5-8.5 and SSP2-4.5. EERIE's smaller ratio is only its smaller warming.
4. **It is the low tail, not the mean.** Mean wind CF changes by -0.3 to -4 %, the 5th percentile of the 48 h mean falls 13-17 %. Removing the mean change takes away only 5-28 % of the increase.
5. **It happens within weather regimes.** Model regime frequencies match ERA5. Scandinavian blocking is the main Dunkelflaute regime. The change splits into a small negative regime-frequency term and a large within-regime term (9 km +1.36 events/winter = -0.20 + 1.55); events per regime day rise in all four regimes in almost every run (exceptions: 31 km NAO-, EERIE hist AR).
6. **Literature check** (section "Literature positioning"): the direction (more NH low-wind events) is published (Qu et al. 2025, Brayshaw et al. 2026, Heinrich & Geyer 2026). Not found anywhere: global coupled km-scale models for this question, a resolution test of the change signal, the per-K collapse, the regime decomposition, the quantified mean-vs-tail split, the sampling lesson. Three papers still to read in full (access).
7. **CMIP6 on Levante** (`/pool/data/CMIP6/data/ScenarioMIP`): about 20 models with daily `sfcWind` and `rsds` under SSP5-8.5, 8 with 3-hourly `uas` and `rsds` (AWI-CM-1-1-MR, BCC-CSM2-MR, ACCESS-CM2, MPI-ESM1-2-HR/LR, MRI-ESM2-0, CNRM-CM6-1, NESM3). HighResMIP is there too. Earlier CMIP6 studies report per scenario or warming level, region by region, with one fixed threshold and within-run tests, and do not separate mean from tail or dynamics from within-regime change.
8. **CMIP6 through our pipeline** (2026-11-05, Phase 2b): all 20 models show the wind-only increase (median +7.5 %/K, 10-90% +4.7 to +10.2; regression through the origin); our high-resolution runs (+5-6 %/K) sit in the lower half. The daily footing reproduces our 6-hourly per-K values. The solar + wind mix is weaker and model-dependent (+4.4 %/K, 18 of 20 up), and Germany is one of its weakest regions (11 of 20 models up), which fits the published Germany null results.
9. **The daily footing is safe, and the cause is circulation** (2026-11-05, Phase 2b).
   - Three footings of the same 4 CMIP6 models (6-hourly, daily from 3-hourly, daily files) agree within about 0.7 %/K.
   - Geostrophic wind from the msl gradient explains most of the increase: median 91 % across 20 CMIP6 models, about all of it at 9 and 31 km. The surface/geostrophic ratio does not rise in CMIP6 and rises slightly in our runs, so a more stable boundary layer is not the cause.
   - Pressure gradients weaken within every regime in the CMIP6 medians and under Scandinavian blocking and Atlantic ridge in every run (most under Scandinavian blocking); under NAO+/- not in every run (31 km NAO- +3.2 %).

### Decisions

- **Paper 1 is global, organised around one physical result, not a catalogue of countries:** how the winter wind-and-solar drought tail responds to warming, and why. Germany stays as the validated worked example, not the headline.
- **Gridded per-K maps replace the 24 boxes** as the main signal figure (no arbitrary region choice; readers can read off any country). The boxes stay for statistics and for the regime analysis.
- **CMIP6 goes through the same pipeline** (each model's own warming, base-rate-matched thresholds, block nulls), to show whether the per-K agreement holds beyond AWI-CM3 and EERIE, and whether CMIP6's spread shrinks per K. Done 2026-11-05 (Phase 2b): all 20 models show the wind-only increase, median +7.5 %/K; our runs (+5-6 %/K) sit in the lower half.
- **Developing regions only where the physics differs**, not as more of the same. The Southern Hemisphere goes into paper 1 as a mechanism test (Qu et al. find SH wind droughts decrease): wind-only at 9 km (log law on the raw 10 m wind for Apr-Sep), full mix at 31 km, 100 km and EERIE, which have all months.
- **A second paper stays open as a parking lot**, so ideas that would overload paper 1 are kept, not binned (list below).

### Paper 1 outline (working)

1. Signal: per-K maps at grid-cell level, 9 / 31 / 100 km and EERIE; field significance from block nulls that keep spatial correlation.
2. Robustness: CMIP6 daily ensemble through the identical pipeline; resolution (our three plus HighResMIP if cheap); perfect-model test of the reconstruction (TCO95).
3. Mechanism: mean vs tail everywhere; weather regimes for Europe; geostrophic vs 10 m wind everywhere (weaker pressure gradients vs a more stable boundary layer; needs no regime definition, so it works globally); NH vs SH contrast.
4. Europe in depth: Germany validated against SMARD and ERA5; regimes; whether events become more Europe-wide at once (grid balancing).
5. Impact, brief: which regions with large or fast-growing wind fleets sit in the strongest-response areas (Europe, US, China; then India, Brazil, South Africa).

### Paper 2 parking lot (ideas kept for later, not for paper 1)

- South Korea / East Asia with IBS: EPSIS capacity weights and hourly validation, year-round vs Sep-Mar season given the monsoon (already planned, see "Later").
- Monsoon and low-wind regions (India, East China, Sahel): solar- and cloud-led droughts, year-round season, solar-led definitions. Our combined metric fails there (no base-rate match in TR, US_CA, CN_N, CN_E, KR, IN_N).
- "Darker but warmer": temperature during events, a demand model from SMARD load, residual-load-based detection (Biewald et al. 2025).
- Fleet sensitivity: 2018 / 2024 / 2045-target mixes, offshore weight, country-specific mixes.
- Hub height 120 / 140 m from TCO95 model levels.
- Li-style and other definitions as a full intercomparison.
- EERIE ICON-ESM-ER (~10 km, 3/6-hourly wind, no sub-daily radiation): wind-only check, or moved into paper 1 robustness if cheap.
- Z500 blocking index (TCO95 daily, EERIE 6-hourly) as a cross-check of the msl regimes.
- Sub-national patterns, onshore vs offshore during events, icing and snow on panels.

## Where the project stands

- Done: wind (log law, 10 m to 100 m) and solar (v5 ML model, monthly to 3-hourly `ssrd`) reconstructions for TCo1279-DART 1950C and 2080C, validated against ERA5 and SMARD for Germany.
- Done: Mockert-type events for Germany, Sep to Mar, Mockert's 2018 weights (44% solar, 50% onshore, 6% offshore):

| Source | Winters | Events/winter |
|---|---|---|
| ERA5, reconstructed fields | 11 (2016-26) | 5.64 |
| ERA5, real fields | 11 | 6.82 |
| SMARD | 11 | 5.82 |
| DART 1950C | 17 (1953-69) | 3.47 |
| DART 2080C | 10 (2083-92) | 5.50 |

- The 1950C to 2080C increase is ~~significant at the headline threshold~~ significant only by within-run tests (superseded 2026-11-05: against the block null from long controls, Germany alone gives p = 0.13-0.20; Phase 2b). Exact Poisson test p = 0.015, and resampling whole winters gives a ratio of 1.58 with a 90% interval of 1.19 to 2.12. It is a low-tail shift, though:
  - The mean combined CF is unchanged (0.201 to 0.199).
  - The 5th percentile of the 48 h mean drops from 0.070 to 0.061.
  - The ratio falls from 2.3 at a threshold of 0.04 to 1.03 at 0.10.
  - A Li-style definition gives only +16%.
- Not started: Li threshold on DART, South Korea, sub-national patterns, synoptic composites, onshore vs offshore during events, uncertainty quantification.
- The country folders (BE, DK, FI, IE, NL, NO, PL, SE, plus Korea's EPSIS) hold downloads only.

## Why the data is interesting

- Same model at three resolutions with matching forcing windows, in `/work/ab0995/ICCP_AWI_hackthon_2025`:
  - TCO95L91-CORE2 (about 100 km) and TCo319-DART (31 km): `ctl1950d` 1950-69 and `ssp585d` 2080-99. Both have 6-hourly `10u/10v/2t/tcc/lcc/msl` and monthly `ssrd`.
  - TCo1279-DART (9 km): 1950C 1950-69 and 2080C 2080-92, 3-hourly.
  - TCo319 `ssp585d` also has monthly `z500` for 2015-2099.
  - TCo2559-DART-1950C exists but has no atmospheric output in this archive.
- The ocean is eddy-resolving, which is relevant for North Sea and Baltic offshore wind and for storm tracks.
- The 1950-forcing to 2080s SSP5-8.5 contrast is very large, and the domain is global.
- No km-scale Dunkelflaute projection study exists. The closest comparison, Strnad et al. 2026 (Mockert definition, 5 CMIP6 models), finds no significant change under SSP5-8.5. That is for one region (Germany/Europe); our Germany at 31/100 km agrees, but the multi-region result says a single region with 10-20 winters can miss a common signal.

## Method verdict

- Solar is good enough (r = 0.998 vs ERA5 over Germany, same skill as ERA5 against SMARD). ~~For TCo95/TCo319 it needs a reduced model without `tsr/mcc/hcc`, at 6 h.~~ (TCO95 rerun writes 3 h `ssrd`; TCo319 has all four cloud layers, so v7 uses `mcc/hcc` too.)
- Wind needs calibration. On ERA5 2016-26:
  - The reconstructed 48 h combined CF is about 10% too high in the event range, and onshore CF is 46% too high when real onshore CF is below 0.05.
  - The reconstruction counts 17% fewer events than real ERA5.
  - Year-to-year event counts correlate at r = 0.77 with real ERA5.
  - Agreement with SMARD comes from two errors cancelling.
- Statistics are the weak point: 17 vs 10 winters at 9 km, and a strong threshold dependence. More winters (31 and 100 km, a long control, ERA5 back to 1979) are needed.

## Todo

### Phase 1: housekeeping and method

- [x] Repoint `_ROOT = "/work/ab0995/a270321/AWI-dark-doldrums-project"` in all scripts and sbatch files to this folder.
- [ ] Optional: `chgrp -R ab0246` on the copied data, and `chmod a+r solar/data/era5_1h/*`. SLURM jobs still use `--account=ab0995`; change if the budget should come from ab0246.
- [x] Python environment: `analysis-toolbox/python-04.2026`, which the README and older sbatch files load, is not available to a270092. Two interpreters are in use:
  - `python3/2025.01-gcc-13.3.0` module (sklearn 1.6.1): wind calibration and DART jobs.
  - `/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin/python` (sklearn 1.1.1, the login default): all solar model steps. The v1-v5 pickles fail to load in 1.6.1 (`__pyx_unpickle_CyHalfSquaredError`), and the solar sbatch files already put this interpreter first on `PATH`. v6 is trained here too, so v5 and v6 can be evaluated side by side.
- [x] Fix factual errors in `README.md`, `solar/README.md`, `wind/README.md` and `Project_Proposal.md`:
  - [x] The weights are Mockert's 2018 mix (44/50/6), not "2024".
  - [x] Mockert citation: Mockert, Grams, Brown, Neumann (2023), Meteorol. Appl. 30, e2141, doi:10.1002/met.2141.
  - [x] "Lohmann et al. 2025" is Biewald et al. (2025), Environ. Res.: Energy 2, 025007. The 23% optimum uses 24 h, not 48 h. F = 0.14 belongs to Li's method, and F = 0.41 to residual load.
  - [x] Mockert reports about 4.3 events/yr (1979-2018, max 8 d). ERA5 at 5.64 is above that, not "inside 4-6".
  - [x] Drop the "30 and 54 day events in 2022/23" claim, which has no peer-reviewed source.
  - [x] Offshore r against SMARD is 0.86-0.87 once the 2015-17 capacity-step years are excluded, not 0.70.
  - [x] `dunkelflaute/data/germany_era5/README.md` still points at the old `boundaries/land.geojson` paths.
- [x] Note the known data-label issues in the READMEs, without rewriting about 2 TB of output:
  - The 2080C wind files' `description` says "1950C".
  - The solar `ssrd:comment` says v3.
  - The solar DART files lack `valid_training_region`.
- [x] Promote the robustness analysis into `dunkelflaute/scripts/analysis/` with outputs in `dunkelflaute/results/`. It covers significance, threshold and window sensitivity, 2024 weights, the Li-style definition, low tails and seasonality. Done: `dunkelflaute/scripts/analysis/robustness_mockert_events.py`, output `dunkelflaute/results/robustness/robustness_summary.json`.
- [ ] Wind calibration on ERA5 (2015-26 for calibration, plus 1979-2014 winters as an independent test):
  - [x] Target: ERA5's native 100 m wind. Hub height stays fixed at 100 m onshore and offshore with the generic power curves; no fleet realism (hub heights, park ages, per country). SMARD is a sanity check only, not a calibration target.
  - [x] Option A: quantile mapping of log-law 100 m wind onto ERA5 100 m wind (bins: month x 3 h slot x offshore / 3 onshore fsr classes). Scripts in `wind/scripts/methods/calibration/`, table `wind/models/quantile_mapping_v1/`, results `dunkelflaute/results/wind_calibration/`. 2-fold CV on 2016-26: onshore CF 0.252 to 0.207 (real 0.207), calm-hour bias 1.46x to 1.10x, events 5.64 to 6.18 (real 7.00), per-winter r 0.84 to 0.92. Independent 1980-2014 Oct-Feb winters confirm the mean fix (0.281 to 0.236, real 0.233). Applied to DART: events 1950C 3.47 to 3.59, 2080C 5.50 to 6.00, ratio 1.58 to 1.67, threshold curve unchanged in shape. The signal does not depend on the residual wind error.
  - [ ] Optional, only if needed for absolute event counts: Option B, a GBM for the 100 m / 10 m speed ratio using fields every resolution has (`u10`, `2t`, `tcc`, `lcc`, hour, month, land-sea mask, z0).
  - [ ] Pass criterion: event statistics from the calibrated reconstruction match those from real ERA5 100 m wind (6.82/winter at the headline threshold) across thresholds 0.04 to 0.10 and windows 24/48/72 h, and the per-winter correlation improves on 0.77.
- [x] 6-hourly sampling alone is harmless (`dunkelflaute/scripts/analysis/sampling_3h_vs_6h.py`, output `dunkelflaute/results/robustness/sampling_3h_vs_6h.json`). Degrading the 3 h series to 6 h output (wind instantaneous at 00/06/12/18, solar summed over 6 h, 48 h = 8 steps) changes events/winter at 0.06 by 0% for ERA5 real, ERA5 recon and 2080C, +2% for SMARD and -8% for 1950C, with per-winter r 0.93-0.99. The 2080C/1950C ratio at 0.06 goes 1.58 to 1.73.
- [x] Reduced solar model (v6): features `tcc`, `lcc` at both ends of a 6 h window plus solar geometry, trained on ERA5 Sep-Mar 2015-26 (`solar/scripts/core/reconstruct_ssrd_v6.py`, `extract_v6_6h.py`, `train_v6_model.py`, `jobs/train_v6_model.sbatch`, models in `solar/models/model_v6/`).
  - [x] Germany, held-out folds, 6 h windows (`solar/scripts/evaluation/evaluate_v6_germany.py`, output `solar/models/model_v6/evaluation/germany_v6_vs_v5_6h.{npz,json}`): v5 r = 0.996, RMSE 9.1 W/m2; v6 r = 0.982, RMSE 19.8 W/m2, bias +0.2 W/m2. With real ERA5 wind, events/winter at 0.06: real solar 6.73, v5 6.91, v6 6.73; v6 tracks real within 0.5 events at every threshold 0.04-0.10. v6 is good enough for event statistics.
  - [x] Without the monthly rescale v6 is biased +17 W/m2 over Germany in winter and undercounts events by 15% (5.73 at 0.06). The rescale is essential, ~~so TCO95 is unusable for solar unless its monthly `ssrd` turns up~~ (moot: the TCO95 rerun writes native 3 h `ssrd`).
  - Note: the legacy interpreter needs `LD_LIBRARY_PATH=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/lib` for `shapely` (system `libstdc++` too old).
- [x] ~~Decide scope (open, see conversation 2026-10-28): proposed is a 9 km paper at full 3 h resolution, with TCo319 (31 km) as a supporting wind-led check and v6 solar, and TCO95 dropped unless monthly `ssrd` exists.~~ Superseded 2026-11-05: all three resolutions plus EERIE, many regions, ratios normalised per K of warming. The long MR control (open question below) would matter more than the 20-year slices.

### TCO95 rerun (decided 2026-10-28: cheap, rerun with the output we need)

Full design spec: [`docs/tco95_rerun_spec.md`](docs/tco95_rerun_spec.md). Adds winds at 120 and 140 m (interpolated from model levels, validated against native 100 m) for a hub-height sensitivity at 100 km only.

- [ ] Output, matched to the 9 km run:
  - 3-hourly: `10u`, `10v`, `100u`, `100v`, `ssrd`, `tsr`, `tcc`, `lcc`, `mcc`, `hcc`, `2t`, `msl`, plus `sshf` and `blh` (or `skt`) for stability.
  - Daily: `z500`.
  - Monthly: `ssrd`, `fsr`.
  - Check that `100u`/`100v` can be enabled in the XIOS output for this OpenIFS version.
- [ ] Runs: a 1950 control (ideally 100+ years, as the natural-variability baseline) and SSP5-8.5 over at least 2080-2099 (ideally the full 1950-2100 transient, for time evolution).
- [ ] Perfect-model test: apply the log law (plus calibration) to TCO95's own `10u/10v` and v5 to its own clouds, then score against its own `100u/100v` and `ssrd`, in both 1950C and 2080C. This tests whether reconstruction error changes with climate, the main caveat on the 9 km 2080C result.
- [ ] Compare events from TCO95's native 100 m wind and solar with events from the reconstruction. If the 1950C to 2080C ratio agrees, the 9 km method is defended.

- [x] Runs done 2026-11-05: `ctl1950g` (E1) and `hi1950g` (E2), both 1950-2099, in `/work/ab0995/a270092/runtime/awicm3-frontiers-xios/`. E2 2100 still missing (1-year extra leg).
- [x] Perfect-model test (`dunkelflaute/scripts/tco95/perfect_model_series.py` + `perfect_model_events.py`, job `dunkelflaute/jobs/perfect_model_tco95.sbatch`, series in `dunkelflaute/data/tco95/`, summary `dunkelflaute/results/tco95/perfect_model_summary.json`). Germany on TCO95 = 55 land + 10 EEZ cells.
  - Reconstruction error barely changes with climate. Onshore CF from the log law vs native 100 m wind is 1.239 (control) and 1.264 (2081-99); calibrated 1.035 and 1.055. Offshore 0.889 and 0.885 (calibrated 0.929 and 0.925). v5 solar after rescale is exact by construction, r = 0.996 in both. Before the rescale, v5 shifts from 0.974 to 1.019, so the rescale also absorbs a climate-dependent solar error.
  - The reconstruction undercounts events as on ERA5 (control at 0.06: native 7.08, calibrated 6.48, uncalibrated 6.00) but reproduces the native future/control ratio: 2081-99 vs 1953-2099 at 0.06 is native 0.99, calibrated 1.03; 0.95-1.10 natively across thresholds 0.05-0.10.
  - ~~TCO95 itself shows NO Dunkelflaute increase under SSP5-8.5~~ TCO95 shows no increase **for Germany** with the Mockert mix (ratio about 1.0 at all thresholds; decadal means along E2 6.3-9.2/winter inside the control's 4.9-8.9). Superseded 2026-11-05: over 24 NH regions TCO95 does show the increase (wind-only 1.30, 16 of 24 regions significant, 2081-99 vs the 147-winter control), and in Germany wind-only 1.18. Germany is a weak-signal region at 100 km (Phase 2b).
  - Natural variability from the 147-winter control: a 10-winter block exceeds a disjoint 17-winter block by 1.58x in 0.3% of placements at 0.06 (7 events/winter), but in 2.3-5.3% at thresholds where the control has a 9 km-like base rate (2.9-5.9 events/winter). The 9 km result is therefore marginal at about the 95% level; spread is wider than Poisson, so within-run Poisson/bootstrap tests overstate significance.
  - ~~Open: is the 9 km increase resolution (blocking, sub-grid wind), tuning (TCO95 runs differ, see Notes for the paper), or a ~1-in-25 chance? The 31 km d runs (same tuning as 9 km, 20 + 20 years) are now the deciding comparison.~~ Answered 2026-11-05 by the multi-region test (Phase 2b): neither resolution-only nor chance; a broad signal present at all resolutions, with the Germany value sampling-dependent.
- [x] 31 km done 2026-11-05 (`dunkelflaute/scripts/tco319/`, summary `dunkelflaute/results/tco319/tco319_events_summary.json`). v7 solar (four cloud layers): r = 0.989, RMSE 15.6 W/m2 over Germany (v6 0.982 / 19.8, v5 0.996 / 9.1); events with native wind within 0.2 of native solar at every threshold.
  - All on a 6-hourly footing, 2080s/control at 0.06: 9 km 1.73 (uncalibrated) and 1.76 (calibrated, p = 0.003); 31 km 1.02 on the 9 km spans (p = 0.93) and 1.08 for 2081-99 (p = 0.66); 100 km native 0.96. At 0.04 the 31 km ratio is 1.3-1.5, not significant.
  - The 31 km run has the same tuning as 9 km, so the 9 km increase is not a tuning effect~~, and it is not shared by the coarser resolution~~ (superseded 2026-11-05: true for Germany only; over 24 regions 31 km gives 1.31 wind-only against 1.37 at 9 km). Mean onshore CF (calibrated) changes by 0% at 9 km (0.228 to 0.228) and +6% at 31 km (0.211 to 0.224). Offshore -4% at 9 km (0.527 to 0.508), +1% at 31 km.
  - 9 km and 31 km winters are independent samples: per-winter event counts correlate at r = -0.35 (control) and 0.04 (future), winter-mean CF at about 0 and +0.4-0.5 (n = 10, not significant). So no shared ocean-driven variability explains the difference.
  - ~~Remaining explanations for 9 km: a genuine resolution effect, or chance (TCO95 block null at a 9 km-like base rate: 2-5%).~~ Superseded 2026-11-05 (Phase 2b). Still worth doing: test the 9 km control for drift (within-run trend -1.69 events/decade, p = 0.06) and look at what differs in the low-wind tail (offshore CF falls, onshore tail).

### Phase 2: resolution comparison

- [x] Extract Germany land and EEZ from the TCo319 and TCo95 regular grids. Build z0 for those grids from ERA5 `fsr` (land) and Charnock (ocean).
- [x] Run wind (log law plus calibration) and solar (~~v6~~ v7 at 31 km, native at 100 km) at 100, 31 and 9 km, all at 6 h. The 48 h window is 8 steps.
- [x] Use matching windows: 1953-69 for the control, and both 2083-92 and ~~2083-99~~ 2081-99 for the future, since 9 km stops at 2092.
- [ ] Events per winter with intervals from resampling whole winters, full threshold and window curves, durations, and flagged time.
- [x] Natural-variability baseline: how often do two 10-winter windows differ by about 60%? Done with the TCO95 147-winter control and the EERIE 79-winter control (block nulls, Sections 10 and 10c of the methods notes). ERA5 1979-2026 not yet.

### Phase 2b: is the 9 km increase real? (proposed 2026-11-05)

- [x] Multi-region test, done 2026-11-05 (`dunkelflaute/scripts/regions/`, results `dunkelflaute/results/regions/region_events_summary.json`, methods in docs/methods_notes.md 10c). 24 NH land regions, base rate matched to 5 events/winter in each control, uncalibrated log law everywhere, 6-hourly.
  - **The increase is not a 9 km feature. It is a broad NH increase under SSP5-8.5 at all three AWI-CM3 resolutions.** Mean ratio over regions (9 km spans, 17 vs 10 winters):

    | metric | 9 km | 31 km | 100 km | EERIE SSP2-4.5 2041-70 vs control |
    |---|---|---|---|---|
    | wind only (24 regions) | 1.37 (22 >1, 14 sig.) | 1.31 (20, 11) | 1.25 (20, 6) | 1.14 (21, 6) |
    | 0.44 solar + 0.56 wind (18 regions*) | 1.26 (15, 5) | 1.18 (15, 4) | 1.13 (13, 1) | 1.07 (15, 3) |

    *Six sunny, low-wind regions (TR, US_CA, CN_N, CN_E, KR, IN_N) never reach 5 separate combined-CF events per winter, so they are left out of the combined field test.
  - Null for the mean over regions (same block placements in every region): 5-95% 0.92-1.08 (combined) and 0.96-1.04 (wind), TCO95 control; the EERIE control is similar. All three AWI-CM3 resolutions are far outside it, p < 0.002.
  - Full periods: 31 km 2081-99 1.29 (wind) / 1.23 (combined); 100 km 2081-99 against the 147-winter control 1.30 / 1.19, with 16 of 24 regions significant.
  - 9 km minus 31 km: ratio of ratios 1.05 (wind) / 1.06 (combined), with 9 km higher in 15/24 and 12/18 regions; p = 0.11-0.29. The gradient 1.25 / 1.31 / 1.37 is monotonic but not significant. The across-region patterns agree (r = 0.52-0.57 between 9 and 31 km).
  - Germany is a weak-signal region at 31 and 100 km: combined 1.04-1.09, wind only 1.17-1.20. At 9 km it is 1.23 (combined, matched rate) and 1.64 (Mockert mix at fixed 0.06, full EEZ). Per-region block null for Germany at 9 km: p = 0.13-0.20 (combined), 0.06-0.10 (wind).
  - So "do we need 9 km to see the change?" No. The signal is there at every resolution in most regions. The 9 km Germany value is the high end of a common signal, not a resolution-only effect.
- [x] EERIE IFS-FESOM2-SR extracted (hist, control, ssp245; 9-15 min per run). Weaker forcing gives a smaller but same-sign increase: wind 1.14 (2041-70 vs control 1972-2050), 21/24 regions > 1. Control drift (1972-2000 vs 2021-2050): 1.03-1.04, 0 significant. Germany, Mockert mix: 1.00. Rough scaling with warming (to be checked with global 2t): ln(1.37)/~4.5 K at 9 km vs ln(1.14)/~2 K in EERIE, about 0.07 per K in both.
- [x] Per K of global warming, done 2026-11-05 (`dunkelflaute/scripts/regions/global_t2m.py`, `region_per_K.py`, output `dunkelflaute/results/regions/region_per_K_summary.json`). Global annual-mean 2t between the compared periods: 9 km +5.50 K, 31 km +5.52 / +5.68 K, 100 km +4.67 / +4.90 K, EERIE +2.30 K (control) / +2.02 K (hist). Field mean, 5-95% from a winter bootstrap that uses the same winters in all regions:

    | comparison | wind only, %/K (24 regions) | solar + wind, %/K (18 regions) |
    |---|---|---|
    | 9 km | +5.9 [4.7, 7.1] | +4.2 [2.9, 5.4] |
    | 31 km, 9 km spans | +5.0 [3.9, 6.1] | +3.1 [1.2, 4.8] |
    | 31 km, 2081-99 | +4.6 [3.7, 5.6] | +3.8 [2.5, 5.1] |
    | 100 km, 9 km spans | +4.8 [3.3, 6.1] | +2.6 [0.4, 4.6] |
    | 100 km, 2081-99 vs 147-winter control | +5.5 [4.5, 6.3] | +3.5 [2.2, 4.6] |
    | EERIE SSP2-4.5 2041-70 vs control | +5.7 [3.9, 7.5] | +3.2 [0.8, 5.3] |
    | EERIE SSP2-4.5 2041-70 vs hist 1971-2000 | +5.9 [3.6, 8.3] | +4.4 [1.3, 7.6] |

  - **Per K the sensitivity is the same in every model, resolution and pathway: about +5-6% more low-wind events per K, +3-4% for the solar + wind mix.** EERIE's smaller ratio is entirely its smaller warming. The 9 km value is at the top of the range, with intervals overlapping all others.
  - The 9 km 1950C control drifts cold: global mean 13.25 C in 1950-54, 12.49 C in 1965-69 (-0.76 K in 20 years), consistent with the falling event counts within that control. dT uses the means of the compared years, so the drift is inside dT.
- [x] Mean shift or tail (same script). Future CF rescaled to the control mean and events recounted at the same threshold. Wind only: mean CF changes by only -0.3% to -4.2%, but the 5th percentile of the 48 h mean falls by 13-17% (AWI-CM3) and 7% (EERIE). Removing the mean change takes away 5-28% of the increase (9 km 1.37 -> 1.26, 31 km 1.29 -> 1.26, 100 km 1.30 -> 1.24, EERIE 1.14 -> 1.13). **It is a lower and/or more persistent low tail, not mean stilling.** Solar + wind: 15-45% from the mean.
- [x] Weather regimes, done 2026-11-05 (`dunkelflaute/scripts/regimes/`, output `dunkelflaute/results/regimes/regime_summary.json`, methods in docs/methods_notes.md 13). Four ERA5 Euro-Atlantic regimes from daily msl (NAO+, NAO-, Scandinavian blocking, Atlantic ridge; ERA5 Sep-Mar frequencies 27/24/26/23%), model days assigned to the ERA5 centroids.
  - Model control frequencies match ERA5 within 2-3 points in every model.
  - Under warming AWI-CM3 loses NAO- days (0.23-0.26 to 0.15-0.17) and gains NAO+ or AR; Scandinavian blocking is unchanged (0.27-0.28 to 0.24-0.28). With each period's own climatology the shifts nearly vanish (NAO- 0.22-0.23), so they are mostly the forced mean msl change, not a change in the weather sequence. EERIE: small shifts (NAO- 0.23 to 0.21).
  - Events per 100 regime days, control (mean over 8 European regions, wind only): ScBL 3.4-3.8, NAO+ 1.8-2.2, NAO- 2.0-2.7, AR 1.3-1.5. Blocking is the main Dunkelflaute regime, as in Mockert et al. and Grams et al.
  - **Decomposition: the increase is within-regime.** Mean over 8 European regions, wind only, events per winter: 9 km +1.36 = -0.20 (frequency) + 1.55 (within); 31 km +1.33 = -0.13 + 1.46; 100 km +1.08 = -0.07 + 1.15; EERIE +0.31 = -0.06 + 0.38. Events per regime day rise in all four regimes except 31 km NAO- (2.66 to 2.13) and EERIE hist AR (flat), most in ScBL (9 km 3.66 to 4.89, 31 km 3.36 to 4.91). Regime frequency changes slightly reduce events.
  - So the same circulation types produce more low-wind spells in the warmer climate. Candidate causes: weaker pressure gradients within regimes (reduced baroclinicity) or a more stable boundary layer at the same pressure gradient. ~~Next test: geostrophic wind from the msl gradient against 10 m wind, per regime, control vs future.~~ Done below: weaker pressure gradients, not a more stable boundary layer.
- [x] CMIP6 through the same pipeline, done 2026-11-05 (`dunkelflaute/scripts/cmip6/`: `inventory.py`, `extract_daily.py`, `daily_events.py`; output `dunkelflaute/results/cmip6/daily_events_summary.json`; methods in docs/methods_notes.md 14). 20 models with daily `sfcWind`, `rsds` and `Amon/tas` in historical and SSP5-8.5 (first member with all three). Daily footing (48 h = 2 days), historical 1951-1980 against 2071-2100 (and 2041-2070), base rate matched to 5 events/winter, 24 regions.
  - **Bridge: the daily footing gives the same per K as the 6-hourly one.** Wind only: EERIE daily +5.6 / +5.9 %/K (6-hourly +5.7 / +5.9), TCO95 daily +5.4 %/K (6-hourly +5.5), and +6.0 %/K with CMIP6-style periods (hist1950g 1951-80 vs 2071-2100). Solar + wind: EERIE +4.6 / +5.1, TCO95 +4.2 / +4.7 %/K. Null for the field mean (TCO95 daily control, 30 vs 30 winters): 0.98-1.02 (wind), 0.96-1.04 (combined).
  - **CMIP6, wind only: every model shows the increase.** Field-mean ratio 1.10-1.65 (median 1.31), 17-24 of 24 regions up per model. Per K: median +7.5 %/K, 10-90% +4.7 to +10.2. The ensemble regression of ln(ratio) on dT gives +7.8 %/K with intercept 0.98 (passes through the origin); r(dT, ln ratio) = 0.67; dividing by each model's warming cuts the spread modestly (CV 0.39 to 0.30). Mid-century (2041-2070) gives +8.6 %/K, so roughly linear in warming.
  - **Our runs sit in the lower half of CMIP6**: AWI-CM3 and EERIE +5-6 %/K against the CMIP6 median +7.5 %/K. The high-resolution runs are not more extreme than CMIP6.
  - **Mean vs tail differs:** in CMIP6 the winter mean wind CF falls more (-2 to -19 %, median about -9 %) and the mean explains about 40 % of the increase (20-58 %); in AWI-CM3 and EERIE the mean falls -1 to -3 % and explains 16-20 % on the daily footing. In both, the larger part is the tail.
  - **Solar + wind: weaker and less consistent.** CMIP6 median +4.4 %/K (10-90% +1.4 to +6.5), 18 of 20 models up (ACCESS-CM2 0.93 and MRI-ESM2-0 0.98 down). Normalising by dT barely helps (CV 0.63 to 0.56; r(dT, ln ratio) = 0.27), so the spread of the combined metric is model physics (the solar/cloud response), not warming.
  - **Germany in CMIP6:** wind only up in 19 of 20 models (median 1.29, +6.8 %/K); solar + wind up in only 11 of 20 (median 1.13, +2.9 %/K), one of the weakest regions in the combined metric (with BK, UA, IT, IB). That is consistent with Strnad et al. (5 CMIP6 models, Germany, no significant change): Germany's combined signal is small and model-dependent, while the wind signal is robust.
  - Regions, CMIP6 multi-model median %/K, wind only: from +3.0 (CN_E) and +3.9 (CA_PR) to +10.8 (IT) and +12.9 (IN_N); 20 of 20 models up in 13 of 24 regions.
  - On the daily footing all 24 regions reach the base rate for the combined metric too (daily means merge fewer low spells than the 6-hourly series).
- [x] Sub-daily check, done 2026-11-05 (`dunkelflaute/scripts/cmip6/extract_3hr.py`, `subdaily_events.py`; output `dunkelflaute/results/cmip6/subdaily_summary.json`; methods 14b). The 4 CMIP6 models with 3-hourly `uas`/`vas`/`rsds` (MPI-ESM1-2-HR, MPI-ESM1-2-LR, IPSL-CM6A-LR, MRI-ESM2-0). Winters 1961-1990 vs 2071-2100, because MRI 3hr output starts in 1960. Each model on three footings: 6-hourly as our runs, daily from the 3-hourly data, and the daily files.
  - **The footing does not matter.** Wind only, %/K on 6h / daily-from-3h / daily files: IPSL +7.5 / +7.9 / +7.6, MPI-HR +7.9 / +7.7 / +7.7, MPI-LR +7.1 / +7.0 / +7.2, MRI +4.9 / +5.5 / +5.6.
  - The regional patterns correlate r = 0.91-0.99 between footings. Solar + wind is within 0.6 %/K across footings, MRI negative on all three.
- [x] Geostrophic vs surface wind, done 2026-11-05 (`dunkelflaute/scripts/cmip6/extract_psl.py`, `dunkelflaute/scripts/regimes/geostrophic.py`; output `dunkelflaute/results/regimes/geostrophic_summary.json`; methods 13b).
  - Method: daily geostrophic speed G from the 2.5 deg msl gradient per region, and a pseudo-surface speed s* = (control surface/geostrophic ratio per calendar month) x G. Events in s* give the circulation part of the change.
  - Data: our runs for the 8 European regions; 20 CMIP6 models (`day/psl`) for the 22 regions north of 29 N.
  - **The increase is circulation: weaker pressure gradients, not a more stable boundary layer.**
    - CMIP6: actual ratio 1.34, circulation-only 1.29, median share from circulation 0.91 (10-90% 0.60-1.22). Mean G -3.6 %, transfer ratio unchanged (median 0.0 % on 22 regions; -0.9 % on the 8 European regions).
    - Our runs: share 1.32 (9 km), 1.06 (31 km), 0.54 (100 km), 1.10 / 1.55 (EERIE ctl / hist). Mean G -1.4 to -4.6 %.
  - **In our runs the transfer goes up** (+1.6 to +4.1 %), partly offsetting the gradient loss at the surface. This is a candidate reason why AWI-CM3 and EERIE lose less mean wind than CMIP6.
  - **Within regimes, G weakens in every regime in the CMIP6 medians** (not in every run under NAO+/-: 31 km NAO- +3.2 %). It weakens most under Scandinavian blocking in our runs (9 km -11.5 %, 31 km -7.6 %). In CMIP6 the largest weakening is under Atlantic ridge and blocking (-5.4 %, -4.5 %). Low-s* days rise within every regime. This is the mechanism behind the within-regime result above.
  - Still open: whether the 100 km run's larger transfer part (1.11) is real or the coarse msl, and a stability diagnostic (`sshf`/`blh` from the TCO95 rerun) to test the transfer term directly.
- [x] Gridded per-K maps, done 2026-11-05 (5 deg land boxes, 20-70 N; `REGION_SET=grid5` in `dunkelflaute/scripts/regions/regions.py`, `grid_per_K.py`, `dunkelflaute/scripts/plotting/plot_per_K_maps.py`; figures in `dunkelflaute/results/maps/`; methods 15).
  - **The increase is broad: 70-80 % of land boxes go up in every comparison.** Wind only, area mean per K: CMIP6 median map +3.7, 9 km +3.5, 31 km +2.8, 100 km +2.2, EERIE +4.7 / +4.8. At least 80 % of CMIP6 models agree on an increase in half the boxes; the Sahara and Arabia go down.
  - The CMIP6 signal peaks at 35-40 N (+7.1 %/K). Pattern correlation with the CMIP6 median is 0.31-0.63 for our runs, against 0.60 for a single CMIP6 model against the rest.
  - **In CMIP6 the change grows with aggregation scale:** +4 %/K for 5 deg boxes, +7 for 20 deg tiles, +9 for 40 deg tiles. The mean-removed part stays at +2 to +4, so the growth is the mean wind decline. Our runs stay flat from 5 to 20 deg (smaller mean decline). This matters for grid-scale balancing, because Europe-wide low-wind spells grow faster than local ones. In CMIP6, Europe as a whole is +9.6 %/K.
  - [x] Done 2026-11-06 (methods 15a, 15b). Field significance: every run's area mean lies above all null placements; 99-252 of 420 boxes FDR-significant. Variance explanation holds (mean-only part = elasticity x mean decline, r = 0.99). Europe-wide spells (>= 4 of 8 regions) +8 to +16 %/K in our runs, +25 %/K in CMIP6, mostly from more events per region; coherence itself +2 to +6 %/K.
- [x] Why CMIP6 gives more than AWI-CM3, 2026-11-05 (methods 14, last two bullets).
  - **Across CMIP6, per K follows the mean wind decline:** r = -0.90, per K = 3.0 - 2.07 x (mean CF change, %/K). AWI-CM3 and EERIE lose little mean wind (-0.3 to -0.75 %/K against a CMIP6 median of -1.9), and they sit above that line. So their increase is mostly a broader low tail, while in CMIP6 a large part is the mean shift.
  - **No resolution effect in CMIP6 pairs:** the finer member is slightly lower in all three pairs.
  - **It is not an IFS trait:** EC-Earth3 (IFS cy36r4, wind only, outside the 20-model ensemble) is the most extreme CMIP6 model. It gives +12.6 %/K, mean CF -16.9 %, mean G -7.9 % and a flat transfer.
  - Paper framing: the increase does not need a mean wind decline. AWI-CM3 is the conservative case.
  - Still open: why our runs' transfer rises. Candidates are newer IFS physics, FESOM2, or the AWI setup. Tests: the stability diagnostic, and the TCO95 historical transfer trend against ERA5 and station data.
- [ ] Turning the rising transfer into a finding (proposed 2026-11-06). Do the steps in order; step A can kill it.
  - [x] A done 2026-11-06 (methods 13c): **not an artifact.**
    - The rise survives every check, at or above the CMIP6 90th percentile in every variant (EERIE hist on 22 regions is the exception).
    - On the 22 regions: +0.8 to +2.7 % (CMIP6 0.0 %), rising in 18-20 of 22 regions (CMIP6 median 10). Per K: +0.4 to +0.6 %/K against 0.0.
    - At native 100 m it is about a third smaller (+1.7 % at 100 km).
    - Side finding: G over IB, IT and BK depends on how msl is put onto the 2.5 deg grid.
  - [x] B done 2026-11-06 (methods 13d, TCO95). **The cause is weaker gradients plus a less stable boundary layer.**
    - Weaker G gives about +1.0 %, the same as CMIP6.
    - A thermodynamically less stable boundary layer gives about +1.6 %: inversion -0.3 K and blh +5 % at fixed wind. Joint G x stability classes explain all of it (within-class part -0.1 %).
    - Snow loss amplifies it (+3.4 % in snowy region-months, +1.4 % elsewhere), but snow-free IB and FR rise too.
    - Ruled out: roughness change (wrong sign), the 10 m diagnostic (+0.35 %), control drift.
    - CMIP6 at fixed G: -2.0 % (EC-Earth3 -4.3 %) against TCO95 +1.5 %. That term is where the models differ.
    - No signal in 1951-1980 against 1991-2020, so D has to use the interannual sensitivity, not the trend.
  - [x] C done 2026-11-06 (methods 13e). **Not stability; mostly land-use roughening.**
    - CMIP6 destabilises like our runs (22 regions dtheta -0.42 K; TCO95 -0.55, EERIE -0.24). Stability, hfss, snow and warming ratio do not predict the at-fixed-G term (|r| <= 0.24, wrong sign for stability).
    - Crop fraction change predicts it (r = +0.89, n = 10); tree fraction and LAI too. INM (no land-cover change) sits at -0.1 to -0.3 %.
    - At zero land-cover change, our runs are still 1 to 1.5 points above the CMIP6 fit. So about 2 points of the gap is SSP land-use forcing that IFS does not have, and 1 to 1.5 points is physics (IFS transfer responds to destabilisation).
    - Paper consequence: the CMIP6 10 m wind decline over Europe partly reflects prescribed land use, which matters less at hub height. Worth a sentence and maybe a check of sfcWind against ua at the lowest level.
    - Gaps: KACE and EC-Earth3 have no land-cover output; no model writes roughness.
  - [x] D done 2026-11-06 (methods 13f). **No observational constraint.**
    - Interannual sensitivity at fixed G does not predict the forced change across models (r = -0.27 over CMIP6, +0.02 over all 25). Not an emergent constraint.
    - HadISD (423 stations, 1980-2025): -0.14 %/K [-0.94, +0.60], against our runs +0.6 to +0.9 and ERA5 +0.76. Stations minus ERA5 at the same points: -0.62 [-1.30, -0.02]. So IFS-type boundary layers (ours and ERA5) overstate the warm-winter coupling, but that does not transfer to the forced change.
  - Outcome of A-D for the paper:
    - The rise is real (A) and has a mechanism (B: weaker G plus destabilisation).
    - About 1.4 points of the CMIP6 gap is prescribed SSP land-use roughening (C, measured directly in 13g; C's cross-model estimate was 2). That is sellable: a forcing difference, and it acts on 10 m wind, less at hub height.
    - The remaining physics difference is unconstrained by observations (D), with a mild warning sign that IFS couples too strongly in warm winters. Present it as structural uncertainty, not as AWI-CM3 being right.
  - [x] LUMIP direct test done 2026-11-06 (methods 13g). **Land use confirmed, but smaller: about -1.4 % of the CMIP6 at-fixed-G term (range -1.0 to -1.9), not 2 points.**
    - Crop to forest (CanESM5, MPI-LR, BCC): 10 m wind -2.2 to -3.3 % at fixed circulation, 3-5 standard errors. Crop to grass (UKESM): no drop.
    - Per pp: tree -0.33 %/pp matches the cross-model slope; crop +0.25 %/pp is half of it.
    - So CMIP6 at zero land-cover change is about -1.0 %, and the physics gap to our runs grows to 2.4-3.0 points.
  - [x] Hub height done 2026-11-06 (methods 13h). **The CMIP6 10 m decline about doubles the hub-height decline over Europe.**
    - CMIP6 median, 8 regions: sfcWind -5.1 %, V1000 (typically 100-140 m on valid days) -2.2 %, model-level 100 m -2.7 % (4 models).
    - The difference is near-surface shear (s10/V1000 -2.4 %), which tracks crop loss (r = +0.83) and tree gain (r = -0.86). V1000/G is +2.0 %, close to TCO95's V100/G +1.7 %, with no land-cover signal.
    - Crude: CMIP6 per K at hub height about +5 %/K instead of +7.5; TCO95 about +6.1 instead of +5.4. They roughly meet.
    - Weak spots: 1000 hPa is above ground on only 10-47 % of European winter days (calmer days); model-level V100/G is only +0.2 % in the 4 models that have it.
  - Revised outcome: the model difference sits mostly in the surface layer, where CMIP6 roughens with land cover and IFS does not. At about 100 m the transfer change at fixed G agrees. For the paper: compute wind CF at hub height for CMIP6 (from V1000 or model levels where available) rather than log-law from sfcWind, or at least show the sensitivity. Next decision: is that worth rebuilding the CMIP6 event pipeline on 1000 hPa / CFday winds for the subset of models that have them?
  - A. Is it real?
    - The same regions and G method as CMIP6. CMIP6 transfer over just the 8 European regions; our G from bilinear msl as for CMIP6.
    - Our runs over all 22 regions: global daily msl from the TCO95 rerun (1d), TCo319 and EERIE (6h).
    - The raw 10 m speed instead of the log-law 100 m one.
    - For TCO95, the native 100u/100v. Does the rise survive above the 10 m diagnostic?
  - B. What causes it? Use the TCO95 rerun (ctl1950g/hi1950g, 1950-2099); its 3-hourly output has 10u/10v, 100u/100v, 2t, skt, sshf, blh, msl and model levels 85-91.
    - Stratify the transfer by stability (bulk Richardson number from the lowest levels, the sign of sshf, blh). Is the change a shift between stability classes or a change within a class?
    - Snow and roughness: is the change in transfer largest in the regions and months losing snow? Does `fsr` change?
    - The 100 m / 10 m speed ratio.
  - C. Why does CMIP6 differ? For 21 models, set the transfer change against the change in lower-tropospheric stability (Amon `ta` 850 hPa minus `tas`, all 21) and in snow cover (LImon `snc`, 15 of 21). Either our runs sit at the end of a common relationship (a larger stability loss), or they lie off it (a different boundary-layer sensitivity).
  - D. Which is right? Observations.
    - The interannual sensitivity of the transfer to winter temperature or snow, detrended: station winds against G from ERA5 msl, compared with each model's control-period sensitivity. Test first, in the models, whether the interannual sensitivity predicts the forced change. That would make it a constraint.
    - The 1950-2024 transfer trend in hi1950g against stations, with the stilling and roughness caveats.
    - ERA5 land 10 m wind is IFS cy41r2 physics, not an observation, so stations are needed (HadISD, DWD CDC; not on Levante).
- [ ] Next: regional 2t during events ("darker but warmer"); a figure set for the per-K result.
- [x] EERIE IFS-FESOM2-SR (TCo1279 + FESOM2 NG5, AWI), a second 9 km model with other winters. Checked on Levante 2026-11-05:
  - Runs: hist-1950 1970-2014 (r1, plus r2/r3), eerie-control-1950 1971-2050, highres-future-ssp245 2015-2070.
  - 6-hourly on 0.25 deg (same grid as ERA5, so the QM calibration applies unchanged): `10u`, `10v`, `tcc`, `lcc`, `mcc`, `hcc`, `msl`, `2t`, and 6 h accumulated `ssrd` and `tisr`. Native-grid 6 h exists too (ssp245: instant fields only).
  - Native `ssrd` means no solar reconstruction. Only the log law for 100 m.
  - Raw data: FDB GRIB under `/work/bm1344/a270228/fdb_hist1950/REGULARLL_grids/...` (readable, regular_ll 1440 x 721). Catalogue: `/work/bm1344/DKRZ/intake_catalogues/dkrz/disk/model-output/ifs-fesom2-sr/`; the kerchunk parquet refs did not open with zarr 3 (module python3/2025.01), so read the GRIB directly or use an older zarr.
  - Limitation: SSP2-4.5 to 2070 is a weaker forcing than SSP5-8.5 in the 2080s. Compare at matched warming, or use control 1971-2050 against ssp245 2041-2070.
- [ ] EERIE ICON-ESM-ER (~10 km): 6 h and 3 h `uas`/`vas`/`psl`, but no sub-daily radiation or cloud. Wind-only calm-spell check.
- Not usable: DestinE storyline copy at `/work/bk1415/a270229/awi_hackathon` (2017-2024 only; sub-daily wind and cloud folders empty, only hourly 2 m T; nudged storylines fix the circulation). EERIE HadGEM3-GC5 on DKRZ has daily atmosphere only. EERIE IFS AMIP TCo399/1279/2559/3999 are historical only: useful for resolution vs climatology, not for change.
- In the hackathon catalogue but not yet checked: `TCo1279-DART-2000` (atmos entry commented out, ocean only) and `TCo319-DART-hi1950d` (remapped daily and monthly only).

### Phase 3: mechanism and impact

- [x] Weather regimes ~~or blocking from 6-hourly `msl`~~ from daily `msl` at all resolutions and EERIE, checked against ERA5. Split the change in event frequency into regime frequency versus within each regime. Done 2026-11-05, see Phase 2b: the increase is within-regime. A Z500 blocking index is still open (TCO95 daily z500, EERIE 3D 6-hourly, TCo319 monthly only).
- [ ] Composites of `msl` and `2t` during events, control vs future, at each resolution.
- [ ] "Darker but warmer?": compare temperature anomalies during events across climates. Fit a temperature-to-demand model to SMARD load (to be downloaded), then do residual-load-based detection as recommended by Biewald et al. 2025.
- [ ] Onshore vs offshore CF during events.

### Phase 4: fleet sensitivity and paper

- [ ] Recompute with the 2018, 2024 and 2045-target mixes, and compare the size of the fleet effect with the climate effect.
- [ ] Li definition at all resolutions, as a robustness check. (Fleet mixes and the full definition intercomparison are parked for paper 2; see "Paper 2 parking lot".)
- [ ] Paper draft 2026-11-06: paper repo as submodule `paper/`, local branch `draft-outline` (not pushed). Target Nature Climate Change, then Nature Communications. All sections drafted, 9 figures, builds clean. Open: Jan's review, main vs Extended Data split, SH contrast, author affiliations.
  - 2026-11-07, after Jan's questions.docx: solar-only analysis (methods_notes 15c), SH test (13i), figure revision (13: per-K decomposition, centroid check). NCC layout: Intro, Results, Discussion, Methods, 7 Extended Data figures, Supplement. Main figures: maps, perK (wind | solar | mix | mean), mechanism, southern, transfer, europe. Style pass (no fragments, no colons in flow text, no far references) done. Affiliations done. Open: main text 4200 words against the 3000 limit; aerosol question (do AWI-CM3 and EERIE prescribe time-varying aerosols?); Qu2025 season check; title still says wind-and-solar.
- Title 2026-11-08 (Jan): "Warming brings more winter wind-and-solar droughts as blocked flow weakens" (descriptive term instead of Dunkelflauten in the title; no similar title in Crossref/OpenAlex).
- [ ] Paper (earlier title): ~~"Resolution dependence of projected Dunkelflaute change in Germany" (AWI-CM3 at 100, 31 and 9 km).~~ Working title (2026-11-05): "More winter wind-and-solar droughts in a warmer climate: a Northern-Hemisphere signal robust across model resolution" (AWI-CM3 at 100, 31 and 9 km plus EERIE IFS-FESOM2-SR), with Germany as the worked example and its weak signal at 31/100 km as a sampling lesson. Scope and outline: see "Current direction" at the top (global, per-K maps, CMIP6 through the same pipeline, NH vs SH as a mechanism test).

### Later

- [ ] Pan-European coherence of events, using the country downloads for validation. Moved into paper 1, section 4 (whether events become more Europe-wide at once).
- [ ] South Korea (second paper, with IBS; see "Paper 2 parking lot"): capacity weights from EPSIS, hourly validation data, and a Sep-Mar vs year-round season choice given the monsoon.

## Notes for the paper

Detailed methods log, updated with every step: [`docs/methods_notes.md`](docs/methods_notes.md).

- Limitation, one sentence: real output can fall below the meteorological index during cold Dunkelflaute because of turbine icing (mostly low-mountain hilltops, rare on the northern lowlands and offshore) and snow on PV panels (small, since winter solar CF is only ~0.03-0.08). Both make real events slightly deeper and become rarer under warming, the opposite direction to the projected increase. Not modelled; to quantify, use the share of event hours with `2t` below 0 °C in 1950C vs 2080C from the Phase 3 temperature analysis.
- TCO95 is a different model configuration, not just a third resolution. The 9 and 31 km d runs (Aleph) used the aleph_production TCO319 tuning (RPRCON 1e-3, ENTRORG 1.25e-3, GGAUSSB -0.45). At TCO95 that tuning did not work at all (Jan, 2026-11-02): ctl1950df drifted, 2t 14.36 to 12.99 C over 1950-68. The production runs `ctl1950g` and `hi1950g` therefore use release-style AWI-CM3 v3.1 parameters (`namelists/fesom_own_gm`, aerosol-scaling OIFS branch, no TCO319 tuning), spun up 20 years from `own1950g`. In the paper, present TCO95 as the perfect-model and natural-variability testbed, and flag that any 100 vs 31/9 km difference mixes resolution with tuning. Details in `/work/ab0995/a270092/runtime/awicm3-frontiers-xios/` (drift figures `ctl1950df_drift.png`, `own1950a_drift.png`, `own1950g_drift.png`).
- Curtailment (grid congestion, negative prices) does not affect events: it happens under oversupply, not when wind and solar are both low. It only lowers SMARD CFs in high-output hours.

## Literature positioning (checked 2026-11-05)

- Already published, direction only: NH mid-latitude wind droughts get longer (Qu et al. 2025, Nat. Clim. Change 15, 842-849; CMIP6, 20th-percentile droughts, about +10/15/20% by SSP, attributed to fewer cyclones and weaker meridional gradients); European Dunkelflaute +5-25% at 2 C (Brayshaw et al. 2026, Met. Appl. e70192, C3S/ECEM); Germany slightly more Dunkelflaute at 3 and 12 km, attributed to the mean wind decline (Heinrich & Geyer 2026, ERL 21 124022; NUKLEUS, SSP3-7.0, 2/3 K). Germany flat: Strnad et al. 2025/26 (CMIP6 plus generative downscaling), DWD.
- Resolution: Liu et al. 2024 (GRL 51, HighResMIP ~25 vs ~100 km, W. North America/Texas) find that resolution changes the projected change. Heinrich & Geyer compare 3 and 12 km but do not test the change signal.
- Not found anywhere: (1) a Dunkelflaute/wind-drought projection with global coupled km-scale models; (2) a resolution test of the change signal at 9/31/100 km in one model family; (3) the per-K collapse across resolutions, models and pathways; (4) the regime decomposition showing a within-regime increase; (5) the quantified mean-vs-tail split (closest: heavier-tailed low-solar-low-wind distributions in China, Natl Sci Rev 2025); (6) the sampling lesson that a single region with 10-20 winters misses a common signal, which reconciles the conflicting Germany studies.
- Not read in full (access): "Global warming intensifies compound renewable energy droughts" (Commun. Earth Environ. 2026, warming levels 1.5/2/3 C, solar as the main driver; probably the work of Aminimehr, Hellinckx & Tabari, EGU26-13163); "Widespread increase in frequency and duration of European wind droughts based on CMIP6 projections" (iScience 2026, fewer weak-wind events in winter in some regions); Qu et al. 2025 full text.

## Open questions

- [ ] Where is the long MR control (the proposal mentions 184 years) and the MR SSP5-8.5 transient (2015-2099), with sub-daily wind?
- [ ] Does TCo2559-DART-1950C have atmospheric output anywhere?
- [x] Is Amanda continuing, and is she a co-author? (Jan 2026-11-08: no longer our intern, but she is a co-author of paper 1 and will read and edit the draft.)
- [ ] When does Korea need to happen for the IBS collaboration?
