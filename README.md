# 🌙 AWI Dark Doldrums Project

> Changes in Dunkelflaute (dark-doldrum) risk under climate change: a 9 km-resolution pipeline that reconstructs the missing wind and solar fields in TCo1279-DART, validates the result against ERA5 reanalysis and real German grid data (SMARD), and compares detected events against four published definitions from the literature.

---

## 📌 Project Purpose

A *Dunkelflaute* is a period when wind and solar output are both low at once — the case that matters most for grid resilience, since the two resources are normally expected to compensate for each other. Germany is where the term originates; under the Mockert et al. (2023) definition used here, ERA5 1979–2018 gives about 4.3 events per year with a median length of 3.2 days and a longest event of 8 days. This project builds the pipeline needed to detect these events in a high-resolution (~9 km) climate simulation, so that how Dunkelflaute risk changes under warming can actually be assessed.

The catch: **TCo1279-DART doesn't carry everything needed directly.** It writes 100 m wind only as 10 m (`u10`/`v10`), and solar irradiance (`ssrd`) only as a monthly mean, not 3-hourly. Both gaps have to be filled by a reconstruction trained and validated on ERA5 — which does have the real fields — before being trusted on DART, where there's no ground truth to check against directly. That reconstruction-and-validation work is most of what this repo contains; detecting and characterizing Dunkelflaute events is what it's *for*.

This is an extension of the proposal in [`Project_Proposal.md`](Project_Proposal.md), a 4-week study originally scoped to cover Germany **and** South Korea, both the 1950-control and 2080s-SSP5-8.5 climates, sub-national spatial patterns, synoptic drivers, and formal uncertainty quantification. The **Status** section below says plainly which of that is done, and which is still ahead.

---

## ✅ Status at a Glance

| Proposal item | Status |
|---|---|
| Wind reconstruction (10 m → 100 m), validated vs. ERA5 + SMARD | **Done** — log law, see [`wind/README.md`](wind/README.md) |
| Solar reconstruction (monthly → 3-hourly `ssrd`), validated vs. ERA5 + SMARD | **Done** — ML model v5, see [`solar/README.md`](solar/README.md) |
| Real Germany boundary (land + EEZ), not a bounding box | **Done** — [`boundaries/`](boundaries/) |
| Pipeline validation against ERA5 (proposal §5.8) | **Done, with a gap** — Mockert threshold gives 5.64 events/winter on ERA5 2016–2026 (6.82 with real ERA5 100 m wind and `ssrd`, 5.82 from SMARD), about 30% above Mockert's published ~4.3/year for 1979–2018; Li/Kaspar/Biewald comparisons also run, see Literature Comparison below |
| Baseline (1950C) and future (2080C) Dunkelflaute climatology for Germany, Mockert threshold | **Done** — see Key Findings below |
| Same comparison under the Li et al. (2021) dual-technology threshold, applied to DART | **Not yet** — Li/Kaspar/Biewald so far only validated on ERA5, not re-run on 1950C/2080C |
| South Korea domain | **Not started** — only a `dunkelflaute/context/` visualization of the proposal's pattern-scaling coefficients exists, not an actual Dunkelflaute detection run |
| Sub-national spatial analysis (coast vs. inland vs. mountains) | **Not started** |
| Synoptic composite analysis (mean sea-level pressure during events) | **Not started** |
| Onshore vs. offshore wind behaviour specifically *during events* | **Not started** (onshore/offshore wind CF is validated generally; event-conditional breakdown isn't built) |
| Threshold/window sensitivity tests (CF 0.04/0.08/0.10, 24h/72h) | **Not started** |
| Uncertainty quantification (year-block bootstrap, internal consistency) | **Not started** |

In short: the hard validation work — proving the reconstructed wind and solar fields can be trusted — is done and holds up well. The comparative climate-change science the proposal was ultimately for (Germany vs. Korea, synoptic mechanisms, uncertainty bounds) is the work that's still ahead.

---

## 🔑 Key Findings

**Solar reconstruction is essentially at the ERA5 ceiling, and the raw model is already doing the work before any correction.** The v5 model reproduces real ERA5 `ssrd` at r = 0.998 over Germany, and real SMARD generation at r = 0.931 — matching ERA5's *own* ability to predict SMARD. Even before the exact monthly rescale that forces each month's total to match truth, the model alone reaches R² = 0.994. Its one real weak spot is terrain it's never trained on — Antarctica, high mountains — where raw R² drops to 0.36, now quantified for the first time ([v5 Raw vs Rescaled](https://claude.ai/artifact/Yb8k4qxYzs8wEBebFfU9cr)). See [solar/README.md](solar/README.md#-key-findings).

**Wind reconstruction is solid onshore, noisier offshore — for reasons unrelated to the method.** Log-law 100 m wind matches real ERA5 at RMSE ≈ 1.19 m/s onshore / 0.60 m/s offshore (r = 0.97 / 0.997). Against real SMARD generation, onshore capacity factor reaches r = 0.91 (vs. real ERA5's own r = 0.95); offshore pools at r ≈ 0.70–0.71 for *both* the reconstruction and real ERA5, but r = 0.86 / 0.87 once the 2015–17 offshore capacity-step years are excluded, so the pooled figure mostly reflects SMARD capacity data rather than the reconstruction. See [wind/README.md](wind/README.md#-key-findings).

**The pipeline reproduces realistic Dunkelflaute statistics, and detects a real climate-change signal for Germany.** Using the Mockert et al. (2023) threshold (48h rolling combined CF < 6%, Mockert's own weights: Germany's 2018 mix, 44% solar / 50% onshore / 6% offshore):

| Period | Winters | Events/winter | Mean combined CF |
|---|---|---|---|
| ERA5 (validation, 2016–2026) | 11 | 5.64 | 0.186 |
| **DART 1950C (control)** | 17 (1953–1969) | **3.47** | 0.201 |
| **DART 2080C (SSP5-8.5 future)** | 10 (2083–2092) | **5.50** | 0.199 |

The ERA5 validation figure (5.64/winter) is about 30% above Mockert's published ~4.3/year (1979–2018). Their events fall almost entirely in SON and DJF, so counting Sep–Mar only does not explain the gap; their capacity factors are bias-corrected against observed generation, this pipeline's are not. The reconstruction also counts 17% fewer events than real ERA5 fields do (5.64 vs 6.82), which [`PLAN.md`](PLAN.md) addresses. Comparing the two DART endpoints: event frequency rises from 3.47/winter under 1950-control forcing to 5.50/winter under the 2080s SSP5-8.5 future — a **~59% increase**, and the future figure lands close to where present-day ERA5 already sits. This is a genuine, if single-threshold, answer to proposal questions 1 and 2 (baseline climatology and the climate-change signal) for Germany.

How robust it is (`dunkelflaute/scripts/analysis/robustness_mockert_events.py`, output in `dunkelflaute/results/robustness/robustness_summary.json`):
- ~~Significant at the headline threshold~~ Significant only by within-run tests (2026-11-05: against block nulls from long control runs Germany alone gives p = 0.13-0.20, and the increase turns out to be a broad Northern-Hemisphere signal at all three resolutions; see PLAN.md Phase 2b): exact Poisson test p = 0.015, and resampling whole winters puts the 2080C/1950C ratio at 1.58 with a 90% interval of 1.19 to 2.12. Germany's 2024 weights give a similar ratio (3.88 to 6.00).
- It is a shift in the low tail, not the mean. Mean combined CF is unchanged (0.201 to 0.199), while the 5th percentile of the 48 h mean drops from 0.070 to 0.061, close to the threshold itself. Onshore CF below 0.05 goes from 14% to 17% of the time.
- The ratio depends strongly on the threshold: at 48 h it is 2.32 at CF < 0.04, 1.58 at 0.06, 1.18 at 0.08 and 1.03 at 0.10. A Li-style count (wind and solar CF both below 0.2 for more than 24 h) rises only 16%, from 9.35 to 10.8 per winter.
- DART has fewer events than ERA5 at every threshold, so 2080C landing near today's ERA5 is not meaningful on its own.
- It does not depend on the wind reconstruction's error. The log law overestimates 100 m wind against ERA5's native 100 m wind (onshore CF 0.252 vs 0.207), so it undercounts events. Quantile-mapping it onto ERA5 100 m wind (`wind/scripts/methods/calibration/`, results in `dunkelflaute/results/wind_calibration/`) fixes the mean and most of the calm-hour bias. Applied to DART, it moves events from 3.47 to 3.59 (1950C) and from 5.50 to 6.00 (2080C), a ratio of 1.67 instead of 1.58, with the same threshold dependence.

Source data: `dunkelflaute/results/{1950c,2080c,validation_era5_2015_2026}/summary.json`; write-up: `dunkelflaute/results/report_dunkelflaute.html`.

---

## 📂 Layout

- **[`solar/`](solar/)** — reconstructs 3-hourly `ssrd` from the fields DART *does* carry 3-hourly (cloud fractions + `tsr`), via a `HistGradientBoostingRegressor` predicting clearness index, then an exact monthly rescale to match DART's real monthly total. `solar/models/` holds five training iterations, v1→v5; **use v5**. `solar/reference/` documents the OpenIFS solar-geometry physics this reimplements.
- **[`wind/`](wind/)** — reconstructs 100 m wind (`u100`/`v100`) from 10 m wind via the log law, using ERA5 land roughness (year-cycled) over land and the Charnock relation over ocean. `wind/reports/FINDINGS.md` has the log-law-vs-regression validation writeup.
- **[`dunkelflaute/`](dunkelflaute/)** — combines the solar + wind capacity factors into Germany's combined CF and detects Dunkelflaute events. This is where the actual science question lives:
  - `scripts/core/capacity_factor.py`, `domain.py` — shared CF formulas and the real Germany land+EEZ boundary mask (from [`boundaries/`](boundaries/), not a bounding box).
  - `scripts/era5/` — the ERA5-side pipeline. Start with `compute_germany_dunkelflaute_2015_2026_mockert.py`, the current, validated version.
  - `scripts/dart/` — the DART-side pipeline (1950C/2080C), Mockert-weighted; `compute_germany_dunkelflaute_dart_mockert.py` + `make_report.py`.
  - `smard_validation/` — SMARD real-generation data and CF computation, used as ground truth throughout the ERA5 validation.
  - `results/` — finished 1950C/2080C/ERA5-validation summaries and figures (the Key Findings table above comes from here). `reports/` — the earlier Mockert/Li/Kaspar/Biewald literature-comparison write-ups.
  - `context/` — not part of the reconstruction/validation pipeline itself, but the motivation for the 1950C-vs-2080C comparison: TCo1279-DART's climate pattern-scaling coefficients (local change in wind speed and cloud cover — the two Dunkelflaute drivers — per 1K of global warming), visualized from an externally-supplied dataset. Also shows the South Korea domain from the proposal, alongside Germany, though no Dunkelflaute detection has been run there yet.
- **[`boundaries/`](boundaries/)** — Germany's real land (Natural Earth) and EEZ (Marine Regions) polygons, used by every Germany-domain script instead of a lat/lon bounding box.

Everything else at the repo root is either a fixed project file (`LICENSE`, `Project_Proposal.md`, this `README.md`) or data/tooling excluded from version control — see Repo Hygiene below.

---

## 🧹 Repo Hygiene

- **`scratch/`** (gitignored) — a working directory for one-off, exploratory analysis that isn't part of the documented pipeline: draft versions of HTML reports before they're published and intermediate JSON summaries not referenced anywhere else. The rule is firm: **if a scratch check's result ends up quoted in a README as a finding, the script and its output get promoted out of `scratch/`** into the real structure (tracked script in `scripts/evaluation/`, output under the relevant `models/model_vN/evaluation/`) — not left behind as the uncited source of a cited number. Treat `scratch/` as a notebook margin, not a fourth project folder.
- **Superseded outputs are archived with a README, not scattered.** `dunkelflaute/data/germany_era5/archive/` is the model for this: old combined-CF runs that predate the current methodology are kept (for traceability) in their own `archive/` subfolder with a README explaining exactly why each one is superseded and that nothing should be built on them. Follow this pattern rather than leaving stale output next to current output with no explanation.
- **2026-10-01 cleanup**: removed three directories that didn't follow either rule above — `solar/data/dart_eval_1950_check/`, `dart_eval_1950/`, and `dart_eval_1950_winter1/` were pilot/diagnostic runs from the earliest DART evaluation attempts (Sep 21–22, scoring the old v1/v2 models, before the v5 model and the production `dart_reconstruct_*` pipeline existed), undocumented and unreferenced by anything except the job scripts that produced them. Removed those job scripts (`solar/jobs/dart_eval_year.sbatch`, `dart_recon_year.sbatch`, `dart_winter_1950.sbatch`) along with their logs. Also removed `dunkelflaute/scripts/compute_germany_combined_cf.py` and its output `dunkelflaute/data/germany/` — the script was already documented as superseded (see the solar README's Recommendations table) and its single output directory, confusingly similarly named to `dunkelflaute/data/germany_era5/`, had no archive note of its own.
- **2026-10-02 follow-up**: applied the `scratch/` promotion rule above to two checks that were already being cited as findings — `heldout_germany_v4v5.py` and `prerescale_vs_postrescale_germany.py` moved into `solar/scripts/evaluation/`, their results into `solar/models/model_v4/` and `model_v5/evaluation/`, and the duplicates removed from `scratch/`.

---

## 📖 Literature Comparison

Four published Dunkelflaute definitions are checked against this pipeline's ERA5 output — algorithm faithfully reproduced from each paper's own text, not a secondhand summary:

- **Mockert et al. (2023)** — 48h rolling-mean combined CF < 6%, their exact window-expansion event-construction rule, their capacity weights (Germany 2018: 44% solar, 50% onshore, 6% offshore). ERA5 1979–2018: about 4.3 events/year, median 3.2 days, longest 8 days, mostly under blocked weather regimes. The only definition also run on DART 1950C/2080C so far (see Key Findings).
- **Li et al. (2021)** — instantaneous wind CF < 20% *and* solar CF < 20%, sustained > 24h, no smoothing. Renewables.ninja (MERRA-2) CFs, 1985–2016; about 5–10 events/year longer than one day in Germany.
- **Kaspar et al. (2019)** — the paper Mockert calibrated their threshold against; CF ≤ 10%, ≥ 48h, no smoothing. COSMO-REA6 wind and SARAH-2 PV, 1995–2015.
- **Biewald et al. (2025)** — cited as "Lohmann et al." in earlier versions of this repo and the proposal. Not a new definition but a report card on the others, scored against shortage days from the ERAA 2023 adequacy simulation (PECD v3.1, 1982–2016, 2033 capacities): Li's per-technology method peaks at F = 0.14 (23% threshold, 24h), high residual load reaches F = 0.41.

References:
- Mockert, F., Grams, C. M., Brown, T., and Neumann, F.: Meteorological conditions during periods of low wind speed and insolation in Germany: the role of weather regimes, Meteorol. Appl., 30, e2141, https://doi.org/10.1002/met.2141, 2023.
- Li, B., Basu, S., Watson, S. J., and Russchenberg, H. W. J.: A brief climatology of Dunkelflaute events over and surrounding the North and Baltic Sea areas, Energies, 14, 6508, https://doi.org/10.3390/en14206508, 2021.
- Kaspar, F., Borsche, M., Pfeifroth, U., Trentmann, J., Drücke, J., and Becker, P.: A climatological assessment of balancing effects and shortfall risks of photovoltaics and wind energy in Germany and Europe, Adv. Sci. Res., 16, 119–128, https://doi.org/10.5194/asr-16-119-2019, 2019.
- Biewald, Cozian, Dubus, Zappa, and Stoop: Evaluation of 'Dunkelflaute' event detection methods considering grid operators' needs, Environ. Res.: Energy, 2, 025007, https://doi.org/10.1088/2753-3751/adcf29, 2025.

Visual reports: [Mockert Comparison](https://claude.ai/artifact/54xCd9Uaf2QZiD9jRb4yQm), [Li Comparison](https://claude.ai/artifact/3W6FQWUhVoSF8rKKxXR9pk), [Kaspar and Lohmann Comparison](https://claude.ai/artifact/5j1r6dcDGJukME3RfFsLPL) (the last one predates the Biewald correction). Write-ups: `dunkelflaute/reports/`.

---

## 🗺️ What's Next

Directly from the Status table above, roughly in the order the original proposal sequenced them:

1. **Re-run the Li et al. threshold on DART 1950C/2080C**, not just ERA5 — gives the dual-threshold robustness check the proposal calls for (§5.4), using machinery that already exists.
2. **South Korea domain** — extract the Korean TCo1279-DART fields, build the equivalent boundary mask, and run the same validated pipeline. The hard part (reconstruction + validation) transfers directly; what's missing is just running it on a second domain.
3. **Sub-national spatial patterns** — North Sea vs. Baltic vs. inland Germany; the data already carries this resolution, it isn't aggregated out yet.
4. **Synoptic composite analysis** — mean sea-level pressure during Dunkelflaute events vs. all-time mean, control vs. future, for the physical "why" behind the frequency change above.
5. **Uncertainty quantification** — year-block bootstrap and the internal-consistency check described in proposal §5.9; right now the 59% frequency increase is a point estimate with no reported confidence interval.
6. **Threshold/window sensitivity** — CF thresholds of 0.04/0.08/0.10 and 24h/72h windows, to see how much the headline finding depends on the specific Mockert calibration.

---

## 🖥️ Running on the HPC

`python` is not on `PATH` by default. Load the environment module first:

```bash
module load analysis-toolbox/python-04.2026
```

Many scripts have a matching `.sh` wrapper that does this for you — prefer the wrapper where one exists.

Downloads use the [CDS API](https://cds.climate.copernicus.eu/how-to-api) and need a personal access token in `~/.cdsapirc`; no credentials are stored here.

## 💾 Data

Input/output NetCDF and most intermediate `.npz`/`.csv` files are **not tracked** in git (see `.gitignore`) — they run to hundreds of GB. Recreate them by running the relevant `download_*` / `compute_*` / `dart_reconstruct_*` scripts. Small, genuinely load-bearing summary files (report data, training configs) are tracked.

## 🌐 ERA5 Grid Note

ERA5 here is 0.25° regular lat/lon. Raw percentiles over that grid over-count polar points — 48% of land points sit poleward of 60°N/S but only 21% of land *area* does. All spatial averages in this repo are area-weighted by cos(lat) unless a script says otherwise.

## 📜 License

MIT — see [LICENSE](LICENSE).
