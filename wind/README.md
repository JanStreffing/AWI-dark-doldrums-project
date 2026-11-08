# 💨 Wind Capacity Factor Reconstruction

> Reconstructs 100 m wind (`u100`/`v100`) and capacity factor for TCo1279-DART — which only carries 10 m wind — validated against real ERA5 100 m wind and real German grid data (SMARD), for use in Dunkelflaute (dark-doldrum) risk analysis.

---

## 📌 Project Purpose

TCo1279-DART carries `u10`/`v10` only, and no usable roughness field to extrapolate from — but driving a wind capacity-factor calculation needs wind at hub height (100 m). This project closes that gap with a physics-based extrapolation, and answers the questions that matter before the result can be trusted:

- Which extrapolation method — log law, power law, or an ERA5-fitted regression — best predicts real 100 m wind from 10 m wind?
- Does DART carry a usable roughness-length field of its own, and if not, what should stand in for it?
- Does the chosen method hold up against **real German grid generation** (SMARD) — an independent check that doesn't depend on ERA5 being right?
- What's still unresolved, and what should be checked before applying this to a future-climate scenario (2080C) with no ground truth at all?

---

## 📂 Data Sources

| Source | Period | Role |
|--------|--------|------|
| **ERA5** reanalysis | `training_data/` spans 1979–2026; method comparison and SMARD-aligned validation use full-year 2015–2026 | Validation testbed (real `u100`/`v100`) + production land-z0 source (`fsr` climatology, 2016–2025, year-cycled) |
| **TCo1279-DART 1950C** | 1950–1969 (20 years) | Present-day-analogue scenario the method is applied to — reconstruction complete, 140/140 months |
| **TCo1279-DART 2080C** | 2080–2092 (13 years) | High-emission future scenario the method is applied to — reconstruction complete, 91/91 months |
| **SMARD** (Bundesnetzagentur) | 2015–2026 | Independent ground truth, unrelated to ERA5 |
| **Boundaries** ([`boundaries/`](../boundaries/)) | Static — Natural Earth land polygon + Marine Regions EEZ polygon | Real Germany land+EEZ mask; not read directly here — applied one step downstream in `dunkelflaute/` to turn the global reconstruction into the Germany-only aggregate SMARD is compared against |

### ERA5 fields

| Field | Description | Role |
|-------|-------------|------|
| `u10`/`v10` | 10 m wind components | Input to every method |
| `u100`/`v100` | 100 m wind components, **real** | Validation target (held out, never fitted on for log law/power law) |
| `fsr` | Forecast surface roughness | Production land-z0 source for DART (DART's own is unusable — see Key Findings) |
| `lsm` | Land-sea mask | Onshore/offshore split, same convention used throughout |

### TCo1279-DART fields

| Field | Description | Role |
|-------|-------------|------|
| `u10`/`v10` | 10 m wind components | Only wind field DART carries — the gap this project fills |
| `sr` (param 173) | Nominal roughness length | Present but unusable — a flat 0.0001 m placeholder in every init file checked |
| `sdor`/`sdfor` | Orography (terrain relief) | Not a substitute for roughness — ~29× apart in magnitude, and exactly zero over ocean |
| `lsm` | Land-sea mask | Native grid, no regridding needed; 28.2% land |

### SMARD fields

| Field | Description | Unit | Timestep |
|-------|-------------|------|----------|
| Onshore/offshore wind generation | Realisierte Erzeugung | MW | Hourly (native); resampled to 3-hourly to match ERA5/DART |
| Installed capacity | Onshore/offshore wind | MW | Yearly |

### Derived quantity

| Field | Description |
|-------|-------------|
| `z0` | Roughness length — DART's own is unusable; land borrows ERA5's `fsr` climatology, ocean is solved from the Charnock relation using that month's own DART wind |

---

## 🛠️ Technical Stack

| Layer | Tool |
|-------|------|
| **Method** | Pure physics formulas (log law, power law) + one ERA5-fitted linear regression kept only for comparison, never used in production |
| **Data handling** | xarray, netCDF4, numpy |
| **Compute** | AWI HPC (SLURM), `analysis-toolbox/python-04.2026` module |
| **Persistence** | NetCDF4 (reconstructed `u100`/`v100`), json/npz (comparison statistics) |

---

## 🏗️ Pipeline

```
ERA5 (real u10/v10/u100/v100, fsr, lsm)
        │
        ▼
Phase 1: Download + Preprocess
  - Pull ERA5 winds/z0/fsr/lsm (scripts/download/)
  - QC: orography vs. roughness, z0 distribution, time coverage (scripts/preprocessing/)
        │
        ▼
Phase 2: Validate three methods against real ERA5 u100/v100
  - Log law, power law, ERA5-fitted regression (scripts/methods/{loglaw,powerlaw,comparison}/)
  - Multiple timestamps, full year 2015–2025, and 10/20-winter pooled runs,
    area-weighted by cos(lat); regression fit on 80% of 10°×10° spatial
    blocks, all three scored on the held-out 20%
        │
        ▼
Phase 3: Decide the z0 source for DART
  - DART's own roughness field is unusable (see Key Findings)
  - Settled on ERA5's fsr climatology, after a sensitivity check showed the
    25 km→9 km coarsening cost is small (scripts/preprocessing/z0_scale_sensitivity.py)
        │
        ▼
Phase 4: Production reconstruction (log law, the only method actually run on DART)
  - scripts/methods/dart_loglaw/reconstruct_month.py, chunked streaming
  - Land z0: ERA5 fsr, year-cycled onto DART's longer run
  - Ocean z0: Charnock relation, solved iteratively from that month's own DART u10/v10
  - Finished: 1950C 140/140 months, 2080C 91/91 months
        │
        ▼
Phase 5: Validate the resulting wind capacity factor against SMARD real generation
  - dunkelflaute/smard_validation/, independent of the ERA5-only checks above
        │
        ▼
Phase 6: Literature validation (NOT in wind/ — the gate before DART is trusted)
  - Wind CF alone isn't what gets checked against the literature — it's
    weighted-combined with solar CF first, in dunkelflaute/, on the same
    2015-2026 ERA5 period, THEN checked against four published Dunkelflaute
    definitions (Mockert, Li, Kaspar, Biewald). See solar/README.md's Pipeline
    section for the linked visual reports — the same combined pipeline covers
    both wind and solar.
```

**Visual reports — Phase 2, method comparison**: log law vs. power law vs. regression, scored against real ERA5 `u100`/`v100`, onshore vs. offshore, overall and in the low-wind band that defines a dark doldrum.

- [All-Season Windshear](https://claude.ai/artifact/G4tqL2B4wBcygXRx5S9nhk) — full calendar year, 2015–2025 (132 months), plus a check of how all three methods hold up in the regions solar's model had to exclude from training (south of 60°S, above 3000 m).

---

## 🗺️ Method Roadmap

| Method | Needs fitted params? | Needs a roughness field? | Status |
|--------|----------------------|---------------------------|--------|
| **Power law** (IEC/Moemken et al. 2018 exponents) | No | No | Rejected — no clear accuracy advantage to justify dropping the log law's physical grounding |
| **Regression** (ERA5-fitted `u10→u100`, per region) | Yes | No | Rejected — better *overall* RMSE onshore, but loses in the low-wind dark-doldrum regime that matters here, and its fitted coefficients are untested on a 9 km grid quite unlike the 25 km one they were fit on |
| **Log law** | No | Yes (`z0`) | **Production** — `v100 = v10 · ln(100/z0)/ln(10/z0)`. DART's own `z0` is unusable, so land `z0` borrows ERA5's `fsr` climatology and ocean `z0` is solved from the Charnock relation using DART's own wind — no coefficients fitted on one grid and transferred to another |

All three methods were validated against real ERA5 `u100`/`v100` before this choice was made (Phase 2 above); only the log law was ever run on DART (Phase 4).

---

## 📐 Method Specification

**Production method (log law)**: `v100 = v10 · ln(100/z0)/ln(10/z0)`, applied independently to `u` and `v` components.

| Parameter | Value |
|-----------|-------|
| Land `z0` | ERA5 `fsr` monthly climatology, 2016–2025, year-cycled onto the DART run (`tco1279_land_z0_2016_2025.npz`) |
| Ocean `z0` | Charnock relation, solved iteratively per month from that month's own DART `u10`/`v10`: `κ, g, ν, α_Charnock = 0.4, 9.81, 1.5e-5, 0.018` |
| Fitted coefficients | None |
| Grid | TCo1279-DART native (6,599,680 cells), no regridding |

**Comparison-only regression fit** (full-year 2015–2025, ERA5-only, never applied to DART):

| Region | Slope | Intercept |
|--------|-------|-----------|
| Onshore | 1.387 | 0.590 |
| Offshore | 1.269 | −0.713 |

**Capacity factor** (`dunkelflaute/scripts/core/capacity_factor.py`): simplified cubic turbine power curve, same 100 m hub height onshore and offshore, only the cut-in/rated/cut-out speeds differ:

| | `v_in` | `v_rated` | `v_out` |
|---|---|---|---|
| Onshore | 3.0 m/s | 12.0 m/s | 25.0 m/s |
| Offshore | 3.0 m/s | 13.0 m/s | 25.0 m/s |

---

## 🔑 Key Findings

### 1. Log law vs. real ERA5 u100/v100

RMSE ≈ 1.19 m/s onshore / 0.60 m/s offshore, r = 0.97 / 0.997. Stable across seasons, no drift. Offshore is ~2× better than onshore; onshore error concentrates over mountains and forest.

### 2. Log law vs. regression, low-wind (<3 m/s) regime — the dark-doldrum regime this whole project cares about

The two methods disagree by metric: regression wins on RMSE (its fitted intercept helps at low wind), log law wins on MAPE (28–31% vs 35–39% onshore, 15% vs 18–24% offshore) because that same intercept is a large *fraction* of a 1–2 m/s true value. **Log law chosen for production** on that basis, plus it needs no fitted coefficients to transfer to DART's very different grid. Neither method is precise in this regime (MAPE 15–39%) — stated as a limitation, not hidden.

### 3. TCo1279's own roughness field is unusable

`sr` (param 173) is a flat placeholder (0.0001 m everywhere) in every DART init file checked, not a real field — confirmed by map, with `sdor` as a control that renders correctly on the same grid, so the flatness is the data, not a reading error. Orography (`sdor`/`sdfor`) isn't a substitute either: it's terrain relief (land median 8.7 m), not aerodynamic roughness (ERA5 land median 0.30 m) — ~29× apart in magnitude — and is exactly zero over ocean, where the log law works best.

### 4. Borrowing ERA5's z0 is defensible

The log-law ratio `ln(100/z0)/ln(10/z0)` is only logarithmically sensitive to z0: a 10× z0 error costs ~19% in wind at the smooth end but ~76% at the rough end. Coarsening ERA5 z0 by 4× (25→111 km) costs only ~0.29–0.32 m/s RMSE at 100 m; the actual mismatch (25→9 km, 2.8×) is estimated at ~0.2–0.25 m/s — about 3% degradation on top of the existing 1.19 m/s onshore error, small next to the method's own baseline error.

### 5. Two caveats bigger than the z0 transfer itself, not yet addressed

- **Fetch**: an internal boundary layer grows ~1:100, so a 100 m wind needs ~10 km of upwind fetch to equilibrate. ERA5 at 25 km exceeds that fetch (local log law defensible); DART's 9 km cells are *at or below* it, so the local log law is on weaker footing there regardless of which z0 is supplied.
- **Displacement height**: ignored entirely. `z0 ≈ 0.1h`, `d ≈ 0.7h` for canopy height `h`; for a 20 m forest, `d ≈ 14 m`, which puts the 10 m level *below* `d`, where `ln((10−d)/z0)` is undefined. Likely a first-order cause of onshore error, concentrated in the same tall-canopy regions where performance is worst.

### 6. Against SMARD real generation (pooled, full year 2015–2026, Germany)

| | RMSE | bias | r | mean CF |
|---|---|---|---|---|
| onshore, reconstruction | 0.098 | +0.005 | 0.911 | 0.205 |
| onshore, real ERA5 | 0.072 | −0.044 | 0.948 | 0.205 |
| offshore, reconstruction | 0.274 | −0.033 | 0.704 | 0.433 |
| offshore, real ERA5 | 0.270 | −0.008 | 0.714 | 0.433 |

Onshore reconstruction tracks real ERA5 closely (r 0.911 vs 0.948). Offshore pools at r ~0.70–0.71 for *both* reconstruction and real ERA5, but that figure is dominated by SMARD's offshore capacity data: offshore capacity roughly tripled within 2015 (993 to 3283 MW), and excluding the 2015–17 capacity-step years gives r = 0.855 (reconstruction) and 0.871 (real ERA5), RMSE 0.158 / 0.162 (`dunkelflaute/smard_validation/data/wind_cf_comparison_fullyear_2015_2026.json`, flag set in `compare_wind_cf_fullyear.py`).

**Visual reports — Phase 5, SMARD validation**: reconstructed and real-ERA5 wind capacity factor against SMARD's actual metered generation, Germany, real land/EEZ boundary.

- [Wind CF Against SMARD](https://claude.ai/artifact/9pqGyEX1KsWmPuMfQudDGv) — the full head-to-head this table is drawn from: Jan 2015–Aug 2026, all months, pooled and by year, with an error decomposition (recon-vs-real, real-vs-SMARD, recon-vs-SMARD) that separates the log law's own error from what a wind-speed-only power curve misses about the real grid.
- [Wind CF, Full Year](https://claude.ai/artifact/JmNk29ZKBqXyWhhdtb7Ezb) — the same Jan 2015–Aug 2026 period as a single continuous daily-mean timeline (not just pooled/by-year stats), plus the installed-capacity chart that's the denominator behind every SMARD CF value and a walkthrough of the capacity-factor formula itself.
- [Wind CF, All Winters](https://claude.ai/artifact/2fvVfhVGFcKR1PquRWw2jw) — the same comparison restricted to Sep–Mar, the season this project actually uses, season by season across all 11 winters 2015-16–2025-26.

### 7. DART reconstruction status

Complete for both 1950C (140/140 months, 1950–1969) and 2080C (91/91 months, 2080–2092) — `data/dart_1950c_100m/`, `data/dart_2080c_100m/`.

---

## 💡 Recommendations Before Trusting DART-1950C / 2080C Further

| Scenario | Status | Recommendation |
|----------|--------|-----------------|
| **1950C** | Reconstruction finished (140/140 months) | Closer to the ERA5 validation climate — lower distribution-shift risk, but still an out-of-sample test with no direct ground truth. |
| **2080C** | Reconstruction finished (91/91 months) | No ground truth exists for a future climate. Frame any result as *"under the assumption that today's log-law/roughness physics holds"*, not as validated. |
| **Both** | Not yet run | Fetch and displacement-height corrections remain open (Key Finding 5) — both are plausible first-order sources of the onshore error and neither has been quantified or corrected yet. |
| **Both** | Done | Germany combined CF and Mockert events for 1950C and 2080C (real boundary, Mockert weights) are in `dunkelflaute/results/` — see `report_dunkelflaute.html`, built from this method's output combined with the solar side. |

---

## 📁 Project Structure

```
wind/
├── scripts/
│   ├── download/       ERA5 winds/z0/fsr/lsm downloads, single timestamps or full winters
│   ├── preprocessing/  QC: orography vs. roughness, z0 distribution, time coverage
│   ├── methods/
│   │   ├── loglaw/       log-law reconstruction vs. real ERA5 (validation)
│   │   ├── powerlaw/     power-law alternative (validation)
│   │   ├── dart_loglaw/  production DART reconstruction (reconstruct_month.py; sbatch launchers in jobs/)
│   │   └── comparison/   method-vs-method scoring, including the only regression fit in the repo
│   └── analysis/        exploratory extensions beyond the three methods (e.g. testing
│                         solar radiation/cloud cover as a log-law stability-error proxy)
├── data/
│   ├── era5/             raw ERA5 downloads (single timestamps + training_data/ winters)
│   ├── dart_extract/     TCo1279-DART single-field extracts
│   ├── dart_1950c_100m/, dart_2080c_100m/   production reconstruction output
│   └── results/          derived json/npz from the method-comparison scripts
├── figures/               era5/ (raw-field diagnostics), methods/{loglaw,powerlaw,comparison}/
├── jobs/                  SLURM sbatch scripts (downloads + DART array reconstruction)
├── logs/                  run logs
└── reports/               FINDINGS.md and other write-ups
```

Every script `chdir`s to this `wind/` root on import, so they run from anywhere; SLURM jobs in `jobs/` invoke them by their `scripts/<category>/` path.

---

## ⚙️ Technical Notes

### Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| Log law over regression for production | No fitted coefficients needed to transfer to DART's very different 9 km grid; wins MAPE in the low-wind dark-doldrum regime that matters most, even though regression wins RMSE |
| Borrow ERA5 `fsr` for land z0, year-cycled | DART carries no usable roughness field at all (flat placeholder); ERA5's own climatology (2016–2025) is cycled onto DART's much longer run |
| Charnock relation for ocean z0, solved per month from DART's own wind | No fitted parameters, physically grounded, computed natively from each month's own DART winds rather than borrowed from ERA5 |
| Chunked streaming in `reconstruct_month.py` | A full-month load was OOM-killed at ~50 GB resident; chunking (default 40 steps, ~5 days) keeps peak memory to a few GB |
| Same 100 m hub height onshore and offshore (no 150 m offshore extrapolation) | Reuses the existing 100 m log-law reconstruction directly instead of adding a second, unvalidated extrapolation step |

### Known Data Quality Issues

- **TCo1279's own roughness field is unusable** (Key Finding 3) — not flagged anywhere in the DART init files themselves; discovered by direct inspection.
- **Fetch and displacement height are both unaddressed** (Key Finding 5) — stated as open limitations, not corrected for.
- **Neither method is precise in the low-wind regime** (MAPE 15–39% at <3 m/s, Key Finding 2) — the regime this whole project's dark-doldrum analysis actually cares about most.
- **Offshore SMARD comparison is noisier for both reconstruction and real ERA5** (Key Finding 6, r ~0.70–0.71 pooled) — mostly the 2015–17 offshore capacity-step artifact in SMARD; without those years r = 0.855 / 0.871 for reconstruction / real ERA5.
- **The 2080C output files are mislabelled.** Every file in `data/dart_2080c_100m/` has the global attribute `description = "...TCo1279-DART-1950C"`, from a hardcoded string in `reconstruct_month.py` (fixed 2026-10-28 for future runs). The data are 2080C; `title` and `land_z0_source` in the same files carry the correct 2080s years.
- **Untested cross-resolution transfer, beyond the z0 coarsening check already run**: the fetch and displacement-height caveats (Key Finding 5) are both resolution-dependent effects that the z0 sensitivity test (Key Finding 4) doesn't cover.

---

## 📊 Key Formulas

```
Log law (production):   v100 = v10 · ln(100/z0) / ln(10/z0)
Power law (comparison):  v100 = v10 · (100/10)^α            α = 0.20 land, 0.14 sea (IEC / Moemken et al. 2018)
Regression (comparison): v100 = slope · v10 + intercept     (ERA5-fitted per region, never applied to DART)

Ocean z0 (Charnock, solved iteratively):
  z0 = α_Charnock · u*² / g + 0.11 · ν / u*
  u* = κ · v10 / ln(10 / z0)
  κ, g, ν, α_Charnock = 0.4, 9.81, 1.5e-5, 0.018

Capacity factor (simplified cubic turbine power curve):
  CF(v) = 0                                            v < v_in
  CF(v) = (v³ − v_in³) / (v_rated³ − v_in³)             v_in ≤ v < v_rated
  CF(v) = 1                                             v_rated ≤ v ≤ v_out
  (v_in, v_rated, v_out) = (3.0, 12.0, 25.0) m/s onshore, (3.0, 13.0, 25.0) m/s offshore
```

---

## 📜 License

MIT — see the repository [LICENSE](../LICENSE).
