# Research Proposal

# Changes in Dunkelflaute Events in Germany and South Korea Under High-Emission Climate Change: A 9 km Resolution Analysis Using AWI-CM3

**An extension study of the AWI wind and solar resource assessment project**

Author: Amanda
Affiliation: Alfred Wegener Institute, Helmholtz Centre for Polar and Marine Research, Bremerhaven
Date: September 2026
Duration: 4 weeks (short-term internship extension)

> **Errata (2026-10-28).** The proposal is kept as written; these points were found wrong when checked against the papers:
> - Mockert et al. (2023) is Mockert, Grams, Brown, Neumann, Meteorol. Appl., 30, e2141, https://doi.org/10.1002/met.2141, not Meteorol. Z. 32, 181–197. Kaspar is not an author.
> - "Lohmann et al. (2025)" is Biewald et al. (2025), Environ. Res.: Energy, 2, 025007. Its 23% optimum uses a 24 h duration, not 48 h.
> - Mockert's 0.06 threshold gives about 4.3 events per year (1979–2018), and their weights are Germany's 2018 mix (44% solar, 50% onshore, 6% offshore). The 2024 weights in Sec. 5.4 are this project's alternative, not Mockert's.
> - The "30 and 54 day" Dunkelflaute periods of winter 2022/23 have no peer-reviewed source. Mockert's longest event in 40 years is 8 days.
>
> Current status and next steps: [`README.md`](README.md) and [`PLAN.md`](PLAN.md).

---

## 1. Background and motivation

### 1.1 The Dunkelflaute challenge

The global energy transition relies increasingly on wind and solar power. Both resources are inherently variable, and grid stability depends on the assumption that when one source is weak, the other often compensates. Dunkelflaute events violate that assumption: they are sustained periods when both wind and solar resources fall below usable levels simultaneously. During such events, electricity systems must draw on storage, backup generation or imports, making Dunkelflaute critical for infrastructure planning and energy security.

**Germany** is where the term *Dunkelflaute* originates and where the challenge is most acute. The country has the highest installed wind and solar capacity in Europe, with electricity demand peaking in winter when solar contribution is lowest. Between December 2022 and March 2023, two extended Dunkelflaute periods lasting 30 and 54 days were observed, underscoring the severity of the problem. Germany's *Bundesnetzagentur* and ENTSO-E have identified Dunkelflaute as a key challenge for system adequacy.

**South Korea** faces an analogous but climatologically distinct challenge. The country's energy transition targets substantial increases in wind and solar capacity under the 10th Basic Plan for Electricity Supply and Demand. South Korea's climate is shaped by the East Asian monsoon, producing seasonal patterns fundamentally different from those in Germany: hot, humid summers with monsoon cloud cover that reduces solar output, and cold winters when high-pressure systems from Siberia can suppress both wind and solar across the peninsula simultaneously. The complex coastline and mountainous interior (~70% of land area) make South Korea an ideal test case for high-resolution climate modelling.

### 1.2 Dunkelflaute definitions in the literature

There is no single standard definition. The two most widely adopted approaches are:

**Combined CF threshold (Mockert et al., 2023; KIT/DWD):** A Dunkelflaute event occurs when the 48-hour running mean of the combined (capacity-weighted) wind and solar capacity factor falls below 0.06. This threshold was calibrated to reproduce the event frequency reported by Kaspar et al. (2019), and the 48-hour window focuses on multi-day events that pose the greatest challenge for energy storage (Schmidt et al., 2019). This identifies approximately 4 severe events per year in Germany. The methodology has been widely adopted in subsequent studies (e.g. Strnad et al., 2025).

**Individual technology CF threshold (Li et al., 2021):** A Dunkelflaute occurs when both wind CF and solar PV CF individually fall below 20%. This identifies approximately 5–10 events per year lasting longer than one day in Germany, and is the most cited definition for climatological studies.

**Evaluation against grid adequacy (Lohmann et al., 2025):** A systematic evaluation against actual grid stress data found that detection performs best near a 23% individual-technology CF threshold with a 48-hour duration criterion, though no simple meteorological threshold perfectly predicts grid stress (peak F-score of 0.14).

### 1.3 Gap in current research

Most climate projections of renewable resources examine wind and solar changes separately. The parent AWI research project uses AWI-CM3 to assess future wind and solar resources, but its wind–solar complementarity assessment is limited to the monthly scale and explicitly notes that the available radiation data do not support weather-scale joint deficit analysis.

Furthermore, existing Dunkelflaute climatologies focus almost exclusively on Northern Europe. No study has applied the same event-based methodology to compare Dunkelflaute risk across fundamentally different climate regimes using high-resolution climate model output.

### 1.4 Relationship to existing work

This study builds on:

**The parent AWI wind and solar resource project**, which provides the scientific framework and methodological conventions.

**Moon et al. (2025)**, "Earth's future climate and its variability simulated at 9 km global resolution" (Earth Syst. Dynam., 16, 1103–1134), which describes the AWI-CM3 high-resolution simulation protocol. The HR data originate from the AWI–IBS (Institute for Basic Science, Busan, South Korea) collaboration.

**An existing Dunkelflaute detection codebase developed internally at AWI**, which implements the methodology of Mockert et al. (2023; KIT/DWD) — a 48-hour moving average combined CF threshold of 0.06 — with Germany's 2024 installed capacity weights and pre-built configurations for both ERA5 reanalysis and TCo1279 model data.

---

## 2. Research questions

### 2.1 Full research framework

The complete set of questions this programme could address:

1. Baseline Dunkelflaute climatology under early-industrial (1950) forcing
2. Climate change signal under SSP5-8.5
3. Temporal evolution through the 21st century
4. Natural variability context from the long MR control
5. Resolution sensitivity (9 km vs 31 km vs 100 km)
6. Scenario dependence (SSP1-2.6, SSP3-7.0, SSP5-8.5)
7. Onshore vs offshore wind behaviour during events
8. Cross-regional comparison: Germany vs South Korea
9. Synoptic mechanisms driving Dunkelflaute events and their change under warming
10. Pipeline validation against ERA5 reanalysis

### 2.2 Scope of this 4-week study

This study focuses on questions 1, 2, 7, 8, 9 and 10 using high-resolution (TCo1279-DART, ~9 km) data:

- **Validate** the analysis pipeline against ERA5 reanalysis for Germany, confirming it produces realistic Dunkelflaute statistics before applying it to climate projections
- **Compare** Dunkelflaute events between the 1950 control and 2080s SSP5-8.5 at 9 km for Germany and South Korea
- **Characterise** sub-national spatial patterns, seasonal cycles and event duration distributions
- **Assess** onshore vs offshore wind behaviour during events
- **Examine** synoptic pressure patterns associated with Dunkelflaute events and whether they change under warming

Medium-resolution analysis, additional scenarios, temporal evolution and the full natural variability assessment from the 184-year MR control are documented as direct future extensions.

---

## 3. Study domains

### 3.1 Germany

The German domain covers approximately 47–55°N, 6–15°E, including the North Sea and Baltic Sea coastal waters where offshore wind is deployed.

Germany offers meaningful internal contrasts at 9 km resolution:

- **North Sea coast and offshore** — Europe's offshore wind hub; strongly Atlantic-influenced; frequent passage of extratropical cyclones
- **Baltic coast** — more sheltered, different wind regime from the North Sea
- **North German Plain** — large onshore wind capacity; flat terrain; continental influence increases eastward
- **Central and southern Germany** — low mountain ranges (Harz, Thuringian Forest, Black Forest), higher solar share, Alpine fringe in Bavaria

These sub-national contrasts allow the study to show where within Germany Dunkelflaute risk is highest and how it changes under warming, rather than producing only a national average.

### 3.2 South Korea

The Korean domain covers approximately 33–39°N, 124–132°E, encompassing the peninsula and surrounding seas (Yellow Sea, East Sea / Sea of Japan, Korea Strait).

South Korea provides sub-national contrasts of its own:

- **West coast (Yellow Sea)** — shallow waters, growing offshore wind development zone; tidal flats
- **South coast** — warmer, higher solar potential; key offshore wind development area
- **East coast** — steep coastal topography facing the East Sea
- **Interior mountains** — Taebaek range covers much of the eastern peninsula; complex terrain that 9 km resolution resolves far better than 31 km

### 3.3 Why Germany and South Korea

Comparing these two countries with the same methodology, model and thresholds provides:

- **Climate regime contrast:** Germany's Dunkelflaute events are driven by Atlantic blocking and NAO-related circulation patterns. South Korea's are driven by the East Asian monsoon and Siberian High. If the climate change response is similar, that supports a robust signal; if it differs, the contrasting drivers help explain why.
- **Resolution test:** Both have complex coastlines and terrain, but shaped by entirely different geography. 9 km resolution should matter in both, but for different features.
- **Policy relevance:** Both countries are actively planning large-scale renewable deployment and need Dunkelflaute risk information under climate change.
- **Data collaboration:** The HR data originate from the AWI–IBS collaboration. Producing results for South Korea alongside Germany directly serves both partner institutions.

---

## 4. Data and simulations

### 4.1 What "1950 control" means

The IPCC defines "pre-industrial" as the 1850–1900 global average. Pre-industrial CO₂ was approximately 280 ppm. By 1950, atmospheric CO₂ had risen to approximately 310 ppm — about 11% above pre-industrial — and global temperature was roughly 0.2–0.3°C above the 1850–1900 baseline.

The 1950 control therefore represents **early industrial conditions**, not a pristine pre-industrial state. However, the vast majority of anthropogenic warming occurs after 1950: CO₂ has since risen from 310 to over 420 ppm, with an additional ~0.9°C of warming by the 2020s, and approximately 4–5°C more by the 2080s under SSP5-8.5 (Moon et al., 2025).

The 1950-to-2080s comparison captures the bulk of the anthropogenic climate change signal. Throughout this study, "control" refers to the 1950-forcing simulation and "future" to the 2080s SSP5-8.5 simulation.

### 4.2 Available data: TCo1279-DART (~9 km)

The data originate from the AWI–IBS collaboration described in Moon et al. (2025). Both simulations use the TCo1279 atmospheric resolution (~9 km) coupled with the DART variable-resolution ocean mesh (4–25 km). Data are global; German and Korean domains are extracted during processing.

**TCo1279-DART-1950C (control)**

| Grid | Frequency | Variables | Years |
|---|---|---|---|
| Native (6,599,680 cells) | 3-hourly | 10u, 10v, 2t, cp, lsp, tcc, hcc, mcc, lcc, tsr, w850 | 1950–1969 |
| 5136 × 2560 (remapped) | 3-hourly | msl | 1950–1969 |
| 5136 × 2560 (remapped) | Monthly | fal, sst, tcwv, ssrd, ssrc, strc, strd, tsrc, ttr, ttrc | 1950–1969 |
| 5136 × 2560 (remapped) | Monthly | tsr | 1966–1969 only |

**TCo1279-DART-2080C (SSP5-8.5 future)**

| Grid | Frequency | Variables | Years |
|---|---|---|---|
| Native (6,599,680 cells) | 3-hourly | 10u, 10v, 2t, cp, lsp, tcc, hcc, mcc, lcc, tsr, w850 | 2080–2092 |
| 5136 × 2560 (remapped) | 3-hourly | msl | 2080–2092 |
| 5136 × 2560 (remapped) | Monthly | fal, sst, tcwv, ssrd, ssrc, str, strc, strd, tsr, tsrc, ttr, ttrc | 2080–2092 |
| 5136 × 2560 (remapped) | Monthly | ssr | 2082–2092 only |

### 4.3 Spin-up exclusion and usable periods

The HR simulations are time-slice experiments branched from the MR transient. When the model resolution jumps from 31 km to 9 km, the atmosphere and land surface need time to adjust: small-scale features such as coastal wind patterns, orographic effects and soil moisture that the 9 km grid can resolve were not present in the 31 km initial conditions. Following Moon et al. (2025), the first 2 years of each simulation are excluded.

| Simulation | Total period | Spin-up exclusion | Usable period | Usable years |
|---|---|---|---|---|
| 1950 Control | 1950–1969 | 1950–1951 | 1952–1969 | 18 years |
| 2080s SSP5-8.5 | 2080–2092 | 2080–2081 | 2082–2092 | 11 years |

### 4.4 ERA5 reanalysis for validation

**ERA5** (Hersbach et al., 2020) provides a global atmospheric reanalysis at approximately 31 km resolution with hourly output from 1940 to present. It includes all variables needed for Dunkelflaute detection (10 m wind components, surface solar radiation, 2 m temperature) at hourly resolution.

The existing AWI codebase includes a pre-built `ERA5Config`. Running the Dunkelflaute pipeline on ERA5 for Germany over a recent period (e.g. 2000–2020) produces a reference Dunkelflaute climatology that can be compared against published results (Li et al., 2021; Mockert et al., 2023; Kaspar et al., 2019). This serves two purposes:

- **Pipeline validation:** If the pipeline reproduces known Dunkelflaute statistics from ERA5 (approximately 4 Mockert-threshold events per year in Germany), it confirms the implementation is correct before applying it to climate model data.
- **Model evaluation:** Comparing the AWI-CM3 1950 control Dunkelflaute climatology against ERA5 helps assess whether the model produces realistic event characteristics, acknowledging that the 1950 control and ERA5 represent different climate states and that ERA5 itself has biases.

### 4.5 Key variables

| Variable | Code | Available at | Used for |
|---|---|---|---|
| 10 m eastward wind | 10u | 3h native | Wind speed → wind CF |
| 10 m northward wind | 10v | 3h native | Wind speed → wind CF |
| 2 m temperature | 2t | 3h native | Solar PV temperature correction |
| Total cloud cover | tcc | 3h native | Surface solar estimation |
| High cloud cover | hcc | 3h native | Surface solar estimation |
| Medium cloud cover | mcc | 3h native | Surface solar estimation |
| Low cloud cover | lcc | 3h native | Surface solar estimation |
| TOA net solar radiation | tsr | 3h native | Surface solar estimation |
| Surface solar rad. downward | ssrd | Monthly remapped | Validation of 3h solar estimate |
| Mean sea-level pressure | msl | 3h remapped | Synoptic composite analysis |

### 4.6 Why SSP5-8.5 only

HR simulations exist exclusively under SSP5-8.5. Running additional scenarios requires approximately 676,000 core hours per 10-year chunk (Moon et al., 2025), far beyond the scope of this short-term internship. SSP5-8.5 provides the strongest climate change signal, maximising detectability within the short analysis windows.

---

## 5. Methods

### 5.1 Wind capacity factor

Wind capacity factors are computed separately for **onshore** and **offshore** grid cells.

**Step 1 — Wind speed:**

v₁₀ = √(u₁₀² + v₁₀²)

**Step 2 — Hub-height extrapolation (power law):**

v_hub = v₁₀ × (H_hub / 10)^α

| Parameter | Onshore | Offshore |
|---|---|---|
| Hub height H_hub | 100 m | 150 m |
| Wind shear exponent α | 1/7 ≈ 0.143 | 0.11 |

**Step 3 — Simplified turbine power curve:**

| Wind speed range | CF_wind(v) |
|---|---|
| v < v_in (cut-in) | 0 |
| v_in ≤ v < v_rated | (v³ − v_in³) / (v_rated³ − v_in³) |
| v_rated ≤ v ≤ v_out | 1.0 |
| v > v_out (cut-out) | 0 |

| Parameter | Onshore | Offshore |
|---|---|---|
| v_in | 3 m/s | 3 m/s |
| v_rated | 12 m/s | 13 m/s |
| v_out | 25 m/s | 25 m/s |

**Step 4 — Apply per sample, then average:**
The power curve is applied to each 3-hourly sample before temporal averaging. CF(v̄) ≠ mean(CF(v)).

The model's land-sea mask classifies grid cells as onshore or offshore.

### 5.2 Solar PV capacity factor

**Primary (simplified) formula:**

CF_PV(t) ≈ G_est(t) / 1000

where G_est(t) is the estimated surface solar radiation and 1000 W/m² is the standard test condition irradiance.

**With temperature correction (sensitivity test):**

CF_PV(t) = [G_est(t) / 1000] × [1 + γ × (T_2m(t) + c_T × G_est(t) − 25°C)]

where γ ≈ −0.005 per °C and c_T ≈ 0.035 °C·m²/W.

### 5.3 Estimating sub-daily surface solar radiation

Surface solar radiation (ssrd) is available only monthly. Sub-daily estimates use 3-hourly cloud cover and TOA radiation.

**Cloud-attenuation model:**

G_est(t) = S_clear(t) × K_clear × (1 − a_low × lcc − a_mid × mcc − a_high × hcc)

where S_clear(t) is the theoretical clear-sky TOA incoming solar (from solar geometry), K_clear ≈ 0.75 is clear-sky atmospheric transmittance, and a_low ≈ 0.7, a_mid ≈ 0.5, a_high ≈ 0.3 are cloud-type attenuation coefficients.

**Validation:** 3-hourly estimates are aggregated to monthly means and compared with the actual monthly ssrd for both simulations and both domains. This dual-domain validation is important because the cloud-attenuation relationship may perform differently under German (Atlantic-dominated) and Korean (monsoon-dominated) cloud regimes.

### 5.4 Dunkelflaute event definition

Two complementary thresholds are applied:

#### Primary: Mockert et al. (2023)

**Combined capacity factor:**

CF_combined(t) = w_solar × CF_PV(t) + w_onshore × CF_wind,onshore(t) + w_offshore × CF_wind,offshore(t)

A Dunkelflaute event occurs when the **48-hour moving average** of CF_combined falls below **0.06**.

**Germany capacity weights (2024):** w_solar = 0.577, w_onshore = 0.369, w_offshore = 0.054 (source: Bundesnetzagentur).

**South Korea capacity weights:** Weights reflecting South Korea's current or planned capacity mix are applied for the Korean domain. If precise weights are unavailable, equal weighting (1/3 each) serves as a neutral alternative.

#### Secondary: Li et al. (2021)

Both wind CF < 20% **and** solar PV CF < 20% simultaneously. This captures more frequent, less severe events (~5–10/year vs ~4/year), providing better statistics for the short analysis windows.

#### Rationale for dual thresholds

Both are applied to all analyses. Agreement strengthens conclusions; divergence reveals threshold-dependent findings. The Li threshold is particularly important given the short future window (11 years × ~4 Mockert events/year ≈ 44 events total, marginal for seasonal breakdown).

#### Sensitivity tests

- Combined CF thresholds: 0.04, 0.08, 0.10
- Moving-average windows: 24 h, 72 h
- Country-specific vs equal capacity weights for South Korea

#### Polar night handling

Germany (47–55°N) experiences very short winter days but not full polar night. Solar CF is naturally very low in December–January but not zero. No special polar-night handling is required for either study domain.

### 5.5 Event metrics

For each grid cell, sub-national region, season (DJF, MAM, JJA, SON) and simulation:

- **Event frequency:** Distinct events per year or season
- **Mean and median duration:** Days
- **Maximum duration:** Longest event in the analysis period
- **90th-percentile duration:** Duration exceeded by 10% of events
- **Dunkelflaute fraction:** Share of time classified as Dunkelflaute (%)
- **Spatial coherence:** Fraction of country area simultaneously below CF threshold (%)

### 5.6 Quantifying the climate change signal

Δ_metric = metric_2080s − metric_1950control

Reported as absolute differences, percentage changes and per degree of global warming.

### 5.7 Synoptic composite analysis

Mean sea-level pressure (msl) is available at 3-hourly resolution for both simulations. For each country:

- Composite msl patterns during Dunkelflaute events are computed and compared with the all-time mean, revealing the characteristic synoptic situations (e.g. blocking anticyclones over Central Europe, Siberian High extensions over Korea).
- Composites are computed separately for the 1950 control and 2080s future, testing whether the synoptic patterns associated with Dunkelflaute change under warming — whether existing patterns intensify, shift geographically, or whether new patterns emerge.

This provides physical interpretation of any detected changes in event frequency or duration.

### 5.8 ERA5 validation protocol

The analysis pipeline is first applied to ERA5 reanalysis for Germany over 2000–2020 using the existing `ERA5Config`. The resulting Dunkelflaute statistics (event frequency, mean duration, seasonal distribution) are compared with published values from Mockert et al. (2023), Li et al. (2021) and Kaspar et al. (2019).

If the pipeline reproduces approximately 4 severe Mockert-threshold events per year with a winter DJF peak and mean durations of 2–5 days, it is considered validated. Significant discrepancies are investigated and documented before proceeding with climate model analysis.

The AWI-CM3 1950 control climatology for Germany is then compared qualitatively with the ERA5 results. Agreement in spatial patterns and seasonal cycles — acknowledging the different climate states — builds confidence in the model's representation of the meteorological conditions that produce Dunkelflaute.

### 5.9 Uncertainty estimation

**Year-block bootstrap:** Whole years are resampled with replacement. Median and 5th–95th percentile ranges are reported. Events are not joined across resampled year boundaries.

**Internal consistency check:** The 18-year control is split into two non-overlapping segments and event metrics are compared, providing a limited estimate of internal variability.

**Minimum sample size:** Grid cells or seasons with fewer than ~20 events are flagged as having insufficient data for robust statistics.

---

## 6. Limitations

**Only two time periods with usable HR data.** Ideally, this study would analyse multiple time slices across the 21st century — including the HR 2000s historical slice and intermediate periods (2030s, 2060s) — to track how Dunkelflaute risk evolves with progressive warming. However, the HR 2000s slice (TCo1279-DART-2000) contains only daily-mean 2 m temperature, with no wind, radiation or cloud cover at sub-daily resolution, making it unusable for Dunkelflaute detection. The 2030s and 2060s HR slices are not present in the accessible data archive. The study is therefore limited to comparing two endpoints: the 1950 control and the 2080s SSP5-8.5. This captures the full magnitude of the climate change signal but cannot determine whether changes emerge gradually, accelerate after a threshold, or appear only under strong warming. Filling this gap would require locating the missing HR slices or running the analysis on the MR transient (2015–2099), which has 6-hourly wind spanning the full century.

**Short analysis windows.** 18 usable years of control and 11 years of future. Under the Mockert threshold (~4 events/year), the future period yields approximately 44 total events — marginal for seasonal or sub-national breakdown. The Li threshold partially mitigates this. Longer simulation periods would substantially improve statistical robustness, particularly for characterising rare long-duration events and for sub-national analysis where events are further split by region.

**Early-industrial baseline, not pre-industrial.** The 1950 control includes CO₂ at ~310 ppm (~11% above the pre-industrial ~280 ppm). The climate change signal slightly underestimates total change since true pre-industrial times. A control simulation at 1850 conditions would provide a cleaner reference but does not exist in the available HR data.

**Estimated sub-daily solar radiation.** Surface solar is reconstructed from TOA radiation and cloud cover. The cloud-attenuation model introduces uncertainty, and its performance may differ between German and Korean cloud regimes. Dual-domain validation against monthly ssrd addresses this partially, but residual 3-hourly errors may affect event detection near the CF threshold. If ssrd were available at 3-hourly resolution, this estimation step would be unnecessary.

**No full natural variability assessment.** Without the 184-year MR control, only a limited internal consistency check is possible. The full MR-based variability assessment is a priority future extension.

**Single emission scenario.** Only SSP5-8.5. Whether Dunkelflaute changes under lower-emission pathways are proportionally smaller — or whether the response is nonlinear — remains open.

**Single model family.** AWI-CM3 only. Different climate models produce different regional circulation responses to warming. Multi-model comparison would strengthen confidence.

**Fixed capacity weights.** Germany's 2024 capacity mix is the primary reference. Actual capacity distributions will change as the energy transition progresses. South Korean weights are approximate.

**Meteorological potential only.** Actual electricity shortfalls depend on demand patterns, storage capacity, grid interconnection, dispatchable backup and curtailment.

---

## 7. Future extensions

These extensions require no new climate simulations:

- **MR analysis and natural variability:** Apply the pipeline to TCo319 (31 km) control and SSP5-8.5; use the 184-year MR control for robust natural variability estimates.
- **Resolution comparison:** Systematically compare 9 km and 31 km Dunkelflaute results for Germany and South Korea.
- **Additional scenarios:** Apply the pipeline to MR SSP1-2.6 and SSP3-7.0.
- **Temporal evolution:** Use MR SSP5-8.5 transient (2015–2099) to track decadal changes.
- **LR comparison:** Apply the pipeline to TCO95L91 (~100 km) for three-scale resolution analysis.
- **Pan-European extension:** Expand from Germany to all of Europe.
- **Global extension:** Apply to other regions (broader East Asia, North America, Australia).
- **Energy-system coupling:** Combine with demand profiles and storage models.
- **Country-specific capacity weights:** Apply national weights for all European countries.
- **Stability-corrected wind extrapolation:** Replace the power-law assumption with Monin-Obukhov similarity theory, using additional fields (friction velocity, surface heat flux) where available.

---

## 8. Expected contribution

This study provides the first comparative analysis of Dunkelflaute events under climate change at 9 km atmospheric resolution for Germany and South Korea — two countries with high renewable ambitions in fundamentally different climate regimes. It delivers:

- **Validated pipeline:** Confirmation against ERA5 reanalysis that the detection methodology reproduces known Dunkelflaute statistics for Germany, providing confidence before applying it to climate projections.
- **Baseline Dunkelflaute climatologies** at 9 km for both countries, with sub-national spatial detail (coastal vs inland vs mountainous) that coarser assessments cannot provide.
- **Climate change projections** showing how event frequency, duration and seasonality change under approximately 5°C of global warming.
- **Cross-regional comparison** revealing whether the Dunkelflaute response is consistent between Atlantic-influenced and monsoon-influenced climates, or whether different circulation drivers produce fundamentally different outcomes.
- **Physical interpretation** through synoptic composite analysis, connecting event changes to shifts in large-scale pressure patterns.
- **Dual-threshold robustness:** Results reported under both the Mockert et al. (2023) severe-event threshold and the Li et al. (2021) moderate-event threshold.
- **Onshore vs offshore characterisation** of wind behaviour during events in both countries.

---

## References

Hersbach, H., Bell, B., Berrisford, P., et al.: The ERA5 global reanalysis, Q. J. R. Meteorol. Soc., 146, 1999–2049, https://doi.org/10.1002/qj.3803, 2020.

Kaspar, F., Borsche, M., Pfeifroth, U., Trentmann, J., Drücke, J., and Becker, P.: A climatological assessment of balancing effects and shortfall risks of photovoltaics and wind energy in Germany and Europe, Adv. Sci. Res., 16, 119–128, https://doi.org/10.5194/asr-16-119-2019, 2019.

Li, B., Basu, S., Watson, S. J., and Russchenberg, H. W. J.: A Brief Climatology of Dunkelflaute Events over and Surrounding the North and Baltic Sea Areas, Energies, 14, 6508, https://doi.org/10.3390/en14206508, 2021.

Lohmann, J., et al.: Evaluation of 'Dunkelflaute' event detection methods considering grid operators' needs, Environ. Res.: Energy, https://doi.org/10.1088/2753-3751/adcf29, 2025.

Mockert, F., Grams, C. M., Brown, T., and Kaspar, F.: Meteorological conditions during periods of low wind speed and insolation in Germany, Meteorol. Z., 32, 181–197, https://doi.org/10.1127/metz/2023/2141, 2023.

Moon, J.-Y., Streffing, J., Lee, S.-S., et al.: Earth's future climate and its variability simulated at 9 km global resolution, Earth Syst. Dynam., 16, 1103–1134, https://doi.org/10.5194/esd-16-1103-2025, 2025.

Schmidt, O., Melchior, S., Hawkes, A., and Staffell, I.: Projecting the future levelized cost of electricity storage technologies, Joule, 3, 81–100, https://doi.org/10.1016/j.joule.2018.12.008, 2019.

---

## Appendix A: Work plan (4 weeks)

### Week 1 — Validation, pipeline and control-run analysis

**Days 1–2: ERA5 validation**
Run the existing Dunkelflaute pipeline on ERA5 for Germany using the pre-built `ERA5Config`. Compute event statistics for 2000–2020. Compare against published Dunkelflaute climatologies (Mockert et al., Li et al., Kaspar et al.). Confirm pipeline produces realistic results. Document any calibration needed.

**Days 3–4: Data access and pipeline adaptation for AWI-CM3**
Extract German and Korean domains from the TCo1279-DART native-grid output. Verify variables. Obtain land-sea mask. Adapt pipeline for TCo1279 data:
- Compute onshore and offshore wind CF from 3-hourly 10u/10v
- Implement cloud-attenuation model for 3-hourly surface solar estimation
- Validate estimated solar against monthly ssrd for both domains
- Apply both Mockert and Li threshold definitions

Test on small subsets.

**Days 5–6: 1950 control Dunkelflaute climatology**
Run pipeline on 1950 control for both Germany and South Korea. Produce maps of seasonal event frequency, mean duration and Dunkelflaute fraction. Compute sub-national statistics. Compare German model results with ERA5 reference.

**Day 7: Documentation**
Write methods section. Document solar estimation validation and ERA5 comparison.

### Week 2 — Climate change signal and cross-regional comparison

**Days 8–10: 2080s future analysis**
Apply identical thresholds to the 2080s SSP5-8.5 simulation. Compute all metrics. Produce change maps (2080s minus control) for both countries.

**Day 11: Seasonal analysis**
Monthly histograms for both countries. Does Germany's winter peak intensify? Does Korea's winter Siberian-High-driven peak change? Any emerging monsoon-season risk?

**Day 12: Cross-regional comparison**
Side-by-side Germany vs Korea: baseline, changes, seasonal patterns. Identify commonalities and differences.

**Day 13: Onshore vs offshore**
Compare wind CF during events. Does offshore wind provide residual generation? Compare German North Sea / Baltic with Korean Yellow Sea / south coast.

**Day 14: Threshold sensitivity**
Repeat with alternative thresholds. Test Korean capacity weights. Identify robust findings.

### Week 3 — Sub-national analysis, synoptic patterns and robustness

**Days 15–16: Sub-national spatial analysis**
Germany: North Sea coast vs Baltic coast vs North German Plain vs south. Korea: west coast vs east coast vs interior mountains. Map spatial patterns within each country.

**Days 17–18: Synoptic composite analysis**
Composite msl during Dunkelflaute events vs all-time mean for each country. Compare composites between control and future. Are blocking patterns changing? New patterns emerging?

**Day 19: Spatial coherence**
What fraction of each country simultaneously in Dunkelflaute? Does coherence change under warming?

**Day 20: Uncertainty and robustness**
Bootstrap for key results. Internal consistency check. Summarise robust findings.

**Day 21: Buffer**
Catch up. Begin figures.

### Week 4 — Synthesis, figures and write-up

**Days 22–23: Figures**
1. ERA5 validation results for Germany
2. Control Dunkelflaute climatology — Germany (9 km)
3. Control Dunkelflaute climatology — South Korea (9 km)
4. Climate change signal — Germany
5. Climate change signal — South Korea
6. Cross-regional comparison (seasonal cycles, metrics)
7. Sub-national spatial patterns
8. Synoptic composites (control vs future)
9. Onshore vs offshore wind CF during events
10. Threshold sensitivity summary

**Days 24–25: Writing**

**Day 26: Review and cross-check**

**Day 27: Revision**

**Day 28: Final delivery**

**Deliverables:**
- Complete written report with publication-quality figures
- Slide summary for colleague presentation
- Reproducibility archive (Python code, environment, processed data)