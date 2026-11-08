# Germany ERA5 combined-CF files — which one to use

**Use `combined_cf_2015_2026_sepmar_mockertweights.npz`.** It's the current,
validated version: real Germany land+EEZ boundary (`boundaries/germany/germany_land.geojson` +
`boundaries/germany/germany_eez.geojson`, not a bounding box), solar reconstructed with
`model_v5` (trained on genuine 3-hour-sum ERA5 ssrd/tsr), Sep–Mar extended
winter 2015–2026, Mockert et al. (2023)'s own capacity weights (44% solar /
50% onshore / 6% offshore), and their exact event-construction rule (a
below-threshold 48h rolling-mean sample pulls its whole contributing window
into the event, not just that instant). Built by
`dunkelflaute/scripts/era5/compute_germany_dunkelflaute_2015_2026_mockert.py`.
Has `cf_combined_{recon,real,smard}`, the `roll_*` 48h rolling means, and the
`results`/`per_winter` event-frequency summary already computed.

`wind_cf_era5_fullyear_2015_2026.npz` — a supporting input to the file above:
onshore/offshore wind CF only (no solar), real boundary, full calendar year
(not Sep-Mar restricted). Built by `compute_germany_wind_cf_fullyear.py`.

`era5_report_data.json` — pre-formatted data for one of the published HTML
comparison reports (see `dunkelflaute/reports/`); not a general-purpose CF
dataset.

## `archive/` — superseded, kept for reference only

Three earlier iterations of the combined-CF computation, all predating the
corrected event-construction rule and the Sep-Mar window used above, and all
generated the same afternoon (2026-09-23) while the methodology was still
being worked out:

- `combined_cf_era5_2007_2025.npz`, `combined_cf_era5_2007_2025_mockertweights.npz`,
  `combined_cf_era5_2007_2025_realboundary.npz` — three distinct runs (verified
  different content via checksum, not copies of each other) from an earlier
  version of `compute_germany_combined_cf_era5.py`/`_mockertweights.py`, most
  likely differing in boundary definition (bbox vs real) and/or capacity
  weights per their filenames, but not independently re-verified file-by-file
  against the current script's exact settings — treat the filenames as a
  rough guide, not a precise spec.

None of the three implement Mockert's actual window-expansion event rule
(Section 2.4) or the Sep-Mar window. Don't build new analysis on these —
they're kept so the evolution of the methodology stays traceable, not as a
source of truth. Use `combined_cf_2015_2026_sepmar_mockertweights.npz` above
instead.
