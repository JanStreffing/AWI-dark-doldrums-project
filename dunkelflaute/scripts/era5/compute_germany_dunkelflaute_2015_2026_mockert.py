"""
Germany combined capacity factor + Mockert et al. (2023) Dunkelflaute events,
ERA5 2015-2026, Sep-Mar extended winter, real land+EEZ boundary.

Sources (all already verified to use the real Germany boundary, not a bbox):
  - solar: dunkelflaute/smard_validation/data/germany_solar_cf_v5_sepmar_2015_2026.npz
           (model_v5, trained on genuine 3-hour-sum ERA5 ssrd/tsr -- NOT the
           earlier 1h-accum-sampled-every-3h data model_v3 used)
  - wind:  dunkelflaute/data/germany_era5/wind_cf_era5_fullyear_2015_2026.npz
           (full year; restricted to Sep-Mar here to match solar's coverage)
  - SMARD: dunkelflaute/smard_validation/data/smard/realisierte_erzeugung/*, resampled to
           the same 3-hourly grid (label="right", closed="right", matching
           ERA5's own backward-looking accumulation convention)

Combined CF uses Mockert et al. (2023)'s own weights (44% solar, 50% onshore,
6% offshore -- Germany's 2018 mix, IRENA 2019), not this project's own 2024
default mix -- so this is a literal, weight-matched comparison to Mockert's
published numbers.

Event construction follows Mockert Section 2.4 exactly: 48-hour running mean
(16 samples at 3-hourly resolution) of combined CF, threshold 0.06, and --
critically -- "all time steps contributing to a running mean below the
threshold are considered part of the Dunkelflaute" (window-expansion, not a
naive threshold-crossing run), which guarantees the minimum event duration is
48h by construction, matching the paper's own stated property.

SMARD offshore CF here uses a capacity denominator LINEARLY INTERPOLATED
between each year's Jan-1 installed-capacity anchor, not a flat per-year
value -- 2015's offshore capacity roughly tripled within the year itself
(993 MW -> 3283 MW anchors), and a flat denominator both overstates CF for
most of 2015 (up to 297%) and, worse for event detection, leaves large NaN
gaps if those points are simply excluded instead.
"""
import os

import numpy as np
import pandas as pd

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
os.chdir(_ROOT)

SOLAR_PATH = "dunkelflaute/smard_validation/data/germany_solar_cf_v5_sepmar_2015_2026.npz"
WIND_PATH = "dunkelflaute/data/germany_era5/wind_cf_era5_fullyear_2015_2026.npz"
SMARD_GEN_DIR = "dunkelflaute/smard_validation/data/smard/realisierte_erzeugung"
SMARD_CAP_DIR = "dunkelflaute/smard_validation/data/smard/installierte_erzeugungsleistung_history"
OUT_PATH = "dunkelflaute/data/germany_era5/combined_cf_2015_2026_sepmar_mockertweights.npz"

WEIGHTS = dict(solar=0.44, onshore=0.50, offshore=0.06)
WIN = 16          # 16 * 3h = 48h
THRESH = 0.06
MONTHS_SEPMAR = {1, 2, 3, 9, 10, 11, 12}


def winter_label(ts):
    return ts.year + 1 if ts.month >= 9 else ts.year


def load_solar():
    d = np.load(SOLAR_PATH, allow_pickle=True)
    df = pd.DataFrame({"cf_solar_recon": d["cf_recon"], "cf_solar_real": d["cf_real"]},
                       index=pd.to_datetime(d["time"], utc=True))
    return df


def load_wind():
    d = np.load(WIND_PATH, allow_pickle=True)
    df = pd.DataFrame({
        "cf_on_recon": d["cf_on_recon"], "cf_off_recon": d["cf_off_recon"],
        "cf_on_real": d["cf_on_real"], "cf_off_real": d["cf_off_real"],
    }, index=pd.to_datetime(d["time"], utc=True))
    return df[df.index.month.isin(MONTHS_SEPMAR)]


def resample_3h(df):
    return df.set_index("timestamp_utc")["value_mw"].resample("3h", label="right", closed="right").mean()


def interpolated_capacity(cap_year_csv, index):
    """Linear interpolation between each year's Jan-1 anchor, not a flat per-year value."""
    c = pd.read_csv(cap_year_csv).sort_values("year")
    anchors = pd.Series(
        c["value_mw"].values,
        index=pd.to_datetime(c["year"].astype(str) + "-01-01", utc=True),
    )
    last_year = c["year"].max()
    anchors[pd.Timestamp(f"{last_year+1}-01-01", tz="UTC")] = c["value_mw"].iloc[-1]
    full_index = anchors.index.union(index).sort_values()
    return anchors.reindex(full_index).interpolate(method="time").reindex(index)


def load_smard():
    solar_gen = pd.read_csv(f"{SMARD_GEN_DIR}/smard_solar_hour_2015_2026.csv", parse_dates=["timestamp_utc"])
    on_gen = pd.read_csv(f"{SMARD_GEN_DIR}/smard_wind_onshore_hour_2015_2026.csv", parse_dates=["timestamp_utc"])
    off_gen = pd.read_csv(f"{SMARD_GEN_DIR}/smard_wind_offshore_hour_2015_2026.csv", parse_dates=["timestamp_utc"])

    solar_3h = resample_3h(solar_gen)
    on_3h = resample_3h(on_gen)
    off_3h = resample_3h(off_gen)

    solar_cap = interpolated_capacity(f"{SMARD_CAP_DIR}/smard_installed_solar_year_history.csv", solar_3h.index)
    on_cap = interpolated_capacity(f"{SMARD_CAP_DIR}/smard_installed_wind_onshore_year_history.csv", on_3h.index)
    off_cap = interpolated_capacity(f"{SMARD_CAP_DIR}/smard_installed_wind_offshore_year_history.csv", off_3h.index)

    df = pd.DataFrame({
        "cf_solar_smard": (solar_3h / solar_cap).clip(upper=1.02),
        "cf_onshore_smard": (on_3h / on_cap).clip(upper=1.02),
        "cf_offshore_smard": (off_3h / off_cap).clip(upper=1.02),
    })
    df.index = df.index.tz_convert("UTC")
    return df[df.index.month.isin(MONTHS_SEPMAR)]


def rolling_mean_per_winter(series, winter_labels, winters):
    out = pd.Series(index=series.index, dtype=float)
    for w in winters:
        mask = winter_labels == w
        out[mask] = series[mask].rolling(WIN, min_periods=WIN).mean()
    return out


def mockert_events(roll, winter_labels, winters, dt_h=3.0):
    """Mockert Section 2.4: expand every below-threshold rolling-mean sample to
    its full contributing window, merge overlapping/adjacent windows -> events.
    Guarantees minimum event duration 48h by construction."""
    events = []
    for w in winters:
        sub = roll[winter_labels == w]
        n = len(sub)
        below = (sub < THRESH).fillna(False).values
        flagged = np.zeros(n, dtype=bool)
        for i in np.where(below)[0]:
            lo = max(0, i - WIN + 1)
            flagged[lo:i + 1] = True
        start = None
        for i in range(n):
            if flagged[i] and start is None:
                start = i
            elif not flagged[i] and start is not None:
                events.append((w, (i - start) * dt_h))
                start = None
        if start is not None:
            events.append((w, (n - start) * dt_h))
    return events


def main():
    solar = load_solar()
    wind = load_wind()
    smard = load_smard()

    df = solar.join(wind, how="inner").join(smard, how="inner")
    df = df.sort_index()
    print(f"merged: {len(df)} rows, {df.index.min()} .. {df.index.max()}")
    print(f"NaN counts:\n{df.isna().sum()}")

    df["winter"] = [winter_label(ts) for ts in df.index]
    winters = sorted(w for w in df["winter"].unique())
    counts = df.groupby("winter").size()
    complete_winters = [w for w in winters if counts.get(w, 0) >= 1600]  # full Sep-Mar block
    print(f"complete Sep-Mar winters: {complete_winters}")
    df = df[df["winter"].isin(complete_winters)].copy()

    for src in ["recon", "real"]:
        df[f"cf_combined_{src}"] = (WEIGHTS["solar"] * df[f"cf_solar_{src}"] +
                                     WEIGHTS["onshore"] * df[f"cf_on_{src}"] +
                                     WEIGHTS["offshore"] * df[f"cf_off_{src}"])
    df["cf_combined_smard"] = (WEIGHTS["solar"] * df["cf_solar_smard"] +
                                WEIGHTS["onshore"] * df["cf_onshore_smard"] +
                                WEIGHTS["offshore"] * df["cf_offshore_smard"])

    results = {}
    per_winter = {}
    roll_cols = {}
    for src in ["recon", "real", "smard"]:
        roll = rolling_mean_per_winter(df[f"cf_combined_{src}"], df["winter"], complete_winters)
        roll_cols[src] = roll
        evs = mockert_events(roll, df["winter"], complete_winters)
        durs_days = np.array([d / 24.0 for _, d in evs])
        n = len(evs)
        pw = {}
        for w, d in evs:
            pw[w] = pw.get(w, 0) + 1
        results[src] = dict(
            n_events=n, events_per_winter=round(n / len(complete_winters), 3),
            median_dur=round(float(np.median(durs_days)), 2) if n else 0.0,
            mean_dur=round(float(np.mean(durs_days)), 2) if n else 0.0,
            max_dur=round(float(np.max(durs_days)), 2) if n else 0.0,
        )
        per_winter[src] = {w: pw.get(w, 0) for w in complete_winters}
        print(f"{src}: {results[src]}")

    np.savez(
        OUT_PATH,
        time=df.index.values, winter=df["winter"].values,
        cf_solar_recon=df["cf_solar_recon"].values, cf_solar_real=df["cf_solar_real"].values,
        cf_solar_smard=df["cf_solar_smard"].values,
        cf_on_recon=df["cf_on_recon"].values, cf_off_recon=df["cf_off_recon"].values,
        cf_on_real=df["cf_on_real"].values, cf_off_real=df["cf_off_real"].values,
        cf_onshore_smard=df["cf_onshore_smard"].values, cf_offshore_smard=df["cf_offshore_smard"].values,
        cf_combined_recon=df["cf_combined_recon"].values, cf_combined_real=df["cf_combined_real"].values,
        cf_combined_smard=df["cf_combined_smard"].values,
        roll_recon=roll_cols["recon"].values, roll_real=roll_cols["real"].values, roll_smard=roll_cols["smard"].values,
        weights_solar=WEIGHTS["solar"], weights_onshore=WEIGHTS["onshore"], weights_offshore=WEIGHTS["offshore"],
        winters=np.array(complete_winters),
        results=results, per_winter=per_winter,
    )
    print(f"\nSaved: {OUT_PATH}")


if __name__ == "__main__":
    main()
