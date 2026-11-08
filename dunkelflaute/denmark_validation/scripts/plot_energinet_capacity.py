"""Denmark (Energinet): installed offshore wind, onshore wind and solar capacity, monthly, and capacity added per year."""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = "/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/denmark_validation"
INK, SEC, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7", "#fcfcfb"
COL = {"OffshoreWindCapacity": "#2a78d6", "OnshoreWindCapacity": "#1baf7a", "SolarPowerCapacity": "#eda100"}
LAB = {"OffshoreWindCapacity": "Offshore wind", "OnshoreWindCapacity": "Onshore wind", "SolarPowerCapacity": "Solar"}

c = pd.read_csv(f"{HERE}/data/energinet_capacity_denmark_monthly.csv", parse_dates=["Month"]).set_index("Month") / 1000  # GW
cols = list(COL)
year_end = c.groupby(c.index.year).last()
added = year_end.diff()
added.iloc[0] = year_end.iloc[0] - c.iloc[0]  # 2016: Jan to Dec only
partial_year = c.index[-1].year if c.index[-1].month != 12 else None


def style(ax, title, ylab):
    ax.set_facecolor(SURF)
    ax.set_title(title, loc="left", color=INK, fontsize=12, fontweight="bold")
    ax.set_ylabel(ylab, color=SEC)
    ax.tick_params(colors=SEC, length=0)
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)


fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.8), facecolor=SURF, gridspec_kw={"width_ratios": [1.2, 1]})
a1.stackplot(c.index, [c[k] for k in cols], colors=[COL[k] for k in cols], labels=[LAB[k] for k in cols], edgecolor=SURF, linewidth=0.6)
style(a1, "Installed capacity (monthly)", "GW")
a1.legend(frameon=False, loc="upper left", labelcolor=SEC)
a1.set_xlim(c.index[0], c.index[-1])
tot = c[cols].sum(axis=1)
a1.annotate(f"{tot.iloc[-1]:.1f} GW", (c.index[-1], tot.iloc[-1]), ha="right", va="bottom", color=INK, fontsize=9, xytext=(0, 3), textcoords="offset points")

w = 0.27
x = added.index.values
for i, k in enumerate(cols):
    a2.bar(x + (i - 1) * w, added[k], width=w, color=COL[k], label=LAB[k])
style(a2, "Capacity added per year", "GW")
a2.axhline(0, color=AXIS, lw=0.8)
a2.set_xticks(x)
a2.set_xticklabels([f"{y}*" if y == partial_year else ("2016**" if y == 2016 else f"{y}") for y in x], rotation=45, color=SEC)
fig.suptitle("Denmark (Energinet): installed wind and solar capacity", x=0.01, ha="left", color=INK, fontsize=13)
fig.text(0.01, 0.01, f"* Dec {partial_year - 1} to {c.index[-1]:%b %Y} only.  ** Jan to Dec 2016 only (series starts Jan 2016).  "
         "Added = change in year-end capacity; negative = decommissioned. Source: Energinet CapacityPerMunicipality, summed over municipalities.",
         color=SEC, fontsize=8)
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
out = f"{HERE}/figures/energinet_capacity_wind_solar.png"
fig.savefig(out, dpi=150)
print(out); print(added.round(2).to_string())
