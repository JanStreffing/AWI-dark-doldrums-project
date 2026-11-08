"""EPSIS (South Korea): installed capacity of solar and wind vs. generation output of New & Renewable energy.
Output is only available for the combined New & Renewable category, so the capacity-factor panel uses the capacity of
the whole category (fuel cell, IGCC, solar, wind, small hydro, marine, bio, waste). Output assumed to be in MWh."""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = "/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/epsis_validation"
INK, SEC, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7", "#fcfcfb"
C_SOLAR, C_WIND, C_OTHER, C_OUT = "#eda100", "#2a78d6", "#a8a69c", "#1baf7a"

cap = pd.read_csv(f"{HERE}/data/generation_capacity/HOME>GenerationCapacity>ByFuel.csv", skiprows=2, header=None)
cap["year"] = cap[0].str[:4].astype(int)
cap = cap.sort_values(0).groupby("year").tail(1).set_index("year")  # year-end (latest month of the year)
c = pd.DataFrame({"solar": cap[10], "wind": cap[11], "nre": cap[[8, 9, 10, 11, 12, 13, 14, 15]].sum(axis=1)}) / 1000  # GW
c["other"] = c.nre - c.solar - c.wind

out = pd.read_csv(f"{HERE}/data/generation_output/HOME>GenerationOutputandRetailSales>GenerationOutput>bySource.csv",
                  skiprows=2, header=None)
out = pd.Series((out[16] / 1e6).values, index=out[0].values).sort_index()  # TWh

years = [y for y in c.index if y in out.index]  # overlap 2012-2025 (drops the partial 2026)
c, out = c.loc[years], out.loc[years]
avg_cap = (c.nre + c.nre.shift(1)) / 2                       # GW, mean of start- and end-of-year capacity
cf = out * 1000 / (avg_cap * 8760) * 100                     # TWh -> GWh; % of capacity running all year

def style(ax, title, ylab):
    ax.set_facecolor(SURF)
    ax.set_title(title, loc="left", color=INK, fontsize=12, fontweight="bold")
    ax.set_ylabel(ylab, color=SEC)
    ax.set_xticks(years); ax.set_xticklabels(years, rotation=60, color=SEC)
    ax.tick_params(colors=SEC, length=0)
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)

fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), facecolor=SURF)
ax = axes[0]
ax.bar(years, c.solar, width=0.7, color=C_SOLAR, label="Solar", edgecolor=SURF, linewidth=1)
ax.bar(years, c.wind, width=0.7, bottom=c.solar, color=C_WIND, label="Wind", edgecolor=SURF, linewidth=1)
ax.bar(years, c.other, width=0.7, bottom=c.solar + c.wind, color=C_OTHER, label="Other New & Renewable", edgecolor=SURF, linewidth=1)
style(ax, "Installed capacity (year-end)", "GW")
ax.legend(frameon=False, loc="upper left", labelcolor=SEC, fontsize=9)
ax = axes[1]
ax.bar(years, out, width=0.7, color=C_OUT)
style(ax, "New & Renewable generation output", "TWh")
ax = axes[2]
ax.bar(years[1:], cf.iloc[1:], width=0.7, color=C_OUT)
style(ax, "Implied capacity factor, New & Renewable", "% of capacity running all year")
ax.set_xticks(years[1:]); ax.set_xticklabels(years[1:], rotation=60, color=SEC)
for a, v, f in [(axes[0], c.nre, "{:.1f}"), (axes[1], out, "{:.1f}"), (axes[2], cf.iloc[1:], "{:.1f}%")]:
    a.annotate(f.format(v.iloc[-1]), (v.index[-1], v.iloc[-1]), ha="center", va="bottom", color=INK, fontsize=9,
               xytext=(0, 2), textcoords="offset points")
fig.suptitle("South Korea (EPSIS): installed capacity vs. generation output, 2012-2025", x=0.01, ha="left", color=INK, fontsize=13)
fig.text(0.01, 0.01, "Output exists only for the combined New & Renewable category (no solar/wind split), so the capacity factor "
         "uses the capacity of the whole category. Capacity factor = output / (mean of start and end-of-year capacity x 8760 h).",
         color=SEC, fontsize=8)
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
path = f"{HERE}/figures/epsis_capacity_vs_output.png"
fig.savefig(path, dpi=150)
print(path)
print(pd.DataFrame({"cap_solar": c.solar, "cap_wind": c.wind, "cap_nre": c.nre, "out_TWh": out, "CF_%": cf}).round(2).to_string())
