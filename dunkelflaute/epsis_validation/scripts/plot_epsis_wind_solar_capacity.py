"""EPSIS (South Korea) installed capacity: solar and wind, year-end values, as separate bar charts."""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = "/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/epsis_validation"
INK, SEC, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7", "#fcfcfb"
COLS = {"Solar Power": "#eda100", "Wind Power": "#2a78d6"}

d = pd.read_csv(f"{HERE}/data/generation_capacity/HOME>GenerationCapacity>ByFuel.csv", skiprows=2, header=None)
d = pd.DataFrame({"period": d[0], "Solar Power": d[10], "Wind Power": d[11]})
d["year"] = d.period.str[:4].astype(int)
d["month"] = d.period.str[5:].astype(int)
# last available month of each year (December for complete years, latest month for the current one)
last = d.sort_values("period").groupby("year").tail(1).set_index("year")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), facecolor=SURF)
for ax, name in zip(axes, COLS):
    v = last[name] / 1000  # MW -> GW
    ax.bar(v.index, v.values, width=0.7, color=COLS[name])
    ax.set_facecolor(SURF)
    ax.set_title(f"{name.split()[0]} capacity", loc="left", color=INK, fontsize=12, fontweight="bold")
    ax.set_ylabel("Installed capacity (GW)", color=SEC)
    ax.set_xticks(v.index)
    ax.set_xticklabels([f"{y}" if last.month[y] == 12 else f"{y}*" for y in v.index], rotation=45, color=SEC)
    ax.tick_params(colors=SEC, length=0)
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)
    ax.annotate(f"{v.iloc[-1]:.1f}", (v.index[-1], v.iloc[-1]), ha="center", va="bottom", color=INK, fontsize=9, xytext=(0, 2), textcoords="offset points")
fig.suptitle("South Korea (EPSIS): installed solar and wind capacity", x=0.01, ha="left", color=INK, fontsize=13)
fig.text(0.01, 0.01, f"* latest month available ({last.period.iloc[-1]}), not year-end. Source: EPSIS GenerationCapacity ByFuel.", color=SEC, fontsize=8)
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
out = f"{HERE}/figures/epsis_wind_solar_capacity.png"
fig.savefig(out, dpi=150)
print(out); print(last[["period", "Solar Power", "Wind Power"]])

# ---- added capacity per year (difference of consecutive year-end values) ----
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), facecolor=SURF)
for ax, name in zip(axes, COLS):
    v = (last[name].diff().dropna()) / 1000  # GW added during each year
    ax.bar(v.index, v.values, width=0.7, color=COLS[name])
    ax.set_facecolor(SURF)
    ax.set_title(f"{name.split()[0]} capacity added per year", loc="left", color=INK, fontsize=12, fontweight="bold")
    ax.set_ylabel("Added capacity (GW)", color=SEC)
    ax.set_xticks(v.index)
    ax.set_xticklabels([f"{y}" if last.month[y] == 12 else f"{y}*" for y in v.index], rotation=45, color=SEC)
    ax.tick_params(colors=SEC, length=0)
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)
    ax.annotate(f"{v.iloc[-1]:.2f}", (v.index[-1], v.iloc[-1]), ha="center", va="bottom", color=INK, fontsize=9, xytext=(0, 2), textcoords="offset points")
fig.suptitle("South Korea (EPSIS): solar and wind capacity added per year", x=0.01, ha="left", color=INK, fontsize=13)
fig.text(0.01, 0.01, f"* Dec 2025 to {last.period.iloc[-1]} only (partial year). Added = change in year-end installed capacity. Source: EPSIS GenerationCapacity ByFuel.", color=SEC, fontsize=8)
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
out = f"{HERE}/figures/epsis_wind_solar_capacity_added.png"
fig.savefig(out, dpi=150)
print(out)
