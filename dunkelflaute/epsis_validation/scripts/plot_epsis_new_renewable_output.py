"""EPSIS (South Korea) generation output of New & Renewable energy (all renewables other than large hydro, one combined
column in the source; solar and wind are not separated). Assumes the table is in MWh."""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = "/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/epsis_validation"
INK, SEC, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7", "#fcfcfb"
GREEN = "#1baf7a"

d = pd.read_csv(f"{HERE}/data/generation_output/HOME>GenerationOutputandRetailSales>GenerationOutput>bySource.csv",
                skiprows=2, header=None)
# column 0 = year, 16 = New & Renewable, 25 = total generation (utilities + non-utility in common use)
d = pd.DataFrame({"year": d[0], "nre": d[16] / 1e6, "total": d[25] / 1e6})  # MWh -> TWh
d = d[d.nre > 0].sort_values("year").set_index("year")
d["share"] = d.nre / d.total * 100

fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), facecolor=SURF)
for ax, col, title, ylab, fmt in [
    (axes[0], "nre", "New & Renewable generation output", "Generation (TWh)", "{:.1f}"),
    (axes[1], "share", "Share of total generation", "Share of total generation (%)", "{:.1f}"),
]:
    v = d[col]
    ax.bar(v.index, v.values, width=0.7, color=GREEN)
    ax.set_facecolor(SURF)
    ax.set_title(title, loc="left", color=INK, fontsize=12, fontweight="bold")
    ax.set_ylabel(ylab, color=SEC)
    ax.set_xticks(v.index); ax.set_xticklabels(v.index, rotation=60, color=SEC)
    ax.tick_params(colors=SEC, length=0)
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)
    ax.annotate(fmt.format(v.iloc[-1]), (v.index[-1], v.iloc[-1]), ha="center", va="bottom", color=INK, fontsize=9,
                xytext=(0, 2), textcoords="offset points")
fig.suptitle("South Korea (EPSIS): New & Renewable energy generation output", x=0.01, ha="left", color=INK, fontsize=13)
fig.text(0.01, 0.01, "Solar, wind, bio, waste, fuel cell etc. combined (EPSIS does not split them in this table). "
         "Share = of total generation incl. non-utility in common use. Source: EPSIS Generation Output by Source.",
         color=SEC, fontsize=8)
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
out = f"{HERE}/figures/epsis_new_renewable_output.png"
fig.savefig(out, dpi=150)
print(out); print(d.round(2).tail(6))
