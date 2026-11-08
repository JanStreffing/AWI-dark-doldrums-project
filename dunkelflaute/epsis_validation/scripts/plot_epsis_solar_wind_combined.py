"""EPSIS (South Korea): solar and wind installed and added capacity in one stacked view, with absolute and percent shares.
Installed = latest month of each year (December, except the partial latest year). Added = change in that value."""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = "/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/epsis_validation"
INK, SEC, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7", "#fcfcfb"
C_SOLAR, C_WIND = "#eda100", "#2a78d6"

d = pd.read_csv(f"{HERE}/data/generation_capacity/HOME>GenerationCapacity>ByFuel.csv", skiprows=2, header=None)
d["year"] = d[0].str[:4].astype(int)
last = d.sort_values(0).groupby("year").tail(1).set_index("year")
inst = pd.DataFrame({"solar": last[10], "wind": last[11]}) / 1000  # GW
added = inst.diff().dropna()
latest = last[0].iloc[-1]
partial = last[0].str[5:] != "12"


def pct(df):
    return df.div(df.sum(axis=1), axis=0) * 100


def draw(ax, df, title, ylab, percent=False):
    x = df.index
    ax.bar(x, df.solar, width=0.7, color=C_SOLAR, label="Solar", edgecolor=SURF, linewidth=1)
    ax.bar(x, df.wind, width=0.7, bottom=df.solar, color=C_WIND, label="Wind", edgecolor=SURF, linewidth=1)
    ax.set_facecolor(SURF)
    ax.set_title(title, loc="left", color=INK, fontsize=12, fontweight="bold")
    ax.set_ylabel(ylab, color=SEC)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{y}*" if partial[y] else f"{y}" for y in x], rotation=60, color=SEC)
    ax.tick_params(colors=SEC, length=0)
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)
    if percent:
        ax.set_ylim(0, 100)
        for y in x:  # label the wind share on every bar, the solar share on the bar below it
            ax.text(y, df.solar[y] + df.wind[y] / 2, f"{df.wind[y]:.0f}", ha="center", va="center", color="white", fontsize=7.5, fontweight="bold")
            ax.text(y, df.solar[y] / 2, f"{df.solar[y]:.0f}", ha="center", va="center", color=INK, fontsize=7.5)
    else:
        tot = df.sum(axis=1)
        ax.annotate(f"{tot.iloc[-1]:.1f}", (x[-1], tot.iloc[-1]), ha="center", va="bottom", color=INK, fontsize=9,
                    xytext=(0, 2), textcoords="offset points")


fig, axes = plt.subplots(2, 2, figsize=(13, 9), facecolor=SURF)
draw(axes[0, 0], inst, "Installed capacity", "GW")
axes[0, 0].legend(frameon=False, loc="upper left", labelcolor=SEC)
draw(axes[0, 1], added, "Capacity added per year", "GW")
draw(axes[1, 0], pct(inst), "Installed capacity: solar vs. wind share", "% of solar + wind capacity", percent=True)
draw(axes[1, 1], pct(added), "Capacity added: solar vs. wind share", "% of solar + wind capacity added", percent=True)
fig.suptitle("South Korea (EPSIS): solar and wind capacity combined", x=0.01, ha="left", color=INK, fontsize=14)
fig.text(0.01, 0.005, f"* latest month available ({latest}), not year-end; its addition covers Dec {last.index[-2]} to {latest} only. "
         "Percentages are shares of solar + wind only (other fuels not included). Source: EPSIS GenerationCapacity ByFuel.",
         color=SEC, fontsize=8)
fig.tight_layout(rect=(0, 0.02, 1, 0.96))
out = f"{HERE}/figures/epsis_solar_wind_combined.png"
fig.savefig(out, dpi=150)
print(out)
p = pct(inst); print(p.round(1).iloc[[0, -2, -1]]); print(pct(added).round(1).tail(3))
