"""Germany solar: reconstructed (v4/v5) vs real ERA5-direct ssrd, and both vs
real-world SMARD-observed solar CF, Sep-Mar 2015-2026."""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
SMARD_PATH = f"{ROOT}/dunkelflaute/smard_validation/data/cf_smard_solar_fullyear_2015_2026.csv"
INK, SEC, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
COL = {"v4": "#1baf7a", "v5": "#eda100", "real": "#2a78d6", "smard": "#eb6834"}


def load_variant(variant):
    path = f"{ROOT}/dunkelflaute/smard_validation/data/germany_solar_cf_{variant}_sepmar_2015_2026.npz"
    if not os.path.exists(path):
        return None
    d = np.load(path, allow_pickle=True)
    t = pd.to_datetime([str(x) for x in d["time"]]).tz_localize("UTC")
    return pd.DataFrame({f"cf_recon_{variant}": d["cf_recon"], "cf_real": d["cf_real"]}, index=t)


def style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=9, length=3)
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)


def main():
    variants = [v for v in ("v4", "v5") if load_variant(v) is not None]
    frames = [load_variant(v) for v in variants]
    df = frames[0]
    for f in frames[1:]:
        df = df.join(f.drop(columns=["cf_real"]), how="outer")
    smard = pd.read_csv(SMARD_PATH, parse_dates=["timestamp_utc"]).set_index("timestamp_utc")
    df = df.join(smard[["cf_solar_smard"]], how="inner")

    fig = plt.figure(figsize=(16, 11), facecolor=SURF)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.1, 1.1, 1], hspace=0.42, wspace=0.28)

    # Panel A: representative single winter time series (2024-09 .. 2025-03)
    ax = fig.add_subplot(gs[0, :])
    win = df.loc["2024-09-01":"2025-03-31"]
    ax.plot(win.index, win["cf_real"], color=COL["real"], lw=1.1, label="ERA5-direct (true 3h, real ssrd)")
    for v in variants:
        ax.plot(win.index, win[f"cf_recon_{v}"], color=COL[v], lw=1.0, alpha=0.85, label=f"{v} reconstruction")
    ax.plot(win.index, win["cf_solar_smard"], color=COL["smard"], lw=1.1, label="SMARD (real observed generation)")
    ax.set_ylabel("Capacity factor", fontsize=10, color=SEC)
    ax.set_title("Germany solar CF, winter 2024-2025 (Sep-Mar)", fontsize=12, fontweight="bold", color=INK, loc="left")
    ax.legend(frameon=False, fontsize=8.5, labelcolor=SEC, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    style(ax)

    # Panel B: recon vs real ssrd scatter (pooled, all years) -- the "compare recon and real ssrd" ask
    for i, v in enumerate(variants):
        ax = fig.add_subplot(gs[1, i])
        x_all, y_all = df["cf_real"].values * 1000, df[f"cf_recon_{v}"].values * 1000
        day = x_all > 1.0   # drop the night-zero cluster, which otherwise dominates the plot
        x, y = x_all[day], y_all[day]
        ax.hexbin(x, y, gridsize=45, cmap="Greens" if v == "v4" else "Oranges", mincnt=1, bins="log")
        lim = max(x.max(), y.max())
        ax.plot([0, lim], [0, lim], color=AXIS, lw=1, ls=(0, (3, 2)))
        err = y - x
        rmse = np.sqrt((err ** 2).mean()); bias = err.mean()
        r2 = 1 - (err ** 2).sum() / ((x - x.mean()) ** 2).sum()
        ax.text(0.03, 0.95, f"daylight only:\nR²={r2:.3f}\nRMSE={rmse:.1f} W/m²\nbias={bias:+.2f} W/m²",
                transform=ax.transAxes, va="top", fontsize=9, color=SEC)
        ax.set_xlabel("real ssrd (ERA5-direct, W/m²)", fontsize=9.5, color=SEC)
        ax.set_ylabel(f"{v} reconstructed ssrd (W/m²)", fontsize=9.5, color=SEC)
        ax.set_title(f"Germany domain-mean: {v} recon vs real ssrd", fontsize=10.5, color=INK, loc="left")
        style(ax)

    # Panel C: bar summary -- mean CF, model vs SMARD
    ax = fig.add_subplot(gs[2, 0])
    labels = ["SMARD\n(real observed)", "ERA5-direct\n(true 3h)"] + [f"{v} recon" for v in variants]
    means = [df["cf_solar_smard"].mean(), df["cf_real"].mean()] + [df[f"cf_recon_{v}"].mean() for v in variants]
    cols = [COL["smard"], COL["real"]] + [COL[v] for v in variants]
    bars = ax.bar(labels, means, color=cols, width=0.6)
    for b, m in zip(bars, means):
        ax.text(b.get_x() + b.get_width() / 2, m + 0.001, f"{m:.3f}", ha="center", fontsize=9, color=SEC)
    ax.set_ylabel("Mean CF, Sep-Mar 2015-2026", fontsize=9.5, color=SEC)
    style(ax)

    # Panel D: correlation with SMARD (the real-world check)
    ax = fig.add_subplot(gs[2, 1])
    corrs = [np.corrcoef(df["cf_real"], df["cf_solar_smard"])[0, 1]] + \
            [np.corrcoef(df[f"cf_recon_{v}"], df["cf_solar_smard"])[0, 1] for v in variants]
    clabels = ["ERA5-direct\nvs SMARD"] + [f"{v} recon\nvs SMARD" for v in variants]
    ccols = [COL["real"]] + [COL[v] for v in variants]
    bars = ax.bar(clabels, corrs, color=ccols, width=0.6)
    for b, c in zip(bars, corrs):
        ax.text(b.get_x() + b.get_width() / 2, c + 0.005, f"{c:.3f}", ha="center", fontsize=9, color=SEC)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Correlation with real-world SMARD CF", fontsize=9.5, color=SEC)
    style(ax)

    fig.text(0.02, 0.975, "Germany solar: model reconstruction vs real ERA5 vs real-world SMARD generation",
              fontsize=16, fontweight="bold", color=INK)
    fig.text(0.02, 0.955, "Sep-Mar, 2015-2026 (2026 Sep-Dec not yet available); domain-mean over Germany's land area",
              fontsize=10, color=SEC)
    fig.tight_layout(rect=(0.02, 0.02, 1, 0.94))
    out = f"{ROOT}/dunkelflaute/figures/germany_smard_compare.png"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=140, facecolor=SURF)
    plt.close(fig)
    print("saved", out)
    print("\nvariants included:", variants)
    print(f"mean CF -- SMARD: {df['cf_solar_smard'].mean():.4f}  ERA5-direct: {df['cf_real'].mean():.4f}")
    for v in variants:
        print(f"  {v} recon: {df[f'cf_recon_{v}'].mean():.4f}  corr-with-SMARD: {np.corrcoef(df[f'cf_recon_{v}'], df['cf_solar_smard'])[0,1]:.4f}")


if __name__ == "__main__":
    main()
