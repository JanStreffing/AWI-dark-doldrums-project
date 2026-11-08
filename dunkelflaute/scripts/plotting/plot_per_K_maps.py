"""
Per-K maps of the change in winter low-output events (5 deg land boxes, 20-70 N), from
results/maps/grid_per_K.npz (scripts/regions/grid_per_K.py).

Figure 1 (per_K_maps_wind.png): wind only, %/K.
  a  CMIP6 multi-model median (daily footing); dots: >= 80 % of models increase, x: >= 80 % decrease
  b-e  9 km, 31 km, 100 km (AWI-CM3 TCo1279/319/95, SSP5-8.5) and EERIE IFS-FESOM2-SR (ctl -> SSP2-4.5),
     6-hourly footing; dots: ratio above the per-box 95 % chance range of the TCO95 control
  f  CMIP6 median, wind plus solar (0.44 solar + 0.56 wind)
Figure 2 (per_K_maps_tail.png): CMIP6 median and 9 km, all events against mean removed (the tail part).

Usage: python3 plot_per_K_maps.py      (module python3/2025.01: cartopy)
"""
import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
import numpy as np                                  # noqa: E402
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap   # noqa: E402

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
IN = f"{_ROOT}/dunkelflaute/results/maps/grid_per_K.npz"
OUT = f"{_ROOT}/dunkelflaute/results/maps"
LEVELS = np.array([-15, -12, -9, -6, -3, -1, 1, 3, 6, 9, 12, 15])   # neutral grey bin -1..1
# two hues and a neutral grey midpoint
CMAP = LinearSegmentedColormap.from_list("div", ["#1f4e8c", "#6f9fd0", "#f0efec", "#e58a6e", "#a3221b"])
NORM = BoundaryNorm(LEVELS, CMAP.N, extend="both")
LATS = np.arange(20, 70.1, 5)
LONS = np.arange(-180, 180.1, 5)
INK = "#2b2b2b"


def to_grid(z, values):
    g = np.full((len(LATS) - 1, len(LONS) - 1), np.nan)
    i = ((z["lat0"] - 20) / 5).round().astype(int)
    j = ((z["lon0"] + 180) / 5).round().astype(int)
    g[i, j] = values
    return g


def panel(ax, z, values, title, dots=None, crosses=None):
    ax.set_extent([-180, 180, 20, 70], ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="#e6e6e3", zorder=0)
    m = ax.pcolormesh(LONS, LATS, to_grid(z, values), cmap=CMAP, norm=NORM, transform=ccrs.PlateCarree(), zorder=1)
    ax.coastlines(linewidth=0.4, color="#555555", zorder=3)
    for mask, mk in ((dots, "."), (crosses, "x")):
        if mask is not None and np.any(mask):
            ax.scatter(z["lon0"][mask] + 2.5, z["lat0"][mask] + 2.5, s=3 if mk == "." else 6, marker=mk,
                       color=INK, linewidths=0.4, transform=ccrs.PlateCarree(), zorder=4)
    ax.set_title(title, loc="left", fontsize=8.5, color=INK)
    return m


def area_mean(z, v):
    w = np.cos(np.deg2rad(z["lat0"] + 2.5))
    ok = np.isfinite(v)
    return (v[ok] * w[ok]).sum() / w[ok].sum()


def main():
    z = np.load(IN)
    k = set(z.files)
    rows = [("cmip6|wind|perK_median", "a  CMIP6 median (20 models), wind", "cmip6|wind|agree_up")]
    for name, lab in (("9km", "b  9 km (TCo1279)"), ("31km", "c  31 km (TCo319)"), ("100km", "d  100 km (TCO95)"),
                      ("eerie_ctl", "e  EERIE 9 km, SSP2-4.5")):
        rows.append((f"6h|{name}|wind|perK", lab, f"6h|{name}|wind"))
    rows.append(("cmip6|comb|perK_median", "f  CMIP6 median, wind + solar", "cmip6|comb|agree_up"))
    rows = [r for r in rows if r[0] in k]
    fig, axes = plt.subplots(len(rows), 1, figsize=(7.2, 1.45 * len(rows) + 0.6),
                             subplot_kw=dict(projection=ccrs.PlateCarree()))
    for ax, (key, lab, sig) in zip(np.atleast_1d(axes), rows):
        v = z[key]
        if sig.endswith("agree_up"):
            a = z[sig]
            dots, crosses = a >= 0.8, a <= 0.2
        else:
            dots = (np.log(z[f"{sig}|ratio"]) > z[f"{sig}|null_p95"]) if f"{sig}|null_p95" in k else None
            crosses = (np.log(z[f"{sig}|ratio"]) < z[f"{sig}|null_p05"]) if f"{sig}|null_p05" in k else None
        m = panel(ax, z, v, f"{lab}   area mean {area_mean(z, v):+.1f} %/K", dots, crosses)
    cb = fig.colorbar(m, ax=axes, orientation="horizontal", fraction=0.035, pad=0.03, aspect=45, ticks=LEVELS)
    cb.set_label("change in low-output events per K of global warming (%/K)", fontsize=8, color=INK)
    cb.ax.tick_params(labelsize=7)
    fn = f"{OUT}/per_K_maps_wind.png"
    fig.savefig(fn, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {fn}")

    pairs = [("cmip6|wind|perK_median", "cmip6|wind|perK_mr_median", "CMIP6 median"), ("6h|9km|wind|perK", "6h|9km|wind|perK_mr", "9 km")]
    pairs = [p for p in pairs if p[0] in k]
    fig, axes = plt.subplots(len(pairs), 2, figsize=(11, 1.05 * len(pairs) + 0.7), subplot_kw=dict(projection=ccrs.PlateCarree()),
                             gridspec_kw=dict(hspace=0.35, wspace=0.04))
    axes = np.atleast_2d(axes)
    for r, (a, b, lab) in enumerate(pairs):
        panel(axes[r, 0], z, z[a], f"{'abcd'[2 * r]}  {lab}: all events   {area_mean(z, z[a]):+.1f} %/K")
        m = panel(axes[r, 1], z, z[b], f"{'abcd'[2 * r + 1]}  {lab}: mean change removed   {area_mean(z, z[b]):+.1f} %/K")
    cb = fig.colorbar(m, ax=axes, orientation="horizontal", fraction=0.04, pad=0.04, aspect=50, ticks=LEVELS)
    cb.set_label("change in low-wind events per K (%/K)", fontsize=8, color=INK)
    cb.ax.tick_params(labelsize=7)
    fn = f"{OUT}/per_K_maps_tail.png"
    fig.savefig(fn, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {fn}")
    plot_scale()


SCALES = ["5deg", "10deg", "20deg", "40deg", "named"]
SCALE_LAB = ["5°", "10°", "20°", "40°", "named\nregions"]
# fixed categorical order (reference palette slots 1-5); CMIP6 in ink with its 10-90 % band
RUNS = [("6h|9km", "9 km", "#2a78d6"), ("6h|31km", "31 km", "#eb6834"), ("6h|100km", "100 km", "#1baf7a"),
        ("6h|eerie_ctl", "EERIE ctl", "#eda100"), ("6h|eerie_hist", "EERIE hist", "#e87ba4")]


def plot_scale():
    """Figure 3 (per_K_vs_scale.png): area-mean per K against the scale the series is aggregated to."""
    import json
    s = json.load(open(f"{OUT}/grid_per_K_summary.json"))
    sc, ens = s["scales"], s["cmip6_scales"]
    x = np.arange(len(SCALES))
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.0), sharey=True)
    for ax, (kind, field, title) in zip(axes, (("wind", "area_mean_perK", "a  wind, all events"),
                                               ("wind", "area_mean_perK_meanremoved", "b  wind, mean change removed"),
                                               ("comb", "area_mean_perK", "c  wind + solar, all events"))):
        e = ens[kind]
        if field == "area_mean_perK":
            lo, hi = np.array([e[c]["p10_p90"] for c in SCALES]).T
            ax.fill_between(x, lo, hi, color="#d6d5d0", lw=0, zorder=0, label="CMIP6 10–90 %")
        med = [e[c]["median" if field == "area_mean_perK" else "median_meanremoved"] for c in SCALES]
        ax.plot(x, med, color=INK, lw=2, marker="o", ms=4, label="CMIP6 median", zorder=3)
        for key, lab, col in RUNS:
            r = sc.get(f"{key}|{kind}")
            if r:
                ax.plot(x, [r[c][field] for c in SCALES], color=col, lw=1.5, marker="o", ms=3.5, label=lab, zorder=2)
        ax.axhline(0, color="#999999", lw=0.6, zorder=1)
        ax.set_xticks(x)
        ax.set_xticklabels(SCALE_LAB, fontsize=7.5)
        ax.tick_params(labelsize=7.5, colors=INK)
        ax.grid(axis="y", color="#e6e6e3", lw=0.6)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.set_title(title, loc="left", fontsize=8.5, color=INK)
    axes[0].set_ylabel("area-mean change (%/K)", fontsize=8, color=INK)
    # one legend entry for the ensemble: median line drawn over its 10-90 % band
    h, lab = axes[0].get_legend_handles_labels()
    d = dict(zip(lab, h))
    hand = [(d.pop("CMIP6 10–90 %"), d.pop("CMIP6 median"))] + list(d.values())
    labs = ["CMIP6 median, 10–90 % range"] + list(d.keys())
    axes[-1].legend(hand, labs, fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1))
    fn = f"{OUT}/per_K_vs_scale.png"
    fig.savefig(fn, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {fn}")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    main()
