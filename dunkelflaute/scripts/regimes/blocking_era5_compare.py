"""
Compare the blocking schematic (blocking_schematic.py, present-day panel) with ERA5: one real Scandinavian-blocking
Dunkelflaute day and the ERA5 composite of all ScBL days, in the same style. For the internal report
reports/regime_checks (Jan 2026-11-10).

ERA5 data and regime labels exactly as regime_analysis.py: daily-mean msl on the 2.5 deg grid 20-80 N, 80 W-40 E,
Sep-Mar 1979/80-2024/25, anomalies from the stored mean + 3 harmonics, sqrt(cos lat) weights, projection on the stored
ERA5 EOFs, nearest stored centroid. G = |grad p| / (rho f) with centred differences on that grid, rho = 1.25 kg m-3.
The real day is 2024-11-06 (labelled ScBL; Germany-box G 2.7 m s-1, below the 5th percentile of 4.3 over all days),
during the German Dunkelflaute of early November 2024. The composite shows the mean msl of all ScBL days and the mean
of their daily G (not G of the mean field, which would be smoother).

Environment: levante, module load python3/2025.01-gcc-13.3.0
Usage: python3 blocking_era5_compare.py
Output: dunkelflaute/reports/regime_checks/figs/blocking_era5_compare.png
"""
import sys

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute"
sys.path.insert(0, f"{_ROOT}/scripts/regimes")
sys.path.insert(0, f"{_ROOT}/scripts/regions")
import blocking_schematic as bs  # noqa: E402
import regime_analysis as ra  # noqa: E402
import regions as rg  # noqa: E402

OUT = f"{_ROOT}/reports/regime_checks/figs/blocking_era5_compare.png"
DAY = np.datetime64("2024-11-06")
R_EARTH, OMEGA, RHO = 6.371e6, 7.292e-5, 1.25


def era5():
    P = np.load(f"{_ROOT}/results/regimes/regime_patterns.npz")
    dates, x, lat, lon = ra.load_msl("era5_-")
    a = ra.anomalies(dates, x, P["coef"]) * P["wts"].reshape(1, -1)
    lab = ra.assign(a, P["eofs"], P["cent"])
    names = [str(n) for n in P["names"]]
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    phi = np.deg2rad(LA)
    dpdy = np.gradient(x, np.deg2rad(lat) * R_EARTH, axis=1)
    dpdx = np.gradient(x, np.deg2rad(lon), axis=2) / (R_EARTH * np.cos(phi))
    G = np.hypot(dpdx, dpdy) / (RHO * 2 * OMEGA * np.sin(phi))
    return dates, x / 100.0, G, lat, lon, LA, LO, lab, names


def main():
    dates, p, G, lat, lon, LA, LO, lab, names = era5()
    sc = names.index("ScBL")
    i = int(np.where(dates == DAY)[0][0])
    assert lab[i] == sc, "the example day is not labelled ScBL"
    la0, la1, lo0, lo1 = rg.DE_BOX
    de = (LA >= la0) & (LA <= la1) & (LO >= lo0) & (LO <= lo1)
    print(f"{DAY}: {names[lab[i]]}, Germany-box G {G[i][de].mean():.1f} m/s; ScBL days {int((lab == sc).sum())}")

    panels = [("(a) ERA5, 6 November 2024 (ScBL)", LO, LA, p[i], G[i]),
              (f"(b) ERA5, mean of all {int((lab == sc).sum())} ScBL days", LO, LA, p[lab == sc].mean(0),
               G[lab == sc].mean(0)),
              ("(c) Schematic, present day", bs.LON, bs.LAT, bs.pressure(1.0), bs.geostrophic(bs.pressure(1.0)))]

    proj = ccrs.PlateCarree()   # ERA5 grid ends at 40 E, so no conic projection (it would leave a blank wedge)
    fig = plt.figure(figsize=(17.0 / 2.54, 11.5 / 2.54), dpi=200)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.75], hspace=0.3, wspace=0.04,
                          left=0.07, right=0.9, top=0.94, bottom=0.11)
    levels_g = np.arange(0, 22.5, 2.5)
    plt.rcParams["hatch.linewidth"] = 0.6
    for k, (title, X, Y, P, S) in enumerate(panels):
        ax = fig.add_subplot(gs[0, k], projection=proj)
        ax.set_extent([-30, 40, 35, 74], crs=ccrs.PlateCarree())
        cf = ax.contourf(X, Y, S, levels=levels_g, cmap="Blues", extend="max", transform=ccrs.PlateCarree())
        plt.rcParams["hatch.color"] = "0.25"
        ax.contourf(X, Y, S, levels=[0, bs.G_LOW], colors="none", hatches=["//////"], transform=ccrs.PlateCarree())
        cs = ax.contour(X, Y, P, levels=np.arange(972, 1052, 4), colors="0.15", linewidths=0.5,
                        transform=ccrs.PlateCarree())
        ax.clabel(cs, fmt="%d", fontsize=5, inline_spacing=1)
        ax.add_feature(cfeature.COASTLINE.with_scale("110m"), linewidth=0.4, edgecolor="0.45")
        ax.plot([lo0, lo1, lo1, lo0, lo0], [la0, la0, la1, la1, la0], color="k", lw=0.8, transform=ccrs.PlateCarree())
        ax.plot(np.full(200, bs.LON_SECTION), np.linspace(36, 72, 200), color="k", lw=0.7, ls="--",
                transform=ccrs.PlateCarree())
        ax.set_title(title, fontsize=7)
    cax = fig.add_axes([0.915, 0.53, 0.012, 0.38])
    cb = fig.colorbar(cf, cax=cax)
    cb.set_label("G (m s$^{-1}$)", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    bs.section_panel(fig.add_subplot(gs[1, :]),
                     [("ERA5 6 Nov 2024", "k", "-", lat, lon, panels[0][3], panels[0][4]),
                      ("ERA5 ScBL mean", "0.55", "-", lat, lon, panels[1][3], panels[1][4]),
                      ("schematic, present", bs.PRESENT, "-", bs.lat, bs.lon, panels[2][3], panels[2][4])],
                     f"(d) Side view along the dashed line ({bs.LON_SECTION:.0f}$^\\circ$E)",
                     legend=dict(loc="center left", bbox_to_anchor=(0.0, 0.42)))
    fig.text(0.5, 0.005, f"Hatching: G below {bs.G_LOW:.0f} m s$^{{-1}}$. Box: Germany region. ERA5 on the 2.5$^\\circ$ "
             "grid of the regime analysis.", ha="center", fontsize=6, color="0.35")
    fig.savefig(OUT)
    print(OUT)


if __name__ == "__main__":
    main()
