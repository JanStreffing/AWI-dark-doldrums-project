"""
Schematic (not data) of a Scandinavian blocking day in the present and in a warmer climate, for the internal
report reports/regime_checks (Jan 2026-11-10). It draws the finding of regime_G_attribution.py as a cartoon.

What was measured: under ScBL the background (seasonal-cycle) part of the G change is close to zero, and the regime
high hardly changes in peak, extent, or position. G of the day-to-day deviations from the regime composite falls by
5-9 %, while their rms over the European box changes by only 0 to -4 %, i.e. a longer length scale at nearly
unchanged amplitude.

What is drawn: the background and the high are the same in both climates. The disturbances (lows) keep their
positions, get WIDENING_FACTOR wider, and their depths are scaled so that the area-weighted rms of the disturbance
field over the map area (RMS_BOX) equals the present one. So the lows are wider and somewhat shallower at the same
overall amplitude, and their slopes are gentler. The change is exaggerated for clarity; the real ScBL weakening of G
is about 1-2 % per K.

Anchoring (blocking_era5_compare.py, 2026-11-10): the layout follows the ERA5 ScBL day 2024-11-06 (German
Dunkelflaute): a broad high over the North Sea and Scandinavia (about 1028 hPa at 60 N), a deep Atlantic low to the
west, lower pressure to the north-east, and pressure at 40 N about equal to that at 62 N (ERA5 -1.2 hPa, here
-1.3 hPa). A weak low over northern Italy (Genoa-type) makes it a borderline day: Germany-box G about 5 m s-1 today,
just above the illustrative low-wind level, and below it in (b). Which way G moves at one spot depends on where the
lows sit: over 90 tested placements of the southern low, the map-mean G falls in all 90 and the Germany-box G in 54.
The data show a shift of the distribution of days; one day can only illustrate it.

Earlier versions (same day, git history and methods notes): widening the lows in place at the same depth raised the
rms by 25 % and G by 9 %; stretching their pattern outward raised the rms by 16 %. Both contradicted the measurement.
A symmetric omega layout with a strong westerly background put higher pressure south of the high, unlike ERA5.

Top row: sea-level pressure (isobars every 4 hPa) and geostrophic speed G = |grad p| / (rho f) (shading), with the
area below an illustrative low-wind level of G hatched and the Germany box outlined. Bottom: north-south section along 10 E through the Italian low, Germany, and the high, G and pressure.

Two versions: main() draws the exaggerated one (WIDENING_FACTOR 1.25, blocking_schematic.png, no longer in the
report), main_realistic() the one the report uses (WIDEN_REALISTIC 1.13: slope of the lows -6 %, map-mean G -5 %,
inside the measured ranges; made visible with a difference map and zoomed side views, blocking_schematic_realistic.png).
With 1.13 the Germany-box G falls in 51 of the 90 placements of the southern low, the map-mean G in all 90.

Environment: levante, module load python3/2025.01-gcc-13.3.0
Usage: python3 blocking_schematic.py
Output: dunkelflaute/reports/regime_checks/figs/blocking_schematic.png, blocking_schematic_realistic.png
"""
import sys

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute"
sys.path.insert(0, f"{_ROOT}/scripts/regions")
import regions as rg  # noqa: E402

OUT = f"{_ROOT}/reports/regime_checks/figs/blocking_schematic.png"
WIDENING_FACTOR = 1.25
RHO, OMEGA, R_EARTH = 1.25, 7.292e-5, 6.371e6
G_LOW = 4.0                     # m/s, the illustrative low-wind level hatched in the maps
LON_SECTION = 10.0              # north-south side view through Italy, Germany and the high
PRESENT, FUTURE = "#2a78d6", "#eb6834"   # dfroles colours of 9 km and 31 km, reused as present / future
BACKGROUND_SLOPE = 0.006        # hPa per km northward, a weak westerly background (seasonal cycle)

lon = np.linspace(-70, 80, 601)
lat = np.linspace(25, 85, 241)
LON, LAT = np.meshgrid(lon, lat)

# (lon, lat, amplitude hPa, e-folding half width km)
HIGH = (8.0, 60.0, 20.0, 1500.0)                   # the regime pattern, unchanged
DISTURBANCES = [(-25.0, 58.0, -24.0, 600.0),       # deep Atlantic low west of the block
                (38.0, 70.0, -10.0, 600.0),        # low to the north-east
                (10.0, 46.0, -8.0, 500.0)]         # weak Genoa-type low south of the block

RMS_BOX = (LAT > 35) & (LAT < 72) & (LON > -30) & (LON < 45)   # the map area
W_AREA = np.cos(np.deg2rad(LAT))
_LA0, _LA1, _LO0, _LO1 = rg.DE_BOX
DE = (LAT >= _LA0) & (LAT <= _LA1) & (LON >= _LO0) & (LON <= _LO1)


def blob(lo, la, amp, width_km):
    dx = np.deg2rad(LON - lo) * R_EARTH * np.cos(np.deg2rad(0.5 * (LAT + la))) / 1e3
    dy = np.deg2rad(LAT - la) * R_EARTH / 1e3
    return amp * np.exp(-(dx ** 2 + dy ** 2) / (2 * width_km ** 2))


def rms(f):
    return float(np.sqrt(np.sum(W_AREA[RMS_BOX] * f[RMS_BOX] ** 2) / W_AREA[RMS_BOX].sum()))


def area_mean(f, mask=RMS_BOX):
    return float(np.sum(W_AREA[mask] * f[mask]) / W_AREA[mask].sum())


def disturbance_field(dist):
    return sum(blob(*x) for x in dist)


def widened(widen):
    """Disturbances at fixed positions, widths scaled by widen, depths scaled so that the rms of the disturbance
    field over RMS_BOX stays at its present value."""
    d = [(lo, la, amp, w * widen) for lo, la, amp, w in DISTURBANCES]
    if widen != 1.0:
        s = rms(disturbance_field(DISTURBANCES)) / rms(disturbance_field(d))
        d = [(lo, la, amp * s, w) for lo, la, amp, w in d]
    return d


def pressure(widen):
    p = 1013.0 - BACKGROUND_SLOPE * (LAT - 50.0) * 111.0
    return p + blob(*HIGH) + disturbance_field(widened(widen))


def geostrophic(p):
    phi = np.deg2rad(LAT)
    dpdy = np.gradient(p * 100.0, np.deg2rad(lat) * R_EARTH, axis=0)
    dpdx = np.gradient(p * 100.0, np.deg2rad(lon), axis=1) / (R_EARTH * np.cos(phi))
    f = 2 * OMEGA * np.sin(phi)
    return np.hypot(dpdx, dpdy) / (RHO * f)


def section_panel(ax, lines, title, xlim=(36, 72), legend=None):
    """North-south section along LON_SECTION: G (top axis use) as lines, pressure as thin lines on a twin axis,
    the Germany band shaded, the low-wind level drawn. lines: (label, colour, style, lat, lon, p, G), p and G (lat, lon)."""
    ax.axvspan(_LA0, _LA1, color="0.9", lw=0, zorder=0)
    ax.text(0.5 * (_LA0 + _LA1), 27.5, "Germany", ha="center", va="center", fontsize=6.5, color="0.3")
    ax.axhline(G_LOW, color="0.3", lw=0.7, ls="--")
    ax.text(xlim[1] - 0.3, G_LOW + 0.4, "low-wind level", ha="right", va="bottom", fontsize=6, color="0.3")
    axp = ax.twinx()
    for lab, col, ls, la, lo, p, g in lines:
        i = int(np.argmin(np.abs(lo - LON_SECTION)))
        sel = (la >= xlim[0]) & (la <= xlim[1])
        ax.plot(la[sel], g[sel, i], color=col, lw=1.8, ls=ls, label=f"G, {lab}")
        axp.plot(la[sel], p[sel, i], color=col, lw=0.8, ls=":", label=f"pressure, {lab}")
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 30)
    axp.set_ylim(990, 1040)
    ax.set_xlabel(f"latitude along {LON_SECTION:.0f}$^\\circ$E, south to north (dashed line in the maps)", fontsize=7)
    ax.set_ylabel("geostrophic speed G (m s$^{-1}$), thick", fontsize=7)
    axp.set_ylabel("sea-level pressure (hPa), dotted", fontsize=7)
    for a in (ax, axp):
        a.tick_params(labelsize=6)
    ax.set_title(title, fontsize=8)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = axp.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=6, ncol=2, frameon=False,
              **(legend or dict(loc="upper right")))


def main():
    fields = {"present": pressure(1.0), "future": pressure(WIDENING_FACTOR)}
    speeds = {k: geostrophic(v) for k, v in fields.items()}

    proj = ccrs.LambertConformal(central_longitude=10, central_latitude=55, standard_parallels=(45, 65))
    fig = plt.figure(figsize=(17.0 / 2.54, 13.0 / 2.54), dpi=200)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.85], hspace=0.22, wspace=0.05,
                          left=0.08, right=0.9, top=0.95, bottom=0.12)
    levels_g = np.arange(0, 22.5, 2.5)
    for i, (key, title) in enumerate((("present", "(a) Present-day blocking day"),
                                      ("future", "(b) Same blocking day, warmer climate"))):
        ax = fig.add_subplot(gs[0, i], projection=proj)
        ax.set_extent([-30, 45, 35, 74], crs=ccrs.PlateCarree())
        cf = ax.contourf(LON, LAT, speeds[key], levels=levels_g, cmap="Blues", extend="max",
                         transform=ccrs.PlateCarree())
        cs = ax.contour(LON, LAT, fields[key], levels=np.arange(972, 1052, 4), colors="0.15", linewidths=0.6,
                        transform=ccrs.PlateCarree())
        ax.clabel(cs, fmt="%d", fontsize=5.5, inline_spacing=1)
        plt.rcParams["hatch.color"] = PRESENT if key == "present" else FUTURE
        plt.rcParams["hatch.linewidth"] = 0.6
        ax.contourf(LON, LAT, speeds[key], levels=[0, G_LOW], colors="none", hatches=["//////"],
                    transform=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE.with_scale("110m"), linewidth=0.4, edgecolor="0.45")
        ax.plot([_LO0, _LO1, _LO1, _LO0, _LO0], [_LA0, _LA0, _LA1, _LA1, _LA0], color="k", lw=0.8,
                transform=ccrs.PlateCarree())
        ax.plot(np.full(200, LON_SECTION), np.linspace(36, 72, 200), color="k", lw=0.8, ls="--",
                transform=ccrs.PlateCarree())
        ax.text(HIGH[0], HIGH[1] + 2, "H", ha="center", va="center", fontsize=11, weight="bold",
                transform=ccrs.PlateCarree())
        for lo, la, _, _ in DISTURBANCES:
            ax.text(lo, la, "L", ha="center", va="center", fontsize=8, weight="bold", color="0.15",
                    transform=ccrs.PlateCarree())
        ax.set_title(title, fontsize=8)
    cax = fig.add_axes([0.915, 0.53, 0.015, 0.36])
    cb = fig.colorbar(cf, cax=cax)
    cb.set_label("geostrophic speed G (m s$^{-1}$)", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    section_panel(fig.add_subplot(gs[1, :]),
                  [("present", PRESENT, "-", lat, lon, fields["present"], speeds["present"]),
                   ("warmer climate", FUTURE, "-", lat, lon, fields["future"], speeds["future"])],
                  "(c) Side view along the dashed line: the wind over Germany weakens toward the low-wind level")
    fig.text(0.5, 0.005, f"Schematic, not model output, laid out after ERA5 on 6 November 2024.\nIn (b) the lows are "
             f"{WIDENING_FACTOR:.2f} times wider at the same overall amplitude, exaggerated for clarity.",
             ha="center", fontsize=6, color="0.35", linespacing=1.2)
    fig.savefig(OUT)

    d0, d1 = disturbance_field(widened(1.0)), disturbance_field(widened(WIDENING_FACTOR))
    g0, g1 = area_mean(geostrophic(d0)), area_mean(geostrophic(d1))
    print(f"disturbance rms {rms(d0):.2f} -> {rms(d1):.2f} hPa, disturbance G {g0:.2f} -> {g1:.2f} m/s "
          f"({100 * (g1 / g0 - 1):+.0f} %), depth factor {widened(WIDENING_FACTOR)[0][2] / DISTURBANCES[0][2]:.2f}")
    for k in fields:
        print(f"{k}: map-mean G {area_mean(speeds[k]):.1f}, Germany-box G {area_mean(speeds[k], DE):.1f} m/s, "
              f"Germany box below G_LOW {float(np.mean(speeds[k][DE] < G_LOW)):.2f}")
    print(OUT)


WIDEN_REALISTIC = 1.13          # slope of the disturbances -6 %, map-mean G -5 %: inside the measured ranges
OUT_REALISTIC = f"{_ROOT}/reports/regime_checks/figs/blocking_schematic_realistic.png"


def main_realistic():
    """Same day with a realistic widening (WIDEN_REALISTIC), made visible through a difference map and zoomed
    side views instead of an exaggerated change. Measured under ScBL: G of the day-to-day deviations -4 to -9 % at
    0 to -4 % rms, region-mean G -3 to -12 % (100 km to 9 km), CMIP6 median about -5 %."""
    p0, p1 = pressure(1.0), pressure(WIDEN_REALISTIC)
    g0, g1 = geostrophic(p0), geostrophic(p1)
    dg = g1 - g0
    proj = ccrs.LambertConformal(central_longitude=10, central_latitude=55, standard_parallels=(45, 65))
    fig = plt.figure(figsize=(17.0 / 2.54, 12.5 / 2.54), dpi=200)
    gs = fig.add_gridspec(4, 3, height_ratios=[0.5, 0.03, 0.1, 0.62], hspace=0.1, wspace=0.12,
                          left=0.07, right=0.98, top=0.96, bottom=0.09)
    levels_g = [0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20]
    norm_g = matplotlib.colors.BoundaryNorm(levels_g, 256, extend="max")
    levels_d = np.arange(-2.0, 2.01, 0.25)
    for k, (title, g, p, col) in enumerate((("(a) Present-day blocking day", g0, p0, PRESENT),
                                            ("(b) Same day, warmer climate", g1, p1, FUTURE),
                                            ("(c) Change in G, (b) minus (a)", dg, p1, None))):
        ax = fig.add_subplot(gs[0, k], projection=proj)
        ax.set_extent([-30, 45, 35, 74], crs=ccrs.PlateCarree())
        if col is not None:
            cf = ax.contourf(LON, LAT, g, levels=levels_g, cmap="Blues", norm=norm_g, extend="max",
                             transform=ccrs.PlateCarree())
            plt.rcParams["hatch.color"] = "0.3"
            plt.rcParams["hatch.linewidth"] = 0.5
            ax.contourf(LON, LAT, g, levels=[0, G_LOW], colors="none", hatches=["//////"],
                        transform=ccrs.PlateCarree())
            cs = ax.contour(LON, LAT, p, levels=np.arange(972, 1052, 4), colors="0.15", linewidths=0.5,
                            transform=ccrs.PlateCarree())
            ax.clabel(cs, fmt="%d", fontsize=5, inline_spacing=1)
        else:
            # pcolormesh with discrete steps: cartopy 0.24 contourf fills nested polygons of this field wrongly
            norm_d = matplotlib.colors.BoundaryNorm(levels_d, 256, extend="both")
            cd = ax.pcolormesh(LON, LAT, dg, cmap="PuOr", norm=norm_d, shading="auto", transform=ccrs.PlateCarree(),
                               rasterized=True)
            ax.contour(LON, LAT, g0, levels=[G_LOW], colors=[PRESENT], linewidths=1.0, transform=ccrs.PlateCarree())
            ax.contour(LON, LAT, g1, levels=[G_LOW], colors=[FUTURE], linewidths=1.0, transform=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE.with_scale("110m"), linewidth=0.4, edgecolor="0.45")
        ax.plot([_LO0, _LO1, _LO1, _LO0, _LO0], [_LA0, _LA0, _LA1, _LA1, _LA0], color="k", lw=0.8,
                transform=ccrs.PlateCarree())
        ax.plot(np.full(200, LON_SECTION), np.linspace(36, 72, 200), color="k", lw=0.7, ls="--",
                transform=ccrs.PlateCarree())
        ax.set_title(title, fontsize=7.5)
    cb = fig.colorbar(cf, cax=fig.add_subplot(gs[1, 0:2]), orientation="horizontal", ticks=levels_g)
    cb.set_label(f"geostrophic speed G (m s$^{{-1}}$); hatched below {G_LOW:.0f}", fontsize=6.5, labelpad=1)
    cb.ax.tick_params(labelsize=6)
    cb2 = fig.colorbar(cd, cax=fig.add_subplot(gs[1, 2]), orientation="horizontal", ticks=[-2, -1, 0, 1, 2])
    cb2.set_label("change in G (m s$^{-1}$)", fontsize=6.5, labelpad=1)
    cb2.ax.tick_params(labelsize=6)

    i = int(np.argmin(np.abs(lon - LON_SECTION)))
    sel = (lat >= 36) & (lat <= 72)
    ax1 = fig.add_subplot(gs[3, 0:2])
    ax1.axvspan(_LA0, _LA1, color="0.9", lw=0, zorder=0)
    ax1.tick_params(labelsize=6)
    ax1.axhline(G_LOW, color="0.3", lw=0.7, ls="--")
    ax1.text(71.7, G_LOW + 0.15, "low-wind level", ha="right", va="bottom", fontsize=6, color="0.3")
    ax1.text(0.5 * (_LA0 + _LA1), 11.3, "Germany", ha="center", va="center", fontsize=6.5, color="0.3")
    ax1.plot(lat[sel], g0[sel, i], color=PRESENT, lw=1.8, label="present")
    ax1.plot(lat[sel], g1[sel, i], color=FUTURE, lw=1.8, label="warmer climate")
    ax1.fill_between(lat[sel], g1[sel, i], g0[sel, i], where=g1[sel, i] < g0[sel, i], color=FUTURE, alpha=0.25, lw=0)
    ax1.set_ylim(0, 12)
    ax1.set_xlim(36, 72)
    ax1.set_ylabel("G (m s$^{-1}$)", fontsize=7)
    ax1.set_xlabel("latitude", fontsize=7)
    ax1.legend(fontsize=6, loc="upper left", frameon=False)
    ax1.set_title(f"(d) G along {LON_SECTION:.0f}$^\\circ$E, south to north (dashed line in the maps)", fontsize=7.5)

    # (e) idealised day-to-day disturbance: same height, slope weaker by the measured amount
    slope_drop = 0.06                                  # map disturbances here -6 %; measured -4 to -9 %
    wl0 = 2000.0                                       # km
    wl1 = wl0 / (1.0 - slope_drop)
    amp = 8.0                                          # hPa
    x = np.linspace(-0.25 * wl0, 1.75 * wl0, 600)
    ax3 = fig.add_subplot(gs[3, 2])
    ax3.axhline(0, color="0.6", lw=0.5)
    for wl, col, lab in ((wl0, PRESENT, "present"), (wl1, FUTURE, "warmer")):
        ax3.plot(x, amp * np.sin(2 * np.pi * x / wl), color=col, lw=1.6, label=lab)
    ax3.annotate("same height", xy=(1.3 * wl0, amp), xytext=(1.3 * wl0, amp + 4.5), fontsize=6, color="0.25",
                 ha="center", arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))
    ax3.annotate(f"gentler slope, {-100 * slope_drop:.0f} %", xy=(0.5 * wl0, 0.0), xytext=(0.5 * wl0, -amp - 4.5),
                 fontsize=6, color="0.25", ha="center", arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))
    ax3.text(0.02, 0.97, "measured: height 0 to $-$4 %,\nslope $-$4 to $-$9 %", transform=ax3.transAxes,
             ha="left", va="top", fontsize=5.5, color="0.3")
    ax3.set_ylim(-amp - 7, amp + 7)
    ax3.set_xlim(x[0], x[-1])
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.set_xlabel("distance", fontsize=7)
    ax3.set_ylabel("pressure departure", fontsize=7)
    ax3.set_title("(e) Idealised disturbance", fontsize=7.5)
    ax3.legend(fontsize=6, loc="lower right", frameon=False)
    fig.text(0.5, 0.005, f"Schematic, not model output, laid out after ERA5 on 6 November 2024. The lows are "
             f"{100 * (WIDEN_REALISTIC - 1):.0f} % wider at the same overall amplitude, within the measured range.",
             ha="center", fontsize=6, color="0.35")
    fig.savefig(OUT_REALISTIC)
    d0, d1 = disturbance_field(widened(1.0)), disturbance_field(widened(WIDEN_REALISTIC))
    print(f"realistic: disturbance G {100 * (area_mean(geostrophic(d1)) / area_mean(geostrophic(d0)) - 1):+.1f} %, "
          f"map G {100 * (area_mean(g1) / area_mean(g0) - 1):+.1f} %, Germany G {area_mean(g0, DE):.2f} -> "
          f"{area_mean(g1, DE):.2f}, Germany below level {area_mean((g0 < G_LOW).astype(float), DE):.2f} -> "
          f"{area_mean((g1 < G_LOW).astype(float), DE):.2f}")
    print(OUT_REALISTIC)


if __name__ == "__main__":
    main()
    main_realistic()
