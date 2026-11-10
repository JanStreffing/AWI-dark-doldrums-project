"""
Figures for the internal report dunkelflaute/reports/regime_checks/regime_checks.tex (how Fig. 3c and 3d of the
paper are obtained and how robust they are). Plotting only; every number comes from

  results/regimes/regime_patterns.npz                  (regime_analysis.py)       ERA5 centroids
  results/regimes/regime_decomposition_robustness.json (regime_decomposition_robustness.py) Fig. 3c terms
  results/regimes/regime_decomposition_robustness.npz  bootstrap and null draws
  results/regimes/regime_G_robustness.json             (regime_G_robustness.py)   Fig. 3d bootstrap and null
  winter-mean G of the TCO95 control: recomputed with regime_G_robustness.comparison / winter_means (no new
  definitions), the series of its tco95_control_winter_G.png

The regime_G_robustness figures are redrawn at report size (7 pt) from its json, in % per K as in Fig. 3d
(percentiles are converted with the monotone per-K map 100 (exp(ln(1 + x/100) / dT) - 1)).

Environment: levante, module load python3/2025.01-gcc-13.3.0
Usage: python3 regime_checks_report_figs.py
Output: dunkelflaute/reports/regime_checks/figs/*.png and tables/*.tex (table bodies, input by the report)
"""
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regimes")
import regime_analysis as ra               # noqa: E402
import region_events as ev                 # noqa: E402
WIN_6H = ev.WIN
import regime_G_robustness as rgr          # noqa: E402
import regions as rg                       # noqa: E402
ev.WIN = WIN_6H

import cartopy.crs as ccrs                 # noqa: E402
import matplotlib                          # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt            # noqa: E402
from matplotlib.colors import BoundaryNorm # noqa: E402
from matplotlib.lines import Line2D        # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

RES = f"{_ROOT}/dunkelflaute/results/regimes"
FIG = f"{_ROOT}/dunkelflaute/reports/regime_checks/figs"
RUNS = ["9km", "31km", "100km", "eerie_ctl"]
# run colours of the paper (paper/scripts/common/dfroles.py RUN)
COL = {"9km": "#2a78d6", "31km": "#eb6834", "100km": "#1baf7a", "eerie_ctl": "#eda100"}
LAB = {"9km": "AWI-CM3 9 km", "31km": "AWI-CM3 31 km", "100km": "AWI-CM3 100 km", "eerie_ctl": "IFS-FESOM 9 km"}
REG = ["NAO+", "NAO-", "ScBL", "AR"]
REGL = {"NAO+": "NAO+", "NAO-": "NAO−", "ScBL": "ScBL", "AR": "AR"}
EU = ra.EU
W = 6.5                                    # inches, LaTeX text width
NULLBOX = "0.86"

plt.rcParams.update({"font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 7,
                     "ytick.labelsize": 7, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
                     "ytick.major.size": 2.5, "savefig.dpi": 300, "font.family": "sans-serif"})


def per_k(x, dT):
    return 100 * (np.exp(np.log1p(np.asarray(x, float) / 100) / dT) - 1)


def clean(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def save(fig, name):
    fig.savefig(f"{FIG}/{name}", dpi=300)
    plt.close(fig)
    print(f"wrote {FIG}/{name}")


# ------------------------------------------------------------------------------------------------ 1 centroids
def fig_centroids():
    z = np.load(f"{RES}/regime_patterns.npz")
    rs = json.load(open(f"{RES}/regime_summary.json"))
    lat, lon, C = z["lat"], z["lon"], z["centroid_hPa"]
    LO, LA = np.meshgrid(lon, lat)
    PC = ccrs.PlateCarree()
    lev = np.arange(-16, 17, 2)
    cmap = plt.get_cmap("RdBu_r", len(lev) + 1)
    norm = BoundaryNorm(lev, cmap.N, extend="both")
    fig, axs = plt.subplots(1, 4, figsize=(W, 1.95), subplot_kw=dict(projection=PC))
    boxes = dict(rg.BOXES, DE=rg.DE_BOX)
    for ax, k in zip(axs, range(4)):
        ax.set_extent([-80, 40, 20, 80], PC)
        ax.set_aspect("auto")
        m = ax.contourf(LO, LA, C[k], levels=lev, cmap=cmap, norm=norm, extend="both", transform=PC)
        ax.contour(LO, LA, C[k], levels=[l for l in lev if l != 0], colors="0.3", linewidths=0.25, transform=PC)
        ax.coastlines(resolution="110m", linewidth=0.35, color="0.15")
        for r in EU:
            la0, la1, lo0, lo1 = boxes[r]
            ax.add_patch(Rectangle((lo0, la0), lo1 - lo0, la1 - la0, fill=False, ec="k", lw=0.5, transform=PC))
        nm = str(z["names"][k])
        ax.set_title(f"{REGL[nm]}  (ERA5 {100 * rs['era5']['freq'][nm]:.0f} %)", pad=2)
        ax.spines["geo"].set_linewidth(0.5)
    cax = fig.add_axes([0.3, 0.2, 0.4, 0.035])
    cb = fig.colorbar(m, cax=cax, orientation="horizontal", ticks=lev[::2])
    cb.set_label("centroid msl anomaly (hPa); black boxes: the 8 European regions", labelpad=1)
    cb.ax.tick_params(length=1.5, pad=1)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.91, bottom=0.3, wspace=0.05)
    save(fig, "centroids.png")


# ------------------------------------------------------------------------------------------------ 2 Fig. 3c robustness
def fig_decomp_boot_null(D, Z):
    """Row 1: total, frequency and within term per K; row 2: within term per regime. Dot = observed,
    bar = winter bootstrap 5-95 %, grey box = TCO95 block null 5-95 %."""
    fig, axs = plt.subplots(2, 4, figsize=(W, 3.6), sharey="row")
    for i, name in enumerate(RUNS):
        r = D["runs"][name]
        bs = Z[f"{name}__boot"]
        nv = list(r["block_null"])[0]
        nl = Z[f"{name}__null__{nv}"]
        for row, cols, labels in ((0, [0, 1, 2], ["total", "freq.", "within"]),
                                  (1, [7, 8, 9, 10], [REGL[x] for x in REG])):
            ax = axs[row, i]
            for j, c in enumerate(cols):
                lo, hi = np.percentile(nl[:, c], [5, 95])
                ax.add_patch(Rectangle((j - 0.32, lo), 0.64, hi - lo, color=NULLBOX, lw=0, zorder=1))
                b = np.percentile(bs[:, c], [5, 95])
                obs = (r["observed"]["total"], r["observed"]["frequency"], r["observed"]["within"])[c] if c < 3 \
                    else r["observed"]["within_k"][c - 7]
                ax.plot([j, j], b, color=COL[name], lw=1.6, solid_capstyle="butt", zorder=2)
                ax.plot(j, obs, "o", ms=4, mfc="white" if False else COL[name], mec=COL[name],
                        mew=1.0, zorder=3)
            ax.axhline(0, color="k", lw=0.5, zorder=0)
            ax.set_xticks(range(len(cols)), labels)
            ax.set_xlim(-0.6, len(cols) - 0.4)
            clean(ax)
            if row == 0:
                ax.set_title(f"{LAB[name]}\n{r['n_winters'][0]} vs {r['n_winters'][1]} winters", pad=2)
    axs[0, 0].set_ylabel("events per winter per K")
    axs[1, 0].set_ylabel("within-regime part\nper regime (per K)")
    h = [Line2D([], [], ls="", marker="o", ms=4, color="0.3", label="observed"),
         Line2D([], [], color="0.3", lw=1.6, label="winter bootstrap 5–95 %"),
         Patch(color=NULLBOX, label="TCO95 control block null 5–95 %")]
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.1, right=0.995, top=0.9, bottom=0.12, hspace=0.35, wspace=0.12)
    save(fig, "decomp_boot_null.png")


def fig_decomp_split(D):
    """Rows: within-regime term, frequency term (per K). Columns: (a) per region, summed over regimes;
    (b) per regime, mean over the 8 regions. Bars: the four runs."""
    fig, axs = plt.subplots(2, 2, figsize=(W, 3.5), gridspec_kw=dict(width_ratios=[2.1, 1]), sharey=True)
    w = 0.2
    off = (np.arange(4) - 1.5) * w
    for i, name in enumerate(RUNS):
        r = D["runs"][name]
        for row, key, kk in ((0, "within_perK", "within_k"), (1, "frequency_perK", "freq_k")):
            axs[row, 0].bar(np.arange(len(EU)) + off[i], [r["regions"][q][key] for q in EU], w * 0.9, color=COL[name], lw=0)
            axs[row, 1].bar(np.arange(4) + off[i], r["observed"][kk], w * 0.9, color=COL[name], lw=0)
    for row in range(2):
        for col, lab in ((0, EU), (1, [REGL[x] for x in REG])):
            ax = axs[row, col]
            ax.axhline(0, color="k", lw=0.5)
            ax.set_xticks(range(len(lab)), lab)
            clean(ax)
        axs[row, 0].set_ylabel("events per winter per K")
    axs[0, 0].set_title("(a) within-regime term per region, sum over regimes", loc="left", pad=2)
    axs[0, 1].set_title("(b) within-regime term per regime", loc="left", pad=2)
    axs[1, 0].set_title("(c) frequency term per region", loc="left", pad=2)
    axs[1, 1].set_title("(d) frequency term per regime", loc="left", pad=2)
    h = [Patch(color=COL[n], label=LAB[n]) for n in RUNS]
    fig.legend(handles=h, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0), handlelength=1.4)
    fig.subplots_adjust(left=0.08, right=0.995, top=0.95, bottom=0.12, wspace=0.06, hspace=0.35)
    save(fig, "decomp_split.png")


# ------------------------------------------------------------------------------------------------ 3 Fig. 3d robustness
def fig_G_boot_null(G):
    cols = G["columns"]
    fig, axs = plt.subplots(1, 4, figsize=(W, 2.2), sharey=True)
    for ax, name in zip(axs, RUNS):
        r = G["runs"][name]
        dT = r["dT"]
        nv = list(r["block_null"])[0]
        for j, c in enumerate(cols):
            lo, hi = per_k(r["block_null"][nv]["p05_p95"][c], dT)
            ax.add_patch(Rectangle((j - 0.32, lo), 0.64, hi - lo, color=NULLBOX, lw=0, zorder=1))
            ax.plot([j, j], per_k(r["boot_p05_p95"][c], dT), color=COL[name], lw=1.6, solid_capstyle="butt", zorder=2)
            ax.plot(j, r["change_per_K"][c], "o", ms=4, color=COL[name], zorder=3)
        ax.axhline(0, color="k", lw=0.5, zorder=0)
        ax.set_xticks(range(len(cols)), [REGL.get(c, c) for c in cols], rotation=45, ha="right", rotation_mode="anchor")
        ax.set_title(f"{LAB[name]}\n{r['n_winters'][0]} vs {r['n_winters'][1]} winters", pad=2)
        clean(ax)
    axs[0].set_ylabel("change in mean G (% per K)")
    h = [Line2D([], [], ls="", marker="o", ms=4, color="0.3", label="observed"),
         Line2D([], [], color="0.3", lw=1.6, label="winter bootstrap 5–95 %"),
         Patch(color=NULLBOX, label="TCO95 control block null 5–95 %")]
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.08, right=0.995, top=0.86, bottom=0.25, wspace=0.1)
    save(fig, "G_boot_null.png")


def fig_G_cmip6(G):
    cols = G["columns"]
    M = G["cmip6"]
    ms = sorted(M, key=lambda m: M[m]["change_per_K"]["ScBL"])
    fig, axs = plt.subplots(1, len(cols), figsize=(W, 3.3), sharey=True)
    for ax, c in zip(axs, cols):
        for i, m in enumerate(ms):
            v = M[m]
            lo, hi = per_k(v["boot_p05_p95"][c], v["dT"])
            col = "#c0392b" if hi < 0 else ("#2a78d6" if lo > 0 else "0.45")
            ax.plot([lo, hi], [i, i], color=col, lw=1.1)
            ax.plot(v["change_per_K"][c], i, "o", color=col, ms=2.6)
        ax.axvline(0, color="k", lw=0.5)
        e = G["cmip6_ensemble"][c]
        ax.set_title(f"{REGL.get(c, c)}: 90 % range\nbelow 0 in {e['n_boot90_below_zero']}/{e['n_models']}", pad=2)
        ax.set_xlabel("% per K")
        clean(ax)
    axs[0].set_yticks(range(len(ms)), ms, fontsize=6)
    axs[0].tick_params(axis="y", length=0)
    fig.subplots_adjust(left=0.17, right=0.995, top=0.88, bottom=0.13, wspace=0.12)
    save(fig, "G_cmip6_boot.png")


def fig_null_series(G):
    pat = np.load(f"{RES}/regime_patterns.npz")
    ra.NAMES_ORDER = [str(x) for x in pat["names"]]
    P = rgr.comparison("100km", pat)
    x = rgr.winter_means(*P["ctl"]["SN"])
    yrs = P["ctl"]["years"]
    sd = G["null_control"]["interannual_sd_pct"]
    fig, ax = plt.subplots(figsize=(W, 1.7))
    ax.plot(yrs, x[:, 4], color="0.25", lw=0.9, label=f"all days (sd {sd['all']:.1f} %)")
    ax.plot(yrs, x[:, 2], color="#7b3294", lw=0.9, label=f"ScBL days (sd {sd['ScBL']:.1f} %)")
    ax.set_ylabel("winter-mean G,\n8 regions (m s$^{-1}$)")
    ax.set_xlim(yrs[0] - 1, yrs[-1] + 1)
    ax.set_xticks(range(1960, 2100, 20))
    ax.legend(frameon=False, ncol=2, loc="upper left", bbox_to_anchor=(0, 1.12))
    clean(ax)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.9, bottom=0.14)
    save(fig, "null_series.png")


def fig_lam_vs_G(D):
    fig, ax = plt.subplots(figsize=(3.3, 2.5))
    mk = {"NAO+": "o", "NAO-": "s", "ScBL": "D", "AR": "^"}
    xs, ys = [], []
    for name in RUNS:
        r = D["runs"][name]
        lc = np.mean([r["regions"][q]["lam100_ctl"] for q in EU], 0)
        lf = np.mean([r["regions"][q]["lam100_fut"] for q in EU], 0)
        dl = per_k(100 * (lf / lc - 1), r["dT"])
        dg = [D["G_footing"][name]["per_K"]["box"][k] for k in REG]
        for k in range(4):
            ax.plot(dg[k], dl[k], mk[REG[k]], color=COL[name], ms=4.5, mec="white", mew=0.4)
        xs += dg
        ys += list(dl)
    ax.axhline(0, color="k", lw=0.5)
    ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("change in mean G within regime (% per K)")
    ax.set_ylabel("change in events per regime day,\n8-region mean (% per K)")
    h = [Line2D([], [], ls="", marker=mk[k], color="0.35", ms=4, label=REGL[k]) for k in REG] + \
        [Line2D([], [], ls="", marker="o", color=COL[n], ms=4, label=LAB[n].replace("AWI-CM3 ", "")) for n in RUNS]
    ax.legend(handles=h, frameon=False, fontsize=6, ncol=2, loc="lower left", handletextpad=0.2, columnspacing=0.6)
    print(f"lam_vs_G: r = {np.corrcoef(xs, ys)[0, 1]:.3f}, fit {np.polyfit(xs, ys, 1)}")
    clean(ax)
    fig.subplots_adjust(left=0.2, right=0.98, top=0.97, bottom=0.17)
    save(fig, "lam_vs_G.png")


# ------------------------------------------------------------------------------------------------ tables
TAB = f"{_ROOT}/dunkelflaute/reports/regime_checks/tables"
SHORT = {"9km": "9 km", "31km": "31 km", "100km": "100 km", "eerie_ctl": "IFS-FESOM"}


def fp(p, n):
    """p-value of a null with n placements; 0 means none of the n placements reached the observed value."""
    if p == 0:
        return f"$<${1 / n:.0e}".replace("e-0", "e-").replace("1e-", "10^{-") + "}$" if n >= 1000 else f"$<${1 / n:.3f}"
    return f"{p:.3f}"


def fpp(p, n):
    if p == 0:
        return "$<$0.001" if n >= 1000 else f"$<${np.ceil(1000 / n) / 1000:.3f}"
    return f"{p:.3f}"


def rng_(a, nd=3):
    return "{}" + f"[{a[0]:+.{nd}f}, {a[1]:+.{nd}f}]".replace("+", "$+$").replace("-", "$-$")   # {} guards \\[


def sg(x, nd=3):
    return f"{x:+.{nd}f}".replace("+", "$+$").replace("-", "$-$")


def tables(D, G):
    os.makedirs(TAB, exist_ok=True)
    # T1: worked example, DE at 9 km
    e = D["runs"]["9km"]["regions"]["DE"]
    L = []
    for k, nm in enumerate(REG):
        L.append(f"{REGL[nm].replace(chr(0x2212), '$-$')} & {e['f_ctl'][k]:.1f} & {e['f_fut'][k]:.1f} & {e['e_ctl'][k]:.2f} & {e['e_fut'][k]:.2f} & "
                 f"{e['lam100_ctl'][k]:.2f} & {e['lam100_fut'][k]:.2f} & {sg(e['freq_k_perK'][k])} & {sg(e['within_k_perK'][k])} \\\\")
    L.append("\\midrule")
    L.append(f"sum & {sum(e['f_ctl']):.1f} & {sum(e['f_fut']):.1f} & {sum(e['e_ctl']):.2f} & {sum(e['e_fut']):.2f} & & & "
             f"{sg(e['frequency_perK'])} & {sg(e['within_perK'])} \\\\")
    open(f"{TAB}/example_DE_9km.tex", "w").write("\n".join(L) + "\n")
    # T2: terms per K, bootstrap, null
    L = []
    for name in RUNS:
        r = D["runs"][name]
        nv = list(r["block_null"])[0]
        b = r["block_null"][nv]
        n = b["n_placements"]
        cells = []
        for t in ("total", "frequency", "within"):
            cells.append(f"\\makecell{{{sg(r['observed'][t])}\\\\{rng_(r['boot']['p05_p95'][t])}}} & {fpp(b['p_one_sided'][t], n)}")
        L.append(f"{SHORT[name]} & {r['dT']:.2f} & {r['n_winters'][0]}/{r['n_winters'][1]} & " + " & ".join(cells)
                 + f" & {r['within_positive_regions']} \\\\")
    open(f"{TAB}/terms.tex", "w").write("\n".join(L) + "\n")
    # T3: within-regime part per regime
    L = []
    for k, nm in enumerate(REG):
        cells = []
        for name in RUNS:
            r = D["runs"][name]
            nv = list(r["block_null"])[0]
            b = r["block_null"][nv]
            cells.append(f"\\makecell{{{sg(r['observed']['within_k'][k])}\\\\{rng_(r['boot']['within_k_p05_p95'][k])}}} & {fpp(b['within_k_p_one_sided'][k], b['n_placements'])}")
        L.append(f"{REGL[nm].replace(chr(0x2212), '$-$')} & " + " & ".join(cells) + " \\\\")
    L.append("\\midrule")
    L.append("ScBL share & " + " & ".join(
        f"\\multicolumn{{2}}{{c}}{{{100 * D['runs'][n]['observed']['within_k'][2] / D['runs'][n]['observed']['within']:.0f}\\,\\%}}" for n in RUNS) + " \\\\")
    open(f"{TAB}/within_per_regime.tex", "w").write("\n".join(L) + "\n")
    # T4: sensitivities
    L = []
    for name in RUNS:
        r = D["runs"][name]
        ev_ = r["events"]
        tie = 100 * (ev_["ctl"]["tied"] + ev_["fut"]["tied"]) / (ev_["ctl"]["n"] + ev_["fut"]["n"])
        L.append(f"{SHORT[name]} & {tie:.0f}\\,\\% & {sg(r['observed']['frequency'])} & {sg(r['observed']['within'])} & "
                 f"{sg(r['symmetric']['frequency'])} & {sg(r['symmetric']['within'])} & "
                 f"{sg(r['ties_split']['frequency'])} & {sg(r['ties_split']['within'])} & {sg(r['ties_split']['within_k'][2])} & {sg(r['observed']['within_k'][2])} \\\\")
    open(f"{TAB}/sensitivity.tex", "w").write("\n".join(L) + "\n")
    # T5: G per regime per K: ours box [boot] p, bilinear; CMIP6
    L = []
    cols = G["columns"]
    for c in cols:
        cells = []
        for name in RUNS:
            r = G["runs"][name]
            nv = list(r["block_null"])[0]
            p = r["block_null"][nv]["p_one_sided"][c]
            n = r["block_null"][nv]["n_placements"]
            cells.append(f"\\makecell{{{sg(r['change_per_K'][c], 2)}\\\\{rng_(per_k(r['boot_p05_p95'][c], r['dT']), 2)}\\\\"
                         f"{{\\scriptsize $p$ {fpp(p, n)}; bil. {sg(D['G_footing'][name]['per_K']['bilinear'][c], 2)}}}}}")
        cm = np.array([v["change_per_K"][c] for v in G["cmip6"].values()])
        e = G["cmip6_ensemble"][c]
        L.append(f"{REGL.get(c, c).replace(chr(0x2212), '$-$')} & " + " & ".join(cells) +
                 f" & \\makecell{{{sg(np.median(cm), 2)}\\\\{rng_(np.percentile(cm, [10, 90]), 2)}\\\\{{\\scriptsize {e['n_negative']} neg.; {e['n_boot90_below_zero']} sig.}}}} \\\\")
    open(f"{TAB}/G_per_regime.tex", "w").write("\n".join(L) + "\n")
    print(f"wrote {TAB}/*.tex")


def main():
    os.makedirs(FIG, exist_ok=True)
    D = json.load(open(f"{RES}/regime_decomposition_robustness.json"))
    Z = np.load(f"{RES}/regime_decomposition_robustness.npz")
    G = json.load(open(f"{RES}/regime_G_robustness.json"))
    fig_centroids()
    fig_decomp_boot_null(D, Z)
    fig_decomp_split(D)
    fig_G_boot_null(G)
    fig_G_cmip6(G)
    fig_lam_vs_G(D)
    fig_null_series(G)
    tables(D, G)


if __name__ == "__main__":
    main()
