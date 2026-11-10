"""
Diagnostics for the co-authors (not a paper figure): do the free-running models have the ERA5
Euro-Atlantic EOFs and weather regimes that regime_analysis.py projects them onto?

Preprocessing exactly as regime_analysis.py (imported): Sep-Mar daily msl on the 2.5 deg grid
20-80 N, 80 W-40 E, anomalies from mean + 3 annual harmonics, weights sqrt(cos lat), EOFs by SVD,
10 PCs, k-means k = 4 with 50 initialisations. ERA5 EOFs, centroids and weights from
results/regimes/regime_patterns.npz.

Per simulation (run and winter span as in region_events.COMPARISONS; future runs use the control
seasonal cycle of their comparison for the paper's assignment):
  1. own EOFs 1-10 from the run's own anomalies (own seasonal cycle, i.e. exactly the ERA5 recipe),
     explained variance, sign aligned to ERA5, pattern r with the same-rank ERA5 EOF, best match
     among model EOF 1-10 for each ERA5 EOF 1-4, and the congruence of each ERA5 EOF with the span of
     the model's EOF 1-4 (|projection|, 1 = the ERA5 pattern lies in the model's leading subspace,
     so a low same-rank r with a high congruence means reordering or mixing, not a different mode)
  2. own k-means centroids (10 PCs of its own EOFs) matched to the 4 ERA5 centroids by the
     permutation with the largest summed r; frequencies
  3. composites of the days the paper assigns to each ERA5 regime (anomalies from the control
     seasonal cycle, projection on ERA5 EOFs, nearest ERA5 centroid); r and rms amplitude ratio
     against the ERA5 centroid. Also with the run's own seasonal cycle (sensitivity).
  4. fraction of the model's (own-anomaly, centred) variance captured by ERA5 EOF 1-10 vs its own
     EOF 1-10.
Pattern r: uncentred, area-weighted (cos lat) correlation of the hPa maps.
Sampling reference: the same diagnostics for ERA5 sub-periods of 17 and 10 winters (the lengths of
the 9 km control and future).
CMIP6 (summary only): day/psl from extract_psl.py, historical 1951-1980 and ssp585 2071-2100
(control = historical), the windows geostrophic.py uses.

Environment: levante, module load python3/2025.01-gcc-13.3.0
Usage: python3 eof_diagnostics.py      (--summary-only: redraw the two summary figures from the json)
Output: dunkelflaute/results/regimes/eof_diagnostics/{eofs_<run>,regimes_<run>,summary_ours,
        summary_cmip6}.png and eof_diagnostics.json
"""
import glob
import itertools
import json
import os
import sys

import cartopy.crs as ccrs
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                               # noqa: E402
import numpy as np                                                            # noqa: E402
from matplotlib.colors import BoundaryNorm                                    # noqa: E402

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regimes")
import regime_analysis as ra                                                  # noqa: E402

OUT = f"{_ROOT}/dunkelflaute/results/regimes/eof_diagnostics"
DC = f"{_ROOT}/dunkelflaute/data/cmip6"
NPC, K, NE = ra.NPC, ra.K, 4

# (label, run, w0, w1, control run, c0, c1)  spans as region_events.COMPARISONS
OURS = [
    ("9 km ctl 1953-69", "dart9km_1950c", 1953, 1969, "dart9km_1950c", 1953, 1969),
    ("9 km fut 2083-92", "dart9km_2080c", 2083, 2092, "dart9km_1950c", 1953, 1969),
    ("31 km ctl 1953-69", "tco319_ctl1950d", 1953, 1969, "tco319_ctl1950d", 1953, 1969),
    ("31 km fut 2081-99", "tco319_ssp585d", 2081, 2099, "tco319_ctl1950d", 1953, 1969),
    ("100 km ctl 1953-2099", "tco95_ctl1950g", 1953, 2099, "tco95_ctl1950g", 1953, 2099),
    ("100 km fut 2081-99", "tco95_hi1950g", 2081, 2099, "tco95_ctl1950g", 1953, 2099),
    ("IFS-FESOM ctl 1972-2050", "eerie_ctl", 1972, 2050, "eerie_ctl", 1972, 2050),
    ("IFS-FESOM hist 1971-2000", "eerie_hist", 1971, 2000, "eerie_hist", 1971, 2000),
    ("IFS-FESOM ssp245 2041-70", "eerie_ssp245", 2041, 2070, "eerie_ctl", 1972, 2050),
]
ERA5_REF = [("ERA5 1980-96 (17 w)", 1980, 1996), ("ERA5 2009-25 (17 w)", 2009, 2025), ("ERA5 2016-25 (10 w)", 2016, 2025)]

P = np.load(f"{ra.OUT}/regime_patterns.npz")
LAT, LON = P["lat"], P["lon"]
W = P["wts"].ravel()                     # sqrt(cos lat)
C = W ** 2                               # cos lat area weight
E_EOF = P["eofs"]                        # (10, 1225), weighted space, unit norm
E_CENT = P["cent"]                       # (4, 10) in ERA5 PC space
E_CMAP = P["centroid_hPa"].reshape(K, -1)
NAMES = [str(n) for n in P["names"]]


def prc(a, b):
    """Uncentred cos(lat)-weighted pattern correlation of hPa maps (last axis = grid)."""
    return float((C * a * b).sum() / np.sqrt((C * a * a).sum() * (C * b * b).sum()))


def rms(a):
    return float(np.sqrt((C * a * a).sum() / C.sum()))


def eofs_of(a_w):
    """a_w: weighted anomalies (days, grid). Returns unit EOFs (10, grid), explained variance (10),
    hPa-per-sd maps (10, grid), total centred variance."""
    x = a_w - a_w.mean(0)
    _, s, Vt = np.linalg.svd(x, full_matrices=False)
    ev = s ** 2 / (s ** 2).sum()
    maps = Vt[:NPC] * (s[:NPC, None] / np.sqrt(len(x))) / W / 100.0
    return Vt[:NPC], ev[:NPC], maps, float((s ** 2).sum())


def era5_eof_maps():
    d, x, *_ = ra.load_msl("era5_-")
    w = ra.winter_of_day(d)
    k = (w >= 1980) & (w <= 2025)
    d, x = d[k], x[k]
    a = ra.anomalies(d, x, P["coef"]) * W
    V, ev, maps, _ = eofs_of(a)
    sgn = np.sign((V * E_EOF).sum(1))
    assert np.all(np.abs((V * E_EOF).sum(1)) > 0.999), "ERA5 EOFs not reproduced"
    return maps * sgn[:, None], ev, d, x


def diagnose(dates, x, dates_c, x_c, kmeans_seed=0):
    """dates/x: the period analysed; dates_c/x_c: its comparison's control period (for the paper's
    anomalies). Returns a dict of numbers and a dict of maps."""
    own_coef = ra.seasonal_fit(dates, x)
    ctl_coef = ra.seasonal_fit(dates_c, x_c)
    a_own = ra.anomalies(dates, x, own_coef) * W
    a_pap = ra.anomalies(dates, x, ctl_coef) * W
    V, ev, maps, tot = eofs_of(a_own)
    # 1. EOFs: same rank, sign aligned to ERA5
    R = V @ E_EOF.T                                      # (model 10, ERA5 10), unit vectors -> pattern r
    sgn = np.where(np.diag(R) < 0, -1.0, 1.0)
    V, maps, R = V * sgn[:, None], maps * sgn[:, None], R * sgn[:, None]
    same = [float(R[k, k]) for k in range(NE)]
    best_for_era5 = [dict(model_eof=int(np.argmax(np.abs(R[:, k]))) + 1, r=float(np.abs(R[:, k]).max())) for k in range(NE)]
    best_for_model = [dict(era5_eof=int(np.argmax(np.abs(R[k, :]))) + 1, r=float(np.abs(R[k, :]).max())) for k in range(NE)]
    congr4 = [float(np.linalg.norm(V[:NE] @ E_EOF[k])) for k in range(NE)]
    # 4. variance captured
    xc = a_own - a_own.mean(0)
    var_era5_10 = float(((xc @ E_EOF.T) ** 2).sum() / tot)
    var_own_10 = float(ev.sum())
    # 2. own k-means
    pcs = a_own @ V.T
    _, lab = ra.kmeans(pcs, K, seed=kmeans_seed)
    own_cmaps = np.stack([a_own[lab == j].mean(0) / W / 100.0 for j in range(K)])
    own_freq = np.array([(lab == j).mean() for j in range(K)])
    rmat = np.array([[prc(own_cmaps[j], E_CMAP[k]) for k in range(K)] for j in range(K)])
    perm = max(itertools.permutations(range(K)), key=lambda p: sum(rmat[p[k], k] for k in range(K)))
    own_names, _ = ra.name_clusters(own_cmaps.reshape(K, len(LAT), len(LON)), LAT, LON)
    # 3. composites of the paper's assignment
    comp = {}
    for tag, a in (("ctlclim", a_pap), ("ownclim", a_own)):
        lb = ra.assign(a, E_EOF, E_CENT)
        cm = np.stack([a[lb == k].mean(0) / W / 100.0 if (lb == k).any() else np.zeros(a.shape[1]) for k in range(K)])
        comp[tag] = dict(maps=cm, freq=[float((lb == k).mean()) for k in range(K)],
                         r=[prc(cm[k], E_CMAP[k]) for k in range(K)],
                         amp=[rms(cm[k]) / rms(E_CMAP[k]) for k in range(K)])
    num = dict(
        n_days=int(len(dates)), n_winters=int(len(np.unique(ra.winter_of_day(dates)))),
        explained_variance=[float(v) for v in ev], eof_r_same_rank=same,
        eof_best_model_match_for_era5=best_for_era5, eof_best_era5_match_for_model=best_for_model,
        eof_congruence_with_model_eof1_4=congr4,
        eof_r_matrix_model1_4_x_era5_1_4=R[:NE, :NE].round(3).tolist(),
        var_frac_era5_eof1_10=var_era5_10, var_frac_own_eof1_10=var_own_10,
        kmeans=dict(matched_to=[NAMES[k] for k in range(K)],
                    r=[float(rmat[perm[k], k]) for k in range(K)],
                    freq=[float(own_freq[perm[k]]) for k in range(K)],
                    amp=[rms(own_cmaps[perm[k]]) / rms(E_CMAP[k]) for k in range(K)],
                    own_names_by_rule=[own_names[perm[k]] for k in range(K)],
                    r_matrix_own_x_era5=rmat.round(3).tolist()),
        assigned={t: {kk: v for kk, v in c.items() if kk != "maps"} for t, c in comp.items()},
    )
    mp = dict(eof=maps[:NE], own_cent=own_cmaps[list(perm)], comp=comp["ctlclim"]["maps"])
    return num, mp


def sel_winters(d, x, w0, w1):
    w = ra.winter_of_day(d)
    k = (w >= w0) & (w <= w1)
    return d[k], x[k]


# ---------------------------------------------------------------------------------------------- plots
PC = ccrs.PlateCarree()
LONS, LATS = np.meshgrid(LON, LAT)
plt.rcParams.update({"font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 7,
                     "ytick.labelsize": 7, "legend.fontsize": 7})


def sym_levels(vmax, n=8):
    steps = np.array([0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5])
    st = steps[np.searchsorted(steps, vmax / n)] if vmax / n <= steps[-1] else np.ceil(vmax / n)
    m = int(np.ceil(vmax / st))
    return np.arange(-m, m + 1) * st


def mapax(ax, field, levels, cmap):
    ax.set_extent([-80, 40, 20, 80], PC)
    ax.set_aspect("auto")
    norm = BoundaryNorm(levels, plt.get_cmap(cmap).N, extend="both")
    m = ax.contourf(LONS, LATS, field.reshape(len(LAT), len(LON)), levels=levels, cmap=cmap, norm=norm,
                    extend="both", transform=PC)
    ax.contour(LONS, LATS, field.reshape(len(LAT), len(LON)), levels=[l for l in levels if l != 0], colors="0.25",
               linewidths=0.25, transform=PC)
    ax.coastlines(resolution="110m", linewidth=0.35, color="0.2")
    ax.spines["geo"].set_linewidth(0.5)
    return m


def fig_eofs(label, run, num, mp, era_maps, era_ev):
    lev = sym_levels(max(np.abs(era_maps[:NE]).max(), np.abs(mp["eof"]).max()))
    fig, axs = plt.subplots(2, NE, figsize=(7.2, 3.3), subplot_kw=dict(projection=PC))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.86, bottom=0.17, wspace=0.05, hspace=0.42)
    for k in range(NE):
        m = mapax(axs[0, k], era_maps[k], lev, "PuOr_r")
        axs[0, k].set_title(f"ERA5 EOF{k + 1} ({100 * era_ev[k]:.1f} %)", pad=2)
        mapax(axs[1, k], mp["eof"][k], lev, "PuOr_r")
        b = num["eof_best_era5_match_for_model"][k]
        axs[1, k].set_title(f"model EOF{k + 1} ({100 * num['explained_variance'][k]:.1f} %), r = {num['eof_r_same_rank'][k]:.2f}\n"
                            f"best: ERA5 EOF{b['era5_eof']} |r| = {b['r']:.2f}", pad=2)
    fig.suptitle(f"{label} ({run}, {num['n_winters']} winters): own EOFs vs ERA5; variance in EOF1-10: "
                 f"own {100 * num['var_frac_own_eof1_10']:.0f} %, ERA5 EOFs {100 * num['var_frac_era5_eof1_10']:.0f} %",
                 fontsize=7, y=0.985)
    cb = fig.colorbar(m, ax=axs, orientation="horizontal", fraction=0.05, aspect=45, pad=0.04, shrink=0.7)
    cb.set_label("msl regression on standardised PC (hPa per sd); sign aligned to the same-rank ERA5 EOF")
    fig.savefig(f"{OUT}/eofs_{run}.png", dpi=200)
    plt.close(fig)


def fig_regimes(label, run, num, mp, era_freq):
    allm = np.concatenate([E_CMAP, mp["own_cent"], mp["comp"]])
    lev = sym_levels(np.abs(allm).max())
    fig, axs = plt.subplots(K, 3, figsize=(6.4, 6.0), subplot_kw=dict(projection=PC))
    fig.subplots_adjust(left=0.05, right=0.98, top=0.855, bottom=0.10, wspace=0.05, hspace=0.42)
    km, asg = num["kmeans"], num["assigned"]["ctlclim"]
    for k in range(K):
        m = mapax(axs[k, 0], E_CMAP[k], lev, "RdBu_r")
        axs[k, 0].set_title(f"ERA5 centroid\nfreq. {100 * era_freq[k]:.0f} %", pad=2)
        mapax(axs[k, 1], mp["own_cent"][k], lev, "RdBu_r")
        axs[k, 1].set_title(f"own k-means cluster (rule name: {km['own_names_by_rule'][k]})\n"
                            f"freq. {100 * km['freq'][k]:.0f} %, r = {km['r'][k]:.2f}, amp {km['amp'][k]:.2f}", pad=2)
        mapax(axs[k, 2], mp["comp"][k], lev, "RdBu_r")
        axs[k, 2].set_title(f"days assigned to ERA5 {NAMES[k]}\n"
                            f"freq. {100 * asg['freq'][k]:.0f} %, r = {asg['r'][k]:.2f}, amp {asg['amp'][k]:.2f}", pad=2)
        axs[k, 0].text(-0.07, 0.5, NAMES[k], transform=axs[k, 0].transAxes, rotation=90, ha="center", va="center",
                       fontsize=8, fontweight="bold")
    fig.suptitle(f"{label} ({run}, {num['n_winters']} winters): weather regimes\n"
                 "middle: the model's own k = 4 clusters (own EOFs), matched to ERA5; right: composite of the days the paper\n"
                 "assigns to the ERA5 regime (control-climate anomalies); amp = rms amplitude ratio to the ERA5 centroid",
                 fontsize=7, y=0.995)
    cb = fig.colorbar(m, ax=axs, orientation="horizontal", fraction=0.04, aspect=45, pad=0.03, shrink=0.7)
    cb.set_label("msl anomaly (hPa)")
    fig.savefig(f"{OUT}/regimes_{run}.png", dpi=200)
    plt.close(fig)


def fig_summary(rows, fname, title, row_h=0.26):
    """rows: list of (label, num). Heatmap of r: EOF same rank, EOF best match, ERA5 EOF in model
    EOF1-4 span, own k-means, assigned composite; text columns: explained variance and variance
    fractions. Layout in inches so that it holds for any number of rows."""
    blocks = [("EOF\nsame rank", lambda n: n["eof_r_same_rank"], None),
              ("EOF best match\nin model 1-10", lambda n: [b["r"] for b in n["eof_best_model_match_for_era5"]],
               lambda n: [b["model_eof"] for b in n["eof_best_model_match_for_era5"]]),
              ("ERA5 EOF in\nmodel 1-4 span", lambda n: n["eof_congruence_with_model_eof1_4"], None),
              ("own k-means\ncentroid", lambda n: n["kmeans"]["r"], None),
              ("assigned-day\ncomposite", lambda n: n["assigned"]["ctlclim"]["r"], None)]
    nb = len(blocks)
    M = np.array([sum((list(f(n)) for _, f, _ in blocks), []) for _, n in rows])
    Wi, left, hm_w, tab_w = 7.2, 1.55, 3.95, 1.6       # inches
    top, bot = 0.75, 1.1
    H = top + bot + row_h * len(rows)
    fig = plt.figure(figsize=(Wi, H))
    ax = fig.add_axes([left / Wi, bot / H, hm_w / Wi, row_h * len(rows) / H])
    cmap = plt.get_cmap("viridis", 9)
    norm = BoundaryNorm([0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98, 1.0], cmap.N, extend="min")
    im = ax.imshow(M, cmap=cmap, norm=norm, aspect="auto")
    for i, (_, n) in enumerate(rows):
        col = 0
        for _, f, idx in blocks:
            vals = f(n)
            ids = idx(n) if idx else None
            for k, v in enumerate(vals):
                s = f"{v:.2f}".replace("0.", ".") if v < 0.995 else "1.0"
                if ids and ids[k] != k + 1:
                    s += f"\n#{ids[k]}"
                ax.text(col, i, s, ha="center", va="center", fontsize=5.5 if ids and ids[k] != k + 1 else 6.5,
                        color="white" if v < 0.85 else "black", linespacing=0.85)
                col += 1
    xt = sum(([f"EOF{k + 1}" for k in range(NE)] if b < 3 else NAMES for b in range(nb)), [])
    ax.set_xticks(range(len(xt)), xt, rotation=90)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows])
    ax.tick_params(length=0)
    for b in range(1, nb):
        ax.axvline(4 * b - 0.5, color="white", lw=2)
    for b, (bt, _, _) in enumerate(blocks):
        ax.text(4 * b + 1.5, -0.65, bt, ha="center", va="bottom", linespacing=1.0)
    # text table on the right
    tx = fig.add_axes([(left + hm_w + 0.05) / Wi, bot / H, tab_w / Wi, row_h * len(rows) / H])
    tx.set_xlim(0, 1)
    tx.set_ylim(len(rows) - 0.5, -0.5)
    tx.axis("off")
    tx.text(0.3, -0.65, "expl. var.\nEOF1-4 (%)", ha="center", va="bottom", linespacing=1.0)
    tx.text(0.77, -0.65, "10 EOFs (%)\nown / ERA5", ha="center", va="bottom", linespacing=1.0)
    for i, (_, n) in enumerate(rows):
        tx.text(0.3, i, " ".join(f"{100 * v:.0f}" for v in n["explained_variance"][:NE]), ha="center", va="center",
                fontsize=6.5, family="monospace")
        tx.text(0.77, i, f"{100 * n['var_frac_own_eof1_10']:.0f}/{100 * n['var_frac_era5_eof1_10']:.0f}",
                ha="center", va="center", fontsize=6.5, family="monospace")
    cax = fig.add_axes([(left + 0.5) / Wi, 0.36 / H, (hm_w - 1.0) / Wi, 0.07 / H])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("pattern r with ERA5 (uncentred, cos-lat weighted); ERA5 EOF1-4: 20, 15, 14, 8 %", labelpad=1)
    fig.text(0.5, 1 - 0.08 / H, title, ha="center", va="top")
    fig.savefig(f"{OUT}/{fname}", dpi=200)
    plt.close(fig)


def summaries(res):
    rows = list(res["era5_subperiods"].items()) + list(res["ours"].items())
    fig_summary(rows, "summary_ours.png",
                "Pattern r with ERA5: our runs, and ERA5 sub-periods of 17 and 10 winters as a sampling reference")
    cm = [(k.replace(" 1951-80", "").replace(" 2071-2100", ""), v) for k, v in res["cmip6"].items()]
    fig_summary(cm, "summary_cmip6.png",
                "Pattern r with ERA5: CMIP6 day/psl, historical 1951-80 and ssp585 2071-2100 (control climate: historical)",
                row_h=0.2)


# ---------------------------------------------------------------------------------------------- CMIP6
def cmip6_dates(ymd):
    """ymd -> datetime64[D] for the harmonics (360-day calendars: day clipped), as geostrophic.dates_of."""
    ml = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    return np.array([np.datetime64(f"{a:04d}-{b:02d}-{min(c, ml[b - 1]):02d}") for a, b, c in ymd])


def cmip6_runs():
    for fn in sorted(glob.glob(f"{DC}/psl_*.npz")):
        model = os.path.basename(fn)[4:-4]
        Z = np.load(fn)
        ii = np.where(Z["lat"] >= 20.0 - 1e-6)[0]
        jj = np.where((Z["lon"] >= -80.0 - 1e-6) & (Z["lon"] <= 40.0 + 1e-6))[0]
        assert np.allclose(Z["lat"][ii], LAT) and np.allclose(Z["lon"][jj], LON)
        S = {}
        for exp, (y0, y1) in (("historical", (1951, 1980)), ("ssp585", (2071, 2100))):
            if f"{exp}|psl" not in Z:
                break
            ymd = Z[f"{exp}|ymd"]
            w = np.where(ymd[:, 1] >= 9, ymd[:, 0] + 1, ymd[:, 0])
            k = (w >= y0) & (w <= y1) & np.isin(ymd[:, 1], (9, 10, 11, 12, 1, 2, 3))
            x = Z[f"{exp}|psl"][k][:, ii][:, :, jj].astype(np.float64).reshape(k.sum(), -1)
            S[exp] = (cmip6_dates(ymd[k]), x)
        if len(S) == 2 and np.isfinite(S["historical"][1]).all() and np.isfinite(S["ssp585"][1]).all():
            yield model, S


# ---------------------------------------------------------------------------------------------- main
def main():
    os.makedirs(OUT, exist_ok=True)
    era_maps, era_ev, ed, ex = era5_eof_maps()
    S = json.load(open(f"{ra.OUT}/regime_summary.json"))
    era_freq = [S["era5"]["freq"][n] for n in NAMES]
    res = dict(method=__doc__.strip().split("\n\n")[1:4], names=NAMES,
               era5=dict(explained_variance=[float(v) for v in era_ev], freq=era_freq),
               era5_subperiods={}, ours={}, cmip6={})

    for label, w0, w1 in ERA5_REF:
        d, x = sel_winters(ed, ex, w0, w1)
        num, _ = diagnose(d, x, d, x)
        res["era5_subperiods"][label] = num
        print(f"{label}: EOF r {np.round(num['eof_r_same_rank'], 2)}, kmeans r {np.round(num['kmeans']['r'], 2)}", flush=True)

    cache = {}
    for label, run, w0, w1, crun, c0, c1 in OURS:
        if run not in cache:
            cache[run] = ra.load_msl(run)[:2]
        if crun not in cache:
            cache[crun] = ra.load_msl(crun)[:2]
        d, x = sel_winters(*cache[run], w0, w1)
        dc, xc = sel_winters(*cache[crun], c0, c1)
        x, xc = x.reshape(len(d), -1), xc.reshape(len(dc), -1)
        num, mp = diagnose(d, x, dc, xc)
        num.update(run=run, winters=[w0, w1], control=[crun, c0, c1])
        res["ours"][label] = num
        fig_eofs(label, run, num, mp, era_maps, era_ev)
        fig_regimes(label, run, num, mp, era_freq)
        print(f"{label}: EOF r {np.round(num['eof_r_same_rank'], 2)} ev {np.round(100 * np.array(num['explained_variance'][:4]), 1)}"
              f" | kmeans r {np.round(num['kmeans']['r'], 2)} | assigned r {np.round(num['assigned']['ctlclim']['r'], 2)}"
              f" | var own/ERA5 {num['var_frac_own_eof1_10']:.2f}/{num['var_frac_era5_eof1_10']:.2f}", flush=True)

    for model, Sx in cmip6_runs():
        dc, xc = Sx["historical"]
        for exp, tag in (("historical", "hist 1951-80"), ("ssp585", "ssp585 2071-2100")):
            d, x = Sx[exp]
            num, _ = diagnose(d, x, dc, xc)
            res["cmip6"][f"{model} {tag}"] = num
        print(f"{model}: hist EOF r {np.round(res['cmip6'][f'{model} hist 1951-80']['eof_r_same_rank'], 2)}", flush=True)

    json.dump(res, open(f"{OUT}/eof_diagnostics.json", "w"), indent=1)
    summaries(res)
    print(f"Saved {OUT}")


if __name__ == "__main__":
    if "--summary-only" in sys.argv:          # redraw the summary figures from the json
        summaries(json.load(open(f"{OUT}/eof_diagnostics.json")))
    else:
        main()
