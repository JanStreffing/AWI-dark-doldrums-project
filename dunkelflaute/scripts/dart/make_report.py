"""
Visual report: Germany Dunkelflaute (Mockert et al. 2023 definition) over time for
  validation_era5_2015_2026 (ERA5 reconstruction), 1950c and 2080c (TCo1279-DART).

Reads dunkelflaute/results/<period>/combined_cf_*.npz, writes
  dunkelflaute/results/<period>/figures/overview.png       (per-period development)
  dunkelflaute/results/comparison/*.png                    (cross-period)
  dunkelflaute/results/<period>/summary.json
  dunkelflaute/results/report_dunkelflaute.html            (self-contained)
"""
import base64
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
os.chdir(_ROOT)
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/dart")
from compute_germany_dunkelflaute_dart_mockert import DT_H, mockert_events  # noqa: E402

R = "dunkelflaute/results"
PERIODS = {
    "validation_era5_2015_2026": dict(label="ERA5 2015–2026 (validation)", short="ERA5", color="#1baf7a"),
    "1950c": dict(label="DART 1950C", short="1950C", color="#2a78d6"),
    "2080c": dict(label="DART 2080C", short="2080C", color="#eb6834"),
}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
MONTH_ORDER = [9, 10, 11, 12, 1, 2, 3]
MONTH_NAMES = ["Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar"]
MOCKERT_PUBLISHED = 4.0   # events per winter, Mockert et al. (2023), Germany 1979-2018

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": MUTED, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "figure.facecolor": "white", "axes.facecolor": "white", "text.color": INK,
})


def load(period):
    """-> dict(time, winter, solar, on, off, combined, roll) + events DataFrame-like list."""
    if period == "validation_era5_2015_2026":
        d = np.load(f"{R}/{period}/combined_cf_2015_2026_sepmar_mockertweights.npz", allow_pickle=True)
        out = dict(time=d["time"].astype("datetime64[s]"), winter=d["winter"],
                   solar=d["cf_solar_recon"], on=d["cf_on_recon"], off=d["cf_off_recon"],
                   combined=d["cf_combined_recon"], roll=d["roll_recon"])
        extra = {"ERA5 real": d["roll_real"], "SMARD": d["roll_smard"]}
    else:
        d = np.load(f"{R}/{period}/combined_cf_{period}_sepmar_mockertweights.npz")
        out = dict(time=d["time"], winter=d["winter"], solar=d["solar"], on=d["on"], off=d["off"],
                   combined=d["combined"], roll=d["roll"])
        extra = {}
    return out, extra


def events_of(winter, time, roll):
    rows = []
    for w in np.unique(winter):
        m = winter == w
        _, ev = mockert_events(roll[m])
        t = time[m]
        for s, e in ev:
            rows.append((int(w), t[s], t[e - 1], (e - s) * DT_H))
    return rows


def winter_table(o, ev):
    ws = np.unique(o["winter"])
    tab = []
    for w in ws:
        m = o["winter"] == w
        evw = [r for r in ev if r[0] == w]
        flagged_h = sum(r[3] for r in evw)
        tab.append(dict(winter=int(w), n_events=len(evw), flagged_days=flagged_h / 24.0,
                        flagged_pct=100 * flagged_h / (m.sum() * DT_H),
                        longest_days=max([r[3] for r in evw], default=0) / 24.0,
                        mean_cf=float(o["combined"][m].mean()),
                        min_roll=float(np.nanmin(o["roll"][m]))))
    return tab


def summarize(o, ev, tab):
    dur = np.array([r[3] for r in ev]) / 24.0
    n_w = len(tab)
    ev_pw = np.array([t["n_events"] for t in tab], float)
    ws = np.array([t["winter"] for t in tab], float)
    sl = stats.linregress(ws, ev_pw) if n_w > 2 else None
    sl_cf = stats.linregress(ws, [t["mean_cf"] for t in tab]) if n_w > 2 else None
    return dict(
        n_winters=n_w, winters=f"{int(ws[0])}–{int(ws[-1])}", n_events=len(ev),
        events_per_winter=round(len(ev) / n_w, 2), events_per_winter_sd=round(float(ev_pw.std(ddof=1)), 2),
        median_dur_days=round(float(np.median(dur)), 2), mean_dur_days=round(float(dur.mean()), 2),
        max_dur_days=round(float(dur.max()), 2),
        flagged_pct=round(float(np.mean([t["flagged_pct"] for t in tab])), 2),
        mean_cf_combined=round(float(o["combined"].mean()), 3), mean_cf_solar=round(float(o["solar"].mean()), 3),
        mean_cf_onshore=round(float(o["on"].mean()), 3), mean_cf_offshore=round(float(o["off"].mean()), 3),
        events_trend_per_decade=None if sl is None else round(sl.slope * 10, 2),
        events_trend_p=None if sl is None else round(sl.pvalue, 2),
        cf_trend_per_decade=None if sl_cf is None else round(sl_cf.slope * 10, 4),
        cf_trend_p=None if sl_cf is None else round(sl_cf.pvalue, 2),
        winters_without_event=int((ev_pw == 0).sum()),
    )


def day_of_winter(t, w):
    return (t - np.datetime64(f"{w-1}-09-01")) / np.timedelta64(1, "D")


def overview_figure(period, o, ev, tab, path):
    cfg = PERIODS[period]
    c = cfg["color"]
    ws = [t["winter"] for t in tab]
    fig = plt.figure(figsize=(12, 7.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.25])
    for k, (key, title, unit) in enumerate([("n_events", "Events per winter", "events"),
                                            ("flagged_days", "Days inside an event", "days"),
                                            ("mean_cf", "Mean combined CF", "CF")]):
        ax = fig.add_subplot(gs[0, k])
        y = [t[key] for t in tab]
        if key == "mean_cf":
            ax.plot(ws, y, color=c, lw=2, marker="o", ms=4)
            ax.set_ylim(min(y) * 0.9, max(y) * 1.08)
        else:
            ax.bar(ws, y, color=c, width=0.7)
        if len(ws) > 2:
            sl = stats.linregress(ws, y)
            ax.plot(ws, sl.intercept + sl.slope * np.array(ws), color=MUTED, lw=1.2, ls="--")
            ax.text(0.02, 0.97, f"trend {sl.slope*10:+.2f}/decade (p={sl.pvalue:.2f})", transform=ax.transAxes,
                    va="top", fontsize=8.5, color=MUTED,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.5))
        if key == "n_events":
            ax.axhline(MOCKERT_PUBLISHED, color=MUTED, lw=1, ls=":")
            ax.text(ws[0] - 0.4, MOCKERT_PUBLISHED, "Mockert published (4) ", va="bottom", ha="left", fontsize=8, color=MUTED,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1))
        ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
        ax.set_ylabel(unit); ax.set_xticks(ws[::2]); ax.tick_params(axis="x", labelsize=8.5)
    ax = fig.add_subplot(gs[1, :])
    for i, w in enumerate(ws):
        segs = [(day_of_winter(r[1], w), (r[2] - r[1]) / np.timedelta64(1, "D") + DT_H / 24.0)
                for r in ev if r[0] == w]
        ax.broken_barh(segs, (i - 0.36, 0.72), facecolors=c, edgecolor="white", linewidth=0.6)
    month_starts = [0, 30, 61, 91, 122, 153, 181]
    ax.set_xticks(month_starts); ax.set_xticklabels(MONTH_NAMES, ha="left")
    ax.set_xlim(0, 212); ax.set_yticks(range(len(ws))); ax.set_yticklabels(ws); ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_title("When events happen inside each winter (Sep → Mar); one bar = one Dunkelflaute event",
                 loc="left", fontsize=11, fontweight="bold")
    fig.suptitle(f"{cfg['label']} — Germany Dunkelflaute, Mockert et al. definition", x=0.01, ha="left",
                 fontsize=13, fontweight="bold")
    fig.savefig(path, dpi=130); plt.close(fig)


def comparison_figures(data, outdir):
    os.makedirs(outdir, exist_ok=True)
    paths = {}
    # 1. events per winter, one dot per winter + mean
    fig, ax = plt.subplots(figsize=(7.5, 4.2), constrained_layout=True)
    rng = np.random.default_rng(0)
    for i, (p, d) in enumerate(data.items()):
        y = np.array([t["n_events"] for t in d["tab"]])
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(y)), y, s=34, color=PERIODS[p]["color"], alpha=0.8,
                   edgecolor="white", linewidth=0.8)
        ax.hlines(y.mean(), i - 0.3, i + 0.3, color=INK, lw=2)
        ax.text(i + 0.33, y.mean(), f"{y.mean():.1f}", va="center", fontsize=10, fontweight="bold")
    ax.axhline(MOCKERT_PUBLISHED, color=MUTED, ls=":", lw=1)
    ax.text(len(data) - 0.5, MOCKERT_PUBLISHED, "Mockert published (4)", ha="right", va="bottom", fontsize=8.5, color=MUTED)
    ax.set_xticks(range(len(data))); ax.set_xticklabels([PERIODS[p]["short"] for p in data])
    ax.set_ylabel("events per winter"); ax.set_xlim(-0.5, len(data) - 0.4)
    ax.set_title("Events per winter — each dot is one winter, bar = mean", loc="left", fontsize=11, fontweight="bold")
    paths["events_per_winter"] = f"{outdir}/events_per_winter.png"; fig.savefig(paths["events_per_winter"], dpi=130); plt.close(fig)
    # 2. duration ECDF
    fig, ax = plt.subplots(figsize=(7.5, 4.2), constrained_layout=True)
    for p, d in data.items():
        x = np.sort([r[3] / 24 for r in d["ev"]]); y = np.arange(1, len(x) + 1) / len(x)
        ax.step(x, y, where="post", color=PERIODS[p]["color"], lw=2, label=PERIODS[p]["short"])
    ax.set_xlabel("event duration (days)"); ax.set_ylabel("share of events ≤ duration")
    ax.legend(frameon=False, loc="lower right"); ax.set_xlim(left=2)
    ax.set_title("Event duration distribution", loc="left", fontsize=11, fontweight="bold")
    paths["duration"] = f"{outdir}/duration_ecdf.png"; fig.savefig(paths["duration"], dpi=130); plt.close(fig)
    # 3. share of time flagged by month
    fig, ax = plt.subplots(figsize=(7.5, 4.2), constrained_layout=True)
    wd = 0.27
    for i, (p, d) in enumerate(data.items()):
        o = d["o"]; flagged = np.zeros(len(o["winter"]), bool)
        for w in np.unique(o["winter"]):
            m = np.where(o["winter"] == w)[0]
            f, _ = mockert_events(o["roll"][m]); flagged[m] = f
        mon = o["time"].astype("datetime64[M]").astype(int) % 12 + 1
        share = [100 * flagged[mon == m].mean() for m in MONTH_ORDER]
        d["month_share"] = share
        ax.bar(np.arange(7) + (i - 1) * wd, share, wd, color=PERIODS[p]["color"], label=PERIODS[p]["short"], edgecolor="white", linewidth=1)
    ax.set_xticks(range(7)); ax.set_xticklabels(MONTH_NAMES); ax.set_ylabel("% of time inside an event")
    ax.legend(frameon=False)
    ax.set_title("Seasonal distribution", loc="left", fontsize=11, fontweight="bold")
    paths["monthly"] = f"{outdir}/monthly_share.png"; fig.savefig(paths["monthly"], dpi=130); plt.close(fig)
    # 4. mean CF by component
    fig, ax = plt.subplots(figsize=(7.5, 4.2), constrained_layout=True)
    comps = [("solar", "Solar"), ("on", "Wind onshore"), ("off", "Wind offshore"), ("combined", "Combined")]
    for i, (p, d) in enumerate(data.items()):
        vals = [float(d["o"][k].mean()) for k, _ in comps]
        ax.bar(np.arange(4) + (i - 1) * wd, vals, wd, color=PERIODS[p]["color"], label=PERIODS[p]["short"], edgecolor="white", linewidth=1)
    ax.set_xticks(range(4)); ax.set_xticklabels([n for _, n in comps]); ax.set_ylabel("mean capacity factor, Sep–Mar")
    ax.legend(frameon=False)
    ax.set_title("What drives the difference: mean CF by component", loc="left", fontsize=11, fontweight="bold")
    paths["components"] = f"{outdir}/cf_components.png"; fig.savefig(paths["components"], dpi=130); plt.close(fig)
    return paths


def b64(path):
    return base64.b64encode(open(path, "rb").read()).decode()


def main():
    data = {}
    for p in PERIODS:
        if not os.path.exists(f"{R}/{p}") or not any(f.endswith(".npz") for f in os.listdir(f"{R}/{p}")):
            print(f"skip {p}: no data yet"); continue
        o, extra = load(p)
        ev = events_of(o["winter"], o["time"], o["roll"])
        tab = winter_table(o, ev)
        s = summarize(o, ev, tab)
        for name, roll in extra.items():           # ERA5 real / SMARD event rate, for validation context
            e2 = events_of(o["winter"], o["time"], roll)
            s[f"events_per_winter_{name.replace(' ', '_').lower()}"] = round(len(e2) / len(tab), 2)
        os.makedirs(f"{R}/{p}/figures", exist_ok=True)
        overview_figure(p, o, ev, tab, f"{R}/{p}/figures/overview.png")
        json.dump(dict(summary=s, per_winter=tab), open(f"{R}/{p}/summary.json", "w"), indent=1)
        data[p] = dict(o=o, ev=ev, tab=tab, s=s)
        print(p, s)
    paths = comparison_figures(data, f"{R}/comparison")
    html = build_html(data, paths)
    open(f"{R}/report_dunkelflaute.html", "w").write(html)
    print("wrote", f"{R}/report_dunkelflaute.html")


def build_html(data, paths):
    def row(label, key, fmt="{}"):
        return "<tr><th>%s</th>%s</tr>" % (label, "".join(
            f"<td>{fmt.format(d['s'][key]) if d['s'][key] is not None else '–'}</td>" for d in data.values()))
    head = "".join(f"<th><span class='dot' style='background:{PERIODS[p]['color']}'></span>{PERIODS[p]['short']}</th>" for p in data)
    table = f"""<table><thead><tr><th></th>{head}</tr></thead><tbody>
{row('Winters (Sep–Mar)', 'winters')}{row('Winter count', 'n_winters')}{row('Events, total', 'n_events')}
{row('<b>Events per winter</b>', 'events_per_winter')}{row('  spread between winters (SD)', 'events_per_winter_sd')}
{row('Winters with no event', 'winters_without_event')}
{row('Median duration (days)', 'median_dur_days')}{row('Longest event (days)', 'max_dur_days')}
{row('Time inside events (%)', 'flagged_pct')}
{row('Mean CF combined', 'mean_cf_combined')}{row('  solar', 'mean_cf_solar')}{row('  wind onshore', 'mean_cf_onshore')}{row('  wind offshore', 'mean_cf_offshore')}
{row('Trend in events per decade', 'events_trend_per_decade')}{row('  p-value', 'events_trend_p')}
</tbody></table>"""
    sections = ""
    for p, d in data.items():
        sections += (f"<h3>{PERIODS[p]['label']}</h3><img alt='{p} overview' src='data:image/png;base64,{b64(f'{R}/{p}/figures/overview.png')}'>")
    comp = "".join(f"<img alt='{k}' src='data:image/png;base64,{b64(v)}'>" for k, v in paths.items())
    era = data.get("validation_era5_2015_2026", {}).get("s", {})
    era_note = ""
    if era:
        era_note = (f"ERA5 check: reconstruction {era['events_per_winter']} events/winter vs real ERA5 "
                    f"{era.get('events_per_winter_era5_real', '–')} vs SMARD {era.get('events_per_winter_smard', '–')}.")
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Dunkelflaute Development</title><style>
:root{{--bg:#fcfcfb;--ink:#0b0b0b;--muted:#52514e;--line:#e4e3df}}
@media (prefers-color-scheme:dark){{:root:not([data-theme=light]){{--bg:#1a1a19;--ink:#fff;--muted:#c3c2b7;--line:#3a3a38}}}}
body{{background:var(--bg);color:var(--ink);font:15px/1.55 system-ui,sans-serif;margin:0 auto;max-width:980px;padding:24px 16px}}
h1{{font-size:26px;margin:0 0 4px}}h2{{margin-top:36px;border-top:1px solid var(--line);padding-top:18px}}h3{{margin:26px 0 8px}}
.sub{{color:var(--muted)}}img{{max-width:100%;background:#fff;border-radius:6px;margin:6px 0}}
table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{padding:5px 10px;border-bottom:1px solid var(--line);text-align:right}}
th:first-child{{text-align:left;font-weight:500;white-space:pre}}thead th{{font-weight:700}}.dot{{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:6px}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:8px}}.note{{background:var(--line);padding:10px 14px;border-radius:6px}}
</style></head><body>
<h1>How Dunkelflaute develops over time</h1>
<p class="sub">Germany combined solar + wind capacity factor, Mockert et al. (2023) definition: 48 h running mean &lt; 0.06, weights 44 % solar / 50 % onshore / 6 % offshore, Sep–Mar winters. ERA5 2015–2026 is the validation period; TCo1279-DART 1950C and 2080C are the reconstructed model runs.</p>
<h2>Summary</h2>{table}
<p class="sub">{era_note} Published Mockert reference: about 4 events per winter (1979–2018).</p>
<h2>Comparison</h2><div class="grid">{comp}</div>
<h2>Development within each period</h2>{sections}
<h2>How to read this</h2>
<div class="note"><b>Trends inside a run are not climate change.</b> 1950C and 2080C are fixed-forcing time slices: a trend across their winters is internal variability, and the p-values show that. The meaningful comparisons are between runs and against ERA5.<br><br>
<b>Caveats.</b> Solar is a reconstruction (model v5, exact monthly rescale) whose skill is validated on ERA5 and SMARD only; no ground truth exists for DART. The 2080C result assumes today's cloud-to-irradiance physics holds. Wind uses a 100 m log-law reconstruction. 1952 Feb/Mar solar inputs are missing from the archive but fall in the spin-up-excluded winter. Winters used: 1953–1969 (1950C), 2083–2092 (2080C), spin-up years excluded.</div>
</body></html>"""


if __name__ == "__main__":
    main()
