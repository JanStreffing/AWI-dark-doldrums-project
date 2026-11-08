"""Visual report (self-contained HTML) for the EPSIS solar/wind installed-capacity data.
Run plot_epsis_wind_solar_capacity.py first; writes reports/epsis_wind_solar_capacity_report.html."""
import base64
import pandas as pd

HERE = "/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/epsis_validation"

d = pd.read_csv(f"{HERE}/data/generation_capacity/HOME>GenerationCapacity>ByFuel.csv", skiprows=2, header=None)
d = pd.DataFrame({"period": d[0], "solar": d[10] / 1000, "wind": d[11] / 1000})
d["year"] = d.period.str[:4].astype(int)
last = d.sort_values("period").groupby("year").tail(1).set_index("year")
added = last[["solar", "wind"]].diff()
latest = last.period.iloc[-1]
y0, y1 = last.index[0], last.index[-1]


def img(name):
    with open(f"{HERE}/figures/{name}", "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode()


def cagr(col):
    full = last.iloc[:-1][col]  # complete years only (drop partial latest year)
    n = full.index[-1] - full.index[0]
    return ((full.iloc[-1] / full.iloc[0]) ** (1 / n) - 1) * 100


rows = "".join(
    f"<tr><td>{y}{'*' if last.period[y][5:] != '12' else ''}</td>"
    f"<td>{last.solar[y]:.2f}</td><td>{'' if pd.isna(added.solar[y]) else f'{added.solar[y]:+.2f}'}</td>"
    f"<td>{last.wind[y]:.2f}</td><td>{'' if pd.isna(added.wind[y]) else f'{added.wind[y]:+.2f}'}</td></tr>"
    for y in last.index
)
peak_s = added.solar.idxmax(); peak_w = added.wind.idxmax()

html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>EPSIS solar and wind capacity</title>
<style>
body{{font:15px/1.55 system-ui,sans-serif;color:#0b0b0b;background:#fcfcfb;max-width:980px;margin:0 auto;padding:24px 16px}}
h1{{font-size:24px;margin:0 0 4px}} h2{{font-size:18px;margin-top:32px}} .sub{{color:#52514e}}
.tiles{{display:flex;gap:12px;flex-wrap:wrap;margin:20px 0}}
.tile{{flex:1 1 200px;border:1px solid #e1e0d9;border-radius:8px;padding:12px 14px;background:#fff}}
.tile b{{display:block;font-size:26px}} .tile span{{color:#52514e;font-size:13px}}
img{{width:100%;height:auto;border:1px solid #e1e0d9;border-radius:8px;background:#fff}}
table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}}
th,td{{padding:5px 10px;text-align:right;border-bottom:1px solid #e1e0d9}} th:first-child,td:first-child{{text-align:left}}
th{{color:#52514e;font-weight:600}} .note{{color:#52514e;font-size:13px}}
</style></head><body>
<h1>South Korea: installed solar and wind capacity</h1>
<p class="sub">EPSIS, Generation Capacity by Fuel, monthly, {y0}–{latest}. Year-end (December) values; the latest year uses {latest}.</p>
<div class="tiles">
<div class="tile"><span>Solar, {latest}</span><b>{last.solar.iloc[-1]:.1f} GW</b><span>{last.solar.iloc[0]:.1f} GW in {y0}</span></div>
<div class="tile"><span>Wind, {latest}</span><b>{last.wind.iloc[-1]:.2f} GW</b><span>{last.wind.iloc[0]:.2f} GW in {y0}</span></div>
<div class="tile"><span>Peak solar addition</span><b>{added.solar.max():.1f} GW</b><span>in {peak_s}</span></div>
<div class="tile"><span>Peak wind addition</span><b>{added.wind.max():.2f} GW</b><span>in {peak_w}</span></div>
<div class="tile"><span>Annual growth (complete years)</span><b>{cagr('solar'):.0f}% / {cagr('wind'):.0f}%</b><span>solar / wind, compound</span></div>
</div>
<h2>Installed capacity</h2>
<img src="{img('epsis_wind_solar_capacity.png')}" alt="Bar charts of installed solar and wind capacity per year">
<h2>Capacity added per year</h2>
<img src="{img('epsis_wind_solar_capacity_added.png')}" alt="Bar charts of solar and wind capacity added per year">
<h2>Data table (GW)</h2>
<table><tr><th>Year</th><th>Solar installed</th><th>Solar added</th><th>Wind installed</th><th>Wind added</th></tr>{rows}</table>
<p class="note">* Latest available month ({latest}), not year-end; its addition covers Dec {y1-1} to {latest} only.
Added capacity is the change in year-end installed capacity. Solar and wind panels use separate y-axes.</p>
</body></html>"""

out = f"{HERE}/reports/epsis_wind_solar_capacity_report.html"
with open(out, "w") as f:
    f.write(html)
print(out)
