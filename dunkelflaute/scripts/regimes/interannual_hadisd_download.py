"""
Download the HadISD (v3.4.3.2025f) station files needed for the observed interannual transfer test.

Selection from the station full-info list (hadisd_station_fullinfo_v343_2025f.txt):
  inside one of the 8 European regions (regime_analysis.EU; DE = the Germany land polygon, the others
  the region boxes of regions.py), record starting by 1994-09-01 and running to at least 2024-04-01
  (so that 30 Sep-Mar winters up to 2025 are possible), station elevation below 800 m (mountain
  summits sample the free troposphere, not the boundary layer).
File names are taken from the per-station download pages (station_download_<k>.html).

Raw files stay gzipped in dunkelflaute/data/obs/hadisd_raw/ (gitignored); nothing in them is run.

Usage: python3 -I interannual_hadisd_download.py      (login node: needs internet)
"""
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/regions")
import regions as rg                        # noqa: E402

VER = "v343_2025f"
BASE = f"https://www.metoffice.gov.uk/hadobs/hadisd/{VER}"
RAW = f"{_ROOT}/dunkelflaute/data/obs/hadisd_raw"
EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
START_BY, END_FROM, MAX_ELEV = "1994-09-01", "2024-04-01", 800.0


def station_table():
    rows = [l.split() for l in open(f"{RAW}/hadisd_station_fullinfo_{VER}.txt")]
    ids = [r[0] for r in rows]
    lat = np.array([float(r[-5]) for r in rows])
    lon = np.array([float(r[-4]) for r in rows])
    el = np.array([float(r[-3]) for r in rows])
    s = np.array([r[-2] for r in rows])
    e = np.array([r[-1] for r in rows])
    return ids, lat, lon, el, s, e


def select():
    ids, lat, lon, el, s, e = station_table()
    ok = (s <= START_BY) & (e >= END_FROM) & (el < MAX_ELEV)
    cells = rg.region_cells(lat, lon, np.ones(len(lat), bool))
    sel = {r: [ids[i] for i in cells[r] if ok[i]] for r in EU}
    return sel, dict(zip(ids, zip(lat.tolist(), lon.tolist(), el.tolist(), s.tolist(), e.tolist())))


def file_names():
    names = {}
    for k in range(10):
        html = subprocess.run(["curl", "-s", "-m", "120", f"{BASE}/station_download_{k}.html"],
                              capture_output=True, text=True, check=True).stdout
        for f in re.findall(r'href="data/(hadisd\.[^"]+\.nc\.gz)"', html):
            names[f.split("_")[-1].replace(".nc.gz", "")] = f
    return names


def fetch(fn):
    out = f"{RAW}/{fn}"
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return fn, "exists"
    r = subprocess.run(["curl", "-s", "-f", "-m", "600", "--retry", "3", "-o", out + ".part", f"{BASE}/data/{fn}"])
    if r.returncode == 0:
        os.rename(out + ".part", out)
        return fn, "ok"
    return fn, f"failed {r.returncode}"


def main():
    sel, info = select()
    names = file_names()
    want = sorted({i for v in sel.values() for i in v})
    missing = [i for i in want if i not in names]
    json.dump(dict(version=VER, criteria=dict(start_by=START_BY, end_from=END_FROM, max_elev_m=MAX_ELEV),
                   regions=sel, info={i: info[i] for i in want}, files={i: names[i] for i in want if i in names},
                   missing=missing), open(f"{_ROOT}/dunkelflaute/data/obs/hadisd_selection.json", "w"), indent=1)
    print(f"{len(want)} stations ({', '.join(f'{r} {len(v)}' for r, v in sel.items())}); {len(missing)} without a file", flush=True)
    with ThreadPoolExecutor(4) as ex:
        for n, (fn, st) in enumerate(ex.map(fetch, [names[i] for i in want if i in names])):
            if st != "exists" and (n % 50 == 0 or st.startswith("failed")):
                print(f"  {n} {fn} {st}", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
