"""List the CMIP6 files (not directories) that the extraction scripts read, by category, as NUL-separated lists
for dust --files0-from. Selection: the variable/table/experiment combinations and year ranges of
extract_daily.py, extract_3hr.py, extract_psl.py, interannual_extract_cmip6.py, hubheight_extract.py,
stability_cmip6_extract.py, lumip_extract_daily.py, lumip_extract_monthly.py; files kept if their date range
overlaps the years read. Members and versions: main ensemble from data/cmip6/inventory.json (latest version
directory where a variable is not in the inventory), LUMIP from results/regimes/lumip_inventory.json (all members
listed there, so an upper bound)."""
import glob
import json
import os
import re
import sys

R = "/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute"
OUT = sys.argv[1]
inv = json.load(open(f"{R}/data/cmip6/inventory.json"))
linv = json.load(open(f"{R}/results/regimes/lumip_inventory.json"))
DATE = re.compile(r"_(\d{4})\d{0,8}-(\d{4})\d{0,8}\.nc$")


def files(d, y0, y1):
    out = []
    for f in sorted(glob.glob(f"{d}/*.nc")):
        m = DATE.search(f)
        if not m:                                   # fx and other undated files
            out.append(f)
        elif int(m.group(1)) <= y1 and int(m.group(2)) >= y0:
            out.append(f)
    return out


def vardir(base, table, var):
    """base: a known variable directory .../<member>/<table>/<var>/<grid>/<version>; same member, other table/var."""
    mem = base.split("/")[:-4]
    cands = sorted(glob.glob("/".join(mem + [table, var, "*", "v*"])))
    return cands[-1] if cands else None


cats = {k: [] for k in ("A_daily_wind_solar_tas", "B_3hr", "C_psl", "D_hubheight", "E_stability", "F_lumip")}
for m, v in inv.items():
    P = v.get("paths", {})
    for exp, yrs in (("historical", (1950, 2014)), ("ssp585", (2040, 2100))):
        pp = P.get(exp, {})
        base = pp.get("day/sfcWind")
        if not base:
            continue
        hist = exp == "historical"
        y_main = (1950, 1980) if hist else (2040, 2100)
        for key in ("day/sfcWind", "day/rsds", "Amon/tas"):
            if pp.get(key):
                cats["A_daily_wind_solar_tas"] += files(pp[key], *(yrs if (hist and key == "day/sfcWind") else y_main))
        for key in ("3hr/uas", "3hr/vas", "3hr/rsds"):
            if pp.get(key):
                cats["B_3hr"] += files(pp[key], *((1960, 1990) if hist else (2070, 2100)))
        d = vardir(base, "day", "psl")
        if d:
            cats["C_psl"] += files(d, *yrs)
        yh = (1950, 1980) if hist else (2070, 2100)
        for table, var in (("day", "ua"), ("day", "va"), ("day", "zg"), ("day", "uas"), ("day", "vas"),
                           ("Eday", "ua"), ("Eday", "va"), ("CFday", "ua"), ("CFday", "va"), ("CFday", "ps"), ("fx", "orog")):
            d = vardir(base, table, var)
            if d:
                cats["D_hubheight"] += files(d, *yh)
        for table, var in (("Amon", "ts"), ("Amon", "ps"), ("Amon", "hfss"), ("Amon", "hfls"), ("Amon", "sfcWind"),
                           ("Amon", "ta"), ("LImon", "snc"), ("LImon", "snw"), ("Lmon", "lai"), ("Lmon", "treeFrac"),
                           ("Lmon", "cropFrac"), ("Lmon", "grassFrac"), ("fx", "sftlf")):
            d = vardir(base, table, var)
            if d:
                cats["E_stability"] += files(d, *yh)

LUMIP_EXPS = {"hist-noLu": (1950, 2014), "historical": (1950, 2014), "ssp370": (2015, 2100), "ssp126": (2015, 2100),
              "ssp370-ssp126Lu": (2015, 2100), "ssp126-ssp370Lu": (2015, 2100), "esm-ssp585": (2015, 2100),
              "esm-ssp585-ssp126Lu": (2015, 2100)}
LUMIP_VARS = {"day/sfcWind", "day/psl", "Amon/sfcWind", "Amon/tas", "Amon/psl", "Lmon/cropFrac", "Lmon/treeFrac",
              "Lmon/grassFrac", "Lmon/pastureFrac", "Lmon/lai", "Eyr/cropFrac", "Eyr/treeFrac"}
seen = set()
for m, v in linv.items():
    for exp, yrs in LUMIP_EXPS.items():
        for key, mems in (v.get(exp) or {}).items():
            if key not in LUMIP_VARS or not isinstance(mems, dict):
                continue
            for mem, info in mems.items():
                d = info.get("dir") if isinstance(info, dict) else None
                if d and d not in seen:
                    seen.add(d)
                    y = yrs if exp not in ("hist-noLu", "historical") or key.startswith(("Amon", "Lmon", "Eyr")) else (1985, 2014)
                    cats["F_lumip"] += files(d, *y)

for k, L in cats.items():
    L = sorted(set(L))
    with open(f"{OUT}/{k}.files0", "w") as fh:
        fh.write("\0".join(L))
    print(k, len(L), "files")
allf = sorted(set(f for L in cats.values() for f in L))
with open(f"{OUT}/ALL.files0", "w") as fh:
    fh.write("\0".join(allf))
print("ALL", len(allf), "files")
