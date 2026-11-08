"""
LUMIP land-use test (methods notes 13e follow-up): what is on /pool for the LUMIP experiments that
differ from a standard experiment only in land use, and for their paired controls.

Pairs (land-use effect = control minus LUMIP experiment):
  historical       - hist-noLu              historical land use minus land use fixed at 1850
  ssp370           - ssp370-ssp126Lu        SSP3 minus SSP1 land use under SSP3-7.0 forcing
  ssp126-ssp370Lu  - ssp126                 SSP3 minus SSP1 land use under SSP1-2.6 forcing
  esm-ssp585       - esm-ssp585-ssp126Lu    SSP5 minus SSP1 land use, emission driven SSP5-8.5
  deforest-globe   - piControl              idealised global deforestation (sensitivity only)
Per model, experiment and variable: members (latest version, preferred grid gn > gr > gr1), and the
first and last year over the files of each member.

Usage: python3 lumip_inventory.py        (any python3; reads directory listings only)
Output: dunkelflaute/results/regimes/lumip_inventory.json
"""
import json
import os
import re

POOL = "/pool/data/CMIP6/data"
_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
OUT = f"{_ROOT}/dunkelflaute/results/regimes/lumip_inventory.json"
VARS = ["day/sfcWind", "day/psl", "day/uas", "day/vas", "Amon/sfcWind", "Amon/psl", "Amon/ua", "Amon/va",
        "Amon/uas", "Amon/vas", "Amon/tas", "Lmon/cropFrac", "Lmon/treeFrac", "Lmon/grassFrac", "Lmon/pastureFrac",
        "Lmon/lai", "Eyr/cropFrac", "Eyr/treeFrac", "Eyr/grassFrac", "fx/sftlf"]
EXPS = {"hist-noLu": "LUMIP", "ssp370-ssp126Lu": "LUMIP", "ssp126-ssp370Lu": "LUMIP", "esm-ssp585-ssp126Lu": "LUMIP",
        "deforest-globe": "LUMIP", "historical": "CMIP", "piControl": "CMIP", "ssp370": "ScenarioMIP",
        "ssp126": "ScenarioMIP", "esm-ssp585": "C4MIP"}
GRID_PREF = ["gn", "gr", "gr1", "gr2"]


def ls(p):
    try:
        return sorted(os.listdir(p))
    except (FileNotFoundError, NotADirectoryError, PermissionError):
        return []


def latest(member_dir, var):
    base = f"{member_dir}/{var}"
    grids = ls(base)
    if not grids:
        return None
    g = sorted(grids, key=lambda x: GRID_PREF.index(x) if x in GRID_PREF else 9)[0]
    vs = [v for v in ls(f"{base}/{g}") if v.startswith("v")]
    if not vs:
        return None
    return f"{base}/{g}/{sorted(vs)[-1]}"


def years(d):
    ys = []
    for f in ls(d):
        m = re.search(r"_(\d{4})\d*-(\d{4})\d*\.nc$", f)
        if m:
            ys += [int(m.group(1)), int(m.group(2))]
        elif f.endswith(".nc"):
            ys += [-1]
    return (min(ys), max(ys)) if ys else None


def main():
    models = {}
    for inst in ls(f"{POOL}/LUMIP"):
        for model in ls(f"{POOL}/LUMIP/{inst}"):
            models[model] = inst
    inv = {}
    for model, inst in models.items():
        inv[model] = {"inst": inst}
        for exp, mip in EXPS.items():
            ed = f"{POOL}/{mip}/{inst}/{model}/{exp}"
            mems = ls(ed)
            if not mems:
                continue
            e = {}
            for var in VARS:
                per = {}
                for m in mems:
                    d = latest(f"{ed}/{m}", var)
                    if d:
                        per[m] = dict(dir=d, years=years(d))
                if per:
                    e[var] = per
            inv[model][exp] = e
    json.dump(inv, open(OUT, "w"), indent=1)
    short = {"day/sfcWind": "dW", "day/psl": "dP", "day/uas": "du", "Amon/sfcWind": "mW", "Amon/psl": "mP",
             "Amon/ua": "mU", "Amon/va": "mV", "Amon/tas": "mT", "Lmon/cropFrac": "Cr", "Lmon/treeFrac": "Tr",
             "Lmon/lai": "La", "Eyr/cropFrac": "yC", "Eyr/treeFrac": "yT"}
    for model, ex in inv.items():
        print(f"\n{model} ({ex['inst']})")
        for exp, e in ex.items():
            if exp == "inst":
                continue
            parts = []
            for v, s in short.items():
                if v in e:
                    yy = [x["years"] for x in e[v].values() if x["years"]]
                    y0 = min(a for a, b in yy) if yy else -1
                    y1 = max(b for a, b in yy) if yy else -1
                    parts.append(f"{s}:{len(e[v])}({y0}-{y1})")
            print(f"  {exp:20s} " + " ".join(parts))
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
