"""
Inventory of CMIP6 data on the DKRZ pool for the per-K Dunkelflaute pipeline: for every model with
daily sfcWind under SSP5-8.5, which members have day/sfcWind, day/rsds, Amon/tas in both historical
and ssp585, and 3hr/uas, 3hr/vas, 3hr/rsds (for the sub-daily subset).

Usage: python3 inventory.py
       python3 inventory.py --add-wind-only MODEL    add one model without day/rsds (wind-only checks,
                                                     e.g. EC-Earth3, the IFS-family CMIP6 model); the
                                                     entry gets wind_only: true and is kept out of the
                                                     20-model ensemble statistics
Output: dunkelflaute/data/cmip6/inventory.json  {model: {member, inst, paths: {exp: {table/var: dir}}}}
"""
import glob
import json
import os

P = "/pool/data/CMIP6/data"
OUT = "/work/ab0246/a270092/software/dark-doldrums-project/dunkelflaute/data/cmip6"
NEED = ["day/sfcWind", "day/rsds", "Amon/tas"]
SUB = ["3hr/uas", "3hr/vas", "3hr/rsds"]


def vdir(base, tv):
    t, v = tv.split("/")
    d = glob.glob(f"{base}/{t}/{v}/*/*")
    d = sorted(x for x in d if glob.glob(f"{x}/*.nc"))
    return d[-1] if d else None          # latest version


def add_wind_only(model):
    """Merge one model that lacks day/rsds into the existing inventory."""
    inv = json.load(open(f"{OUT}/inventory.json"))
    sdir = glob.glob(f"{P}/ScenarioMIP/*/{model}/ssp585")[0]
    hdir = glob.glob(f"{P}/CMIP/*/{model}/historical")[0]
    inst = sdir.split("/")[-3]
    for mem in sorted(os.listdir(sdir), key=lambda m: (not m.startswith("r1i1p1"), m)):
        paths = {e: {tv: vdir(b, tv) for tv in NEED + SUB} for e, b in (("historical", f"{hdir}/{mem}"), ("ssp585", f"{sdir}/{mem}"))}
        if all(paths[e][tv] for e in paths for tv in ("day/sfcWind", "Amon/tas")):
            inv[model] = dict(inst=inst, member=mem, subdaily=False, wind_only=True, paths=paths)
            json.dump(inv, open(f"{OUT}/inventory.json", "w"), indent=1)
            print(f"added {model} {mem} (wind only)")
            return
    raise SystemExit(f"{model}: no member with day/sfcWind and Amon/tas in both experiments")


def main():
    os.makedirs(OUT, exist_ok=True)
    inv = {}
    for sdir in sorted(glob.glob(f"{P}/ScenarioMIP/*/*/ssp585")):
        inst, model = sdir.split("/")[-3:-1]
        hdirs = glob.glob(f"{P}/CMIP/*/{model}/historical")
        if not hdirs:
            continue
        hdir = hdirs[0]
        for mem in sorted(os.listdir(sdir), key=lambda m: (not m.startswith("r1i1p1"), m)):
            if not os.path.isdir(f"{hdir}/{mem}"):
                continue
            paths = {"historical": {}, "ssp585": {}}
            ok = True
            for exp, base in (("historical", f"{hdir}/{mem}"), ("ssp585", f"{sdir}/{mem}")):
                for tv in NEED:
                    d = vdir(base, tv)
                    if d is None:
                        ok = False
                    paths[exp][tv] = d
                for tv in SUB:
                    paths[exp][tv] = vdir(base, tv)
            if ok:
                sub = all(paths[e][tv] for e in paths for tv in SUB)
                inv[model] = dict(inst=inst, member=mem, subdaily=sub, paths=paths)
                print(f"{model:18s} {mem:10s} grid {paths['ssp585']['day/sfcWind'].split('/')[-2]}  3hr {'yes' if sub else 'no'}")
                break
    json.dump(inv, open(f"{OUT}/inventory.json", "w"), indent=1)
    print(f"{len(inv)} models, {sum(v['subdaily'] for v in inv.values())} with 3-hourly wind and rsds")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2 and sys.argv[1] == "--add-wind-only":
        add_wind_only(sys.argv[2])
    else:
        main()
