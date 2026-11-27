"""
Data citations for every CMIP6 dataset (model x experiment) the paper read, for the Supplementary table of models
(Jan 2026-11-28: model list with dataset DOIs in the Supplement, as NCC papers with CMIP6 ensembles do).

  input   dunkelflaute/results/cmip6/used_files/{A..F}_*.txt (DRS paths of the files actually read, from used_files.py)
  lookup  DataCite REST API, exact title "<institution> <source> model output prepared for CMIP6 <activity> <experiment>",
          DOI prefix 10.22033 (ESGF/WDCC CMIP6 citation service)
  output  dunkelflaute/results/cmip6/data_citations.json            per dataset: doi, title, year, creators
          paper/cmip6_data.bib                                       one @misc per dataset (key CMIP6<source><experiment>)
          paper/tables/tab_cmip6_models.tex                          Supplementary table, one row per model

Usage: python3 data_citations.py          (any python3 with internet access; about 2 min)
"""
import json
import os
import re
import time
import unicodedata
import urllib.parse
import urllib.request
from collections import OrderedDict, defaultdict

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
UF = f"{_ROOT}/dunkelflaute/results/cmip6/used_files"
OUT_JSON = f"{_ROOT}/dunkelflaute/results/cmip6/data_citations.json"
OUT_BIB = f"{_ROOT}/paper/cmip6_data.bib"
OUT_TAB = f"{_ROOT}/paper/tables/tab_cmip6_models.tex"
LISTS = OrderedDict([
    ("A", ("A_daily_wind_solar_tas.txt", "daily")),
    ("C", ("C_psl.txt", "G")),
    ("B", ("B_3hr.txt", "3h")),
    ("D", ("D_hubheight.txt", "ML")),
    ("E", ("E_stability.txt", "BL")),
    ("F", ("F_lumip.txt", "LU")),
])
# Datasets whose DataCite title does not follow the usual pattern, checked by resolving the DOI (2026-11-28).
MANUAL = {
    "CMIP.NOAA-GFDL.GFDL-CM4.historical": dict(
        title="NOAA-GFDL GFDL-CM4 model output historical", doi="10.22033/ESGF/CMIP6.8594", year=2018,
        creators=["Guo, Huan", "John, Jasmin G", "Blanton, Chris", "McHugh, Colleen", "Nikonov, Serguei",
                  "Radhakrishnan, Aparna", "Rand, Kristopher", "Zadeh, Niki T.", "Balaji, V", "Durachta, Jeff"]),
}
EXP_ORDER = ["historical", "ssp585", "ssp126", "ssp370", "esm-ssp585", "hist-noLu", "ssp370-ssp126Lu",
             "ssp126-ssp370Lu", "esm-ssp585-ssp126Lu"]


_ACC = {"\u0308": '"', "\u0301": "'", "\u0300": "`", "\u0302": "^", "\u0303": "~", "\u0327": "c ", "\u030c": "v ",
        "\u0306": "u ", "\u030a": "r ", "\u0328": "k ", "\u0304": "="}
_SPECIAL = {"\u00f8": "{\\o}", "\u00d8": "{\\O}", "\u00df": "{\\ss}", "\u00e6": "{\\ae}", "\u00c6": "{\\AE}", "\u0142": "{\\l}",
            "\u0141": "{\\L}", "\u0131": "{\\i}"}


def tex(name):
    """Non-ASCII letters to LaTeX accents, so BibTeX can form initials (e.g. Gürses, Özgür)."""
    out = []
    for ch in name:
        if ord(ch) < 128:
            out.append(ch)
        elif ch in _SPECIAL:
            out.append(_SPECIAL[ch])
        else:
            d = unicodedata.normalize("NFD", ch)
            base, marks = d[0], d[1:]
            if len(marks) == 1 and marks in _ACC:
                out.append("{\\" + _ACC[marks] + base + "}")
            else:
                out.append(ch)
    return "".join(out)


def key(src, exp):
    return "CMIP6" + re.sub(r"[^A-Za-z0-9]", "", src) + re.sub(r"[^A-Za-z0-9]", "", exp)


def lookup(act, inst, src, exp):
    title = f"{inst} {src} model output prepared for CMIP6 {act} {exp}"
    q = f'titles.title:"{title}"'
    u = "https://api.datacite.org/dois?" + urllib.parse.urlencode({"query": q, "page[size]": 25})
    for attempt in range(4):
        try:
            d = json.load(urllib.request.urlopen(u, timeout=60))
            break
        except Exception as e:                      # noqa: BLE001
            if attempt == 3:
                return dict(title=title, doi=None, error=str(e))
            time.sleep(3)
    for x in d["data"]:
        a = x["attributes"]
        t = a["titles"][0]["title"].strip()
        if t == title and x["id"].lower().startswith("10.22033/"):
            return dict(title=t, doi=x["id"].replace("esgf/cmip6", "ESGF/CMIP6"), year=a.get("publicationYear"),
                        creators=[c.get("name") for c in a.get("creators", [])])
    # Fallback: some dataset titles spell the source differently from its DRS name (UKESM1.0-LL, MPIESM1.2-LR,
    # MRI-ESM2.0, KACE1.0-G, GISS-E2.1G). Search all titles of the experiment and match the normalised source.
    norm = lambda z: re.sub(r"[^a-z0-9]", "", z.lower())        # noqa: E731
    tail = f"model output prepared for CMIP6 {act} {exp}"
    q = f'titles.title:"{tail}" AND titles.title:"{inst}"'
    u = "https://api.datacite.org/dois?" + urllib.parse.urlencode({"query": q, "page[size]": 200})
    try:
        d = json.load(urllib.request.urlopen(u, timeout=60))
    except Exception as e:                          # noqa: BLE001
        return dict(title=title, doi=None, error=str(e))
    for x in d["data"]:
        a = x["attributes"]
        t = a["titles"][0]["title"].strip()
        if not (t.endswith(tail) and t.startswith(inst + " ") and x["id"].lower().startswith("10.22033/")):
            continue
        if norm(t[len(inst) + 1:-len(tail)]) == norm(src):
            return dict(title=t, doi=x["id"].replace("esgf/cmip6", "ESGF/CMIP6"), year=a.get("publicationYear"),
                        creators=[c.get("name") for c in a.get("creators", [])])
    return dict(title=title, doi=None, error="no 10.22033 match")


def main():
    use = defaultdict(set)            # (act, inst, src, exp) -> analyses
    members = defaultdict(set)        # src -> members
    inst_of = {}
    for tag, (fn, lab) in LISTS.items():
        for line in open(f"{UF}/{fn}"):
            p = line.strip().split("/")
            if len(p) < 5:
                continue
            act, inst, src, exp, mem = p[:5]
            use[(act, inst, src, exp)].add(lab)
            if tag != "F":                          # LUMIP averages all members of each pair; not shown per model
                members[src].add(mem)
            inst_of.setdefault(src, set()).add(inst)
    cache = json.load(open(OUT_JSON)) if os.path.exists(OUT_JSON) else {}
    res = {}
    for (act, inst, src, exp) in sorted(use):
        k = f"{act}.{inst}.{src}.{exp}"
        r = cache.get(k)
        if not r or not r.get("doi"):
            r = lookup(act, inst, src, exp)
            time.sleep(0.3)
        if not r.get("doi") and k in MANUAL:              # resolved by hand on api.datacite.org/dois/<doi>
            r = dict(MANUAL[k])
        r["analyses"] = sorted(use[(act, inst, src, exp)])
        r["key"] = key(src, exp)
        res[k] = r
        print(f"{k:60s} {r.get('doi')}", flush=True)
    json.dump(res, open(OUT_JSON, "w"), indent=1)

    # bib
    with open(OUT_BIB, "w") as f:
        f.write("% CMIP6 data citations, generated by dunkelflaute/scripts/cmip6/data_citations.py; do not edit by hand\n")
        for k, r in res.items():
            if not r.get("doi"):
                continue
            cr = r["creators"]
            auth = " and ".join(tex(c) for c in cr[:5]) + (" and others" if len(cr) > 5 else "")
            f.write(f"@misc{{{r['key']},\n  author = {{{auth}}},\n  title = {{{{{r['title']}}}}},\n"
                    f"  howpublished = {{Earth System Grid Federation}},\n  year = {{{r['year']}}},\n  doi = {{{r['doi']}}},\n"
                    f"  url = {{https://doi.org/{r['doi']}}}\n}}\n")

    # table: one row per model
    flags = OrderedDict([("daily", "D"), ("G", "G"), ("3h", "3"), ("ML", "M"), ("BL", "B"), ("LU", "L")])
    rows = defaultdict(lambda: dict(exps=[], ana=set()))
    for k, r in res.items():
        act, inst, src, exp = k.split(".")
        rows[src]["exps"].append((EXP_ORDER.index(exp) if exp in EXP_ORDER else 99, exp, r))
        rows[src]["ana"] |= set(r["analyses"])
    lines = []
    for src in sorted(rows, key=str.lower):
        ex = sorted(rows[src]["exps"])
        cites = ", ".join(f"{exp}\\citep{{{r['key']}}}" if r.get("doi") else f"{exp} [no DOI]" for _, exp, r in ex)
        ms = sorted(members[src])
        mem = ", ".join(ms) if ms else "--"
        an = "".join(v for a, v in flags.items() if a in rows[src]["ana"])
        lines.append(f"{src} & {mem} & {an} & {cites} \\\\")
    os.makedirs(os.path.dirname(OUT_TAB), exist_ok=True)
    with open(OUT_TAB, "w") as f:
        f.write("%% Generated by dunkelflaute/scripts/cmip6/data_citations.py from the used-file lists; do not edit by hand.\n")
        f.write("\\begin{table*}[p]\n\\caption{CMIP6 models and datasets used. Analyses: "
                "D daily events (wind, solar, and global temperature), G geostrophic speed and regimes (daily sea-level pressure), "
                "3 three-hourly wind and radiation (footing check), M wind at 1000 and 925\\,hPa or on model levels near 100\\,m, "
                "B monthly boundary-layer and land-cover diagnostics, "
                "L LUMIP land-use experiments. The institution of each model is named in its data citations. Member: the ensemble member of all analyses "
                "except LUMIP, a dash for models used only there; for the LUMIP experiments we averaged over all available members of each experiment pair "
                "(Supplementary Note~\\ref{app:lumip}).}\n"
                "\\label{tab:cmip6models}\n\\scriptsize\n\\setlength{\\tabcolsep}{2pt}\\renewcommand{\\arraystretch}{0.85}\n"
                "\\begin{tabular}{@{}l l l >{\\raggedright\\arraybackslash}p{0.535\\textwidth}@{}}\n\\tophline\n"
                "Model & Member & Analyses & Experiments and data citations \\\\\n\\middlehline\n")
        f.write("\n".join(lines) + "\n\\bottomhline\n\\end{tabular}\n\\end{table*}\n")
    n_ok = sum(1 for r in res.values() if r.get("doi"))
    print(f"{n_ok} of {len(res)} datasets with DOI; {len(rows)} models; wrote {OUT_BIB}, {OUT_TAB}")


if __name__ == "__main__":
    main()
