"""
Winter-bootstrap ranges of the 24-region field-mean event ratio for every comparison on the daily
footing (20 CMIP6 models, EC-Earth3 wind only, and the EERIE / TCO95 daily bridges), so that the CMIP6
dots in the per-K figure carry the same kind of within-run range as our runs (region_per_K.py).

Reuses dunkelflaute/scripts/cmip6/daily_events.py unchanged: the same comparisons, winters, common
region sets, and the matched thresholds stored in results/cmip6/daily_events_summary.json (not
re-searched). Per region, events per control and future winter at that threshold; then 5000
bootstrap resamples of whole winters, the SAME winter indices in every region (keeps the spatial
correlation), field mean of ln(ratio) over the common regions, 5-95 %.
The point estimate is checked against daily_events_summary.json.

Usage: python3 fieldsig_cmip6_bootstrap.py      (module python3; daily region files, about 1 min)
Output: dunkelflaute/results/maps/fieldsig_cmip6_bootstrap.json
"""
import json
import os
import sys

import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/cmip6")
import daily_events as de                                              # noqa: E402

OUT = f"{_ROOT}/dunkelflaute/results/maps"
NBOOT = 5000
RNG = np.random.default_rng(3)


def main():
    S = json.load(open(f"{_ROOT}/dunkelflaute/results/cmip6/daily_events_summary.json"))
    C = de.comparisons()
    out = {}
    for name, (fc, ec, c0, c1, ff, ef, f0, f1, Tc, Tf, group) in C.items():
        if name not in S["comparisons"]:
            continue
        dc, df = de.load(fc, ec), de.load(ff, ef)
        dT = S["comparisons"][name]["dT"]
        row = dict(group=group, dT=dT)
        for kind in de.kinds_of(group):
            V = S["common_regions"][kind]
            A, B = [], []
            for r in V:
                i = de.R.index(r)
                thr = np.array([S["regions"][name]["metrics"][kind][r]["threshold"]])
                A.append(de.counts(dc, kind, i, c0, c1, thr)[:, 0])
                B.append(de.counts(df, kind, i, f0, f1, thr)[:, 0])
            a = np.stack(A, axis=1).astype(float)        # (Wc, R)
            b = np.stack(B, axis=1).astype(float)        # (Wf, R)
            with np.errstate(divide="ignore", invalid="ignore"):
                lr = np.log(b.mean(0) / a.mean(0))
            obs = float(np.mean(lr[np.isfinite(lr)]))
            ref = np.log(S["comparisons"][name][kind]["mean_ratio"])
            boots = np.empty(NBOOT)
            for k in range(NBOOT):
                ia = RNG.integers(0, len(a), len(a))
                ib = RNG.integers(0, len(b), len(b))
                with np.errstate(divide="ignore", invalid="ignore"):
                    x = np.log(b[ib].mean(0) / a[ia].mean(0))
                boots[k] = np.mean(x[np.isfinite(x)])
            lo, hi = np.percentile(boots, [5, 95])
            pk = lambda x: float(100 * (np.exp(x / dT) - 1))        # noqa: E731
            row[kind] = dict(n_regions=len(V), n_winters=[len(a), len(b)], lnratio=obs, lnratio_p05_p95=[float(lo), float(hi)],
                             pct_per_K=pk(obs), pct_per_K_p05_p95=[pk(lo), pk(hi)], check_vs_summary=float(obs - ref))
            print(f"  {name:40s} {kind:4s} dT {dT:.2f}  ratio {np.exp(obs):.3f} [{np.exp(lo):.3f}, {np.exp(hi):.3f}]  "
                  f"{pk(obs):+.1f} %/K [{pk(lo):+.1f}, {pk(hi):+.1f}]  check {obs - ref:+.1e}", flush=True)
        out[name] = row
    os.makedirs(OUT, exist_ok=True)
    json.dump(out, open(f"{OUT}/fieldsig_cmip6_bootstrap.json", "w"), indent=1)
    print(f"Saved {OUT}/fieldsig_cmip6_bootstrap.json")


if __name__ == "__main__":
    main()
