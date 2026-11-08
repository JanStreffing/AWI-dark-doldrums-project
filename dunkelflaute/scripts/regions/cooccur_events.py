"""
Europe-wide co-occurrence of winter low-wind events: do low-wind spells hit more of the 8 European
regions at once in a warmer climate?

Regions: DE, FR, IB, UK, PL, SC, IT, BK (regions.py; the same series as region_events.py and
daily_events.py). Wind only (onshore CF of the log-law 100 m wind).

Events per region exactly as in the per-region analysis: trailing 48 h mean, Mockert window
expansion, threshold matched so that the comparison's control has TARGET = 5 events per winter
(lowest level of a 400-step grid between the 0.1st and 40th percentile of the whole control file
reaching it, exactly as region_events.matched_threshold and daily_events.matched). A time step is
"in an event" for a region if it belongs to one of that region's events. N(t) = number of the 8
regions in an event at step t. Statistic: P(N >= k) for k = 1..8, the share of winter time steps
with at least k regions in an event at once (mean over winters; x 212 = days per winter).

Future variants:
  all    future at the control's thresholds (the total change)
  mr     mean removed: each region's future series rescaled to the control's Sep-Mar mean
         (x * mean_ctl / mean_fut), same thresholds (as region_per_K.py)
  dep    dependence only: thresholds re-matched in the future so every region again has 5 events
         per winter; a change in P(N >= k) is then a change in how often regions are low together
Also, per period, P(N >= k) for independent regions with the same marginal shares (Poisson
binomial), and the mean pairwise Jaccard index of the event masks (share of time both regions are in
an event among the time either is).

Footings and comparisons:
  6h    our runs (region_events.COMPARISONS: 9km, 31km, 31km_9kmspans, 100km, 100km_9kmspans,
        eerie_ctl, eerie_hist); 48 h = 8 steps
  day   20 CMIP6 models (+ EC-Earth3 wind only, outside the ensemble), historical 1951-1980 against
        ssp585 2071-2100; bridges EERIE ctl/hist -> ssp245 2041-2070, TCO95 ctl -> hi 2081-2099
        (daily_events.comparisons); 48 h = 2 days
Per K: 100 (exp(ln(P_fut / P_ctl) / dT) - 1), dT as in region_per_K.py / daily_events.py.

Null: the TCO95 constant-1950 control (147 winters) on the same footing. NPL random placements of two
disjoint blocks with the comparison's own winter counts (control block cut to W - n_b - 40 when the
two do not fit, as grid_per_K.py; shorter blocks widen the range); in each placement the first block
plays the control (thresholds matched in it) and the second the future, for all three variants.
5-95 % of ln ratio.

Validation: ERA5 winters 1980-2026 (cooccur_extract_era5.py) on both footings, thresholds matched in
ERA5 itself; P(N >= k), the independence ratio and the Jaccard index against every model control.

Usage: python3 cooccur_events.py      (module python3; about 15 min with 32 workers, compute node)
Output: dunkelflaute/results/regions/cooccur_events_summary.json, cooccur_events.npz
"""
import json
import os
import sys
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

_ROOT = "/work/ab0246/a270092/software/dark-doldrums-project"
for p in ("regions", "analysis", "cmip6"):
    sys.path.insert(0, f"{_ROOT}/dunkelflaute/scripts/{p}")
import regions as rg                                          # noqa: E402
import region_events as ev                                    # noqa: E402
import region_per_K as rpk                                    # noqa: E402
import daily_events as de                                     # noqa: E402
from robustness_mockert_events import trailing_mean           # noqa: E402

OUT = f"{_ROOT}/dunkelflaute/results/regions"
ERA5 = f"{_ROOT}/dunkelflaute/data/regions/cooccur_era5_wind.npz"
EU = ["DE", "FR", "IB", "UK", "PL", "SC", "IT", "BK"]
EI = [rg.REGIONS.index(r) for r in EU]
K = np.arange(1, len(EU) + 1)
TARGET = 5.0
NGRID = 400
NPL = int(os.environ.get("NPL", 400))
RNG = np.random.default_rng(7)
WINDAYS = 212
OURS6 = ["9km", "31km_9kmspans", "31km", "100km_9kmspans", "100km", "eerie_ctl", "eerie_hist"]


# ---------------------------------------------------------------- events
def flags(x, thr, win):
    """Steps of x (T,) that belong to an event at threshold thr (Mockert window expansion)."""
    below = np.nan_to_num(trailing_mean(x, win), nan=np.inf) < thr
    c = np.concatenate([[0], np.cumsum(below)])
    n = len(x)
    hi = np.minimum(np.arange(n) + win, n)
    return (c[hi] - c[:n]) > 0


def counts(x, thr, win):
    """Events of x (T,) for each threshold in thr (G,) (as region_events.counts_grid, explicit win)."""
    roll = trailing_mean(x, win)
    below = (np.nan_to_num(roll, nan=np.inf)[None, :] < thr[:, None]).astype(np.int32)
    c = np.concatenate([np.zeros((len(thr), 1), np.int32), np.cumsum(below, axis=1)], axis=1)
    n = below.shape[1]
    hi = np.minimum(np.arange(n) + win, n)
    f = (c[:, hi] - c[:, :n]) > 0
    return f[:, 0].astype(int) + (f[:, 1:] & ~f[:, :-1]).sum(axis=1)


def match(ws, win, ref=None):
    """ws: list of per-winter (T_w, R) arrays -> thresholds (R,), achieved rate (R,). The search grid
    spans the 0.1-40th percentiles of ref (T, R) (the whole control file), else of ws."""
    R = ws[0].shape[1]
    thr, rate = np.full(R, np.nan), np.zeros(R)
    for j in range(R):
        allx = np.concatenate([w[:, j] for w in ws]) if ref is None else ref[:, j]
        if not np.isfinite(np.concatenate([w[:, j] for w in ws])).all():
            continue
        grid = np.linspace(*np.nanpercentile(allx, [0.1, 40]), NGRID)
        c = np.mean([counts(w[:, j], grid, win) for w in ws], axis=0)
        k = int(np.argmax(c >= TARGET)) if (c >= TARGET).any() else int(np.argmax(c))
        thr[j], rate[j] = grid[k], c[k]
    return thr, rate


def masks(ws, thr, win, scale=None):
    """Per-winter boolean (T_w, R) event masks; scale (R,) multiplies the series first."""
    out = []
    for w in ws:
        x = w if scale is None else w * scale[None, :]
        out.append(np.stack([flags(x[:, j], thr[j], win) for j in range(w.shape[1])], axis=1))
    return out


def shares(ms):
    """Per-winter P(N >= k) (W, 8), per-winter marginal shares (W, R), events per winter."""
    P = np.array([[np.mean(m.sum(1) >= k) for k in K] for m in ms])
    M = np.array([m.mean(0) for m in ms])
    return P, M


def indep(M):
    """P(N >= k) for independent regions with marginal shares M (R,) (Poisson binomial)."""
    p = np.array([1.0])
    for q in M:
        p = np.convolve(p, [1 - q, q])
    return np.array([p[k:].sum() for k in K])


def jaccard(ms):
    m = np.concatenate(ms)
    R = m.shape[1]
    J = np.full((R, R), np.nan)
    for a in range(R):
        for b in range(R):
            if a != b:
                u = (m[:, a] | m[:, b]).sum()
                J[a, b] = (m[:, a] & m[:, b]).sum() / u if u else np.nan
    return J


def describe(ws, win, thr=None, ref=None):
    if thr is None:
        thr, rate = match(ws, win, ref)
    else:
        rate = np.full(len(thr), np.nan)
    ms = masks(ws, thr, win)
    P, M = shares(ms)
    J = jaccard(ms)
    pi = indep(M.mean(0))
    return dict(thr=thr, rate=rate, P=P.mean(0), P_w=P, M=M.mean(0), P_indep=pi, jaccard=J,
                jaccard_mean=float(np.nanmean(J)), ms=ms)


def compare(wc, wf, win, refc=None, reff=None):
    """Control and future winter lists -> dict of P_ctl, P_fut per variant, per-winter arrays.
    refc / reff: whole control / future files for the threshold search grid."""
    c = describe(wc, win, ref=refc)
    mc = np.mean(np.concatenate(wc), axis=0)
    mf = np.mean(np.concatenate(wf), axis=0)
    out = dict(ctl=c, mean_change_pct=100 * (mf / mc - 1))
    out["all"] = describe(wf, win, thr=c["thr"])
    out["mr"] = describe([w * (mc / mf)[None, :] for w in wf], win, thr=c["thr"])
    out["dep"] = describe(wf, win, ref=reff)
    return out


def lnr(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(b / a)


# ---------------------------------------------------------------- data
def winters6(name, y0, y1):
    """Per-winter (T_w, 8) series and the whole file (T, 8) (threshold search range)."""
    d = ev.load(name)
    x = d["wind"][:, EI]
    return [x[m] for _, m in ev.winters(d, y0, y1)], x


def winters_day(fn, exp, y0, y1):
    d = de.load(fn, exp)
    x = d["wind"][:, EI]
    return [x[m] for m in de.winter_masks(d, y0, y1)], x


def era5():
    z = np.load(ERA5)
    names = list(z["regions"])
    ix = [names.index(r) for r in EU]
    t6 = z["time6"].astype("datetime64[s]")
    w6 = ev.winter_of(t6)
    day = z["day"].astype("datetime64[D]")
    m = day.astype("datetime64[M]").astype(int) % 12 + 1
    y = day.astype("datetime64[Y]").astype(int) + 1970
    wd = np.where(m >= 9, y + 1, y)
    X6, XD = z["wind6"][:, ix].astype(float), z["wind_day"][:, ix].astype(float)
    return ([X6[w6 == yy] for yy in range(1980, 2027) if (w6 == yy).sum() >= ev.MIN_STEPS],
            [XD[wd == yy] for yy in range(1980, 2027) if (wd == yy).sum() >= de.MIN_DAYS])


# ---------------------------------------------------------------- null
_NULL = {}


def _placement(args):
    a, na, b, nb, win = args
    ws, ref = _NULL["ws"], _NULL["ref"]
    o = compare(ws[a:a + na], ws[b:b + nb], win, ref, ref)
    return np.stack([lnr(o["ctl"]["P"], o[v]["P"]) for v in ("all", "mr", "dep")])      # (3, 8)


def null_ranges(wsref, win, sizes, nproc):
    """{(na, nb): (3 variants, 2 [p05, p95], 8 k)} from NPL random disjoint placements."""
    ws, ref = wsref
    _NULL["ws"], _NULL["ref"] = ws, ref
    W = len(ws)
    out = {}
    for na0, nb in sizes:
        na = na0 if na0 + nb <= W - 40 else max(W - nb - 40, nb)
        pairs = ev.block_pairs(W, na, nb)
        sel = pairs[RNG.choice(len(pairs), min(NPL, len(pairs)), replace=False)]
        with Pool(nproc) as pool:
            L = np.array(pool.map(_placement, [(int(a), na, int(b), nb, win) for a, b in sel], chunksize=4))
        out[(na0, nb)] = np.stack([np.nanpercentile(L[:, v], [5, 95], axis=0) for v in range(3)])
        print(f"    null {na0}/{nb} (used {na}/{nb}, {len(sel)} placements) done", flush=True)
    return out


# ---------------------------------------------------------------- main
def summarise(o, dT):
    row = dict(dT=float(dT), mean_change_pct=dict(zip(EU, np.round(o["mean_change_pct"], 2).tolist())),
               thr=o["ctl"]["thr"].tolist(), rate=o["ctl"]["rate"].tolist())
    row["ctl"] = dict(P=o["ctl"]["P"].tolist(), P_indep=o["ctl"]["P_indep"].tolist(), M=o["ctl"]["M"].tolist(),
                      jaccard_mean=o["ctl"]["jaccard_mean"])
    for v in ("all", "mr", "dep"):
        L = lnr(o["ctl"]["P"], o[v]["P"])
        row[v] = dict(P=o[v]["P"].tolist(), P_indep=o[v]["P_indep"].tolist(), M=o[v]["M"].tolist(),
                      jaccard_mean=o[v]["jaccard_mean"], ratio=np.exp(L).tolist(),
                      perK=(100 * (np.exp(L / dT) - 1)).tolist())
    return row


def main():
    os.makedirs(OUT, exist_ok=True)
    nproc = int(os.environ.get("NPROC", 32))
    T = rpk.global_t()
    res = dict(regions=EU, k=K.tolist(), target=TARGET, days_per_winter=WINDAYS, six_hourly={}, daily={}, null={},
               era5={}, controls={})
    arrays = {}

    # our runs, 6-hourly
    sizes6 = set()
    for name in OURS6:
        (cn, c0, c1), (fn, f0, f1) = ev.COMPARISONS[name]
        (wc, rc), (wf, rf) = winters6(cn, c0, c1), winters6(fn, f0, f1)
        dT = rpk.mean_years(T, fn, f0, f1)[0] - rpk.mean_years(T, cn, c0, c1)[0]
        o = compare(wc, wf, 8, rc, rf)
        res["six_hourly"][name] = summarise(o, dT) | dict(n_winters=[len(wc), len(wf)])
        sizes6.add((len(wc), len(wf)))
        arrays[f"6h|{name}|ctl_P_w"] = o["ctl"]["P_w"]
        arrays[f"6h|{name}|fut_P_w"] = o["all"]["P_w"]
        res["controls"][f"6h|{name}"] = dict(P=o["ctl"]["P"].tolist(), P_indep=o["ctl"]["P_indep"].tolist(),
                                            jaccard=o["ctl"]["jaccard"].tolist(), jaccard_mean=o["ctl"]["jaccard_mean"])
        print(f"  6h {name}: dT {dT:.2f}  P(N>=4) {o['ctl']['P'][3]:.4f} -> {o['all']['P'][3]:.4f} "
              f"(mr {o['mr']['P'][3]:.4f}, dep {o['dep']['P'][3]:.4f})", flush=True)
        ev._CACHE.clear()

    # daily: CMIP6 and bridges
    C = de.comparisons()
    sizesd = set()
    for name, (fc, ec, c0, c1, ff, ef, f0, f1, Tc, Tf, group) in C.items():
        if group.endswith("_mid") or group == "bridge_cmipstyle":
            continue
        (wc, rc), (wf, rf) = winters_day(fc, ec, c0, c1), winters_day(ff, ef, f0, f1)
        dT = de.mean_T(Tf, f0, f1) - de.mean_T(Tc, c0, c1)
        o = compare(wc, wf, 2, rc, rf)
        res["daily"][name] = summarise(o, dT) | dict(group=group, n_winters=[len(wc), len(wf)])
        sizesd.add((len(wc), len(wf)))
        res["controls"][f"day|{name}"] = dict(P=o["ctl"]["P"].tolist(), P_indep=o["ctl"]["P_indep"].tolist(),
                                             jaccard=o["ctl"]["jaccard"].tolist(), jaccard_mean=o["ctl"]["jaccard_mean"])
        print(f"  day {name}: dT {dT:.2f}  P(N>=4) {o['ctl']['P'][3]:.4f} -> {o['all']['P'][3]:.4f} "
              f"(mr {o['mr']['P'][3]:.4f}, dep {o['dep']['P'][3]:.4f})", flush=True)

    # ERA5
    if os.path.exists(ERA5):
        e6, ed = era5()
        for tag, ws, win in (("6h", e6, 8), ("day", ed, 2)):
            o = describe(ws, win)
            res["era5"][tag] = dict(n_winters=len(ws), thr=o["thr"].tolist(), rate=o["rate"].tolist(), P=o["P"].tolist(),
                                    P_indep=o["P_indep"].tolist(), M=o["M"].tolist(), jaccard=o["jaccard"].tolist(),
                                    jaccard_mean=o["jaccard_mean"],
                                    P_first_half=o["P_w"][:len(ws) // 2].mean(0).tolist(),
                                    P_second_half=o["P_w"][len(ws) // 2:].mean(0).tolist())
            arrays[f"era5|{tag}|P_w"] = o["P_w"]
            print(f"  ERA5 {tag}: {len(ws)} winters, P(N>=k) " + " ".join(f"{x:.4f}" for x in o["P"]), flush=True)

    # nulls from the TCO95 control
    print("null 6h", flush=True)
    n6 = null_ranges(winters6("tco95_ctl1950g", 1953, 2099), 8, sorted(sizes6), nproc)
    ev._CACHE.clear()
    print("null day", flush=True)
    nd = null_ranges(winters_day(f"{de.D}/daily_tco95_ctl1950g.npz", "run", 1953, 2099), 2, sorted(sizesd), nproc)
    for tag, N, R in (("6h", n6, res["six_hourly"]), ("day", nd, res["daily"])):
        res["null"][tag] = {f"{a}_{b}": dict(all=v[0].tolist(), mr=v[1].tolist(), dep=v[2].tolist()) for (a, b), v in N.items()}
        for name, row in R.items():
            key = tuple(row["n_winters"])
            v = N[key]
            for j, var in enumerate(("all", "mr", "dep")):
                L = np.log(np.array(row[var]["ratio"]))
                row[var]["outside_null"] = [int(x) for x in np.where(L > v[j, 1], 1, np.where(L < v[j, 0], -1, 0))]

    # CMIP6 ensemble
    ens = {}
    names = [n for n, r in res["daily"].items() if r["group"] == "cmip6"]
    for var in ("all", "mr", "dep"):
        Pk = np.array([res["daily"][n][var]["perK"] for n in names])
        Lr = np.array([np.log(res["daily"][n][var]["ratio"]) for n in names])
        ens[var] = dict(n_models=len(names), perK_median=np.nanmedian(Pk, 0).tolist(),
                        perK_p10=np.nanpercentile(Pk, 10, 0).tolist(), perK_p90=np.nanpercentile(Pk, 90, 0).tolist(),
                        n_up=[int((Lr[:, j] > 0).sum()) for j in range(len(K))],
                        n_outside_null_up=[int(sum(res["daily"][n][var]["outside_null"][j] == 1 for n in names)) for j in range(len(K))])
    ens["ctl_P_median"] = np.median([res["daily"][n]["ctl"]["P"] for n in names], 0).tolist()
    ens["ctl_P_p10_p90"] = np.percentile([res["daily"][n]["ctl"]["P"] for n in names], [10, 90], 0).tolist()
    ens["ctl_jaccard_mean_median"] = float(np.median([res["daily"][n]["ctl"]["jaccard_mean"] for n in names]))
    ens["dep_jaccard_change_median"] = float(np.median([res["daily"][n]["dep"]["jaccard_mean"] - res["daily"][n]["ctl"]["jaccard_mean"] for n in names]))
    ens["mean_change_pct_median"] = float(np.median([np.mean(list(res["daily"][n]["mean_change_pct"].values())) for n in names]))
    res["cmip6_ensemble"] = ens

    json.dump(res, open(f"{OUT}/cooccur_events_summary.json", "w"), indent=1)
    np.savez(f"{OUT}/cooccur_events.npz", **arrays)

    # print
    def line(tag, n, row):
        f = lambda a: " ".join(f"{x:+6.1f}" for x in a)       # noqa: E731
        print(f"  {tag:4s} {n:40s} dT {row['dT']:.2f}  P_ctl(N>=1,4,8) {row['ctl']['P'][0]:.3f} {row['ctl']['P'][3]:.4f} {row['ctl']['P'][7]:.5f}"
              f"  J {row['ctl']['jaccard_mean']:.3f}->{row['dep']['jaccard_mean']:.3f}")
        for v in ("all", "mr", "dep"):
            print(f"        {v:3s} %/K k=1..8: {f(row[v]['perK'])}   outside null: {row[v].get('outside_null')}")
    print("\n=== co-occurrence change (P(N >= k) of 8 European regions in an event at once)")
    for n, row in res["six_hourly"].items():
        line("6h", n, row)
    for n, row in res["daily"].items():
        if row["group"] != "cmip6":
            line("day", n, row)
    for v in ("all", "mr", "dep"):
        e = ens[v]
        print(f"  CMIP6 median {v:3s} %/K k=1..8: " + " ".join(f"{x:+6.1f}" for x in e["perK_median"]) + f"  up {e['n_up']}  >null95 {e['n_outside_null_up']}")
    print(f"  CMIP6 control P(N>=k) median: " + " ".join(f"{x:.4f}" for x in ens["ctl_P_median"]))
    for tag, e in res["era5"].items():
        print(f"  ERA5 {tag} P(N>=k): " + " ".join(f"{x:.4f}" for x in e["P"]) + f"  indep " + " ".join(f"{x:.4f}" for x in e["P_indep"])
              + f"  J {e['jaccard_mean']:.3f}")
    print(f"Saved {OUT}/cooccur_events_summary.json")


if __name__ == "__main__":
    main()
