"""
Does 6-hourly model output (as in the TCo319/TCO95 runs) change Mockert event statistics?

Degrades every 3-hourly combined-CF source in dunkelflaute/results/ (DART 1950C/2080C, ERA5
recon/real, SMARD) to what a 6-hourly run delivers: wind CF sampled instantaneously at
00/06/12/18 UTC, solar CF averaged over the 6 h window (mean of the two 3 h values, i.e. the
6 h accumulation). Events are recounted with the same 48 h window (8 steps instead of 16).

Output: dunkelflaute/results/robustness/sampling_3h_vs_6h.json
"""
import json
import os

import numpy as np

import robustness_mockert_events as R   # chdirs to the repo root

OUT = "dunkelflaute/results/robustness/sampling_3h_vs_6h.json"
THRESHOLDS = (0.04, 0.05, 0.06, 0.07, 0.08, 0.10)


def to_6h(ds):
    t = ds["time"]
    hour = (t - t.astype("datetime64[D]")) // np.timedelta64(1, "h")
    pos = {v: i for i, v in enumerate(t)}
    h3 = np.timedelta64(3, "h")
    end = np.array([i for i in np.where(hour % 6 == 0)[0]
                    if t[i] - h3 in pos and ds["winter"][pos[t[i] - h3]] == ds["winter"][i]])
    prev = np.array([pos[t[i] - h3] for i in end])
    return dict(time=t[end], winter=ds["winter"][end], on=ds["on"][end], off=ds["off"][end],
                solar=0.5 * (ds["solar"][end] + ds["solar"][prev]))


def main():
    data = R.load()
    out = {}
    for k, ds in data.items():
        d6 = to_6h(ds)
        out[k] = {}
        for thr in THRESHOLDS:
            a = R.counts(R.per_winter(ds, win=16, thr=thr))
            b = R.counts(R.per_winter(d6, win=8, thr=thr))
            r = float(np.corrcoef(a, b)[0, 1]) if a.std() and b.std() else None
            out[k][f"{thr:.2f}"] = dict(ev_3h=float(a.mean()), ev_6h=float(b.mean()), per_winter_r=r)
    out["ratio_2080c_1950c"] = {
        f"{thr:.2f}": dict(h3=out["2080c"][f"{thr:.2f}"]["ev_3h"] / out["1950c"][f"{thr:.2f}"]["ev_3h"],
                           h6=out["2080c"][f"{thr:.2f}"]["ev_6h"] / out["1950c"][f"{thr:.2f}"]["ev_6h"])
        for thr in THRESHOLDS}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(out, open(OUT, "w"), indent=1)
    for k in data:
        print(k, {t: f"{v['ev_3h']:.2f}->{v['ev_6h']:.2f}" for t, v in out[k].items()})
    print("ratio 2080C/1950C:", {t: f"{v['h3']:.2f}->{v['h6']:.2f}" for t, v in out["ratio_2080c_1950c"].items()})
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
