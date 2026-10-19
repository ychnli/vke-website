"""Link production budget terms (budget.py) to the NCL-verified port (ong_compat, corrected).

For the VTE metric (= Ong & Yang's IVKE), terms with identical formulas must agree:
    IVTE  == IVKE                     IVT == IVT
    pe    == IPEKE
    fric  == IKEpart - (IHAKE + IVAKE + IPEKE)
    cond  == IVPpart - (IHAV + IVAV)
and the flux-form transport differs from Ong's advective form only by the continuity defect:
    hflux + vflux  vs  IHAKE + IVAKE + IHAV + IVAV
Usage: python verify/production_vs_compat.py 2023010311 4
"""

import sys
from datetime import datetime, timedelta

import numpy as np

from arbudget import budget, config as C, ong_compat
from verify.ong_era5 import load, stats


def main(start, n_tim):
    t0 = datetime.strptime(start, "%Y%m%d%H")
    hours = [t0 + timedelta(hours=i) for i in range(n_tim + 2)]
    arr = load(hours)
    lat = np.arange(C.FETCH_LAT[1], C.FETCH_LAT[0] - 1e-6, -0.25).astype(np.float32)
    lon = np.arange(C.FETCH_LON[0], C.FETCH_LON[1] + 1e-6, 0.25).astype(np.float32)
    lev = np.array(C.LEVELS, dtype=np.float64)
    cmp = ong_compat.analyze(arr["Q"], arr["U"], arr["V"], arr["Z"], arr["W"], arr["SP"],
                             lev, lat, lon, faithful=False)
    p = (lev[::-1].astype(np.float32) * np.float32(100.0))  # bottom-to-top

    def frame(i):
        d = {k: arr[k][i][::-1] for k in ("Q", "U", "V", "W", "Z")}
        d["SP"] = arr["SP"][i]
        return d

    rows = []
    for j in range(n_tim - 1):  # production centred at hours[j+1]; compat inst index j+1, tend index j
        out = budget.frame_terms(frame(j), frame(j + 1), frame(j + 2), p, lat, lon, cyclic=True)
        v, ii, it = out["VTE"], j + 1, j
        inner = (slice(None), slice(1, -1), slice(None))
        pairs = {
            "IVTE vs IVKE": (out["IVTE"], cmp["IVKE"][ii]),
            "IVT": (out["IVT"], cmp["IVT"][ii]),
            "pe vs IPEKE": (v["pe"], cmp["IPEKE"][ii]),
            "fric vs KE residual": (v["fric"], cmp["IKEpart"][it] - (cmp["IHAKE"][ii] + cmp["IVAKE"][ii] + cmp["IPEKE"][ii])),
            "cond vs vapor residual": (v["cond"], cmp["IVPpart"][it] - (cmp["IHAV"][ii] + cmp["IVAV"][ii])),
            "h+vflux vs advective sum": (v["hflux"] + v["vflux"],
                                         cmp["IHAKE"][ii] + cmp["IVAKE"][ii] + cmp["IHAV"][ii] + cmp["IVAV"][ii]),
        }
        for name, (a, b) in pairs.items():
            s = stats(a[1:-1], b[1:-1])
            rows.append((hours[j + 1], name, s))
    for t, name, s in rows:
        print(f"{t:%d %HZ} {name:26s} rel_rms={s['rel_rms']:.2e} corr={s['corr']:.6f} max|d|={s['max_abs']:.3e}")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
