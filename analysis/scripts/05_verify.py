"""Physics checks on data/budget.nc:
  1. Ong & Yang (2024) Fig. 5 growth rates (Eq. 11a): 3-hour means 18-21Z on 3 and 4 Jan,
     projection region 160W-100W, 20-50N; paper (MERRA-2): +0.91/day and -0.75/day.
  2. Budget closure: RMS(resid)/RMS(tend) inside the AR mask.
  3. Effect of an area-based (instead of 50-point) blob filter on the masks.
"""

import json

import numpy as np
import xarray as xr

from arbudget import budget, config as C, growth


def region(ds):
    lon = ((ds.longitude + 180) % 360) - 180
    return ds.where((lon >= C.GROWTH_LON[0]) & (lon <= C.GROWTH_LON[1])
                    & (ds.latitude >= C.GROWTH_LAT[0]) & (ds.latitude <= C.GROWTH_LAT[1]), drop=True)


def fig5(ds, metric, day):
    key = "IVTE" if metric == "VTE" else "IKEV"
    t = [np.datetime64(f"2023-01-{day:02d}T{h:02d}") for h in (18, 19, 20, 21)]
    r = region(ds.sel(time=t))
    wts = np.array([0.5, 1.0, 1.0, 0.5]) / 3.0
    I = np.tensordot(wts, r[key].values, axes=1)
    lat = r.latitude.values
    res = {"tend_3h": growth.growth_rate((r[key].values[-1] - r[key].values[0]) / 10800.0, I, lat) * 86400}
    for term in budget.TERMS + ["resid"]:
        X = np.tensordot(wts, np.nan_to_num(r[f"{metric}_{term}"].values), axes=1)
        res[term] = growth.growth_rate(X, I, lat) * 86400
    return res


def closure(ds, metric):
    key = "IVTE" if metric == "VTE" else "IKEV"
    m = ds[f"mask_{key}"].values.astype(bool)
    tend = ds[f"{metric}_tend"].values[m]
    res = ds[f"{metric}_resid"].values[m]
    sumv = sum(ds[f"{metric}_{k}"].values[m] for k in budget.TERMS[1:])
    rms = lambda a: float(np.sqrt(np.nanmean(a ** 2)))
    return {"rms_resid_over_rms_tend": rms(res) / rms(tend), "corr_tend_sum": float(np.corrcoef(tend, sumv)[0, 1])}


def area_filter(ds):
    """Ong's 50-grid-point filter (calibrated on MERRA-2's 0.5x0.625 grid) vs an equivalent
    area filter (50 MERRA-2 cells ~ 1.5e5 km^2); masks computed in 02_budget on the full domain."""
    out = {}
    for k in ("IVT", "IVTE"):
        a = ds[f"mask_{k}"].values.astype(bool)
        b = ds[f"mask_area_{k}"].values.astype(bool)
        out[k] = {"count50_cells": int(a.sum()), "area_cells": int(b.sum()), "differing_cells": int((a != b).sum())}
    return out


def main():
    ds = xr.open_dataset(C.DATA / "budget.nc")
    report = {"fig5": {m: {d: fig5(ds, m, d) for d in (3, 4)} for m in budget.METRICS},
              "closure": {m: closure(ds, m) for m in budget.METRICS},
              "area_filter": area_filter(ds)}
    mi, mv = ds.mask_IVT.values.astype(bool), ds.mask_IVTE.values.astype(bool)
    mk = ds.mask_IKEV.values.astype(bool)
    jac = lambda a, b: float((a & b).sum() / max((a | b).sum(), 1))
    report["mask_overlap_jaccard"] = {"IVT_IVTE": jac(mi, mv), "IVT_IKEV": jac(mi, mk)}
    print(json.dumps(report, indent=1))
    (C.DATA / "verify_physics.json").write_text(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
