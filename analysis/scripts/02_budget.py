"""Compute column budgets and AR masks for every display frame -> data/budget.nc."""

import time

import numpy as np
import xarray as xr

from arbudget import budget, config as C, detect
from arbudget.fetch import _cache_path

LAT = np.arange(C.FETCH_LAT[1], C.FETCH_LAT[0] - 1e-6, -0.25).astype(np.float32)
LON = np.arange(C.FETCH_LON[0], C.FETCH_LON[1] + 1e-6, 0.25).astype(np.float32)
P = (np.array(C.LEVELS[::-1], dtype=np.float32) * np.float32(100.0))  # bottom-to-top


def load_hour(t):
    d = {short: np.load(_cache_path(var, t))[::-1] for var, short in C.VARS_3D.items()}
    d["SP"] = np.load(_cache_path("surface_pressure", t))
    return d


def main():
    hours = C.fetch_hours()
    frames = hours[1:-1]
    fields = {"IVT": [], "IVTE": [], "IKEV": []}
    terms = {m: {k: [] for k in budget.TERMS + ["resid"]} for m in budget.METRICS}
    cache = {}
    t_start = time.time()
    for i, t in enumerate(frames):
        for tt in (hours[i], hours[i + 1], hours[i + 2]):
            if tt not in cache:
                cache[tt] = load_hour(tt)
        out = budget.frame_terms(cache[hours[i]], cache[hours[i + 1]], cache[hours[i + 2]], P, LAT, LON)
        cache.pop(hours[i], None)
        for k in fields:
            fields[k].append(out[k])
        for m in budget.METRICS:
            for k in terms[m]:
                terms[m][k].append(out[m][k])
        print(f"{t:%m-%d %HZ} done ({time.time() - t_start:.0f}s)", flush=True)

    ds = xr.Dataset(coords={"time": np.array(frames, dtype="datetime64[ns]"),
                            "latitude": LAT, "longitude": LON})
    dims = ("time", "latitude", "longitude")
    for k, v in fields.items():
        ds[k] = (dims, np.stack(v))
    for m in budget.METRICS:
        for k, v in terms[m].items():
            ds[f"{m}_{k}"] = (dims, np.stack(v))

    # detection on the full fetch domain (the 7-degree Laplacian needs the padding)
    op = detect.laplacian_operator(LAT, LON)
    merra_cell_km2 = (6371.0 ** 2) * np.deg2rad(0.5) * np.deg2rad(0.625) * np.cos(np.deg2rad(40.0))
    for k in ("IVT", "IVTE", "IKEV"):
        mask, _ = detect.detect(ds[k].values, LAT, LON, detect.THRESHOLDS[k], op=op)
        ds[f"mask_{k}"] = (dims, mask.astype(np.int8))
        if k != "IKEV":  # area-based size filter, for the comparison in 05_verify
            m2, _ = detect.detect(ds[k].values, LAT, LON, detect.THRESHOLDS[k], op=op,
                                  min_area_km2=50 * merra_cell_km2)
            ds[f"mask_area_{k}"] = (dims, m2.astype(np.int8))

    # keep only the display region plus a margin (disk space)
    ds = ds.sel(latitude=slice(60.0, 15.0), longitude=slice(190.0, 260.0))
    ds.to_netcdf(C.DATA / "budget.nc", encoding={v: {"zlib": True, "complevel": 4} for v in ds.data_vars})
    print(ds)


if __name__ == "__main__":
    main()
