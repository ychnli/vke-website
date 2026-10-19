"""Verification A: Ong & Yang's ERA5 NCL script vs arbudget.ong_compat on identical input.

Usage:  python verify/ong_era5.py START_HOUR N_TIM      e.g.  2023010311 4
Writes an RDA-layout input file (level ascending, lat descending), runs the minimally
patched NCL script (verify/ARanalysisERA5_patched.ncl), runs the Python port in both
`faithful` and `corrected` modes, and prints per-field comparison statistics.
"""

import subprocess
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import xarray as xr

from arbudget import config as C
from arbudget import ong_compat

HERE = Path(__file__).resolve().parent
WORK = C.DATA / "verify_ong"
FIELDS = ["IVT", "IVKE", "IHAKE", "IVAKE", "IPEKE", "IHAV", "IVAV", "IKEpart", "IVPpart"]


def load(hours):
    """Read cached subset; return arrays in RDA layout (level ascending 150..1000)."""
    from arbudget.fetch import _cache_path
    arr = {}
    for var, short in {**C.VARS_3D, **C.VARS_2D}.items():
        arr[short] = np.stack([np.load(_cache_path(var, t)) for t in hours]).astype(np.float32)
    return arr


def write_rda_like(arr, hours, lat, lon, path):
    lev = np.array(C.LEVELS, dtype=np.float64)  # ascending, as in RDA ds633.0 files
    ds = xr.Dataset(
        {k: (("time", "level", "latitude", "longitude"), arr[k]) for k in ("Q", "U", "V", "W", "Z")}
        | {"SP": (("time", "latitude", "longitude"), arr["SP"])},
        coords={"time": np.array(hours, dtype="datetime64[ns]"), "level": lev,
                "latitude": lat, "longitude": lon},
    )
    # RDA stores time as numeric hours since 1900; NCL cannot write int64/datetime64 back out
    enc = {v: {"zlib": True, "complevel": 1} for v in ds.data_vars}  # disk space
    enc["time"] = {"dtype": "float64", "units": "hours since 1900-01-01 00:00:00"}
    ds.to_netcdf(path, encoding=enc)


def stats(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m].astype(np.float64), b[m].astype(np.float64)
    d = a - b
    rms_b = np.sqrt(np.mean(b ** 2))
    return {
        "max_abs": float(np.abs(d).max()),
        "rel_rms": float(np.sqrt(np.mean(d ** 2)) / rms_b) if rms_b > 0 else 0.0,
        "corr": float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else 1.0,
        "exact_frac": float(np.mean(d == 0)),
    }


def main(start, n_tim):
    t0 = datetime.strptime(start, "%Y%m%d%H")
    hours = [t0 + timedelta(hours=i) for i in range(n_tim + 2)]
    WORK.mkdir(parents=True, exist_ok=True)
    lat = np.arange(C.FETCH_LAT[1], C.FETCH_LAT[0] - 1e-6, -0.25)
    lon = np.arange(C.FETCH_LON[0], C.FETCH_LON[1] + 1e-6, 0.25)
    arr = load(hours)
    inp = WORK / f"rda_like_{start}.nc"
    write_rda_like(arr, hours, lat, lon, inp)

    date = f"{start}"
    out_ncl = WORK / f"analysis_AR_era5_{date}.nc"
    out_ncl.unlink(missing_ok=True)
    args = [str(C.ROOT / "data" / "ncl.sh"), f"date={date}"]
    for f in ("q", "u", "v", "w", "z", "sp"):
        args.append(f'f_{f}=addfile("{inp}","r")')
    args.append(str(HERE / "ARanalysisERA5_patched.ncl"))
    r = subprocess.run(args, cwd=WORK, capture_output=True, text=True)
    print(r.stdout[-2000:], r.stderr[-2000:])
    if not out_ncl.exists():
        sys.exit("NCL produced no output")
    ncl = xr.open_dataset(out_ncl, mask_and_scale=True)

    args_py = dict(Q=arr["Q"], U=arr["U"], V=arr["V"], Z=arr["Z"], W=arr["W"], SP=arr["SP"],
                   level=np.array(C.LEVELS, dtype=np.float64), lat=lat, lon=lon)
    faithful = ong_compat.analyze(**args_py, faithful=True)
    corrected = ong_compat.analyze(**args_py, faithful=False)

    print(f"\n{'field':9s} | NCL vs faithful port: max|d|  rel_rms  exact | faithful vs corrected: rel_rms  corr")
    for k in FIELDS:
        a = ncl[k].values
        s1 = stats(faithful[k], a)
        s2 = stats(faithful[k], corrected[k])
        print(f"{k:9s} | {s1['max_abs']:10.3e} {s1['rel_rms']:9.2e} {s1['exact_frac']:6.3f} |"
              f" {s2['rel_rms']:9.2e} {s2['corr']:6.3f}")
    xr.Dataset({k: (("time", "latitude", "longitude"), v) for k, v in corrected.items()}).to_netcdf(
        WORK / f"python_corrected_{date}.nc")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
