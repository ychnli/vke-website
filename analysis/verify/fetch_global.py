"""Fetch one global ERA5 hour (levels 1000..150) for verification B; lat S->N like MERRA-2."""

import sys
from datetime import datetime

import numpy as np
import xarray as xr

from arbudget import config as C
from arbudget.fetch import _open_store, _time_index


def main(stamp):
    t = datetime.strptime(stamp, "%Y%m%d%H")
    store = _open_store()
    lev = store["level"][:]
    ilev = np.array([int(np.where(lev == L)[0][0]) for L in C.LEVELS[::-1]])  # 1000..150
    lat = store["latitude"][:][::-1]  # S -> N
    lon = store["longitude"][:]
    it = _time_index(store, t)
    ds = xr.Dataset(coords={"time": np.array([t], dtype="datetime64[ns]"),
                            "level": lev[ilev].astype(np.float64), "latitude": lat, "longitude": lon})
    for var, short in C.VARS_3D.items():
        a = np.asarray(store[var][it])[ilev][:, ::-1, :].astype(np.float32)
        ds[short] = (("time", "level", "latitude", "longitude"), a[None])
        print(var, flush=True)
    ds["SP"] = (("time", "latitude", "longitude"), np.asarray(store["surface_pressure"][it])[::-1][None].astype(np.float32))
    out = C.DATA / f"global_{stamp}.nc"
    enc = {v: {"zlib": True, "complevel": 1} for v in ds.data_vars}
    enc["time"] = {"dtype": "float64", "units": "hours since 1900-01-01"}
    ds.to_netcdf(out, encoding=enc)
    print(out)


if __name__ == "__main__":
    main(sys.argv[1])
