"""Stream the regional ERA5 subset from ARCO.

Every ARCO chunk is one hour of the full globe at all 37 levels (~100 MB compressed),
so each (variable, hour) read downloads a whole chunk; we slice it immediately and keep
only the regional subset as a small .npy file. Finished files are skipped on re-runs.
"""

from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import gcsfs
import numpy as np
import xarray as xr
import zarr

from . import config as C


def _open_store():
    fs = gcsfs.GCSFileSystem(token="anon")
    return zarr.open_consolidated(fs.get_mapper(C.ARCO_URL), mode="r")


def _indices(store):
    lat = store["latitude"][:]
    lon = store["longitude"][:]
    lev = store["level"][:]
    ilat = np.where((lat >= C.FETCH_LAT[0]) & (lat <= C.FETCH_LAT[1]))[0]
    ilon = np.where((lon >= C.FETCH_LON[0]) & (lon <= C.FETCH_LON[1]))[0]
    ilev = np.array([int(np.where(lev == L)[0][0]) for L in C.LEVELS])
    return (slice(ilat[0], ilat[-1] + 1), slice(ilon[0], ilon[-1] + 1), ilev,
            lat[ilat], lon[ilon], lev[ilev])


def _time_index(store, t):
    i = int((t - C.ARCO_T0).total_seconds() // 3600)
    # time is "hours since 1900-01-01"; make sure the decoded value agrees.
    tv = int(store["time"][i])
    if tv != i:
        # the store does not start at 1900; fall back to offset from its first value
        i = i - int(store["time"][0])
        assert int(store["time"][i]) == int((t - C.ARCO_T0).total_seconds() // 3600)
    return i


def _cache_path(var, t):
    return C.CACHE / f"{var}_{t:%Y%m%d%H}.npy"


def _fetch_one(store, var, t, sl, retries=4):
    path = _cache_path(var, t)
    if path.exists():
        return path, 0.0
    slat, slon, ilev, *_ = sl
    it = _time_index(store, t)
    arr = store[var]
    for attempt in range(retries):
        try:
            t0 = time.time()
            if arr.ndim == 4:
                full = arr[it, :, slat, slon]  # downloads the whole global chunk
                data = full[ilev]
            else:
                data = arr[it, slat, slon]
            data = np.asarray(data, dtype=np.float32)
            if not np.isfinite(data).all():
                raise ValueError(f"non-finite values in {var} {t}")
            tmp = path.with_suffix(".tmp.npy")
            np.save(tmp, data)
            tmp.rename(path)
            return path, time.time() - t0
        except Exception as e:  # network hiccups: back off and retry
            if attempt == retries - 1:
                raise
            print(f"  retry {var} {t:%m-%d %H}Z after {e!r}", flush=True)
            time.sleep(2 ** attempt * 5)


def fetch(workers: int = 6):
    C.CACHE.mkdir(parents=True, exist_ok=True)
    store = _open_store()
    sl = _indices(store)
    tasks = [(v, t) for t in C.fetch_hours() for v in list(C.VARS_3D) + list(C.VARS_2D)]
    todo = [(v, t) for v, t in tasks if not _cache_path(v, t).exists()]
    print(f"{len(tasks)} chunks total, {len(todo)} to download", flush=True)
    t_start = time.time()
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(_fetch_one, store, v, t, sl): (v, t) for v, t in todo}
        for k, f in enumerate(as_completed(futs), 1):
            v, t = futs[f]
            _, dt = f.result()
            el = time.time() - t_start
            print(f"[{k}/{len(todo)}] {v} {t:%m-%d %H}Z {dt:.0f}s  (elapsed {el/60:.1f} min)", flush=True)
    return sl


def assemble():
    """Combine cached .npy files into one NetCDF with Ong/RDA-style variable names."""
    store = _open_store()
    _, _, _, lat, lon, lev = _indices(store)
    hours = C.fetch_hours()
    ds = xr.Dataset(coords={
        "time": ("time", np.array(hours, dtype="datetime64[ns]")),
        "level": ("level", lev.astype(np.int32), {"units": "hPa"}),
        "latitude": ("latitude", lat.astype(np.float64), {"units": "degrees_north"}),
        "longitude": ("longitude", lon.astype(np.float64), {"units": "degrees_east"}),
    })
    for var, short in {**C.VARS_3D, **C.VARS_2D}.items():
        data = np.stack([np.load(_cache_path(var, t)) for t in hours])
        dims = ("time", "level", "latitude", "longitude") if data.ndim == 4 else ("time", "latitude", "longitude")
        ds[short] = (dims, data, {"long_name": var})
    enc = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
    ds.to_netcdf(C.RAW_NC, encoding=enc)
    return ds
