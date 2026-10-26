"""Write the web data for the AR-frequency comparison figure (public/data/ar-frequency-2010-2019/).

Source: AR frequency climatologies (days/year, 2010-2019, MERRA-2) from detection with IVT,
IVTE and IKEV, each with a threshold tuned to give a similar global AR frequency (Zhang et
al. 2026, Fig. 1). The NetCDF has lat -90..90, lon -180..180 on the native 0.5 x 0.625 deg
grid; the tropics (|lat| < 20) are excluded from detection and are zero.
Display: Pacific-centred (lon 0..360), north up, cropped to DISPLAY_LAT. Fields are linear
Uint8 up to VMAX (above every value in the file), gzip'd. Difference maps are computed in
the browser from these fields.
"""

from __future__ import annotations

import json

import numpy as np
import xarray as xr

from . import config as C, geo
from .export import PUBLIC, _write

SOURCE = C.DATA / "climatology" / "AR_frequency_climatology.nc"
CASE = "ar-frequency-2010-2019"
DISPLAY_LAT = (-75.0, 75.0)  # beyond this every map is < 10 days/year
LON_LEFT = 0.0               # map runs 0E .. 360E, centred on the date line
VMAX = 150.0                 # encoding range, days/year
COLOR_MAX = 120.0            # top of the color scale, as in Zhang et al. Fig. 1
DIFF_MAX = 20.0              # difference maps (IVTE or IKEV minus IVT): +-DIFF_MAX, about the
                             # 99th percentile of |difference| for both metrics
FIELDS = {"IVT": "ivt", "IVTE": "vte", "IKEV": "kev"}  # display name -> NetCDF variable


def export():
    out = PUBLIC / CASE
    out.mkdir(parents=True, exist_ok=True)
    full = xr.open_dataset(SOURCE)
    w = np.cos(np.deg2rad(full.lat))  # area weights for the global means
    ds = full.assign_coords(lon=(full.lon.round(6) - LON_LEFT) % 360 + LON_LEFT).sortby("lon")
    ds = ds.sel(lat=slice(DISPLAY_LAT[0], DISPLAY_LAT[1])).isel(lat=slice(None, None, -1))
    lat, lon = ds.lat.values, ds.lon.values
    man = {
        "case": f"AR frequency, {ds.attrs['years']} ({ds.attrs['source']})",
        "grid": {"lon0": float(lon[0]), "dlon": round(float(lon[1] - lon[0]), 6), "nx": int(lon.size),
                 "lat0": float(lat[0]), "dlat": round(float(lat[1] - lat[0]), 6), "ny": int(lat.size)},
        "units": "days/year",
        "colorMax": COLOR_MAX,
        "diffMax": DIFF_MAX,
        "fields": {},
    }
    for name, var in FIELDS.items():
        a = ds[var].fillna(0).values
        assert a.max() <= VMAX, f"{var} exceeds VMAX"
        man["fields"][name] = {
            "file": _write(out, name.lower(), np.round(255 * np.clip(a / VMAX, 0, 1)), np.uint8),
            "encoding": "u8-linear", "vmax": VMAX,
            "globalMean": round(float(full[var].weighted(w).mean()), 2),
            "threshold": int(ds.attrs[f"threshold_laplacian_{var}"]),
        }
    geo.export_global_lines(out / "lines.json", LON_LEFT, (DISPLAY_LAT[0] - 1, DISPLAY_LAT[1] + 1))
    (out / "manifest.json").write_text(json.dumps(man, indent=1, ensure_ascii=False))
    return man
