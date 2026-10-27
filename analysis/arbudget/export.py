"""Write the web data for the hero figure (public/data/<case>/).

Each source dataset has a profile (PROFILES) describing how to display it:
  era5   0.25-degree hourly budget from scripts/02_budget.py. Gaussian-smoothed
         (sigma 0.5 deg) and decimated to 0.5 deg; budget terms get a 3-hour running mean
         (Ong & Yang show 3-hour means; this also damps the jumps between ERA5 4D-Var
         windows at 09/10Z and 21/22Z). cond and fric are residuals.
  merra2 MERRA-2 budget produced on Sherlock with Zhang et al.'s pipeline (see
         VERIFICATION_MERRA2.md): native 0.5 x 0.625 deg grid, 3-hourly tavg windows, terms
         already 3-hour means; cond (condensation net of re-evaporation) and fric (turbulence
         only) from MERRA-2 physics tendencies. Lightly
         smoothed (sigma 0.3 deg); ana and spte are kept in the data file but not displayed.
AR masks are mask fractions, drawn as the 0.5 contour in the browser.
Encoding: terms Int8 with signed-sqrt companding up to E = 3 x the display scale;
IVT linear Uint8; IVTE/IKEV sqrt Uint8; masks Uint8 fraction. Files are gzip'd.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.ndimage import gaussian_filter

from . import budget, config as C, geo

PUBLIC = C.ROOT.parent / "public" / "data"
TERM_UNITS = "J m⁻² h⁻¹"  # displayed: SI W m-2 x 3600

PROFILES = {
    "era5": {
        "source": C.DATA / "budget.nc",
        "case": "landfall-2023",
        "title": "Atmospheric river landfall on the US West Coast, 3-5 January 2023 (ERA5, hourly)",
        "grid": "decimate-0.5",
        "sigma_deg": 0.5,
        "running_mean": True,
        "fps": 6,
        "frame_label": "",
        "notes": "ERA5 hourly; cond and fric are residuals of the vapor and kinetic-energy equations "
                 "(ERA5 has no physics tendencies); terms are 3-hour running means; 0.5 deg display grid.",
    },
    "merra2": {
        "source": C.DATA / "merra2" / "budget_merra2_2023-01.nc",
        "case": "landfall-2023-merra2",
        "title": "Atmospheric river landfall on the US West Coast, 3-5 January 2023 (MERRA-2, 3-hourly)",
        "grid": "native",
        "sigma_deg": 0.3,
        "running_mean": False,
        "fps": 2.5,
        "frame_label": "3-h mean centred on",
        "notes": "MERRA-2 3-hourly time-averaged budgets (Zhang et al. 2026 framework, computed on Sherlock); "
                 "cond is condensation net of re-evaporation and fric is turbulence only, both from MERRA-2 "
                 "physics tendencies; analysis increments and the surface-pressure "
                 "term are not displayed; native 0.5 x 0.625 deg grid.",
    },
}


def _lon180(ds):
    lon = ((ds.longitude + 180) % 360) - 180
    return ds.assign_coords(longitude=lon).sortby("longitude")


def _select_grid(ds, mode):
    """Return (lat, lon) of the display grid; lat north-to-south, lon in -180..180."""
    if mode == "decimate-0.5":
        lat = np.arange(C.DISPLAY_LAT[1], C.DISPLAY_LAT[0] - 1e-6, -0.5)
        lon = np.arange(C.DISPLAY_LON[0], C.DISPLAY_LON[1] + 1e-6, 0.5)
        return lat, lon
    d = _lon180(ds)
    lat = d.latitude.values[(d.latitude >= C.DISPLAY_LAT[0]) & (d.latitude <= C.DISPLAY_LAT[1])]
    lon = d.longitude.values[(d.longitude >= C.DISPLAY_LON[0]) & (d.longitude <= C.DISPLAY_LON[1])]
    return np.sort(lat)[::-1].astype(np.float64), lon.astype(np.float64)


def _prep(da, lat, lon, sigma_deg):
    """Smooth each frame (NaN -> 0) on the source grid, then pick the display points."""
    dlat = float(abs(da.latitude[1] - da.latitude[0]))
    dlon = float(abs(da.longitude[1] - da.longitude[0]))
    a = np.nan_to_num(da.values.astype(np.float64))
    sig = (sigma_deg / dlat, sigma_deg / dlon)
    a = np.stack([gaussian_filter(f, sig, mode="nearest") for f in a])
    sm = _lon180(xr.DataArray(a, dims=da.dims, coords=da.coords))
    return sm.sel(latitude=lat, longitude=lon, method="nearest").values


def _running3(a):
    out = np.empty_like(a)
    out[1:-1] = (a[:-2] + a[1:-1] + a[2:]) / 3.0
    out[0] = (a[0] + a[1]) / 2.0
    out[-1] = (a[-2] + a[-1]) / 2.0
    return out


def _write(out, name, arr, dtype):
    path = out / f"{name}.bin.gz"
    with open(path, "wb") as fh, gzip.GzipFile(fileobj=fh, mode="wb", mtime=0, compresslevel=9) as gz:
        gz.write(np.ascontiguousarray(arr, dtype=dtype).tobytes())
    return path.name


def _enc_signed_sqrt(x, emax):
    y = np.sign(x) * np.sqrt(np.clip(np.abs(x) / emax, 0, 1))
    return np.round(127 * y).astype(np.int8)


def _enc_sqrt_u8(x, vmax):
    return np.round(255 * np.sqrt(np.clip(x / vmax, 0, 1))).astype(np.uint8)


def _nice(x):
    e = 10 ** np.floor(np.log10(x))
    for m in (1, 2, 2.5, 5, 10):
        if x <= m * e:
            return float(m * e)
    return float(10 * e)


def _time_labels(times):
    out = []
    for t in times.astype("datetime64[m]"):
        s = str(t)  # YYYY-MM-DDTHH:MM
        out.append(s + "Z")
    return out


def export(profile: str = "era5"):
    P = PROFILES[profile]
    out = PUBLIC / P["case"]
    out.mkdir(parents=True, exist_ok=True)
    ds = xr.open_dataset(P["source"])
    lat, lon = _select_grid(ds, P["grid"])
    prep = lambda da: _prep(da, lat, lon, P["sigma_deg"])
    dlon = round(float(lon[1] - lon[0]), 6)
    dlat = round(float(lat[1] - lat[0]), 6)
    man = {
        "case": P["title"],
        "grid": {"lon0": float(lon[0]), "dlon": dlon, "nx": int(lon.size),
                 "lat0": float(lat[0]), "dlat": dlat, "ny": int(lat.size)},
        "times": _time_labels(ds.time.values),
        "frameLabel": P["frame_label"],
        "fps": P["fps"],
        "fields": {}, "masks": {}, "terms": {}, "contours": {},
        "termUnits": TERM_UNITS,
        "notes": P["notes"],
    }

    # masks
    masks = {}
    for k in ("IVT", "IVTE", "IKEV"):
        frac = prep(ds[f"mask_{k}"].astype(np.float32))
        masks[k] = frac > 0.5
        man["masks"][k] = {"file": _write(out, f"mask_{k.lower()}", np.round(255 * np.clip(frac, 0, 1)), np.uint8),
                           "encoding": "u8-fraction", "level": 0.5}

    # fields
    ivt = prep(ds["IVT"])
    man["fields"]["IVT"] = {"file": _write(out, "ivt", np.round(255 * np.clip(ivt / 1200.0, 0, 1)), np.uint8),
                            "encoding": "u8-linear", "vmax": 1200.0, "units": "kg m⁻¹ s⁻¹"}
    for k in ("IVTE", "IKEV"):
        f = prep(ds[k])
        vmax = _nice(float(np.percentile(f, 99.9)))
        man["fields"][k] = {"file": _write(out, k.lower(), _enc_sqrt_u8(f, vmax), np.uint8),
                            "encoding": "u8-sqrt", "vmax": vmax, "units": "J m⁻²"}
        p99 = float(np.percentile(f[masks[k]], 90)) if masks[k].any() else vmax / 2
        step = _nice(p99 / 4)
        man["contours"][k] = [round(step * (i + 0.5), 6) for i in range(12) if step * (i + 0.5) < vmax]

    # budget terms (displayed per hour)
    for m in budget.METRICS:
        key = "IVTE" if m == "VTE" else "IKEV"
        smooth_t = _running3 if P["running_mean"] else (lambda a: a)
        arrs = {t: smooth_t(prep(ds[f"{m}_{t}"])) * 3600.0 for t in budget.TERMS}
        inside = masks[key]
        ref = np.concatenate([np.abs(arrs["tend"][inside]), np.abs(arrs["pe"][inside])])
        scale = _nice(float(np.percentile(ref, 99)))
        emax = 3.0 * scale
        man["terms"][key] = {"scale": scale, "emax": emax, "encoding": "i8-signed-sqrt", "files": {}}
        for t, a in arrs.items():
            man["terms"][key]["files"][t] = _write(out, f"{key.lower()}_{t}", _enc_signed_sqrt(a, emax), np.int8)

    geo.export_lines(out / "lines.json")
    (out / "manifest.json").write_text(json.dumps(man, indent=1, ensure_ascii=False))
    return man

