"""Verification B: Zhang et al. (2026) KEV terms (their NCL lines, global ERA5) vs ours.

Their horizontal derivatives of K and q, and div(qu), are spectral (gradsf, uv2dvF) on the
global grid; ours are NCL-style centred finite differences on the regional grid. PE->KE
uses grad_latlon_cfd in both. Comparison over the display domain at one hour.
Usage: python -m verify.zhang_kev 2023010400   (after verify/fetch_global.py)
"""

import subprocess
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import xarray as xr

from arbudget import budget, config as C
from arbudget.fetch import _cache_path
from verify.ong_era5 import stats

HERE = Path(__file__).resolve().parent
LAT = np.arange(C.FETCH_LAT[1], C.FETCH_LAT[0] - 1e-6, -0.25).astype(np.float32)
LON = np.arange(C.FETCH_LON[0], C.FETCH_LON[1] + 1e-6, 0.25).astype(np.float32)
P = (np.array(C.LEVELS[::-1], dtype=np.float32) * np.float32(100.0))


def load_hour(t):
    d = {short: np.load(_cache_path(var, t))[::-1] for var, short in C.VARS_3D.items()}
    d["SP"] = np.load(_cache_path("surface_pressure", t))
    return d


def main(stamp):
    t = datetime.strptime(stamp, "%Y%m%d%H")
    gfile = C.DATA / f"global_{stamp}.nc"
    out = C.DATA / f"zhang_kev_{stamp}.nc"
    if not out.exists():
        r = subprocess.run([str(C.DATA / "ncl.sh"), f'f=addfile("{gfile}","r")', f'outfile="{out}"',
                            str(HERE / "zhang_kev_era5.ncl")], capture_output=True, text=True)
        if not out.exists():
            print(r.stdout[-3000:], r.stderr[-3000:])
            sys.exit("NCL failed")
    g = xr.open_dataset(gfile)
    raw = xr.open_dataset(out, mask_and_scale=True)  # vibeta output has unnamed dims
    z = xr.Dataset({k: (("latitude", "longitude"), np.asarray(raw[k].values).reshape(g.latitude.size, g.longitude.size))
                    for k in raw.data_vars}, coords={"latitude": g.latitude, "longitude": g.longitude})

    ours = budget.frame_terms(load_hour(t - timedelta(hours=1)), load_hour(t), load_hour(t + timedelta(hours=1)),
                              P, LAT, LON)
    o = xr.Dataset({k: (("latitude", "longitude"), v) for k, v in ours["KEV"].items()}
                   | {"IKEV": (("latitude", "longitude"), ours["IKEV"])},
                   coords={"latitude": LAT.astype(np.float64), "longitude": LON.astype(np.float64)})
    dom = dict(latitude=slice(C.DISPLAY_LAT[1], C.DISPLAY_LAT[0]),
               longitude=slice(C.DISPLAY_LON[0] % 360, C.DISPLAY_LON[1] % 360))
    o = o.sel(**dom)
    zz = z.sel(latitude=o.latitude, longitude=o.longitude)
    pairs = {
        "IKEV vs IVKE1": (o.IKEV, zz.IVKE1),
        "pe vs IPEKE": (o.pe, zz.IPEKE),
        "hflux vs IHAKE+IHCV": (o.hflux, zz.IHAKE + zz.IHCV),
        "vflux vs IVAKE+IVCV": (o.vflux, zz.IVAKE + zz.IVCV),
        "h+vflux vs sum of four": (o.hflux + o.vflux, zz.IHAKE + zz.IHCV + zz.IVAKE + zz.IVCV),
    }
    for name, (a, b) in pairs.items():
        s = stats(a.values, b.values)
        print(f"{name:24s} rel_rms={s['rel_rms']:.2e} corr={s['corr']:.4f} exact={s['exact_frac']:.3f}")


if __name__ == "__main__":
    main(sys.argv[1])
