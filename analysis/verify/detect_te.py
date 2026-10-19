"""Verification C: TempestExtremes DetectBlobs vs arbudget.detect on identical fields.

Runs Ong & Yang's exact DetectBlobs commands (reference/ong_yang_2024/TempestExtremes.sh)
on an IVT/IVKE file and compares the masks with the Python port cell by cell.
Usage: python verify/detect_te.py data/verify_ong/python_corrected_2023010311.nc
"""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import xarray as xr

from arbudget import config as C, detect

ENV = C.ROOT / ".ref-env"


def run_te(infile, var, thresh, outfile):
    cmd = [str(ENV / "bin" / "DetectBlobs"), "--in_data", str(infile), "--out", str(outfile),
           "--thresholdcmd", f"_LAPLACIAN{{9,7.0}}({var}),<,{thresh},0.0",
           "--filtercmd", f"_LAPLACIAN{{9,7.0}}({var}),<,{thresh},50",
           "--minabslat", "20.0", "--latname", "latitude", "--lonname", "longitude",
           "--regional"]  # our subset is regional; Ong ran on global files
    env = dict(os.environ, DYLD_FALLBACK_LIBRARY_PATH=str(ENV / "lib"))
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, stdin=subprocess.DEVNULL)
    if r.returncode != 0:
        print(r.stdout[-3000:], r.stderr[-3000:])
        sys.exit("DetectBlobs failed")
    return xr.open_dataset(outfile)


def main(infile):
    infile = Path(infile)
    ds = xr.open_dataset(infile)
    lat, lon = ds.latitude.values, ds.longitude.values
    op = detect.laplacian_operator(lat, lon)
    for var, thresh in (("IVT", "-40000.0"), ("IVKE", "-2000.0")):
        te = run_te(infile, var, thresh, infile.with_name(f"te_{var}_{infile.stem}.nc"))
        te_mask = te["binary_tag"].values.astype(bool) if "binary_tag" in te else te[list(te.data_vars)[0]].values.astype(bool)
        ours, _ = detect.detect(ds[var].values, lat, lon, float(thresh), op=op)
        diff = te_mask != ours
        print(f"{var}: TE cells={te_mask.sum()}  ours={ours.sum()}  mismatched={diff.sum()}"
              f"  ({diff.sum() / max(te_mask.sum(), 1):.2e} of TE AR cells)")


if __name__ == "__main__":
    main(sys.argv[1])
