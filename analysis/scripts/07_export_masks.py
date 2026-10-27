"""Write public/data/ar-masks-2010-01/ for the AR-mask animation.

Usage:  python scripts/07_export_masks.py
Needs data/masks/AR_masks_2010_01.nc (3-hourly AR masks from IVT, IVTE and IKEV).
"""

import json

from arbudget import climatology

if __name__ == "__main__":
    man = climatology.export_masks()
    print(json.dumps({k: man[k] for k in ("grid", "fps", "bits")}, indent=1, ensure_ascii=False))
    print(len(man["times"]), "frames:", man["times"][0], "..", man["times"][-1])
