"""Write public/data/ar-frequency-2010-2019/ for the AR-frequency comparison figure.

Usage:  python scripts/06_export_climatology.py
Needs data/climatology/AR_frequency_climatology.nc (AR frequency from IVT, IVTE and IKEV).
"""

import json

from arbudget import climatology

if __name__ == "__main__":
    man = climatology.export()
    print(json.dumps({k: man[k] for k in ("grid", "fields")}, indent=1, ensure_ascii=False))
