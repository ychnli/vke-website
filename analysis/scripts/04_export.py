"""Write public/data/<case>/ for one dataset profile.

Usage:  python scripts/04_export.py [era5|merra2]     (default: era5)
"""

import json
import sys

from arbudget import export

if __name__ == "__main__":
    profile = sys.argv[1] if len(sys.argv) > 1 else "era5"
    man = export.export(profile)
    print(json.dumps({k: man[k] for k in ("grid", "terms", "contours", "fps")}, indent=1, ensure_ascii=False))
