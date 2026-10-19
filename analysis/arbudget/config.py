"""Case-study configuration: the January 2023 landfalling AR (Ong & Yang 2024, Fig. 5).

Naming. The two papers use different names for the same quantity:
    VTE  = q^2 K  ("vapor transport energy", Zhang et al. 2026)
         = "VKE"/"IVKE" in Ong & Yang (2024) and in their NCL code
    KEV  = q K    ("kinetic energy of vapor", Zhang et al. 2026)
This package uses the Zhang et al. names except in `ong_compat`, which mirrors Ong's code.
"""

from datetime import datetime, timedelta
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CACHE = DATA / "cache"            # per-(variable, hour) .npy files from ARCO
RAW_NC = DATA / "era5_subset.nc"  # assembled regional subset

# --- ARCO ERA5 (Google Cloud public dataset) -------------------------------
ARCO_URL = "gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
ARCO_T0 = datetime(1900, 1, 1)  # time units: hours since 1900-01-01

VARS_3D = {  # ARCO name -> short name (ERA5/RDA convention, as used by Ong's NCL script)
    "specific_humidity": "Q",
    "u_component_of_wind": "U",
    "vertical_velocity": "W",  # omega, Pa s-1
    "geopotential": "Z",       # Phi, m2 s-2 (already geopotential, not height)
    "v_component_of_wind": "V",
}
VARS_2D = {"surface_pressure": "SP"}  # Pa

# --- Time window ------------------------------------------------------------
# Display frames: 2023-01-03 12Z .. 2023-01-05 12Z (49 hourly frames).
# One extra hour on each side for centered time differences.
FRAME_START = datetime(2023, 1, 3, 12)
FRAME_END = datetime(2023, 1, 5, 12)
FETCH_START = FRAME_START - timedelta(hours=1)
FETCH_END = FRAME_END + timedelta(hours=1)


def fetch_hours():
    n = int((FETCH_END - FETCH_START).total_seconds() // 3600) + 1
    return [FETCH_START + timedelta(hours=i) for i in range(n)]


# --- Domain -----------------------------------------------------------------
# Fetch domain is padded well beyond the display domain so the 7-degree TempestExtremes
# Laplacian and the finite differences are unaffected by edges (all reads are global
# chunks anyway, so padding costs disk, not download).
FETCH_LON = (160.0, 285.0)  # degrees east, 0..360
FETCH_LAT = (5.0, 72.0)
LEVELS = [150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750,
          775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]  # hPa (Ong reads 1000..150)

# Display domain (Ong & Yang Fig. 5 projection region is 160W-100W, 20N-50N).
DISPLAY_LON = (-165.0, -105.0)
DISPLAY_LAT = (20.0, 55.0)
GROWTH_LON = (-160.0, -100.0)  # Ong & Yang projection region for Fig. 5c,d
GROWTH_LAT = (20.0, 50.0)

# --- Physical constants (match Ong's NCL script) ------------------------------
G = 9.81
P_BOT = 1000.0e2  # Pa, vibeta lower limit (clipped by surface pressure)
P_TOP = 200.0e2   # Pa, vibeta upper limit
