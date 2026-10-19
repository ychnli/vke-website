"""Coastlines and borders for the web figure, clipped to the display domain."""

import json

import shapefile  # pyshp
from shapely.geometry import LineString, MultiLineString, box, mapping

from . import config as C

SOURCES = {
    "coast": "ne_50m_coastline",
    "states": "ne_50m_admin_1_states_provinces_lines",
    "countries": "ne_50m_admin_0_boundary_lines_land",
}


def _lines(name):
    path = C.DATA / "naturalearth" / name / f"{name}.shp"
    for shp in shapefile.Reader(str(path)).shapes():
        pts, parts = shp.points, list(shp.parts) + [len(shp.points)]
        for a, b in zip(parts[:-1], parts[1:]):
            if b - a >= 2:
                yield LineString(pts[a:b])


def export_lines(out_path, pad=1.0, tolerance=0.03):
    clip = box(C.DISPLAY_LON[0] - pad, C.DISPLAY_LAT[0] - pad, C.DISPLAY_LON[1] + pad, C.DISPLAY_LAT[1] + pad)
    out = {}
    for key, name in SOURCES.items():
        segs = []
        for ln in _lines(name):
            g = ln.intersection(clip)
            if g.is_empty:
                continue
            g = g.simplify(tolerance)
            geoms = g.geoms if hasattr(g, "geoms") else [g]
            for s in geoms:
                if isinstance(s, LineString) and len(s.coords) >= 2:
                    segs.append([[round(x, 3), round(y, 3)] for x, y in s.coords])
        out[key] = segs
    out_path.write_text(json.dumps(out, separators=(",", ":")))
    return out
