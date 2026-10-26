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


def export_global_lines(out_path, lon_left, lat_range, tolerance=0.15, min_length=1.0):
    """Coastlines for a global map whose longitudes run lon_left .. lon_left + 360.

    Lines are cut at the map's left/right edge and shifted into that range, so no segment
    jumps across the map. Lines shorter than min_length degrees (small islands) are dropped.
    """
    lat0, lat1 = lat_range
    # Natural Earth longitudes are -180..180; the map edge falls inside that range at one place.
    edge = (lon_left + 180) % 360 - 180
    # Each side of the edge is shifted by the multiple of 360 that lands it inside the map.
    pieces = []
    for lo, hi in ((-180, edge), (edge, 180)):
        if hi > lo:
            shift = next(s for s in (-360, 0, 360) if lon_left - 1e-9 <= lo + s < lon_left + 360)
            pieces.append((box(lo, lat0, hi, lat1), shift))
    segs = []
    for ln in _lines(SOURCES["coast"]):
        if ln.length < min_length:
            continue
        for clip, shift in pieces:
            g = ln.intersection(clip)
            if g.is_empty:
                continue
            g = g.simplify(tolerance)
            geoms = g.geoms if hasattr(g, "geoms") else [g]
            for s in geoms:
                if isinstance(s, LineString) and len(s.coords) >= 2:
                    segs.append([[round(x + shift, 2), round(y, 2)] for x, y in s.coords])
    out = {"coast": segs}
    out_path.write_text(json.dumps(out, separators=(",", ":")))
    return out
