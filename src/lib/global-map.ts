/**
 * Shared drawing for the global (Pacific-centred, equirectangular) maps written by
 * analysis/arbudget/climatology.py: FrequencyCompare and MaskAnimation.
 * The overlay uses the classes .grat, .coast and .tick, which each component styles.
 */

/** Grid of a manifest: cell centres at lon0 + i·dlon, lat0 + j·dlat (north first, dlat < 0). */
export type Grid = { lon0: number; dlon: number; nx: number; lat0: number; dlat: number; ny: number };
/** lines.json: coastline segments as [lon, lat] points, lon within the map's range. */
export type Lines = { coast: number[][][] };

/** Width / height of the map, one degree of longitude as wide as one of latitude. */
export const aspectRatio = (g: Grid) => (g.nx * g.dlon) / (g.ny * Math.abs(g.dlat));

export const fmtLat = (v: number) => (v === 0 ? '0°' : `${Math.abs(v)}°${v > 0 ? 'N' : 'S'}`);
export const fmtLon = (v: number) => {
  const e = ((v + 180) % 360 + 360) % 360 - 180; // −180..180
  return e === 0 || Math.abs(e) === 180 ? `${Math.abs(e)}°` : `${Math.abs(e)}°${e > 0 ? 'E' : 'W'}`;
};

/**
 * SVG markup (in W × H pixels) for the graticule, coastlines and graticule labels.
 * Latitude labels whose top would sit above `clearTop` px (under a map name in the upper
 * left) are left out; on small maps the longitude labels would crowd, so they are too.
 */
export function overlay(g: Grid, lines: Lines, W: number, H: number, clearTop = 0): string {
  const lonL = g.lon0 - g.dlon / 2, latT = g.lat0 - g.dlat / 2;
  const sx = W / (g.nx * g.dlon), sy = H / (g.ny * g.dlat);
  const X = (lon: number) => (lon - lonL) * sx, Y = (lat: number) => (lat - latT) * sy;
  const coast = lines.coast
    .map((s) => 'M' + s.map(([lo, la]) => `${X(lo).toFixed(1)},${Y(la).toFixed(1)}`).join('L')).join('');
  const lats = [-60, -30, 0, 30, 60], lons = [60, 120, 180, 240, 300];
  const grat = lats.map((la) => `M0,${Y(la).toFixed(1)}H${W}`).join('') + lons.map((lo) => `M${X(lo).toFixed(1)},0V${H}`).join('');
  const small = W < 480;
  const labels =
    lats.filter((la) => Y(la) + 2 > clearTop).map((la) => `<text class="tick" x="4" y="${Y(la)}" dy="1.1em">${fmtLat(la)}</text>`).join('') +
    (small ? '' : lons.map((lo) => `<text class="tick" x="${X(lo)}" y="${H - 5}" text-anchor="middle">${fmtLon(lo)}</text>`).join(''));
  return `<path class="grat" d="${grat}"/><path class="coast" d="${coast}"/>${labels}`;
}
