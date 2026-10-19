"""AR detection: port of TempestExtremes DetectBlobs with `_LAPLACIAN{N,R}` thresholds.

Ong & Yang (2024) / Zhang et al. (2026) settings (reference/*/TempestExtremes.sh):
    --thresholdcmd "_LAPLACIAN{9,7.0}(VAR),<,THRESH,0.0"
    --filtercmd    "_LAPLACIAN{9,7.0}(VAR),<,THRESH,50"
    --minabslat 20
Laplacian (TempestExtremes DataOp.cpp, BuildLaplacianOperator): for each node i, N points
at great-circle distance R (bearings 2*pi*j/N, starting east, rotating toward north) are
snapped to their nearest grid node k (k != i, no duplicates); then
    Lap f_i = sum_k (4/N) / d_ik^2 * (f_k - f_i)
with d_ik the true great-circle distance (radians) between node i and node k.
Blobs: 4-connected (no --diag_connect), regional grid (no longitude wrap); a blob is kept
if at least `count` of its points satisfy the threshold.

Thresholds (kg s-2 rad-2 for IVTE/IKEV, kg m-1 s-1 rad-2 for IVT):
    IVT  -4e4   (Ong & Yang; McClenny et al.)
    IVTE -2e3   (Ong & Yang, their "IVKE")
    IKEV -5.5e5 (fixed reference value in Zhang et al.'s AR_percentile_v3.py; their final
                 detection used the matching climatological percentile instead)
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree

THRESHOLDS = {"IVT": -4.0e4, "IVTE": -2.0e3, "IKEV": -5.5e5}


def _xyz(lat_rad, lon_rad):
    return np.stack([np.cos(lon_rad) * np.cos(lat_rad),
                     np.sin(lon_rad) * np.cos(lat_rad),
                     np.sin(lat_rad)], axis=-1)


def laplacian_operator(lat, lon, npts=9, dist_deg=7.0):
    """Sparse TempestExtremes Laplacian on a lat-lon grid (lat, lon in degrees, 1-D)."""
    LON, LAT = np.meshgrid(np.deg2rad(lon), np.deg2rad(lat))
    P = _xyz(LAT.ravel(), LON.ravel())
    n = P.shape[0]
    tree = cKDTree(P)

    x0, y0, z0 = P.T
    # eastward reference vector (y0, -x0, 0), scaled to tan(dist)
    e = np.stack([y0, -x0, np.zeros_like(x0)], axis=-1)
    e *= (np.tan(np.deg2rad(dist_deg)) / np.linalg.norm(e, axis=-1))[:, None]
    nvec = np.cross(P, e)  # northward, |n| = |e| since |P| = 1
    scale = 4.0 / npts
    rows, cols, vals = [], [], []
    diag = np.zeros(n)
    ks = []
    for j in range(npts):
        ang = 2.0 * np.pi * j / npts
        q = P + e * np.cos(ang) + nvec * np.sin(ang)
        q /= np.linalg.norm(q, axis=-1)[:, None]
        _, k = tree.query(q)
        ks.append(k)
    ks = np.stack(ks, axis=1)  # (n, npts)
    idx = np.arange(n)
    for i_row in range(npts):
        k = ks[:, i_row]
        # skip self and duplicates of earlier ring points
        dup = (k == idx)
        for prev in range(i_row):
            dup |= (k == ks[:, prev])
        ok = ~dup
        chord2 = np.sum((P[k[ok]] - P[idx[ok]]) ** 2, axis=-1)
        surf2 = (2.0 * np.arcsin(0.5 * np.sqrt(chord2))) ** 2
        w = (scale / surf2).astype(np.float32).astype(np.float64)
        rows.append(idx[ok]); cols.append(k[ok]); vals.append(w)
        np.add.at(diag, idx[ok], scale / surf2)
    rows.append(idx); cols.append(idx); vals.append((-diag).astype(np.float32).astype(np.float64))
    return csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))


def laplacian(field, op):
    """Apply operator to field (..., nlat, nlon)."""
    f = np.asarray(field, dtype=np.float32)
    shp = f.shape
    flat = f.reshape(-1, shp[-2] * shp[-1]).astype(np.float64)
    return (op @ flat.T).T.reshape(shp).astype(np.float32)


def detect(field, lat, lon, threshold, op=None, count=50, min_abs_lat=20.0, min_area_km2=None):
    """AR mask(s) for field (..., nlat, nlon). If min_area_km2 is given, blobs are kept by
    area instead of point count (to transfer a coarse-grid count to a finer grid)."""
    if op is None:
        op = laplacian_operator(lat, lon)
    lap = laplacian(field, op)
    cand = (lap < threshold) & (np.abs(np.asarray(lat))[:, None] >= min_abs_lat)
    area = None
    if min_area_km2 is not None:
        dlat = np.deg2rad(abs(lat[1] - lat[0])); dlon = np.deg2rad(abs(lon[1] - lon[0]))
        area = (6371.0 ** 2 * dlat * dlon * np.cos(np.deg2rad(lat)))[:, None] * np.ones((1, len(lon)))
    out = np.zeros(cand.shape, dtype=bool)
    for idx in np.ndindex(cand.shape[:-2]):
        lbl, nb = ndimage.label(cand[idx])  # 4-connectivity
        if nb == 0:
            continue
        if area is None:
            sizes = ndimage.sum(np.ones_like(lbl), lbl, index=np.arange(1, nb + 1))
            keep = np.where(sizes >= count)[0] + 1
        else:
            sizes = ndimage.sum(area, lbl, index=np.arange(1, nb + 1))
            keep = np.where(sizes >= min_area_km2)[0] + 1
        out[idx] = np.isin(lbl, keep)
    return out, lap
