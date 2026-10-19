"""Numerical building blocks, ported from NCL so results match Ong & Yang's code.

Arrays use NCL's (time, level, lat, lon) layout with levels ordered bottom-to-top
(p decreasing), the order NCL's `vibeta` requires.

- `vibeta`               : NCL vibeta (Fortran DVIBETA, linlog=1), vectorized
- `center_finite_diff_n` : NCL center_finite_diff_n (non-cyclic)
- `grad_latlon_cfd`      : NCL grad_latlon_cfd
- `div_latlon_cfd`       : spherical divergence built from the same finite differences
"""

from __future__ import annotations

import numpy as np

R_EARTH = 6.37122e6  # m; value used by NCL's grad_latlon_cfd (see gradient routines)


# --------------------------------------------------------------------------------------
# Vertical integral
# --------------------------------------------------------------------------------------

def vibeta(p, x, ps, pbot, ptop):
    """Port of NCL `vibeta(p, x, linlog=1, psfc, pbot, ptop)`.

    p    : (nlev,) pressure, Pa, strictly decreasing (bottom to top), must contain ptop
    x    : (..., nlev) integrand, level as the LAST axis (as NCL requires)
    ps   : (...) surface pressure, Pa
    Returns the integral over p (not divided by g), shape (...), computed in float64.

    The DVIBETA algorithm, reduced to the case where ptop is one of the levels and x has
    no missing values (true here):
      * keep the levels with p < ps (contiguous from the first level above ground);
      * "pi" levels alternate data levels and midpoints, with pi(0) = pbot;
      * each kept level k contributes x_k * beta_k * (pi_{k-} - pi_{k+}) where
          - bottom kept level: if pbot >= ps: x * (ps - mid_up)
                               elif p_k < pbot: x * (pbot - mid_up)
                               else (p_k == pbot < ps): 0      <- NCL quirk, reproduced
          - level == ptop    : x * (mid_down - ptop)
          - ptop < p_k       : x * (mid_down - mid_up)
          - p_k < ptop       : 0
    """
    p = np.asarray(p, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    ps = np.asarray(ps, dtype=np.float64)
    nlev = p.size
    if not np.all(np.diff(p) < 0):
        raise ValueError("vibeta: p must be ordered from bottom to top (decreasing)")
    ktop = int(np.where(np.trunc(p) == np.trunc(ptop))[0][0])

    mid = 0.5 * (p[:-1] + p[1:])           # mid[k] between level k and k+1
    mid_up = np.append(mid, 0.5 * p[-1])   # pi above level k
    kept = p < ps[..., None]               # (..., nlev)
    # first kept level index (levels are bottom-to-top, so kept is False...False True...True)
    k0 = np.argmax(kept, axis=-1)
    nkept = kept.sum(axis=-1)

    k = np.arange(nlev)
    mid_down = np.concatenate([[np.nan], mid])  # pi below level k (for k > k0)

    w = np.zeros(x.shape, dtype=np.float64)
    interior = (k[None] > k0[..., None]) & (p < 1e30)  # broadcast helper
    # interior levels strictly between ptop and the bottom kept level
    w = np.where(interior & (p > ptop), mid_down - mid_up, w)
    # ptop level
    w = np.where(interior & (k == ktop), mid_down - ptop, w)
    # bottom kept level
    is_k0 = k[None] == k0[..., None]
    ps_b = ps[..., None]
    bottom = np.where(pbot >= ps_b, ps_b - mid_up,
                      np.where(p < pbot, pbot - mid_up, 0.0))
    w = np.where(is_k0, bottom, w)
    # levels above ptop and below-ground levels contribute nothing
    w = np.where(kept & (k <= ktop), w, 0.0)

    out = (w * x).sum(axis=-1)
    # NCL returns missing if 3 or fewer levels are above the surface
    return np.where(nkept > 3, out, np.nan)


# --------------------------------------------------------------------------------------
# Horizontal / vertical finite differences
# --------------------------------------------------------------------------------------

def center_finite_diff_n(q, r, axis, cyclic=False):
    """Port of NCL `center_finite_diff_n(q, r, cyclic, opt=0, dim)` for 1-D coordinate r.

    Interior: (q[i+1] - q[i-1]) / (r[i+1] - r[i-1]).
    Non-cyclic ends are one-sided: (q[1]-q[0])/(r[1]-r[0]) and (q[n-1]-q[n-2])/(r[n-1]-r[n-2]).
    Cyclic ends wrap: (q[1]-q[n-1]) / (2(r[1]-r[0])) and (q[0]-q[n-2]) / (2(r[n-1]-r[n-2]))
    (verified against NCL 6.6.2). Computed in float64, returned in q's float type.
    """
    q = np.moveaxis(np.asarray(q), axis, -1)
    dtype = q.dtype
    qd = q.astype(np.float64)
    r = np.asarray(r, dtype=np.float64)
    out = np.empty_like(qd)
    out[..., 1:-1] = (qd[..., 2:] - qd[..., :-2]) / (r[2:] - r[:-2])
    if cyclic:
        out[..., 0] = (qd[..., 1] - qd[..., -1]) / (2 * (r[1] - r[0]))
        out[..., -1] = (qd[..., 0] - qd[..., -2]) / (2 * (r[-1] - r[-2]))
    else:
        out[..., 0] = (qd[..., 1] - qd[..., 0]) / (r[1] - r[0])
        out[..., -1] = (qd[..., -1] - qd[..., -2]) / (r[-1] - r[-2])
    return np.moveaxis(out.astype(dtype), -1, axis)


def grad_latlon_cfd(z, lat, lon, cyclic):
    """Port of NCL `grad_latlon_cfd(z, lat, lon, cyclic, opt)` (contributed.ncl).

    z: (..., nlat, nlon). lat/lon in degrees; lon assumed equally spaced. Returns
    (dz/dy, dz/dx) in units of z per metre, like NCL's [/ cgy, cgx /]. Pole rows of the
    zonal gradient are set to NaN as in NCL. Meridional sign follows the lat ordering.
    """
    # NCL does this arithmetic in the precision of `lat` (Ong passes float lat/lon), so
    # coordinates are built in that precision too; this matters at the 1e-5 level.
    ft = np.float32 if np.asarray(lat).dtype == np.float32 else np.float64
    lat = np.asarray(lat, dtype=ft)
    lon = np.asarray(lon, dtype=ft)
    rad = ft(4.0 * np.arctan(1.0) / 180.0)
    re = ft(R_EARTH)
    con = ft(re * rad)
    zdt = np.asarray(z).dtype
    cgy = center_finite_diff_n(z, lat, axis=-2, cyclic=False)
    cgy = (cgy / con).astype(zdt)
    dlon_re = ft(re * ft(lon[2] - lon[1]) * rad)
    lon_re = (np.arange(lon.size).astype(ft) * dlon_re).astype(ft)
    cgx = np.empty_like(cgy)
    for j, phi in enumerate(lat):
        lonx = (lon_re * ft(np.cos(np.float64(ft(phi * rad))))).astype(ft)  # NCL: double cos, cast
        cgx[..., j, :] = center_finite_diff_n(np.asarray(z)[..., j, :], lonx, axis=-1, cyclic=cyclic)
    for j in (0, lat.size - 1):
        if abs(lat[j]) >= 89.999:
            cgx[..., j, :] = np.nan
    return cgy, cgx


def div_latlon_cfd(u, v, lat, lon, cyclic):
    """Spherical divergence  (1/(a cos phi)) [du/dlambda + d(v cos phi)/dphi]  with the
    same centered differences and Earth radius as grad_latlon_cfd."""
    lat = np.asarray(lat, dtype=np.float64)
    coslat = np.cos(np.deg2rad(lat))[:, None]
    _, dudx = grad_latlon_cfd(u, lat, lon, cyclic)
    dvcos_dy, _ = grad_latlon_cfd(np.asarray(v, dtype=np.float64) * coslat, lat, lon, cyclic)
    return (dudx.astype(np.float64) + dvcos_dy.astype(np.float64) / coslat).astype(np.asarray(u).dtype)
