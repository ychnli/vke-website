"""Line-by-line Python reproduction of Ong & Yang's `ARanalysisERA5.ncl`.

Used only for verification: running this and their NCL script on identical input must
give identical output. Inputs follow the NCAR RDA ds633.0 file layout the script was
written for: `level` ascending (1 ... 1000 hPa), latitude descending, longitude ascending.

`faithful=True` reproduces the script exactly as published, including two behaviors
documented in VERIFICATION.md:
  * every integral except IVAV passes its integrand with `level|::-1`, while `p` is
    ordered bottom-to-top (NCL returns level({1000:150}) as 1000 ... 150), so those
    integrands are vertically flipped relative to p;
  * tendency integrals (IKEpart, IVPpart) use surface pressure from the instantaneous
    times (one hour earlier than the tendency's centre).
`faithful=False` fixes both (integrand order, and ps at the tendency's centre time);
see `budget.py` for production code.

Time indexing (with the minimal patch needed for the script to run; see
reference/ong_yang_2024/ARanalysisERA5.patch): input has n_tim + 2 times; instantaneous
outputs use file times 0 .. n_tim-1, tendency outputs are centred on times 1 .. n_tim and
stored at output index 0 .. n_tim-1.
"""

from __future__ import annotations

import numpy as np

from .grid import center_finite_diff_n, grad_latlon_cfd, vibeta

F = np.float32
G = F(9.81)
PB, PT, PTP = 1000.0, 200.0, 150.0


def _vint(p, x_lev_desc, ps, flip):
    """vibeta(p, X(time|:,lat|:,lon|:,level|::-1 or :), 1, ps, pB, pT) / g, in float32."""
    x = np.moveaxis(x_lev_desc, 1, -1)  # (time, lat, lon, level), level as in X
    if flip:
        x = x[..., ::-1]
    out = vibeta(p.astype(np.float64), x, ps, PB * 100.0, PT * 100.0).astype(F)
    return out / G


def analyze(Q, U, V, Z, W, SP, level, lat, lon, faithful=True):
    """Q, U, V, Z(=geopotential), W(=omega): (time, level_asc, lat, lon) float32.
    SP: (time, lat, lon). Returns dict of the script's output fields (float32)."""
    level = np.asarray(level)
    sel = np.where((level >= PTP) & (level <= PB))[0][::-1]  # level({1000:150}) -> 1000..150
    n_all = Q.shape[0]
    n_tim = n_all - 2  # patch: file carries one extra hour on each side
    take = lambda A, t: np.asarray(A[t][:, sel], dtype=F)

    t0 = slice(0, n_tim)
    Qi, Ui, Vi, Hi, Wi = (take(A, t0) for A in (Q, U, V, Z, W))
    ps = np.asarray(SP[t0], dtype=F)
    p = (level[sel].astype(F) * F(100.0))
    LAT = np.asarray(lat, dtype=F)
    LON = np.asarray(lon, dtype=F)
    flip_fix = faithful  # published: level|::-1 for all integrals except IVAV

    out = {}
    QU = Qi * Ui
    QV = Qi * Vi
    K = (Ui * Ui + Vi * Vi) * F(0.5)
    VKE = Qi * Qi * K
    IQU = _vint(p, QU, ps, flip_fix)
    IQV = _vint(p, QV, ps, flip_fix)
    out["IVKE"] = _vint(p, VKE, ps, flip_fix)
    out["IVT"] = np.sqrt(IQU * IQU + IQV * IQV)
    Q2U = Qi * QU
    Q2V = Qi * QV
    KQ2 = Ui * QU + Vi * QV

    dUdp = center_finite_diff_n(Ui, p, axis=1)
    dVdp = center_finite_diff_n(Vi, p, axis=1)
    VAKE = -Q2U * Wi * dUdp - Q2V * Wi * dVdp
    out["IVAKE"] = _vint(p, VAKE, ps, flip_fix)

    dQdp = center_finite_diff_n(Qi, p, axis=1)
    VAV = -KQ2 * Wi * dQdp
    out["IVAV"] = _vint(p, VAV, ps, False)  # published: level|: (no flip)

    nlat = LAT.size
    inner = slice(1, nlat - 1)

    def _inner_vint(X):
        res = np.full(ps.shape, np.nan, dtype=F)
        x = np.moveaxis(X[:, :, inner, :], 1, -1)
        if flip_fix:
            x = x[..., ::-1]
        res[:, inner, :] = vibeta(p.astype(np.float64), x, ps[:, inner, :], PB * 100.0, PT * 100.0).astype(F) / G
        return res

    gy, gx = grad_latlon_cfd(K, LAT, LON, True)
    HAKE = -Q2U * gx - Q2V * gy
    out["IHAKE"] = _inner_vint(HAKE)

    gy, gx = grad_latlon_cfd(Hi, LAT, LON, True)
    PEKE = -Q2U * gx - Q2V * gy
    out["IPEKE"] = _inner_vint(PEKE)

    gy, gx = grad_latlon_cfd(Qi, LAT, LON, True)
    HAV = -KQ2 * Ui * gx - KQ2 * Vi * gy
    out["IHAV"] = _inner_vint(HAV)

    # Tendency parts, centred on file times 1..n_tim (patched indexing: Qall etc.)
    Qa, Ua, Va = (take(A, slice(0, n_all)) for A in (Q, U, V))
    c, fwd, bwd = slice(1, n_tim + 1), slice(2, n_tim + 2), slice(0, n_tim)
    Q2Ut = Qa[c] * Qa[c] * Ua[c]
    Q2Vt = Qa[c] * Qa[c] * Va[c]
    Utend = (Ua[fwd] - Ua[bwd]) / F(7200.0)
    Vtend = (Va[fwd] - Va[bwd]) / F(7200.0)
    KEpart = Q2Ut * Utend + Q2Vt * Vtend
    ps_t = ps if faithful else np.asarray(SP[c], dtype=F)  # published: instantaneous-time ps
    out["IKEpart"] = _vint(p, KEpart, ps_t, flip_fix)

    KQ2t = Qa[c] * (Ua[c] * Ua[c] + Va[c] * Va[c])
    Qtend = (Qa[fwd] - Qa[bwd]) / F(7200.0)
    VPpart = KQ2t * Qtend
    out["IVPpart"] = _vint(p, VPpart, ps_t, flip_fix)
    return out
