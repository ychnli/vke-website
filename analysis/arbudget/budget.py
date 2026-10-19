"""Column-integrated VKE budgets from ERA5 pressure-level data (production code).

For X in {VTE = q^2 K, KEV = q K}, with weights (w, s) = (q^2, 2qK) or (q, K):

    tend  = d<X>/dt                               (centred difference of the column integral)
    hflux = -< div_p(u X) >                        (Zhang et al. 2026, eq. 6/7, first term)
    vflux = transport - hflux,  transport = -< w (u.grad K + omega dK/dp) + s (u.grad q + omega dq/dp) >
            i.e. -< d(omega X)/dp > written via continuity, so that hflux + vflux equals Ong &
            Yang's advective transport exactly. Computing -< d(omega X)/dp > directly leaves a
            continuity defect -< X (div u + d omega/dp) > that is ~30% of the transport with
            ERA5 pressure-level omega (VERIFICATION.md); this form puts that defect in vflux.
    pe    = -< w u . grad_p(Phi) >
    cond  = < s (dq/dt + u.grad q + omega dq/dp) >             vapor residual  (Ong & Yang S12)
    fric  = < w (dK/dt + u.grad K + omega dK/dp + u.grad Phi) > KE residual    (Ong & Yang S11)

Discrete forms follow Ong & Yang's ERA5 code: dK/dt = u.du/dt (centred, 2 h) and
omega dK/dp = omega u.du/dp, which differ from differencing K itself by ~10% on hourly data.

ERA5 has no physics tendencies, so cond and fric are residuals of the moisture and kinetic
energy equations (as in Ong & Yang's ERA5 analysis); they also absorb analysis increments
and numerical error. resid = tend - (hflux + vflux + pe + cond + fric) is the surface-
pressure (Leibniz) term plus the error of the discrete product rule in time
(d<q^n K>/dt vs <s dq/dt + w u.du/dt>).

All finite differences and the vertical integral are the NCL ports in `grid.py`, which
reproduce Ong & Yang's NCL code exactly (see VERIFICATION.md). Integrands are passed
bottom-to-top, matching p (the published ERA5 script reverses them; see ong_compat).
"""

from __future__ import annotations

import numpy as np

from . import config as C
from .grid import center_finite_diff_n, div_latlon_cfd, grad_latlon_cfd, vibeta

F = np.float32
TERMS = ["tend", "hflux", "vflux", "pe", "cond", "fric"]
METRICS = {"VTE": 2, "KEV": 1}  # power of q in the metric


def _column(p, X, ps):
    """<X> = (1/g) * integral from ps (clipped at 1000 hPa) to 200 hPa; X: (lev, lat, lon)."""
    return (vibeta(p, np.moveaxis(X, 0, -1), ps, C.P_BOT, C.P_TOP) / C.G).astype(F)


def frame_terms(prev, cur, nxt, p, lat, lon, cyclic=False):
    """Budget terms at the time of `cur`.

    prev/cur/nxt: dicts with Q, U, V, W, Z of shape (lev, lat, lon) [lev bottom-to-top]
    and SP (lat, lon), one hour apart. Returns {metric: {term: (lat, lon)}} plus
    instantaneous column fields IVT, IVTE, IKEV.
    """
    dt = F(7200.0)
    q, u, v, w_, phi, ps = (cur[k] for k in ("Q", "U", "V", "W", "Z", "SP"))
    K = (u * u + v * v) * F(0.5)

    # tendencies of the ingredients (centred, 2 h); dK/dt as u.du/dt (Ong & Yang)
    dqdt = (nxt["Q"] - prev["Q"]) / dt
    dKdt = u * ((nxt["U"] - prev["U"]) / dt) + v * ((nxt["V"] - prev["V"]) / dt)

    # advection of the ingredients and PE->KE (per unit weight)
    Ky, Kx = grad_latlon_cfd(K, lat, lon, cyclic)
    qy, qx = grad_latlon_cfd(q, lat, lon, cyclic)
    py_, px_ = grad_latlon_cfd(phi, lat, lon, cyclic)
    u_gradK = u * Kx + v * Ky
    u_gradq = u * qx + v * qy
    u_gradPhi = u * px_ + v * py_
    w_dKdp = w_ * (u * center_finite_diff_n(u, p, axis=0) + v * center_finite_diff_n(v, p, axis=0))
    w_dqdp = w_ * center_finite_diff_n(q, p, axis=0)

    out = {}
    # instantaneous column fields for detection / display
    iqu = _column(p, q * u, ps)
    iqv = _column(p, q * v, ps)
    out["IVT"] = np.sqrt(iqu * iqu + iqv * iqv)

    for name, n in METRICS.items():
        wq = q ** n                       # w
        sq = F(n) * q ** (n - 1) * K      # s = dX/dq
        X = wq * K
        Xp = prev["Q"] ** n * (prev["U"] ** 2 + prev["V"] ** 2) * F(0.5)
        Xn = nxt["Q"] ** n * (nxt["U"] ** 2 + nxt["V"] ** 2) * F(0.5)
        col = _column(p, X, ps)
        out["I" + name] = col
        hflux = -_column(p, div_latlon_cfd(u * X, v * X, lat, lon, cyclic), ps)
        transport = -_column(p, wq * (u_gradK + w_dKdp) + sq * (u_gradq + w_dqdp), ps)
        terms = {
            "tend": (_column(p, Xn, nxt["SP"]) - _column(p, Xp, prev["SP"])) / dt,
            "hflux": hflux,
            "vflux": transport - hflux,
            "pe": -_column(p, wq * u_gradPhi, ps),
            "cond": _column(p, sq * (dqdt + u_gradq + w_dqdp), ps),
            "fric": _column(p, wq * (dKdt + u_gradK + w_dKdp + u_gradPhi), ps),
        }
        terms["resid"] = terms["tend"] - sum(terms[k] for k in TERMS[1:])
        out[name] = terms
    return out
