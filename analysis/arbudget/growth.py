"""Projection diagnostics (Ong & Yang 2024, Eq. 11; Zhang et al. 2026, Eqs. 8-10).

growth(X) = ||X I|| / ||I^2||                      (per second; x86400 for per day)
movement(X) = ||X It~|| / ||It~^2||,   It~ = dI/dt - (||I dI/dt|| / ||I^2||) I
where ||.|| is an area integral over the projection region (cos-latitude weights).
"""

import numpy as np


def _area_weights(lat, nlon):
    return np.repeat(np.cos(np.deg2rad(np.asarray(lat, dtype=np.float64)))[:, None], nlon, axis=1)


def growth_rate(X, I, lat):
    w = _area_weights(lat, I.shape[-1])
    return float(np.nansum(w * X * I) / np.nansum(w * I * I))


def movement(X, I, dIdt, lat):
    w = _area_weights(lat, I.shape[-1])
    It = dIdt - (np.nansum(w * I * dIdt) / np.nansum(w * I * I)) * I
    return float(np.nansum(w * X * It) / np.nansum(w * It * It))
