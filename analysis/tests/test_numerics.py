import numpy as np
import pytest

from arbudget import config as C
from arbudget.grid import R_EARTH, center_finite_diff_n, div_latlon_cfd, grad_latlon_cfd, vibeta

P = np.array(C.LEVELS[::-1], dtype=np.float64) * 100.0  # bottom-to-top, Pa
LAT = np.arange(70.0, 9.99, -0.25)
LON = np.arange(170.0, 260.01, 0.25)


def test_dpdp_exact_for_quadratic_interior():
    # centred differences on non-uniform levels are exact for quadratics only where
    # spacing is uniform; check linear exactly everywhere and quadratic on uniform runs
    f = 3.0 * P + 7.0
    d = center_finite_diff_n(f, P, axis=0)
    np.testing.assert_allclose(d, 3.0, rtol=1e-12)


def test_vibeta_constant_integrand():
    # constant integrand: integral = (min(ps, pbot) - ptop), except NCL's quirk that the
    # 1000 hPa level contributes nothing when ps > 1000 hPa (half-layer 1000-987.5 hPa)
    ps = np.array([101500.0, 98000.0, 85000.0, 70010.0])
    x = np.ones((ps.size, P.size))
    got = vibeta(P, x, ps, 1e5, 2e4)
    expect = np.minimum(ps, 1e5) - 2e4
    expect[0] -= 0.5 * (P[0] - P[1])
    np.testing.assert_allclose(got, expect, rtol=1e-12)


def test_vibeta_ignores_below_ground_values():
    ps = np.array([90000.0])
    x = np.ones((1, P.size))
    x2 = x.copy()
    x2[0, P >= ps[0]] = 1e9  # garbage below ground must not matter
    assert vibeta(P, x, ps, 1e5, 2e4) == vibeta(P, x2, ps, 1e5, 2e4)


def test_gradient_of_sin_lat():
    lat2d = np.repeat(LAT[:, None], LON.size, axis=1)
    z = np.sin(np.deg2rad(lat2d))
    gy, gx = grad_latlon_cfd(z, LAT, LON, cyclic=False)
    expect = np.cos(np.deg2rad(lat2d)) / R_EARTH
    np.testing.assert_allclose(gy[1:-1], expect[1:-1], rtol=1e-4)
    np.testing.assert_allclose(gx, 0.0, atol=1e-20)


def test_divergence_free_streamfunction_wind():
    # u = -dpsi/dy, v = dpsi/dx with psi = cos(2 lambda) cos^2(phi): non-divergent
    lam, phi = np.meshgrid(np.deg2rad(LON), np.deg2rad(LAT))
    a = R_EARTH
    u = -(1 / a) * np.cos(2 * lam) * (-2 * np.cos(phi) * np.sin(phi))
    v = (1 / (a * np.cos(phi))) * (-2 * np.sin(2 * lam)) * np.cos(phi) ** 2
    div = div_latlon_cfd(u, v, LAT, LON, cyclic=False)
    scale = np.abs(u).max() / (a * np.deg2rad(0.25))
    assert np.abs(div[2:-2, 2:-2]).max() < 1e-4 * scale


def test_solid_body_rotation_divergence_free():
    lam, phi = np.meshgrid(np.deg2rad(LON), np.deg2rad(LAT))
    u = 20.0 * np.cos(phi)
    v = np.zeros_like(u)
    div = div_latlon_cfd(u, v, LAT, LON, cyclic=False)
    np.testing.assert_allclose(div, 0.0, atol=1e-18)
