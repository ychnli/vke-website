import numpy as np

from arbudget import detect

LAT = np.arange(70.0, 4.99, -0.5)
LON = np.arange(170.0, 290.01, 0.5)


def _filament(lat_c=38.0, lon_c=220.0, length=30.0, width=3.0, amp=800.0):
    lo, la = np.meshgrid(LON, LAT)
    # elongated Gaussian ridge oriented SW-NE
    x = (lo - lon_c) * np.cos(np.deg2rad(35)) + (la - lat_c) * np.sin(np.deg2rad(35))
    y = -(lo - lon_c) * np.sin(np.deg2rad(35)) + (la - lat_c) * np.cos(np.deg2rad(35))
    return amp * np.exp(-(x / length) ** 2 - (y / width) ** 2)


def test_laplacian_of_quadratic_bowl_is_positive_constant_near_center():
    op = detect.laplacian_operator(LAT, LON)
    lo, la = np.meshgrid(np.deg2rad(LON - 230.0), np.deg2rad(LAT - 40.0))
    f = (lo * np.cos(np.deg2rad(40.0))) ** 2 + la ** 2  # ~ r^2 (radians); Laplacian ~ 4
    lap = detect.laplacian(f[None], op)[0]
    j, i = np.argmin(abs(LAT - 40.0)), np.argmin(abs(LON - 230.0))
    assert 3.0 < lap[j, i] < 5.0


def test_detects_filament_and_rejects_small_or_tropical_blobs():
    f = _filament()
    f += _filament(lat_c=50.0, lon_c=180.0, length=1.5, width=1.0)   # too small
    f += _filament(lat_c=12.0, lon_c=260.0)                          # equatorward of 20N
    mask, lap = detect.detect(f[None], LAT, LON, threshold=-4.0e4, count=50)
    mask = mask[0]
    lo, la = np.meshgrid(LON, LAT)
    assert mask[(abs(la - 38) < 1) & (abs(lo - 220) < 1)].all()
    assert not mask[(abs(la - 50) < 2) & (abs(lo - 180) < 2)].any()
    assert not mask[la < 20].any()
