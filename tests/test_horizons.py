import contextlib
import pathlib
import tempfile
import numpy as np
import pytest
import sxs
from .conftest import skip_macOS_GH_actions_downloads

file_name = "SXS:BBH:0004/Lev6/Horizons.h5"


@skip_macOS_GH_actions_downloads
def test_spec_format():
    """Ensure that we can slice the output of spec_horizons_h5 just like the original Horizons.h5"""
    import h5py

    quantities = [
        "ArealMass", "ChristodoulouMass", "CoordCenterInertial", "DimensionfulInertialSpin",
        "DimensionfulInertialSpinMag", "chiInertial", "chiMagInertial"
    ]

    sim = sxs.load("SXS:BBH:4000v3.0")
    horizons = sim.horizons
    cached_path = horizons.__file__

    with h5py.File(cached_path, "r") as f:
        for horizon in "ABC":
            for quantity in quantities:
                name = f"Ah{horizon}.dir/{quantity}.dat"
                a, b = f[name], horizons[name]
                if 'Mag' in quantity:
                    assert np.allclose(a, b, atol=1e-11, rtol=1e-11)
                else:
                    assert np.array_equal(a, b)


@skip_macOS_GH_actions_downloads
def test_xmb_format():
    sim = sxs.load("SXS:BBH:4000v3.0")
    horizons_spec = sim.horizons

    with tempfile.TemporaryDirectory() as temp_dir:
        file = pathlib.Path(temp_dir) / 'horizons.h5'
        sxs.horizons.xor_multishuffle_bzip2.save(horizons_spec, file)
        with pytest.raises(ValueError):
            horizons_error = sxs.horizons.spec_horizons_h5.load(file)
        horizons_xmb = sxs.horizons.xor_multishuffle_bzip2.load(file)

    for horizon_name in "ABC":
        h_spec = horizons_spec[horizon_name]
        h_xmb = horizons_xmb[horizon_name]
        for attr in [a for a in dir(h_spec) if not a.startswith('_')]:
            d_spec = getattr(h_spec, attr)
            d_xmb = getattr(h_xmb, attr)
            assert np.array_equal(d_spec, d_xmb)


@skip_macOS_GH_actions_downloads
def test_horizon_existence():
    bhbh = sxs.load("SXS:BBH:0001", auto_supersede=True)
    bhbh_horizon = bhbh.horizons
    assert bhbh_horizon.A is not None
    assert bhbh_horizon.B is not None
    assert bhbh_horizon.C is not None
    bhns = sxs.load("SXS:BHNS:0001", auto_supersede=True)
    bhns_horizon = bhns.horizons
    assert bhns_horizon.A is not None
    assert bhns_horizon.B is None
    assert bhns_horizon.C is None
    nsns = sxs.load("SXS:NSNS:0001", auto_supersede=True)
    with pytest.raises(ValueError, match="Horizons.h5 not found in any form in files"):
        nsns_horizon = nsns.horizons
