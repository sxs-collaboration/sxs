import contextlib
import pathlib
import tempfile
import numpy as np
import pytest
import sxs

file_name = "SXS:BBH:0004/Lev6/Horizons.h5"


def test_spec_format():
    """Ensure that we can slice the output of spec_horizons_h5 just like the original Horizons.h5"""
    import h5py

    quantities = [
        "ArealMass", "ChristodoulouMass", "CoordCenterInertial", "DimensionfulInertialSpin",
        "DimensionfulInertialSpinMag", "chiInertial", "chiMagInertial"
    ]

    with contextlib.redirect_stdout(None):
        catalog = sxs.load("catalog")
    selected = catalog.select_files(file_name)
    selected_path = sxs.utilities.sxs_path_to_system_path(list(selected.values())[0]["truepath"])

    with contextlib.redirect_stdout(None):
        horizons = sxs.load(file_name)
    cached_path = sxs.sxs_directory("cache") / selected_path

    with h5py.File(cached_path, "r") as f:
        for horizon in "ABC":
            for quantity in quantities:
                name = f"Ah{horizon}.dir/{quantity}.dat"
                a, b = f[name], horizons[name]
                if 'Mag' in quantity:
                    assert np.allclose(a, b, atol=1e-11, rtol=1e-11)
                else:
                    assert np.array_equal(a, b)


def test_xmb_format():
    # horizons_spec = sxs.horizons.spec_horizons_h5.load(file_name)
    with contextlib.redirect_stdout(None):
        horizons_spec = sxs.load(file_name)

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
