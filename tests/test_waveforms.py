import contextlib
import pathlib
import tempfile
import numpy as np
import h5py
import pytest
import sxs

# see test_utilities.py for explanation
shortest_h_com_file = "SXS:BBH:0156v1/Lev5/rhOverM_Asymptotic_GeometricUnits_CoM.h5"


def test_backwards_compatibility():
    path = sxs.sxs_directory("cache") / sxs.utilities.sxs_path_to_system_path(shortest_h_com_file)
    with contextlib.redirect_stdout(None):
        with sxs.loadcontext(shortest_h_com_file, extrapolation_order=...) as h:
            with h5py.File(path, "r") as f:
                for group in f:
                    if group == "VersionHist.ver":
                        continue
                    g = f[group]
                    for d in g:
                        if d == "History.txt":
                            continue
                        assert np.array_equal(g[d], h[group][d])
                for group in f:
                    if group != "VersionHist.ver":
                        continue
                    assert np.array_equal(f[group], h[group])
