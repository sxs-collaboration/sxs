import pickle
import copy
import numpy as np
import pytest
import sxs


def test_creation():
    with  pytest.raises(ValueError):
        sxs.TimeSeries(np.random.rand(13, 7, 2), np.random.rand(13), np.random.rand(2))
    with pytest.raises(ValueError):
        sxs.TimeSeries(np.array(3), np.random.rand(13))
    
    n_times = 57
    shapes = [
        (n_times,),
        (n_times, 7),
        (n_times, 3, 7),
        (3, n_times, 7),
        (3, 7, n_times),
    ]
    axes = [
        0,
        0,
        0,
        1,
        2,
    ]

    for data_type in [float, complex]:
        for shape, axis in zip(shapes, axes):
            t = np.random.rand(n_times)
            a = np.random.rand(shape, dtype=data_type)
            with pytest.raises(ValueError):
                sxs.TimeSeries(a, np.random.rand(n_times) + 1j * np.random.rand(n_times))
            with pytest.raises(ValueError):
                sxs.TimeSeries(a, np.random.rand(n_times-1))
            b = sxs.TimeSeries(a, t)
            assert b.time_axis == axis
            b = sxs.TimeSeries(a, time=t)
            assert b.time_axis == axis

    # # Test successful creation with real data of the right shape
    # a = np.random.rand(3, 7, sf.LM_total_size(ell_min, ell_max)*2)
    # m = sf.Modes(a, spin_weight=s, ell_min=ell_min, ell_max=ell_max)
    # assert m.s == s
    # assert m.ell_min == 0  # NOTE: This is hard coded!!!
    # assert m.ell_max == ell_max
    # assert np.array_equal(a.view(complex), m[..., sf.LM_total_size(0, ell_min-1):])
    # assert np.all(m[..., :sf.LM_total_size(0, abs(s)-1)] == 0.0)
    # m = sf.Modes(a, spin_weight=s, ell_min=ell_min)  # ell_max is deduced!
    # assert m.s == s
    # assert m.ell_min == 0  # NOTE: This is hard coded!!!
    # assert m.ell_max == ell_max
    # assert np.array_equal(a.view(complex), m[..., sf.LM_total_size(0, ell_min-1):])
    # assert np.all(m[..., :sf.LM_total_size(0, abs(s)-1)] == 0.0)

    # # Test successful creation with complex data of the right shape
    # a = a.view(complex)
    # m = sf.Modes(a, spin_weight=s, ell_min=ell_min, ell_max=ell_max)
    # assert m.s == s
    # assert m.ell_min == 0  # NOTE: This is hard coded!!!
    # assert m.ell_max == ell_max
    # assert np.array_equal(a, m[..., sf.LM_total_size(0, ell_min-1):])
    # assert np.all(m[..., :sf.LM_total_size(0, abs(s)-1)] == 0.0)
    # m = sf.Modes(a, spin_weight=s, ell_min=ell_min)  # ell_max is deduced!
    # assert m.s == s
    # assert m.ell_min == 0  # NOTE: This is hard coded!!!
    # assert m.ell_max == ell_max
    # assert np.array_equal(a, m[..., sf.LM_total_size(0, ell_min-1):])
    # assert np.all(m[..., :sf.LM_total_size(0, abs(s)-1)] == 0.0)

    # # Test failed creation with complex data of inconsistent shape
    # if ell_min != 0:
    #     with pytest.raises(ValueError):
    #         m = sf.Modes(a, spin_weight=s)
    # with pytest.raises(ValueError):
    #     m = sf.Modes(a, spin_weight=s, ell_min=ell_min-1, ell_max=ell_max)
    # with pytest.raises(ValueError):
    #     m = sf.Modes(a, spin_weight=s, ell_min=ell_min+1, ell_max=ell_max)
    # with pytest.raises(ValueError):
    #     m = sf.Modes(a, spin_weight=s, ell_min=ell_min, ell_max=ell_max-1)
    # with pytest.raises(ValueError):
    #     m = sf.Modes(a, spin_weight=s, ell_min=ell_min, ell_max=ell_max+1)

    # # Test failed creation with complex data of impossible shape
    # with pytest.raises(ValueError):
    #     m = sf.Modes(a[..., 1:], spin_weight=s, ell_min=ell_min)

    # # Test successful creation with complex data containing extraneous data at ell<abs(s)
    # a = np.random.rand(3, 7, sf.LM_total_size(0, ell_max)*2)
    # a = a.view(complex)
    # m = sf.Modes(a, spin_weight=s)
    # assert m.s == s
    # assert m.ell_min == 0  # NOTE: This is hard coded!!!
    # assert m.ell_max == ell_max
    # assert np.all(m[..., :sf.LM_total_size(0, abs(s)-1)] == 0.0)

