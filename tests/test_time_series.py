import pickle
import copy
import numpy as np
import pytest
import sxs


def test_creation():
    np.random.seed(1234)
    n_times = 57

    # Pass more than one positional argument
    with  pytest.raises(ValueError):
        sxs.data.TimeSeries(np.random.rand(13, 7, 2), np.linspace(0, 10, num=13), np.random.rand(2))

    # Pass 0-dimensional input array
    with pytest.raises(ValueError):
        sxs.data.TimeSeries(np.array(3.0), np.linspace(0, 10))
    with pytest.raises(ValueError):
        sxs.data.TimeSeries(np.array(3.0), np.array(4.0))

    # Pass input array with non-finite values
    for non_finite in [np.inf, -np.inf, np.nan]:
        t = np.linspace(0, n_times)
        a = np.random.rand(n_times, 17)
        a[1, 2] = non_finite
        with pytest.raises(ValueError):
            sxs.data.TimeSeries(a, t)

    # Forget to specify time array
    with pytest.raises(ValueError):
        sxs.data.TimeSeries(np.random.rand(n_times, 17))

    # Pass complex time array
    with pytest.raises(ValueError):
        sxs.data.TimeSeries(np.random.rand(n_times, 17), np.linspace(0, n_times) + 0j)

    # Pass time array with ndim != 1
    with pytest.raises(ValueError):
        sxs.data.TimeSeries(np.random.rand(n_times, 17), np.array(3.0))
    with pytest.raises(ValueError):
        sxs.data.TimeSeries(np.random.rand(n_times, 17), np.linspace(0, n_times, 18).reshape((-1, 2)))

    # Pass time array with non-finite values
    for non_finite in [np.inf, -np.inf, np.nan]:
        t = np.linspace(0, 10, num=n_times)
        a = np.random.rand(n_times, 17)
        t[3] = non_finite
        with pytest.raises(ValueError):
            sxs.data.TimeSeries(a, t)

    # Pass time array with non-monotonic values
    t = np.linspace(0, 10, num=n_times)
    a = np.random.rand(n_times, 17)
    t[3] = t[2]
    with pytest.raises(ValueError):
        sxs.data.TimeSeries(a, t)
    t[3] = (t[1] + t[2]) / 2.0
    with pytest.raises(ValueError):
        sxs.data.TimeSeries(a, t)

    # Impossible time_axis
    for d_size in [-2, -1, 1, 2]:
        t = np.linspace(0, 10, num=n_times)
        a = np.random.rand(n_times+d_size)
        with pytest.raises(ValueError):
            sxs.data.TimeSeries(a, t)
        a = np.random.rand(n_times+d_size, 7)
        with pytest.raises(ValueError):
            sxs.data.TimeSeries(a, t)
        a = np.random.rand(7, n_times+d_size)
        with pytest.raises(ValueError):
            sxs.data.TimeSeries(a, t)

    # Incorrectly specified time_axis
    t = np.linspace(0, 10, num=n_times)
    a = np.random.rand(n_times, 7)
    with pytest.raises(ValueError):
        sxs.data.TimeSeries(a, t, time_axis=1)
    a = np.random.rand(7, n_times)
    with pytest.raises(ValueError):
        sxs.data.TimeSeries(a, t, time_axis=0)


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
    for shape, axis in zip(shapes, axes):
        t = np.linspace(0, 10, num=n_times)
        a = np.random.rand(*shape)
        with pytest.raises(ValueError):
            sxs.data.TimeSeries(a, np.linspace(0, n_times) + 1j * np.linspace(0, n_times))
        with pytest.raises(ValueError):
            sxs.data.TimeSeries(a, np.linspace(0, n_times-1))
        b = sxs.data.TimeSeries(a, t)
        assert b.time_axis == axis
        b = sxs.data.TimeSeries(a, time=t)
        assert b.time_axis == axis
