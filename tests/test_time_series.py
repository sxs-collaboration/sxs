import pickle
import copy
import numpy as np
import pytest
import sxs


def test_creation_failures():
    np.random.seed(1234)
    n_times = 57

    # Ask to store in Fortran order
    with pytest.raises(ValueError):
        sxs.TimeSeries(np.random.rand(13, 7, 2), np.linspace(0, 10, num=13), time_axis=0, order="F")
    with pytest.raises(ValueError):
        sxs.TimeSeries(np.random.rand(13, 7, 2), np.linspace(0, 10, num=13), time_axis=0, dtype=float, order="F")

    # Pass more than one positional argument
    with pytest.raises(ValueError):
        sxs.TimeSeries(np.random.rand(13, 7, 2), np.linspace(0, 10, num=13), np.random.rand(2))

    # Pass 0-dimensional input array
    with pytest.raises(ValueError):
        sxs.TimeSeries(np.array(3.0), np.linspace(0, 10))
    with pytest.raises(ValueError):
        sxs.TimeSeries(np.array(3.0), np.array(4.0))

    # Pass input array with non-finite values
    for non_finite in [np.inf, -np.inf, np.nan]:
        t = np.linspace(0, n_times)
        a = np.random.rand(n_times, 17)
        a[1, 2] = non_finite
        with pytest.raises(ValueError):
            sxs.TimeSeries(a, t)

    # Forget to specify time array
    with pytest.raises(ValueError):
        sxs.TimeSeries(np.random.rand(n_times, 17))

    # Pass complex time array
    with pytest.raises(ValueError):
        sxs.TimeSeries(np.random.rand(n_times, 17), np.linspace(0, n_times) + 0j)

    # Pass bool time array
    with pytest.raises(ValueError):
        sxs.TimeSeries(np.random.rand(n_times, 17), np.random.randint(0, 2, size=n_times, dtype=bool))

    # Pass time array with ndim != 1
    with pytest.raises(ValueError):
        sxs.TimeSeries(np.random.rand(n_times, 17), np.array(3.0))
    with pytest.raises(ValueError):
        sxs.TimeSeries(np.random.rand(n_times, 17), np.linspace(0, n_times, 18).reshape((-1, 2)))

    # Pass time array with non-finite values
    for non_finite in [np.inf, -np.inf, np.nan]:
        t = np.linspace(0, 10, num=n_times)
        a = np.random.rand(n_times, 17)
        t[3] = non_finite
        with pytest.raises(ValueError):
            sxs.TimeSeries(a, t)

    # Pass time array with non-monotonic values
    t = np.linspace(0, 10, num=n_times)
    a = np.random.rand(n_times, 17)
    t[3] = t[2]
    with pytest.raises(ValueError):
        sxs.TimeSeries(a, t)
    t[3] = (t[1] + t[2]) / 2.0
    with pytest.raises(ValueError):
        sxs.TimeSeries(a, t)

    # Mismatched sizes
    for d_size in [-2, -1, 1, 2]:
        t = np.linspace(0, 10, num=n_times)
        a = np.random.rand(n_times+d_size)
        with pytest.raises(ValueError):
            sxs.TimeSeries(a, t)
        a = np.random.rand(n_times+d_size, 7)
        with pytest.raises(ValueError):
            sxs.TimeSeries(a, t)
        a = np.random.rand(7, n_times+d_size)
        with pytest.raises(ValueError):
            sxs.TimeSeries(a, t)

    # Incorrectly specified time_axis
    t = np.linspace(0, 10, num=n_times)
    a = np.random.rand(n_times, 7)
    with pytest.raises(ValueError):
        sxs.TimeSeries(a, t, time_axis=1)
    a = np.random.rand(7, n_times)
    with pytest.raises(ValueError):
        sxs.TimeSeries(a, t, time_axis=0)

def test_creation_successes():
    np.random.seed(1234)
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
    for shape, axis in zip(shapes, axes):
        t = np.linspace(0, 10, num=n_times)
        a = np.random.rand(*shape)
        b = sxs.TimeSeries(a, t)
        assert b.time_axis == axis
        assert np.array_equal(b.ndarray, a)
        assert b.base is a
        i = (0,)*b.ndim
        b[i] = 1.2
        assert a[i] == 1.2
        b = sxs.TimeSeries(a, time=t)
        assert b.time_axis == axis
        assert np.array_equal(b.ndarray, a)
        b = sxs.TimeSeries(a, time=t, time_axis=axis)
        assert b.time_axis == axis
        assert np.array_equal(b.ndarray, a)


def test_view():
    np.random.seed(1234)
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
    for shape, axis in zip(shapes, axes):
        t = np.linspace(0, 10, num=n_times)
        a = np.random.rand(*shape)
        b = a.view(sxs.TimeSeries)


def test_slice():
    np.random.seed(1234)
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
    take_indices = np.array([0, 2, -1])
    for shape, axis in zip(shapes, axes):
        t = np.linspace(0, 10, num=n_times)
        a = np.random.rand(*shape)
        b = sxs.TimeSeries(a, t)
        # indices = tuple(slice(None) if i != b.time_axis else take_indices for i in range(b.ndim))
        # c = b[indices]
        # assert isinstance(c, type(b))
        # assert c.time_axis == b.time_axis
        # shape = list(b.shape)
        # shape[b.time_axis] = len(take_indices)
        # assert c.shape == tuple(shape)
        # assert c.n_times == len(take_indices)
        # assert np.array_equal(c.time, b.time[take_indices])
        # # c = np.take(b, take_indices, axis=b.time_axis)

    a = np.array([
        [0.0, 1.0, 2.0],
        [-1.0, 0.0, 1.0],
        [-2.0, -1.0, 0.0]
    ])
    b = sxs.TimeSeries(a, time=np.array([-1, 0, 1]), time_axis=0)
    with pytest.raises(ValueError):
        b[b < 0.0]


def test_basic_slicing():
    np.random.seed(1234)
    n_times = 57
    t = np.linspace(0, 10, num=n_times)

    for axis in [0, 1, 2, 3]:
        shape = [7, 9, 11, 12]
        shape[axis] = n_times
        a = np.random.rand(*shape)
        b = sxs.TimeSeries(a, time=t, time_axis=axis)

        # Only key is newaxis
        c = b[np.newaxis]
        assert type(c) == type(b)
        assert c.time_axis == b.time_axis + 1
        assert c.shape[0] == 1
        assert c.shape[1:] == b.shape
        assert np.array_equal(c.ndarray[0], b.ndarray)

        # Only key is a single integer
        c = b[2]
        assert type(c) == type(b)
        if axis == 0:
            assert c.time_axis == b.time_axis
            assert c.shape[1:] == b.shape[1:]
            assert c.shape[0] == 1
            assert np.array_equal(c.ndarray[0], b.ndarray[2])
            assert np.array_equal(c.time, b.time[2:3])
        else:
            assert c.time_axis == b.time_axis - 1
            assert c.shape == b.shape[1:]
            assert np.array_equal(c.ndarray, b.ndarray[2])
            assert np.array_equal(c.time, b.time)

        # Only key is a slice
        c = b[1:4]
        assert type(c) == type(b)
        assert c.time_axis == b.time_axis
        assert c.shape[0] == 3
        assert c.shape[1:] == b.shape[1:]
        assert np.array_equal(c.ndarray, b.ndarray[1:4])
        if axis == 0:
            assert np.array_equal(c.time, b.time[1:4])
        else:
            assert np.array_equal(c.time, b.time)

        # Only key is an Ellipsis
        c = b[...]
        assert type(c) == type(b)
        assert c.time_axis == b.time_axis
        assert c.shape == b.shape
        assert np.array_equal(c.ndarray, b.ndarray)
        assert np.array_equal(c.time, b.time)
        
        # Ellipsis + last key is a single integer
        c = b[..., 2]
        assert type(c) == type(b)
        if axis == b.ndim - 1:
            assert c.time_axis == b.time_axis
            assert c.shape[:-1] == b.shape[:-1]
            assert c.shape[-1] == 1
            assert np.array_equal(c.ndarray[..., 0], b.ndarray[..., 2])
            assert np.array_equal(c.time, b.time[2:3])
        else:
            assert c.time_axis == b.time_axis
            assert c.shape == b.shape[:-1], (axis, a.shape, b.shape, c.shape)
            assert np.array_equal(c.ndarray, b.ndarray[..., 2])
            assert np.array_equal(c.time, b.time)

        # Ellipsis + last key is a slice
        c = b[..., 1:4]
        assert type(c) == type(b)
        assert c.time_axis == b.time_axis
        assert c.shape[:-1] == b.shape[:-1]
        assert c.shape[-1] == 3
        assert np.array_equal(c.ndarray, b.ndarray[..., 1:4])
        if axis == b.ndim - 1:
            assert np.array_equal(c.time, b.time[1:4])
        else:
            assert np.array_equal(c.time, b.time)

        # First key is a single integer + Ellipsis
        c = b[2, ...]
        assert type(c) == type(b)
        if axis == 0:
            assert c.time_axis == b.time_axis
            assert c.shape[1:] == b.shape[1:]
            assert c.shape[0] == 1
            assert np.array_equal(c.ndarray[0], b.ndarray[2])
            assert np.array_equal(c.time, b.time[2:3])
        else:
            assert c.time_axis == b.time_axis - 1
            assert c.shape == b.shape[1:]
            assert np.array_equal(c.ndarray, b.ndarray[2])
            assert np.array_equal(c.time, b.time)

        # First key is a slice + Ellipsis
        c = b[1:4, ...]
        assert type(c) == type(b)
        assert c.time_axis == b.time_axis
        assert c.shape[0] == 3
        assert c.shape[1:] == b.shape[1:]
        assert np.array_equal(c.ndarray, b.ndarray[1:4])
        if axis == 0:
            assert np.array_equal(c.time, b.time[1:4])
        else:
            assert np.array_equal(c.time, b.time)

        # Single int + slice
        c = b[1, 1:4]
        assert type(c) == type(b)
        if axis == 0:
            assert c.shape[0] == 1
            assert c.shape[1] == 3
            assert c.shape[2:] == b.shape[2:]
            assert c.time_axis == b.time_axis
            assert np.array_equal(c.ndarray, b.ndarray[1:2, 1:4])
            assert np.array_equal(c.time, b.time[1:2])
        elif axis == 1:
            assert c.shape[0] == 3
            assert c.shape[1:] == b.shape[2:]
            assert c.time_axis == b.time_axis - 1
            assert np.array_equal(c.ndarray, b.ndarray[1, 1:4])
            assert np.array_equal(c.time, b.time[1:4])
        else:
            assert c.shape[0] == 3
            assert c.shape[1:] == b.shape[2:]
            assert c.time_axis == b.time_axis - 1
            assert np.array_equal(c.ndarray, b.ndarray[1, 1:4])
            assert np.array_equal(c.time, b.time)

        # Slice + single int
        c = b[1:4, 2]
        assert type(c) == type(b)
        if axis == 0:
            assert c.shape[0] == 3
            assert c.shape[1:] == b.shape[2:]
            assert c.time_axis == b.time_axis
            assert np.array_equal(c.ndarray, b.ndarray[1:4, 2])
            assert np.array_equal(c.time, b.time[1:4])
        elif axis == 1:
            assert c.shape[0] == 3
            assert c.shape[1] == 1
            assert c.shape[2:] == b.shape[2:]
            assert c.time_axis == b.time_axis
            assert np.array_equal(c.ndarray, b.ndarray[1:4, 2:3])
            assert np.array_equal(c.time, b.time[2:3])
        else:
            assert c.shape[0] == 3
            assert c.shape[1:] == b.shape[2:]
            assert c.time_axis == b.time_axis - 1
            assert np.array_equal(c.ndarray, b.ndarray[1:4, 2])
            assert np.array_equal(c.time, b.time)

        # Add a new axis in front of int-indexed first dimension
        c = b[np.newaxis, 3]
        assert type(c) == type(b)
        if axis == 0:
            assert c.shape[0] == 1
            assert c.shape[1] == 1
            assert c.shape[2:] == b.shape[1:]
            assert c.time_axis == b.time_axis + 1
            assert np.array_equal(c.ndarray, b.ndarray[np.newaxis, 3:4])
            assert np.array_equal(c.time, b.time[3:4])
        else:
            assert c.shape[0] == 1
            assert c.shape[1:] == b.shape[1:]
            assert c.time_axis == b.time_axis
            assert np.array_equal(c.ndarray, b.ndarray[np.newaxis, 3])
            assert np.array_equal(c.time, b.time)

        # Add a new axis between two int-indexed dimensions
        c = b[2, np.newaxis, 3]
        assert type(c) == type(b)
        if axis == 0:
            assert c.shape[0] == 1
            assert c.shape[1] == 1
            assert c.shape[2:] == b.shape[2:]
            assert c.time_axis == b.time_axis
            assert np.array_equal(c.ndarray, b.ndarray[2:3, np.newaxis, 3])
            assert np.array_equal(c.time, b.time[2:3])
        elif axis == 1:
            assert c.shape[0] == 1
            assert c.shape[1] == 1
            assert c.shape[2:] == b.shape[2:]
            assert c.time_axis == b.time_axis
            assert np.array_equal(c.ndarray, b.ndarray[2, np.newaxis, 3:4])
            assert np.array_equal(c.time, b.time[3:4])
        else:
            assert c.shape[0] == 1
            assert c.shape[1:] == b.shape[2:]
            assert c.time_axis == b.time_axis - 1
            assert np.array_equal(c.ndarray, b.ndarray[2, np.newaxis, 3])
            assert np.array_equal(c.time, b.time)

        with pytest.raises(ValueError):
            b[..., 3, ...]


def test_compare():
    np.random.seed(1234)
    n_times = 57
    a = np.array([
        [0.0, 1.0, 2.0],
        [-1.0, 0.0, 1.0],
        [-2.0, -1.0, 0.0]
    ])
    b = sxs.TimeSeries(a, time=np.array([-1, 0, 1]), time_axis=0)
    c = b < 0.0
    assert np.array_equal(b.time, c.time)
    assert b.time_axis == c.time_axis
    assert np.array_equal(a < 0.0, c)
    with pytest.raises(ValueError):
        b[c]


# def test_repr():
#     np.random.seed(1234)
#     n_times = 57
#     t = np.linspace(0, 10, num=n_times)
#     a = np.random.rand(n_times, 3, 5)
#     b = sxs.TimeSeries(a, time=t, time_axis=0)
#     print('repr:', repr(b))
#     print('str:', str(b))
    
