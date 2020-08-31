import numpy as np


class TimeSeries(np.ndarray):
    """Array-like object representing time-series data

    This object wraps the basic numpy array object, but stores (at least a
    reference to) a corresponding array of time values, and provides several
    member functions for interpolating, differentiating, and integrating.

    Parameters
    ----------
    input_array : array_like, shape (..., n, ...)
        Input data representing the dependent variable, in any form that can be
        converted to a numpy array.  This includes scalars, lists, lists of
        tuples, tuples, tuples of tuples, tuples of lists, and numpy ndarrays.
        It can have an arbitrary number of dimensions, but the length along
        `time_axis` (see below) must match the length of `time`.  Values must
        be finite.
    time : array_like, shape (n,)
        1-D array containing values of the independent variable.  Values must
        be real, finite, and in strictly increasing order.
    time_axis : int, optional
        Axis along which `input_array` is assumed to be varying in time,
        meaning that for `time[i]` the corresponding values are
        `np.take(input_array, i, axis=time_axis)`.  Default is the first axis
        of `input_array` that has the same length as `time`.

    """

    def __new__(cls, input_array, *args, **kwargs):
        import copy
        if len(args) > 1:
            raise ValueError("Only one positional argument may be passed")
        elif len(args) == 1:
            kwargs["time"] = args[0]
        metadata = copy.copy(getattr(input_array, "_metadata", {}))
        metadata.update(**kwargs)

        # Interpret input_array as some type of array
        input_array = np.asanyarray(input_array, order='C')
        if input_array.ndim == 0:
            raise ValueError("Input array has 0 dimensions; it must have at least one")
        if not np.all(np.isfinite(input_array)):
            raise ValueError("Input array must contain only finite values.")

        # Get time array
        time = metadata.get("time", None)
        if time is None:
            raise ValueError(
                "Time data must be specified as part of input TimeSeries, as "
                "the second positional argument, or as a keyword argument."
            )
        time = np.asarray(time)
        if np.issubdtype(time.dtype, np.complexfloating):
            raise ValueError("Input `time` must contain real values; it has complex type.")
        time = time.astype(float)
        if time.ndim != 1:
            raise ValueError(f"Input `time` array must have exactly 1 dimension; it has {time.ndim}.")
        if not np.all(np.isfinite(time)):
            raise ValueError("Input `time` must contain only finite values.")
        if np.any(np.diff(time) <= 0):
            raise ValueError("Input `time` must be strictly increasing sequence.")
        metadata["time"] = time

        # Get time_axis
        time_axis = metadata.get("time_axis", None)
        if time_axis is None:
            for i in range(input_array.ndim):
                if input_array.shape[i] == time.size:
                    time_axis = i
                    break
        if time_axis is None:
            raise ValueError(
                f"Cannot find axis of size time.size={time.size} in input_array, "
                f"which has shape input_array.shape={input_array.shape}"
            )
        time_axis = time_axis % input_array.ndim  # Map time_axis into [0, input_array.ndim)
        if input_array.shape[time_axis] != time.shape[0]:
            raise ValueError(
                f"input_array.shape[time_axis]={input_array.shape[time_axis]} (with time_axis={time_axis}) "
                f"does not match time.shape[0]={time.shape[0]}"
            )
        metadata["time_axis"] = time_axis

        obj = input_array.view(cls)
        obj._metadata = metadata
        return obj

    def __array_finalize__(self, obj):
        import copy
        if obj is None:
            return
        # # Since ndarray.__array_finalize__ is None, we skip this in its direct descendents:
        # super().__array_finalize__(obj)
        self._metadata = copy.copy(getattr(self, '_metadata', getattr(obj, '_metadata', {})))
        if "time" not in self._metadata:
            self._metadata["time"] = None  # Placeholder about to be updated

    def __getitem__(self, key):
        """Extract a slice of this object

        Note that slicing this object works slightly differently than slicing
        the underlying ndarray object.

        First, if a single element is requested along the time dimension, that
        dimension will not be removed.  For example, with a 2-d ndarray `arr`,
        taking `arr[3]` will return a 1-d array (the first dimension will be
        removed); with a 2-d TimeSeries having `time_axis=0`, the returned
        TimeSeries will still be 2-d, where the first dimension will have size
        1.  If the requested element is not along the time dimension, the
        requested dimension will be removed as usual.

        Also, taking an irregular slice of this object is not permitted.  For
        example:

            >>> a = np.arange(3*4).reshape(3, 4)
            >>> a[a % 5 == 0]
            array([ 0,  5, 10])

        Even though `a % 5 == 0` is a 2-d array, indexing flattens `a` and the
        indexing set, so that the result is a 1-d array.  This probably does
        not make sense for TimeSeries arrays, so attempting to do something
        like this raises a ValueError.

        """
        from numbers import Integral
        from collections.abc import Iterable

        if key == ():
            raise ValueError(f"Empty index to {type(self)} does not make sense")

        if len(key) == self.ndim and all(isinstance(e, Iterable) for e in key) and len(set(len(e) for e in key)) == 1:
            raise ValueError(
                "Indexing an `n`-dimensional TimeSeries with `n` index arrays is not supported, "
                "as it is hard to do, and more often than not is not what the user wants anyway."
            )

        def basic_indexing(key):
            """Test if basic slicing occurs given this key

            Numpy's indexing documentation says "Basic slicing occurs when
            `obj` is a `slice` object (constructed by `start:stop:step`
            notation inside of brackets), an integer, or a tuple of slice
            objects and integers. Ellipsis and newaxis objects can be
            interspersed with these as well."

            <https://numpy.org/doc/stable/reference/arrays.indexing.html#basic-slicing-and-indexing>

            """
            if isinstance(key, Integral):
                return True
            elif isinstance(key, (slice, type(Ellipsis))):
                return True
            if not isinstance(key, tuple):
                return False
            for e in key:
                if not isinstance(e, (slice, Integral, type(Ellipsis))) and e is not np.newaxis:
                    return False
            return True

        def normalize_basic_indexing_key(key, ndim):
            """Translate key into full-size tuples, without Ellipsis

            Returns a pair of keys, one appropriate for the original array
            (with `newaxis` removed), and one for the new array (with `newaxis`
            replaced by `slice(None)`).

            """
            if not isinstance(key, tuple):
                key = (key,)
            n_missing = ndim - len(tuple(e for e in key if e not in [np.newaxis, None]))
            count = key.count(Ellipsis)
            if count == 0:
                full_key = key + (slice(None),) * n_missing
            elif count == 1:
                i_ellipsis = key.index(Ellipsis)
                full_key = key[:i_ellipsis] + (slice(None),) * n_missing + key[i_ellipsis+1:]
            else:
                raise ValueError(f"Cannot use more than one Ellipsis in key; found {count}.")
            old_key = tuple(e for e in full_key if e not in [np.newaxis, None])
            new_key = tuple(e if e not in [np.newaxis, None] else slice(None) for e in full_key)
            return old_key, full_key, new_key

        if basic_indexing(key):
            old_key, full_key, new_key = normalize_key(key, self.ndim)
            # Find the new time_axis
            i = 0
            new_time_axis = 0
            while i < self.time_axis:
                if full_key[new_time_axis] not in [np.newaxis, None]:
                    i += 1
                new_time_axis += 1
            # Index the old time array
            time_key = old_key[self.time_axis]
            if isinstance(time_key, Integral):
                full_key = list(full_key)
                full_key[new_time_axis] = slice(full_key[new_time_axis], full_key[new_time_axis]+1)
                full_key = tuple(full_key)
                time_key = slice(time_key, time_key+1)
            new_time = self.time[time_key]
            new_data = self.ndarray[full_key]
                

        # Slice the data like it's an instance of the parent class
        sliced = super().__getitem__(key)
        new_ndim = sliced.ndim

        # Now slice the time data, if necessary
        if isinstance(key, Integral):
            if self.time_axis == 0:
                new_time = self.time[key]
                new_time_axis = self.time_axis
                sliced = sliced[np.newaxis]
            else:
                new_time = self.time
                new_time_axis = self.time_axis - 1
        elif isinstance(key, (slice, type(Ellipsis))):
            if self.time_axis == 0:
                new_time = self.time[key]
                new_time_axis = self.time_axis
            else:
                new_time = self.time
                new_time_axis = self.time_axis
        elif basic_indexing_tuple(key):
            old_key, new_key = normalize_key(key, new_ndim)
            if len(key) > self.time_axis:
                new_time = self.time[key[self.time_axis]]
            elif key[0] == Ellipsis and self.time_axis - self.ndim >=  1 - len(key):
                new_time = self.time[key[self.time_axis - self.ndim]]
            else:
                new_time = self.time
        elif key is np.newaxis:
            new_time = self.time
            new_time_axis = self.time_axis + 1
        elif isinstance(key, np.ndarray):
            if key.ndim < self.time_axis:
                new_time = self.time
            else:
                raise NotImplementedError()
        else:
            # Advanced indexing
            raise NotImplementedError()
        # Create the new sliced object along with its metadata
        metadata = self._metadata.copy()
        metadata.update(**getattr(sliced, '_metadata', {}))
        metadata["time"] = new_time
        metadata["time_axis"] = new_time_axis
        try:
            sliced = type(self)(sliced, **metadata)
        except:
            pass
        return sliced

    @property
    def ndarray(self):
        """View this array as a numpy ndarray"""
        return self.view(np.ndarray)

    @property
    def time_axis(self):
        return self._metadata["time_axis"]

    @property
    def time(self):
        """Array of the time steps corresponding to the data"""
        return self._metadata["time"]

    @property
    def n_times(self):
        return self.time.size

    @property
    def time_broadcast(self):
        """Array of the time steps broadcast to same shape as data

        This property returns a new view (usually involving no copying of
        memory) of the `time` array, with additional dimensions to match the
        shape of the data.

        """
        new_shape = [1,] * self.ndim
        new_shape[self.time_axis] = self.n_times
        # This method is from the comment in np.broadcast_arrays, used so that
        # we don't have to disturb the data array (e.g., copying because of
        # weird ordering).
        nditer = np.nditer(
            self.time[tuple(new_shape)],
            flags=['multi_index', 'zerosize_ok'],
            itershape=self.shape,
            order='C'
        )
        return nditer.itviews[0]

    def interpolate(self, new_time, derivative_order=0, out=None):
        from scipy.interpolate import CubicSpline
        if derivative_order > 3:
            raise ValueError(
                f"{type(self)} interpolation uses CubicSpline, which cannot take a derivative of order {derivative_order}"
            )
        new_time = np.asarray(new_time)
        if new_time.ndim != 1:
            raise ValueError(f"New time array must have exactly 1 dimension; it has {new_time.ndim}.")
        new_shape = list(self.shape)
        new_shape[self.time_axis] = new_time.size
        if out is not None:
            out = np.asarray(out)
            if out.shape != new_shape:
                raise ValueError(
                    f"Output array should have shape {new_shape} for consistency with new time array and modes array"
                )
            if out.dtype != self.dtype:
                raise ValueError(
                    f"Output array should have same dtype as this array {self.dtype}; it has dtype {out.dtype}"
                )
        result = out or np.empty(new_shape, dtype=self.dtype)
        spline = CubicSpline(self.time, self.ndarray_float, axis=self.time_axis)
        if derivative_order < 0:
            spline = spline.antiderivative(-derivative_order)
        elif 0 < derivative_order <= 3:
            spline = spline.derivative(derivative_order)
        result[:] = spline(new_time)
        metadata = self._metadata.copy()
        metadata["time"] = new_time
        metadata["time_axis"] = self.time_axis
        return type(self)(result, **metadata)

    def antiderivative(self, antiderivative_order=1):
        """Integrate modes with respect to time"""
        return self.interpolate(self.time, derivative_order=-antiderivative_order)

    def derivative(self, derivative_order=1):
        """Differentiate modes with respect to time"""
        return self.interpolate(self.time, derivative_order=derivative_order)

    @property
    def dot(self):
        """Differentiate modes once with respect to time"""
        return self.derivative()

    @property
    def ddot(self):
        """Differentiate modes twice with respect to time"""
        return self.derivative(2)

    @property
    def int(self):
        """Integrate modes once with respect to time"""
        return self.antiderivative()

    @property
    def iint(self):
        """Integrate modes twice with respect to time"""
        return self.antiderivative(2)
