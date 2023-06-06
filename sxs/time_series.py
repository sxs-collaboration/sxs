"""The `TimeSeries` object provides basic functionality for all time-series-like objects, including
`WaveformModes`, the objects contained in `Horizons`, etc.

"""

import numpy as np


class TimeSeries(np.ndarray):
    # noinspection PyUnresolvedReferences
    """Array-like object representing time-series data

    This object wraps the basic numpy array object, but stores (at least a
    reference to) a corresponding array of time values, and provides several member
    functions for interpolating, differentiating, and integrating.

    Parameters
    ----------
    input_array : (..., N, ...) array_like
        Input data representing the dependent variable, in any form that can be
        converted to a numpy array.  This includes scalars, lists, lists of tuples,
        tuples, tuples of tuples, tuples of lists, and numpy ndarrays.  It can have
        an arbitrary number of dimensions, but the length `N` along `time_axis`
        (see below) must match the length of `time`.  Values must be finite.
    time : (N,) array_like
        1-D array containing values of the independent variable.  Values must be
        real, finite, and in strictly increasing order.
    time_axis : int, optional
        Axis along which `input_array` is assumed to be varying in time, meaning
        that for `time[i]` the corresponding values are `np.take(input_array, i,
        axis=time_axis)`.  If this is not given, the first axis of `input_array`
        that has the same length as `time` is chosen as the time axis — which may
        be prone to errors.

    """

    def __new__(cls, input_array, *args, **kwargs):
        import copy

        dtype = kwargs.pop("dtype", None)
        order = kwargs.pop("order", "C")
        if order == "F":
            raise ValueError(f"Requested array order '{order}' is not supported; it must be 'C'.")
        if len(args) > 1:
            raise ValueError("Only one positional argument may be passed")
        elif len(args) == 1:
            kwargs["time"] = args[0]
        metadata = copy.copy(getattr(input_array, "_metadata", {}))
        metadata.update(**kwargs)

        # Interpret input_array as some type of array
        input_array = np.asanyarray(input_array, dtype=dtype, order=order)
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
        if not np.issubdtype(time.dtype, np.number):
            raise ValueError("Input `time` must contain numbers; its dtype is '{time.dtype}'.")
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

    def _slice(self, key):
        """Slice this object

        This is the core function used by __getitem__, and returns not only the sliced
        result, but also the indices used to slice the time array.  This is useful for
        subclasses when they also need the latter to slice additional time-related
        quantities.

        """
        from numbers import Integral

        if isinstance(key, tuple) and len(key) == 0:
            raise ValueError(f"Empty index to {type(self).__name__} does not make sense")

        def newaxis_type(e):
            return isinstance(e, (type(np.newaxis), type(None)))

        def basic_slicing(key_basic):
            """Test if basic slicing occurs given this key_basic

            Numpy's indexing documentation says "Basic slicing occurs when `obj` is
            a `slice` object (constructed by `start:stop:step` notation inside of
            brackets), an integer, or a tuple of slice objects and
            integers. Ellipsis and newaxis objects can be interspersed with these
            as well."

            <https://numpy.org/doc/stable/reference/arrays.indexing.html#basic-slicing-and-indexing>

            """
            if isinstance(key_basic, (Integral, slice, type(Ellipsis))):
                return True
            if not isinstance(key_basic, tuple):
                return False
            for e in key_basic:
                if not isinstance(e, (Integral, slice, type(Ellipsis), type(np.newaxis), type(None))):
                    return False
            return True

        def normalize_basic_slicing_key(key_normalize, ndim):
            """Translate key_basic into full-size tuples, without Ellipsis

            Returns a pair of keys, one appropriate for the original array (with
            `newaxis` removed), and one for the new array (with `newaxis` replaced
            by `slice(None)`).

            """
            if not isinstance(key_normalize, tuple):
                key_normalize = (key_normalize,)
            n_missing = ndim - len(tuple(e for e in key_normalize if not newaxis_type(e)))
            count = key_normalize.count(Ellipsis)
            if count == 0:
                full_key_normalize = key_normalize + (slice(None),) * n_missing
            elif count == 1:
                i_ellipsis = key_normalize.index(Ellipsis)
                full_key_normalize = (
                        key_normalize[:i_ellipsis]
                        + (slice(None),) * (n_missing + 1)
                        + key_normalize[i_ellipsis + 1:]
                )
            else:
                raise ValueError(f"Cannot use more than one Ellipsis in key_basic; found {count}.")
            old_key_normalize = tuple(e for e in full_key_normalize if not newaxis_type(e))
            return old_key_normalize, full_key_normalize

        time_key = slice(None)

        if key is np.newaxis:
            new_time_axis = self.time_axis + 1
            new_data = super().__getitem__(key)

        elif basic_slicing(key):
            old_key, full_key = normalize_basic_slicing_key(key, self.ndim)
            # Find the new time_axis
            i_old, integral_key_correction = 0, 0
            full_time_axis = 0
            for i in range(len(full_key)):
                full_time_axis = i
                key_i = full_key[i]
                if not newaxis_type(key_i):
                    i_old += 1
                if i_old > self.time_axis:
                    break
                if isinstance(key_i, Integral):
                    integral_key_correction += 1
            new_time_axis = full_time_axis - integral_key_correction
            # Expand the dimensions along the time axis if we're taking just one element
            time_key = old_key[self.time_axis]
            if isinstance(time_key, Integral):
                if time_key == -1:
                    time_key = slice(time_key, None)
                else:
                    time_key = slice(time_key, time_key+1)
                full_key = tuple(e if i != full_time_axis else time_key for i, e in enumerate(full_key))
            # Index the data and time arrays
            new_data = super().__getitem__(full_key)

        elif isinstance(key, np.ndarray) and key.ndim-1 <= self.time_axis:
            # This is okay; it's just indexing the leading dimensions
            new_data = super().__getitem__(key)
            new_time_axis = self.time_axis

        else:  # Try the slow way
            new_data = super().__getitem__(key)
            if new_data.ndim != self.ndim:
                raise ValueError(
                    f"\nAdvanced indexing of this {type(self).__name__} object with the key_basic\n\n"
                    + "    " + "\n    ".join(str(key).split("\n")) + "\n\n" +
                    f"changes the shape of the data from {self.shape} to {new_data.shape}.  It is not clear how this\n"
                    f"should affect the time data.  If you still want to index the data like this, extract the\n"
                    f"underlying data as `a.ndarray` and the time as `a.time`, slice them as desired, and reassemble\n"
                    f"them as\n\n"
                    f"    {type(self).__name__}(sliced_data, time=sliced_time, time_axis=sliced_time_axis)\n"
                )
            try:
                time_key = None
                new_time_slice = tuple(0 if i != self.time_axis else slice(None) for i in range(self.ndim))
                new_time = self.time_broadcast[key][new_time_slice].copy()
            except Exception as e:
                raise ValueError(
                    f"\nAdvanced indexing fails when trying to slice time information from this {type(self).__name__}\n"
                    f"object with the key_basic\n\n"
                    + "    " + "\n    ".join(str(key).split("\n")) + "\n\n" +
                    f"The shape of the {type(self).__name__} is {self.shape} and the time axis is {self.time_axis}."
                ) from e
            new_time_axis = self.time_axis

        # Create the new sliced object along with its metadata
        if time_key is not None:
            new_time = self.time[time_key]
        metadata = self._metadata.copy()
        metadata.update(**getattr(new_data, '_metadata', {}))
        metadata["time"] = new_time
        metadata["time_axis"] = new_time_axis

        return type(self)(new_data, **metadata), time_key

    def __getitem__(self, key):
        """Extract a slice of this object

        Note that slicing this object works slightly differently than slicing the
        underlying ndarray object, basically because we want to ensure that the
        returned object is still a TimeSeries object.

        First, if a single element is requested along the time dimension, that
        dimension will not be removed.  For a 2-d ndarray `arr`, taking `arr[3]` will
        return a 1-d array; the first dimension will be removed because only the third
        element is extracted.  For a 2-d TimeSeries `ts` with `time_axis=0`, `ts[3]`
        will return a 2-d TimeSeries; the first dimension will just have size 1,
        representing the third element.  If the requested element is not along the time
        dimension, the requested dimension will be removed as usual.

        Also, taking an irregular slice of this object is not permitted.  For example:

            >>> a = np.arange(3*4).reshape(3, 4)
            >>> a[a % 5 == 0]
            array([ 0,  5, 10])

        Even though `a % 5 == 0` is a 2-d array, indexing flattens `a` and the indexing
        set, so that the result is a 1-d array.  This probably does not make sense for
        TimeSeries arrays, so attempting to do something like this raises a ValueError.

        """
        return self._slice(key)[0]

    @property
    def ndarray(self):
        """View this array as a numpy ndarray"""
        return self.view(np.ndarray)

    @property
    def time_axis(self):
        """Axis of the array along which time varies

        At the time `time[i]`, the corresponding values of the data are
        `np.take(input_array, i, axis=time_axis)`.

        """
        return self._metadata["time_axis"]

    @property
    def time(self):
        """Array of the time steps corresponding to the data"""
        return self._metadata["time"]

    @time.setter
    def time(self, tprime):
        tprime = np.asarray(tprime)
        if np.issubdtype(tprime.dtype, np.complexfloating):
            raise ValueError("Input `tprime` must contain real values; it has complex type.")
        if not np.issubdtype(tprime.dtype, np.number):
            raise ValueError("Input `tprime` must contain numbers; its dtype is '{tprime.dtype}'.")
        tprime = tprime.astype(float)
        if tprime.ndim != 1:
            raise ValueError(f"Input `tprime` array must have exactly 1 dimension; it has {tprime.ndim}.")
        if not np.all(np.isfinite(tprime)):
            raise ValueError("Input `tprime` must contain only finite values.")
        if np.any(np.diff(tprime) <= 0):
            raise ValueError("Input `tprime` must be strictly increasing sequence.")
        if self._metadata["time"] is not None:
            if tprime.shape != self._metadata["time"].shape:
                raise ValueError(
                    f"Input `tprime` must have same shape as original; they are {tprime.shape} "
                    f"and {self._metadata['time'].shape}"
                )
        self._metadata["time"] = tprime

    t = time

    def register_modification(self, func, **kwargs):
        """Add a record of a modification to the metadata

        Note that this function does not actually run the modification; it simply
        records the function name and arguments in this object's metadata.  You are
        expected to run the function for yourself, with the given keyword arguments.

        Also note that the modifications will most likely be written to JSON, so you
        should adjust them to be in basic formats suitable for JSON.  For example, if
        an argument `arr` is ordinarily passed as a numpy array, you should convert to
        a list, with something like `arr.tolist()`.

        Parameters
        ----------
        func : named function
            The function that will modify (or already has modified) this object.  The
            function must have a `__name__` attribute, as most functions do.

        Because we cannot know whether `func` modifies `self` in place or not, we
        cannot design this function to modify the desired result, which is why you must
        call `func` yourself.

        """
        if "modifications" not in self._metadata:
            self._metadata["modifications"] = {}
        self._metadata["modifications"][func.__name__] = kwargs
        return self

    def index_closest_to(self, t):
        """Time index closest to the given time `t`

        Parameters
        ----------
        t : float

        Returns
        -------
        idx : int
            Index such that abs(self.time[idx]-t) is as small as possible

        """
        idx = min(np.searchsorted(self.t, t), self.n_times-1)
        if idx > 0:
            if abs(self.t[idx]-t) > abs(self.t[idx-1]-t):
                idx -= 1
        return idx

    @property
    def n_times(self):
        """Size of the array along the time_axis"""
        return self.time.size

    @property
    def time_broadcast(self):
        """Array of the time steps broadcast to same shape as data

        This property returns a new view (usually involving no copying of memory)
        of the `time` array, with additional dimensions to match the shape of the
        data.

        """
        new_shape = tuple(np.newaxis if i != self.time_axis else slice(None) for i in range(self.ndim))
        # This method is from the comment in np.broadcast_arrays, used so that
        # we don't have to disturb the data array (e.g., copying because of
        # weird ordering).
        # noinspection SpellCheckingInspection
        nditer = np.nditer(
            self.time[new_shape],
            flags=["multi_index", "zerosize_ok"],
            itershape=self.shape,
            order="C"
        )
        return nditer.itviews[0]

    def __repr__(self):
        r = repr(self.ndarray)
        return f"{type(self).__name__}{r[max(0, r.find('(')):-1]}, time={self.time!r}, time_axis={self.time_axis})"

    def __str__(self):
        return repr(self)

    def interpolate(self, new_time, derivative_order=0, out=None):
        """Interpolate this object to a new set of times

        Parameters
        ----------
        new_time : array_like
            Points to evaluate the interpolant at
        derivative_order : int, optional
            Order of derivative to evaluate.  If negative, the antiderivative is
            returned.  Default value of 0 returns the interpolated data without
            derivatives or antiderivatives.  Must be between -3 and 3, inclusive.

        See Also
        --------
        scipy.interpolate.CubicSpline :
            The function that this function is based on.
        antiderivative :
            Calls this funtion with `new_time=self.time` and
            `derivative_order=-antiderivative_order` (defaulting to a single
            antiderivative).
        derivative :
            Calls this function `new_time=self.time` and
            `derivative_order=derivative_order` (defaulting to a single derivative).
        dot :
            Property calling `self.derivative(1)`.
        ddot :
            Property calling `self.derivative(2)`.
        int :
            Property calling `self.antiderivative(1)`.
        iint :
            Property calling `self.antiderivative(2)`.

        Notes
        -----
        This function is essentially a wrapper around `scipy.interpolate.CubicSpline`

        """
        from scipy.interpolate import CubicSpline
        if derivative_order > 3:
            raise ValueError(
                f"{type(self).__name__} interpolation uses CubicSpline, which cannot take a derivative "
                f"of order {derivative_order}."
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
        spline = CubicSpline(self.time, self.ndarray, axis=self.time_axis)
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
        """Integrate modes with respect to time

        Parameters
        ----------
        antiderivative_order : int, optional
            Order of antiderivative to evaluate.  Default value is 1.  Must be between
            -3 and 3, inclusive.

        See Also
        --------
        scipy.interpolate.CubicSpline :
            The function that this function is based on.
        interpolate :
            This function simply calls `self.interpolate` with appropriate arguments.
        int :
            Property calling `self.antiderivative(1)`.
        iint :
            Property calling `self.antiderivative(2)`.

        """
        return self.interpolate(self.time, derivative_order=-antiderivative_order)

    def derivative(self, derivative_order=1):
        """Differentiate modes with respect to time

        Parameters
        ----------
        derivative_order : int, optional
            Order of derivative to evaluate.  Default value is 1.  Must be between -3
            and 3, inclusive.

        See Also
        --------
        scipy.interpolate.CubicSpline :
            The function that this function is based on.
        interpolate :
            This function simply calls `self.interpolate` with appropriate arguments.
        dot :
            Property returning `self.derivative(1)`.
        ddot :
            Property returning `self.derivative(2)`.

        """
        return self.interpolate(self.time, derivative_order=derivative_order)

    @property
    def dot(self):
        """Differentiate modes once with respect to time

        See Also
        --------
        derivative : This property simply returns `self.derivative(1)`
        ddot : Property returning `self.derivative(2)`.

        """
        return self.derivative()

    @property
    def ddot(self):
        """Differentiate modes twice with respect to time

        See Also
        --------
        derivative : This property simply returns `self.derivative(2)`
        dot : Property returning `self.derivative(1)`.

        """
        return self.derivative(2)

    @property
    def int(self):
        """Integrate modes once with respect to time

        See Also
        --------
        antiderivative : This property simply returns `self.antiderivative(1)`
        iint : Property returning `self.antiderivative(2)`.

        """
        return self.antiderivative()

    @property
    def iint(self):
        """Integrate modes twice with respect to time

        See Also
        --------
        antiderivative : This property simply returns `self.antiderivative(2)`
        int : Property returning `self.antiderivative(1)`.

        """
        return self.antiderivative(2)

    def xor(self, reverse=False, preserve_dtype=False, **kwargs):
        """Progressively XOR data along the time axis

        This function steps through an array, starting with the second element, and
        evaluates bitwise XOR on that element and the preceding one.  This is a useful
        step in compressing reasonably continuous data.

        See the documentation of `sxs.utilities.xor` for a full description of this
        function.  Note that this version sets the `axis` argument automatically to be
        the `time_axis`.

        """
        from .utilities import xor
        kw = kwargs.copy()
        kw.update({"axis": self.time_axis})
        return xor(self.ndarray, reverse=reverse, preserve_dtype=preserve_dtype, **kw)

    def truncate(self, abs_tolerance):
        """Truncate the precision of this object's `data` in place

        This function sets bits in the array data to 0 when they have lower
        significance than the number given as or returned by `abs_tolerance`.  This is
        a useful step in compressing data — though it is obviously lossy.

        Parameters
        ----------
        abs_tolerance : {callable, float, array-like}
            If callable, it is called with this object as the parameter, and the
            returned value is treated as a float or array-like would be.  Floats are
            simply treated as a uniform absolute tolerance to be applied at all times.
            Array-like objects must broadcast against this array, and each element is
            treated as the absolute tolerance for all the elements it broadcasts to.

        Returns
        -------
        None
            This value is returned to serve as a reminder that this function operates
            in place.

        Notes
        -----
        The effect is achieved by multiplying the array's data by the same power of 2
        that would be required to bring the `abs_tolerance` to between 1 and 2.  Thus,
        all digits less significant than 1 are less significant than `abs_tolerance` —
        meaning that we can apply the standard `round` routine to set these digits to
        0.  We then divide by that same power of 2 to bring the array data back to
        nearly its original value.  By working with powers of 2, we ensure that the 0s
        at the intermediate stage are represented as 0 bits in the final result.

        For floats and array-like objects, all values must be strictly positive, or
        `inf` or `nan` will result.

        """
        if callable(abs_tolerance):
            abs_tolerance = abs_tolerance(self)
        power_of_2 = (2.0 ** np.floor(-np.log2(abs_tolerance)))
        ndarray = self.ndarray
        ndarray *= power_of_2
        np.round(ndarray, out=ndarray)
        ndarray /= power_of_2
