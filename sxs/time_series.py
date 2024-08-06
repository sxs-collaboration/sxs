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
        self._metadata = copy.copy(getattr(self, "_metadata", getattr(obj, "_metadata", {})))
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
        metadata.update(**getattr(new_data, "_metadata", {}))
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

    @property
    def abs(self):
        """Absolute value of the data

        Returns
        -------
        absolute : TimeSeries
            Because the absolute values make no sense as mode weights, this is just a
            plain TimeSeries object.

        See Also
        --------
        arg

        """
        return np.abs(self)

    @property
    def arg(self):
        """Complex phase angle of the data

        Note that the result is not "unwrapped", meaning that there may be
        discontinuities as the phase approaches ±π.

        Returns
        -------
        phase : TimeSeries
            Values are in the interval (-π, π].

        See Also
        --------
        numpy.angle
        arg_unwrapped

        """
        return np.angle(self)

    @property
    def arg_unwrapped(self):
        """Complex phase angle of the data, unwrapped along the time axis

        The result is "unwrapped", meaning that discontinuities as the phase approaches
        ±π are removed by adding an appropriate amount to all following data points.

        Returns
        -------
        phase : TimeSeries
            Values at the initial time are in the interval (-π, π], but may evolve to
            arbitrary real values.

        See Also
        --------
        numpy.angle
        numpy.unwrap
        arg

        """
        return TimeSeries(np.unwrap(self.arg, axis=self.time_axis), self.time)

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

    def interpolate(self, new_time, derivative_order=0, out=None, padding_points=10):
        """Interpolate this object to a new set of times

        Parameters
        ----------
        new_time : array_like
            Points to evaluate the interpolant at
        derivative_order : int, optional
            Order of derivative to evaluate.  If negative, the
            antiderivative is returned.  Default value of 0 returns
            the interpolated data without derivatives or
            antiderivatives.  Must be between -3 and 3, inclusive.
        out : array_like, optional
            If provided, the result will be placed into this array.
            It must have the same shape as the result would have been
            if it was not provided.
        padding_points : int, optional
            Number of points to include on either side of `new_time`
            when constructing interpolating spline, to ensure that the
            interpolant has enough data to interpolate the whole
            range.  Default value is 10.  If there are not enough
            points in the input data to supply all of these padding
            points, the function will just come as close as possible.

        See Also
        --------
        scipy.interpolate.CubicSpline :
            The function that this function is based on.
        antiderivative :
            Calls this function with `new_time=self.time` and
            `derivative_order=-antiderivative_order` (defaulting to a
            single antiderivative).
        derivative :
            Calls this function `new_time=self.time` and
            `derivative_order=derivative_order` (defaulting to a
            single derivative).
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
        This function is essentially a wrapper around
        `scipy.interpolate.CubicSpline`

        """
        from scipy.interpolate import CubicSpline

        if derivative_order > 3:
            raise ValueError(
                f"{type(self).__name__} interpolation uses CubicSpline, which cannot take a derivative "
                f"of order {derivative_order}."
            )

        # Sort out shapes and storage for output
        new_time = np.asarray(new_time)
        if new_time.ndim != 1:
            raise ValueError(f"New time array must have exactly 1 dimension; it has {new_time.ndim}.")
        new_shape = list(self.shape)
        new_shape[self.time_axis] = new_time.size
        if out is not None:
            out = np.asarray(out)
            if out.shape != new_shape:
                raise ValueError(
                    f"Output array should have shape {new_shape} for consistency with new time "
                    "array and modes array"
                )
            if out.dtype != self.dtype:
                raise ValueError(
                    f"Output array should have same dtype as this array {self.dtype}; "
                    f"it has dtype {out.dtype}"
                )
        result = out or np.empty(new_shape, dtype=self.dtype)

        # Find a range of indices that includes `new_time`, plus
        # `padding_points` points in either direction, to ensure that
        # the spline has enough data to interpolate the whole range.
        idx_min = max(0, self.index_closest_to(new_time.min()) - padding_points)
        idx_max = min(self.n_times, self.index_closest_to(new_time.max()) + padding_points)
        idx = range(idx_min, idx_max)

        # Construct the spline — possibly including derivatives or antiderivatives
        spline = CubicSpline(
            self.time[idx],
            self.ndarray.take(idx, axis=self.time_axis),
            axis=self.time_axis
        )
        if derivative_order < 0:
            spline = spline.antiderivative(-derivative_order)
        elif 0 < derivative_order <= 3:
            spline = spline.derivative(derivative_order)

        # Evaluate the spline at the new time points
        result[:] = spline(new_time)

        # Copy metadata and return the result
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

        This function sets bits in the array data to 0 when they have
        lower significance than the number given as or returned by
        `abs_tolerance`.  This is a useful step in compressing data —
        though it is obviously lossy.

        Parameters
        ----------
        abs_tolerance : {callable, float, array-like}
            If callable, it is called with this object as the
            parameter, and the returned value is treated as a float or
            array-like would be.  Floats are simply treated as a
            uniform absolute tolerance to be applied at all times.
            Array-like objects must broadcast against this array, and
            each element is treated as the absolute tolerance for all
            the elements it broadcasts to.

        Returns
        -------
        None
            This value is returned to serve as a reminder that this
            function operates in place.

        Notes
        -----
        The effect is achieved by multiplying the array's data by the
        same power of 2 that would be required to bring the
        `abs_tolerance` to between 1 and 2.  Thus, all digits less
        significant than 1 are less significant than `abs_tolerance` —
        meaning that we can apply the standard `round` routine to set
        these digits to 0.  We then divide by that same power of 2 to
        bring the array data back to nearly its original value.  By
        working with powers of 2, we ensure that the 0s at the
        intermediate stage are represented as 0 bits in the final
        result.

        For floats and array-like objects, all values must be strictly
        positive, or `inf` or `nan` will result.

        """
        if callable(abs_tolerance):
            abs_tolerance = abs_tolerance(self)
        power_of_2 = (2.0 ** np.floor(-np.log2(abs_tolerance)))
        ndarray = self.ndarray
        ndarray *= power_of_2
        np.round(ndarray, out=ndarray)
        ndarray /= power_of_2

    def taper(self, t1, t2, y1=0, y2=1, *, transition_type="smooth"):
        """Smoothly taper (or transition) the data

        This is essentially half of a standard window function, using
        the C^∞ "compactified tanh" transition function.  By default,
        this function smoothly transitions the data from 0 before `t1`
        to the input data after `t2`, which is useful to prevent Gibbs
        effects caused by a sudden "turn on" of a signal.  By swapping
        the values of `y1` and `y2`, the transition can be reversed.

        Parameters
        ----------
        t1 : float
            Time before which the output will equal `y1` times the
            data.
        t2 : float
            Time after which the output will equal `y2` times the
            data.
        y1 : float, optional
            Value before `t1`.  Default value is 0.
        y2 : float, optional
            Value after `t2`.  Default value is 1.
        transition_type : str, optional
            Type of transition to apply.  Default value is `"smooth"`,
            which corresponds to the C^∞ bump (compactified tanh)
            function.  The other option is `"cosine"`, which
            corresponds to the Tukey (cosine-taper) function.

        Returns
        -------
        TimeSeries
            New object with transitioned data.

        See Also
        --------
        transition_to_constant : Smoothly transition to a constant
            value
        sxs.utilities.transition_function : Function that does the
            work

        """
        from .utilities import transition_function
        t2 = t2 or self.time[-1]
        if transition_type == "smooth":
            τ = transition_function(self.time, t1, t2, y1, y2)
        elif transition_type == "cosine":
            τ = np.empty_like(self.time)
            τ[self.time <= t1] = y1
            τ[(self.time > t1) & (self.time < t2)] = (y1 + y2) / 2 - (y2 - y1) / 2 * np.cos(np.pi * (self.time[(self.time > t1) & (self.time < t2)] - t1) / (t2 - t1))
            τ[self.time >= t2] = y2
        else:
            raise ValueError(f"Unknown transition type {transition_type=}")
        data = np.moveaxis(self.ndarray.copy(), self.time_axis, -1)
        data = data * τ
        data = np.moveaxis(data, -1, self.time_axis)
        return type(self)(data, **self._metadata.copy())

    def transition_to_constant(self, t1, t2=None):
        """Smoothly transition from the data to a constant

        This function produces a copy of the input, where the data is
        smoothly transitioned from the value at `t1` to a constant
        value at `t2` (which is the last time step if not given).  The
        precise value of this constant will depend on the behavior of
        the data between `t1` and `t2`.

        Parameters
        ----------
        t1 : float
            Time after which the output will equal the input data.
        t2 : float, optional
            Time after which the output will be constant.  Default
            value is the last time step.

        Returns
        -------
        TimeSeries
            New object with transitioned data.

        See Also
        --------
        taper : Smoothly transition the data
        sxs.utilities.transition_to_constant : Function that does
            the work
        
        """
        from .utilities import transition_to_constant_inplace
        t2 = t2 or self.time[-1]
        if self.time_axis != 0:
            data = np.moveaxis(self.ndarray.copy(), self.time_axis, 0)
        else:
            data = self.ndarray.copy()
        data = transition_to_constant_inplace(data, self.time, t1, t2)
        if self.time_axis != 0:
            data = np.moveaxis(data, 0, self.time_axis)
        return type(self)(data, **self._metadata.copy())
    
    def window(self, t1, t2, t3, t4, y1=0, y23=1, y4=None, *, window_type="smooth"):
        """Window the data smoothly

        This function creates a copy of the input and, by default,
        zeroes it outside of the times `t1` and `t4` while reproducing
        it precisely between `t2` and `t3`.

        Parameters
        ----------
        t1 : float
            Time before which the output will be multiplied by `y1`.
        t2 : float
            Time after which the output will be multiplied by `y23`.
        t3 : float
            Time before which the output will be multiplied by `y23`.
        t4 : float
            Time after which the output will be multiplied by `y4`.
        y1 : float, optional
            Multiplicative value before `t1`.  Default value is 0.
        y23 : float, optional
            Multiplicative value between `t2` and `t3`.  Default value
            is 1.
        y4 : float, optional
            Multiplicative value after `t4`.  Default value is `y1`.
        window_type : str, optional
            Type of window to apply.  Default value is `"smooth"`,
            which corresponds to the C^∞ bump (compactified tanh)
            function.  The other option is `"cosine"`, which
            corresponds to the Tukey (cosine-tapered) window.

        Returns
        -------
        TimeSeries
            New object with windowed data.

        See Also
        --------
        sxs.utilities.bump_function : Function that does the work

        """
        from .utilities import bump_function
        y4 = y4 or y1
        if window_type == "smooth":
            τ = bump_function(self.time, t1, t2, t3, t4, y1, y23, y4)
        elif window_type == "cosine":
            t = self.time
            τ = np.empty_like(t)
            τ[t <= t1] = y1
            τ[(t > t1) & (t < t2)] = (y1 + y23) / 2 - (y23 - y1) / 2 * np.cos(np.pi * (t[(t > t1) & (t < t2)] - t1) / (t2 - t1))
            τ[(t >= t2) & (t <= t3)] = y23
            τ[(t > t3) & (t < t4)] = (y23 + y4) / 2 - (y4 - y23) / 2 * np.cos(np.pi * (t[(t > t3) & (t < t4)] - t3) / (t4 - t3))
            τ[t >= t4] = y4
        else:
            raise ValueError(f"Unknown window type {window_type=}")
        data = np.moveaxis(self.ndarray.copy(), self.time_axis, -1)
        data = data * τ
        data = np.moveaxis(data, -1, self.time_axis)
        return type(self)(data, **self._metadata.copy())
    
    def pad(self, pad_length=None, mode="edge", **kwargs):
        """Pad the data values along the time axis
        
        This function is based on `numpy.pad`, but the `pad_length`
        argument is given in units of time, rather than numbers of
        elements, and the `mode` argument defaults to `"edge"`.

        As with `numpy.pad`, the `pad_length` argument may be a single
        number, a tuple with one element, or a tuple with two
        elements.  For just a single number or tuple with one element,
        the padding is applied symmetrically to both ends of the data.
        Unlike `numpy.pad`, there may only be two elements of this
        tuple; padding other dimensions is not allowed.

        Parameters
        ----------
        pad_length : {None, float, tuple}
            Amount of time to pad on both sides (for a single float)
            or the beginning and end (for a tuple) of the data.  By
            default, the full length of the input data is added to
            both the beginning and end of the data; that is, the
            length of the data is tripled.
        mode : str, optional
            Padding mode.  Default value is `"edge"`.
        kwargs : dict
            Additional keyword arguments to pass to `numpy.pad`.

        Returns
        -------
        TimeSeries
            New object with padded data.

        See Also
        --------
        numpy.pad : Similar function with padding in number of
            elements, rather than time, and a different default
            `mode`.
        
        """
        if pad_length is None:
            pad_length = self.time[-1] - self.time[0]
        pad_length = np.asarray(pad_length)
        if pad_length.shape == ():
            pad_length = np.array([pad_length, pad_length])
        elif pad_length.shape == (1,):
            pad_length = np.array([pad_length[0], pad_length[0]])
        elif pad_length.shape != (2,):
            raise ValueError(
                f"Input `pad_length` must be a scalar, a tuple with one element, "
                f"or a tuple with two elements; it has shape {pad_length.shape}."
            )
        dt0 = self.time[1] - self.time[0]
        dt1 = self.time[-1] - self.time[-2]
        pad_width0 = int(np.ceil(pad_length[0] / dt0))
        pad_width1 = int(np.ceil(pad_length[1] / dt1))
        pad_width = list(
            (0, 0) if self.time_axis != i else (pad_width0, pad_width1)
            for i in range(self.ndarray.ndim)
        )
        data = np.pad(self.ndarray, pad_width, mode, **kwargs)
        metadata = self._metadata.copy()
        metadata["time"] = np.concatenate(
            (
                self.time[0] - np.arange(pad_width0, 0, -1) * dt0,
                self.time,
                self.time[-1] + np.arange(1, pad_width1+1) * dt1
            )
        )
        return type(self)(data, **metadata)

    def line_subtraction(self, treat_as_zero="begin"):
        """Subtract a linear function of time from the data

        This is very similar to the `detrend` function found in the
        signal-processing literature and in SciPy and MATLAB, for
        example, except that rather than fitting a line to all of the
        data (which is more appropriate for noisy data), this function
        subtracts the line connecting the two ends (first and last
        time steps) of the data.  This ensures that there is no
        discontinuity in the data at the boundaries, which is not
        guaranteed by the usual `detrend` functions.

        Parameters
        ----------
        treat_as_zero : str, optional
            Whether to treat the data at the boundaries as zero.
            Default value is `"begin"`, meaning that we pretend that
            the data at the beginning of the time array is oscillating
            around zero, so we can treat it as if it *is* zero.  If
            the early data has already been tapered, it will probably
            already be zero, so this shouldn't matter.  In that case,
            `"neither"` would also be correct.  The option `"end"` is
            also accepted, but this will be less likely to be needed.

        Returns
        -------
        TimeSeries
            New object with the linear function of time subtracted.

        Notes
        -----
        There is a minor subtlety in the implementation of this
        function.  It is not quite enough to simply subtract the line
        connecting the first and last data points.  The discrete
        Fourier transform implicitly assumes that the data is
        periodic, so that the line should connect the first data point
        to *one beyond the last data point*.  It is not generally
        clear how to define that, but assuming that the data has
        roughly settled down to be roughly constant by the end, we can
        say that the line should connect `(time[0], data[0])` to
        `(time[-1]+dt, data[-1])`, where `dt` is the time-step size.
        This is the approach taken here.  Even this is only
        approximate, but typically this should be such a tiny effect
        as to not matter at all, though it could improve results for
        very short or highly changing time series.

        """
        N = self.n_times
        data = np.moveaxis(self.ndarray.copy(), self.time_axis, -1)
        data_begin = 0 if treat_as_zero == "begin" else data[..., [0]]
        data_end = 0 if treat_as_zero == "end" else data[..., [-1]]
        line = (
            data_begin
            + (data_end - data_begin)
            * ((self.time - self.time[0]) / (self.time[-1] - self.time[0]))
            * ((N - 1) / N)
        )
        data -= line
        data = np.moveaxis(data, -1, self.time_axis)
        return type(self)(data, **self._metadata.copy())