import numpy as np
from scipy.interpolate import CubicSpline
import spherical_functions


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

        # Get time array
        time = metadata.get("time", None)
        if time is None:
            raise ValueError("Time data must be specified as part of input TimeSeries or as keyword argument")
        time = np.asarray(time).view(float)
        metadata["time"] = time
        if time.ndim != 1:
            raise ValueError(f"Input `time` array must have exactly 1 dimension; it has {time.ndim}.")

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
        if abs(time_axis) >= input_array.ndim:
            raise ValueError(f"time_axis={time_axis} does not exist in input_array with {input_array.ndim} dimensions")
        if time_axis < 0:
            time_axis = input_array.ndim + time_axis
        if input_array.shape[time_axis] != time.size:
            raise ValueError(
                f"input_array.shape[time_axis]={input_array.shape[time_axis]} (with time_axis={time_axis}) "
                f"does not match time.size={time.size}"
            )
        metadata["time_axis"] = time_axis

        obj = input_array.view(cls)
        obj._metadata = metadata
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        super().__array_finalize__(obj)
        self._metadata = copy.copy(getattr(self, '_metadata', getattr(obj, '_metadata', {})))
        if "time" not in self._metadata:
            raise ValueError("Cannot create {type(self)} without `time` data.")

    def __getitem__(self, key):
        # Slice this thing like it's an instance of the parent class
        sliced = super().__getitem__(key)
        # Now slice the time data, if necessary
        if isinstance(tuple, key):
            if len(key) > self.time_axis:
                new_time = self.time[key[self.time_axis]]
            elif key[0] == Ellipsis and self.time_axis - self.ndim >=  1 - len(key):
                new_time = self.time[key[self.time_axis - self.ndim]]
            else:
                new_time = self.time
        elif self.time_axis == 0:
            new_time = self.time[key]
        else:
            new_time = self.time
        # Create the new sliced object along with its metadata
        metadata = self._metadata.copy()
        metadata.update(**getattr(sliced, '_metadata', {}))
        metadata["time"] = new_time
        return type(self)(sliced, **metadata)

    @property
    def ndarray(self):
        """View this array as a numpy ndarray"""
        return self.view(np.ndarray)

    @property
    def ndarray_float(self):
        """View this array as a numpy ndarray of real data

        This function returns a real-valued view (without copying) of the
        underlying data.  If the underlying data is complex, this function
        results in an extra dimension of size 2 (corresponding to the real and
        imaginary parts, respectively).  Because the `time_axis` of this object
        is normalized to a positive number, the `time_axis` remains correct
        even if the extra dimension is added.

        """
        if np.iscomplexobj(self.ndarray):
            return self.view(np.ndarray, dtype=float).reshape((-1, 2))
        return self.view(np.ndarray)

    @property
    def time(self):
        """Return an array of the time steps"""
        return self._metadata["time"]

    # Alias
    t = time

    @property
    def time_axis(self):
        return self._metadata["time_axis"]

    def interpolate(self, new_time, derivative_order=0, out=None):
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
        if np.iscomplexobj(self):
            result = result.view(dtype=float).reshape((-1, 2))
            result[:] = spline(new_time)
            result = result.view(dtype=self.dtype)
        else:
            result[:] = spline(new_time)
        metadata = self._metadata.copy()
        metadata["time"] = new_time
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
