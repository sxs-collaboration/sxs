def minimal_grid(x, y, tol=1e-6, error_scale=1.0):
    """Modified greedy algorithm for building a reduced-order spline

    This function is a modification of the algorithm found in
    `greedy_spline`.  Rather than choosing the single
    worst-approximated point to include in the new spline on each
    iteration (as in the standard greedy algorithm), this function
    chooses a set of points based on peaks in the error function.  To
    do that, it uses the scipy `find_peaks` function.

    Parameters
    ==========
    x: array of float
        An ordered series of sampling points.
    y: array of float
        Function values sampled at points given by `x`.
    tol: float [defaults to 1e-6]
        This function will iteratively increase the number of
        sampling points until the spline formed from those points is
        within this value of the input data at each input sample.
    error_scale: float or array [defaults to 1.0]
        Multiply the error value by this scale before finding peaks
        or comparing to `tol`.  With this can be achieved by
        rescaling `tol` when this is a float, it is also possible to
        pass an array to this function so that each interpolated
        value may be targeted with a different tolerance.  In that
        case, the input array must have the same size as `x` and `y`.

    Returns
    =======
    include_sample: array of bool
        Boolean array of the same length as `x` determining whether
        or not the corresponding point should be included in a
        reduced-order spline.  Note that numpy arrays can be indexed
        with this array, so that `x[include_sample]` will give the
        set of knots and `y[include_sample]` will give the
        corresponding set of data values with which to build the
        reduced-order spline.

    """
    import numpy as np
    from scipy.interpolate import CubicSpline as spline
    from scipy.signal import find_peaks

    deg = 3  # This used to be a parameter, before switching to CubicSpline

    error_scale = np.asarray(error_scale)

    if np.ndim(error_scale) == 0:
        if error_scale != 1.0:
            tol /= error_scale[()]
        def next_sample(y, y_greedy, sign, prominence):
            if sign[0] == 1:
                errors = y - y_greedy
            else:
                errors = y_greedy - y
            peaks = find_peaks(errors, height=tol, prominence=prominence)[0]
            sign[0] *= -1
            if not peaks.size:
                peaks = find_peaks(-errors, height=tol, prominence=prominence)[0]
                sign[0] *= -1
                if not peaks.size:
                    return None
            return peaks
    else:
        def next_sample(y, y_greedy, sign, prominence):
            if sign[0] == 1:
                errors = error_scale * (y - y_greedy)
            else:
                errors = error_scale * (y_greedy - y)
            peaks = find_peaks(errors, height=tol, prominence=prominence)[0]
            sign[0] *= -1
            if not peaks.size:
                peaks = find_peaks(-errors, height=tol, prominence=prominence)[0]
                sign[0] *= -1
                if not peaks.size:
                    return None
            return peaks

    # Create an array describing whether or not to include a given
    # sample in the spline.  Start with all False.
    include_sample = np.zeros(len(x), dtype=bool)

    # Initialize greedy algorithm with evenly spaced indices, as many
    # as needed for the desired spline degree.
    include_sample[np.linspace(0, len(x)-1, num=deg+1, dtype=int)] = True

    # Peak-greed algorithm
    sign = [1]
    for prominence in [tol, None]:
        for _ in range(len(x)):
            # Spline interpolant on current set of knots
            s = spline(x[include_sample], y[include_sample])

            # Evaluate this spline
            i_next = next_sample(y, s(x), sign, prominence)

            # Break out of this loop if `tol` is satisfied
            if i_next is None:
                break

            # Include data point that gives largest interpolation error
            include_sample[i_next] = True
        
    return include_sample
