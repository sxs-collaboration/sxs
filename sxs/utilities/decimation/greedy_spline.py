"""Greedy algorithm"""

def minimal_grid(x, y, tol=1e-6, rel=False):
    """Greedy algorithm for building a reduced-order spline

    This function is a slight generalization of the core function from `romspline`
    -- namely, `romspline.greedy._greedy` -- tailored to the needs of this package.
    This function should behave nearly the same as the `romspline` function (except
    that this function returns fewer quantities), but also accepts a callable
    object as the `tol` parameter, to enable more general evaluations of the error.

    Parameters
    ----------
    x : array of float
        An ordered series of sampling points.
    y : array of float
        Function values sampled at points given by `x`.
    tol : {float, callable}, optional
        If given a float, this function will iteratively increase the number of
        sampling points until the spline formed from those points is within this
        value of the input data at each input sample.  If given a callable object,
        that object must be as described in the Notes section below.  The default
        value is 1e-6.
    rel : bool, optional
        If `tol` is a float and this is True, the tolerance will be taken relative
        to the maximum absolute value in the input `y` data.  If False (the
        default), the tolerance will be treated as an absolute value.  Note that if
        `tol` is callable this parameter is ignored.

    Returns
    -------
    include_sample : array of bool
        Boolean array of the same length as `x` determining whether or not the
        corresponding point should be included in a reduced-order spline.  Note
        that numpy arrays can be indexed with this array, so that
        `x[include_sample]` will give the set of knots and `y[include_sample]` will
        give the corresponding set of data values with which to build the
        reduced-order spline.

    Notes
    -----
    The `tol` callable, if given, is used to decide whether or not the current
    trial spline is acceptable, and if not which sample point to include for the
    next trial spline.  It must be callable with the signature (x, y, y_greedy),
    and is expected to return either None (if the trial is acceptable) or the index
    of the next sample point to include.  Here, `x` and `y` are just the
    corresponding inputs to this function, and `y_greedy` is a set of data points
    representing the values of the trial reduced-order spline evaluated on `x`.
    For example, the following function would behave equivalently to the default
    parameters:

        def tol(x, y, y_greedy):
            errors = np.abs(y-y_greedy)
            i_next = np.argmax(errors)
            if errors[i_next] < 1e-6:
                return None
            return i_next

    """
    import numpy as np
    from scipy.interpolate import CubicSpline as spline

    deg = 3  # This used to be a parameter, before switching to CubicSpline

    if callable(tol):
        next_sample = tol
    else:
        if rel:
            # Rescale the tolerance
            ymax = np.max(np.abs(y))
            if ymax == 0.0:
                raise ValueError("All input `y` data samples are zero.")
            tol *= ymax

        # Construct the function used to evaluate the greedy spline
        def next_sample(x, y, y_greedy):
            errors = np.abs(y-y_greedy)
            i_next = np.argmax(errors)
            if errors[i_next] < tol:
                return None
            return i_next

    # Create an array describing whether or not to include a given
    # sample in the spline.  Start with all False.
    include_sample = np.zeros(len(x), dtype=bool)

    # Initialize greedy algorithm with evenly spaced indices, as many
    # as needed for the desired spline degree.
    include_sample[np.linspace(0, len(x)-1, num=deg+1, dtype=int)] = True

    # Greedy algorithm
    for _ in range(len(x)):
        # Spline interpolant on current set of knots
        s = spline(x[include_sample], y[include_sample])

        # Evaluate this spline
        i_next = next_sample(x, y, s(x))

        # Break out of this loop if `tol` is satisfied
        if i_next is None:
            break

        # Include data point that gives largest interpolation error
        include_sample[i_next] = True

    return include_sample
