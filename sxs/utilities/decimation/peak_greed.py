"""Modified greedy algorithm"""

def minimal_grid(x, y, tol=1e-6, error_scale=1.0, y_reference=None):
    """Modified greedy algorithm for building a reduced-order spline

    This function is a modification of the algorithm found in `greedy_spline`.
    Rather than choosing the single worst-approximated point to include in the new
    spline on each iteration (as in the standard greedy algorithm), this function
    chooses a set of points based on peaks in the error function.  To do that, it
    uses the scipy `find_peaks` function.

    Parameters
    ----------
    x : array_like[float]
        An ordered series of sampling points.
    y : array_like[float]
        Function values sampled at points given by `x`.
    tol : float, optional
        This function will iteratively increase the number of sampling points until
        the spline formed from those points is within this value of the input data
        at each input sample.  Defaults to 1e-6.
    error_scale : {float, array_like[float]}, optional
        Multiply the error value by this scale before finding peaks or comparing to
        `tol`.  While this can be achieved by rescaling `tol` when this is a float,
        it is also possible to pass an array to this function so that each
        interpolated value may be targeted with a different tolerance.  In that
        case, the input array must have the same size as `x` and `y`.  The default
        value is 1.0.
    y_reference : {None, array_like[float]}, optional
        Reference `y` values with which to compare.  By default, this is precisely
        `y`; an array similar to `y` can be passed instead, so that `y` will be
        used to construct the splines, but `y_reference` will be used to evaluate
        the errors.  This can be helpful, for example, if we are truncating the
        precision of `y` in the output, but want to ensure that the result is still
        within the given tolerance of the original (pre-truncation) `y`.

    Returns
    -------
    include_sample : array_like[bool]
        Boolean array of the same length as `x` determining whether or not the
        corresponding point should be included in a reduced-order spline.  Note
        that numpy arrays can be indexed with this array, so that
        `x[include_sample]` will give the set of knots and `y[include_sample]` will
        give the corresponding set of data values with which to build the
        reduced-order spline.

    """
    import numpy as np
    from scipy.interpolate import CubicSpline as spline
    from scipy.signal import find_peaks

    deg = 3  # This used to be a parameter, before switching to CubicSpline

    if y_reference is None:
        y_reference = y

    error_scale = np.asarray(error_scale)

    if np.ndim(error_scale) == 0:
        if error_scale != 1.0:
            tol /= error_scale[()]
        def error(ydiff):
            return ydiff
    else:
        def error(ydiff):
            return error_scale * ydiff

    # Implementation notes
    # --------------------
    #
    # The following function is the core of this algorithm.  There are reasonable
    # different choices for the parameters passed to `find_peaks`.  In particular
    # the "prominence" argument seems like it should help in some situations by
    # ensuring that only peaks that actually stand out above their neighbors by
    # `tol` will be added.  However, tests on our waveforms suggest that this only
    # saves ~4% in file size, while slowing things down by ~20%.  I also worry that
    # this could leave off peaks from the final product that should have been
    # included to actually get within the tolerances.  The solution to that would
    # be to start with prominence=tol, and then switch to prominence=None for
    # safety.  Again, however, tests suggest that this has precisely no effect on
    # file sizes (meaning that no more points get added on the second round), but
    # slows down the whole thing by ~5%.  I suspect that the reason this is not
    # very helpful is specific to our data, because so many of the peaks are *very*
    # sharp and very large (mostly coming from junk radiation, and made worse by
    # extrapolation).  In any case, I think the tradeoff is in favor of ignoring
    # the prominence for our use case.
    #
    # Also, the `find_peaks` function has a `height` argument that could select
    # only peaks outside of the desired tolerance (i.e., doing essentially what
    # this code does in the lines following each call to `find_peaks`).  Tests show
    # that using that argument somehow actually increases file size by ~0.2%, while
    # slowing down the code by ~7%.  Evidently, I don't understand precisely what
    # that argument does.  But I am confident that the code below does what we want
    # -- in addition to being faster.

    def next_sample(y, y_greedy, sign):
        if sign[0] == 1:
            errors = error(y - y_greedy)
        else:
            errors = error(y_greedy - y)
        peaks = find_peaks(errors)[0]
        peaks = peaks[np.abs(errors[peaks])>tol]
        sign[0] *= -1
        if not peaks.size:
            peaks = find_peaks(-errors)[0]
            peaks = peaks[np.abs(errors[peaks])>tol]
            sign[0] *= -1
            if not peaks.size:
                return None
        return peaks

    # Create an array describing whether or not to include a given sample in the
    # spline.  Start with all False.
    include_sample = np.zeros(len(x), dtype=bool)

    # Initialize greedy algorithm with evenly spaced indices, as many
    # as needed for the desired spline degree.
    include_sample[np.linspace(0, len(x)-1, num=deg+1, dtype=int)] = True

    # Peak-greed algorithm
    sign = [1]
    for _ in range(len(x)):
        # Compute the spline interpolant on the current set of knots
        s = spline(x[include_sample], y[include_sample])

        # Evaluate this spline on all input data points and find peak errors
        i_next = next_sample(y_reference, s(x), sign)

        # Break out of this loop if `tol` is satisfied
        if i_next is None:
            break

        # Include data point that gives largest interpolation error
        include_sample[i_next] = True

    return include_sample
