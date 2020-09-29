def suppressor(data, tolerance, pad=0):
    """Identify elements of data close to 0

    Parameters
    ----------
    data : array
        This can be either real or complex data.
    tolerance : float
        Data of absolute value below this tolerance will be set to zero.
    pad : int, optional
        The function will only suppress data that has at least `pad` points between
        it and the nearest value that is above the tolerance.  Defaults to 0.

    """
    import numpy as np
    # Start off suppressing everything below tolerance
    suppressed = np.abs(data) < tolerance
    if pad > 0:
        # Count how many neighbors (within `pad` points to either side) are also suppressed
        counter = np.ones(2*pad+1, dtype=int)
        count_neighboring_suppressions = np.convolve(suppressed, counter, mode='same')
        # Only suppress those whose neighbors are all suppressed
        suppressed = count_neighboring_suppressions == np.sum(counter)
    return suppressed


def suppress(data, tolerance, pad=0, inplace=True):
    """Set data close to 0 to exactly 0

    Parameters
    ----------
    data : array
        This can be either real or complex data.
    tolerance : float
        Data of absolute value below this tolerance will be set to zero.
    pad : int, optional
        The function will only suppress data that has at least `pad` points between
        it and the nearest value that is above the tolerance.  Defaults to 0.
    inplace : bool, optional
        If True (the default), overwrite the `data` array and return it; if False,
        copy it, suppress small numbers, and return the copy.

    """
    import numpy as np
    if inplace:
        data[suppressor(data, tolerance, pad=pad)] = 0.0
        output = data
    else:
        output = np.copy(data)
        output[suppressor(data, tolerance, pad=pad)] = 0.0
    return output
