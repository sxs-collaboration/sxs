"""Functions used to preprocess time-domain signals, including tail-transitioning,
padding, tapering, and detrending
"""

import numpy as np
from .taper_functions import taper_function

# Require the scri package in the dependency! Need to discuss!
from scri.utilities import transition_to_constant


def transition_tail_to_constant(time, data, transition_length=0.0):
    """Transition tails of time-domain signals to a constant. See the function
    scri.utilities.transition_to_constant for details.

    Parameters
    ----------
    time : array_like
        Array of time stamps; must have attribute .ndim (e.g., a numpy array)
    data : array_like, 1 or 2 dimensional
        Time-domain signal; data.shape[-1] must match len(time).
    transition_length : float
        Length of the signal's tail to be transitioned

    Returns
    -------
    data_transitioned : array_like
        Tail-transitioned signals; have the same shape as the data array
    """
    assert (
        transition_length <= time[-1] - time[0]
    ), "transition_length must be shorter than or equal to the length of the data."

    assert data.ndim <= 2, "The data array must be one or two dimensional."

    if transition_length == 0.0:
        return data

    if data.ndim == 1:
        data_transitioned = transition_to_constant(data, time, time[-1] - transition_length, time[-1])
    else:
        data_transitioned = np.array(
            [transition_to_constant(i, time, time[-1] - transition_length, time[-1]) for i in data]
        )

    return data_transitioned


def pad(time, data, pad_type="edge", left_pad_length=0.0, right_pad_length=0.0):
    """Pad time-domain signals.

    Parameters
    ----------
    time : array_like
        Array of time stamps; must have attribute .ndim (e.g., a numpy array)
    data : array_like, 1 or 2 dimensional
        Time-domain signal; data.shape[-1] must match len(time);
        must be equally spaced
    pad_type : str
        'edge' : pad signals by the end values
        'zero' : pad signals by zeros
    left_pad_length : float, optional
        Pad the left sides of signals to an additional length of left_pad_length
    right_pad_length : float, optional
        Pad the right sides of signals to an additional length of right_pad_length

    Returns
    -------
    time_padded : array_like
        Padded time array
    data_padded : array_like
        Padded signals; data_padded.shape[-1] matches len(time_padded);
        data_padded.ndim matches data.ndim.
    """
    if pad_type == "edge":
        pass
    elif pad_type == "zero":
        pad_type = "constant"
    else:
        raise ValueError("pad_type must be either 'edge' or 'zero'.")

    dt = time[1] - time[0]
    left_pad_num = int(np.ceil(left_pad_length / dt))
    right_pad_num = int(np.ceil(right_pad_length / dt))

    # pad data array
    if data.ndim == 1:
        data_padded = np.pad(data, (left_pad_num, right_pad_num), pad_type)
    elif data.ndim == 2:
        data_padded = np.pad(data, ((0, 0), (left_pad_num, right_pad_num)), pad_type)
    else:
        raise ValueError("The data array must be one or two dimensional.")

    # pad time array
    time_padded = np.pad(time, (left_pad_num, right_pad_num), "reflect", reflect_type="odd")

    return time_padded, data_padded


def taper(time, data, taper_type, sigmoid_type="smooth", **kwargs):
    """Taper time-domain signals. Also see taper_function.taper_function.

    Parameters
    ----------
    time : array_like
        Array of time stamps; must have attribute .ndim (e.g., a numpy array)
    data : array_like, 1 or 2 dimensional
        Time-domain signal; data.shape[-1] must match len(time).
    taper_type : str
        'left' : left tapering, i.e., an increasing sigmoid function; must specify
            either (both t0_taper and t1_taper) or left_taper_length
        'right' : right tapering, i.e., a decreasing sigmoid function; must specify
            either (both t2_taper and t3_taper) or right_taper_length
        'both' : tapering both sides, i.e., a window function; must specify the
            arguments required by both the 'left' and 'right' taper_type's
    sigmoid_type : str
        'smooth' : a smooth sigmoid function or a bump window
        'cosine' : a cosine sigmoid function or a Tukey window
    t0_taper : float, optional
        Lower bound of the transition time range for left tapering
    t1_taper : float, optional
        Upper bound of the transition time range for left tapering
    left_taper_length : float, optional
        Length of left tapering; equivalent to the setting of (t0_taper = time[0]
        and t1_taper = time[0]+left_taper_length)
    t2_taper : float, optional
        Lower bound of the transition time range for right tapering
    t3_taper : float, optional
        Upper bound of the transition time range for right tapering
    right_taper_length : float, optional
        Length of right tapering; equivalent to the setting of (t3_taper = time[-1]
        and t2_taper = time[-1]-right_taper_length)

    Returns
    -------
    data_tapered : array_like
        Tapered signals; have the same shape as the data array
    """
    return data * taper_function(time, taper_type, sigmoid_type, **kwargs)


def detrend(time, data, detrend_type):
    """Detrend time-domain signals to a constant, i.e., subtract signals by a
    line.

    Parameters
    ----------
    time : array_like
        Array of time stamps; must have attribute .ndim (e.g., a numpy array);
        must be equally spaced
    data : array_like, 1 or 2 dimensional
        Time-domain signal; data.shape[-1] must match len(time).
    detrend_type : str
        'none' : no detrending
        'zero right end' : bring the right ends of signals to 0; the line to be
            subtracted connects the points (time[0], 0) and (time[-1]+dt, data[:,-1])
        'zero both ends' : bring both ends of signals to 0; the line to be
            subtracted connects the points (time[0], data[:,0]) and
            (time[-1]+dt, data[:,-1])
        where dt := time[1]-time[0] is the time spacing

    Returns
    -------
    data_detrended : array_like
        Detrended signals; have the same shape as the data array

    Notes
    -----
    A detrended signal is assumed to be later Fourier transformed. Discrete Fourier
    transform assumes the signal is periodic and identifies the time stamps time[0]
    and time[-1]+dt (instead of time[-1]). We choose these two time stamps to be
    the anchor points of the subtraction line. We also assume the signal value
    at time[-1]+dt is the same as the value at time[-1], i.e., data[:,-1].
    """
    assert detrend_type in {
        "none",
        "zero right end",
        "zero both ends",
    }, "detrend_type must be 'none', 'zero right end', or 'zero both ends'."

    assert data.ndim <= 2, "The data array must be one or two dimensional."

    if detrend_type == "none":
        return data

    time_common_factor = (time - time[0]) / (time[-1] + time[1] - 2 * time[0])
    if data.ndim == 1:
        offset = 0 if detrend_type == "zero right end" else data[0]
        data_detrended = data - (data[-1] - offset) * time_common_factor - offset
    else:
        offset = np.zeros(data.shape[0]) if detrend_type == "zero right end" else data[:, 0]
        data_detrended = data - np.outer(data[:, -1] - offset, time_common_factor) - offset[:, None]

    return data_detrended


def preprocess_general(
    time,
    data,
    taper_type,
    detrend_type,
    transition_length=0.0,
    pad_type="edge",
    left_pad_length=0.0,
    right_pad_length=0.0,
    sigmoid_type="smooth",
    **kwargs,
):
    """Preprocess time-domain signals with full control. Preprocessing includes
    tail-transitioning, padding, tapering, and detrending.

    Signal Parameters
    -----------------
    time : array_like
        Array of time stamps; must have attribute .ndim (e.g., a numpy array)
    data : array_like, 1 or 2 dimensional
        Time-domain signal; data.shape[-1] must match len(time).

    Tail-transitioning Parameters
    -----------------------------
    transition_length : float
        Length of the signal's tail to be transitioned

    Padding Parameters
    ------------------
    pad_type : str
        'edge' : pad signals by the end values
        'zero' : pad signals by zeros
    left_pad_length : float, optional
        Pad the left sides of signals to an additional length of left_pad_length
    right_pad_length : float, optional
        Pad the right sides of signals to an additional length of right_pad_length

    Tapering Parameters
    -------------------
    taper_type : str
        'left' : left tapering, i.e., an increasing sigmoid function; must specify
            either (both t0_taper and t1_taper) or left_taper_length
        'right' : right tapering, i.e., a decreasing sigmoid function; must specify
            either (both t2_taper and t3_taper) or right_taper_length
        'both' : tapering both sides, i.e., a window function; must specify the
            arguments required by both the 'left' and 'right' taper_type's
    sigmoid_type : str
        'smooth' : a smooth sigmoid function or a bump window
        'cosine' : a cosine sigmoid function or a Tukey window
    t0_taper : float, optional
        Lower bound of the transition time range for left tapering
    t1_taper : float, optional
        Upper bound of the transition time range for left tapering
    left_taper_length : float, optional
        Taper from the left end of the PADDED signal to left_taper_length away
        from the left end. Equivalent to the setting of (t0_taper = time_padded[0]
        and t1_taper = time_padded[0]+left_taper_length)
    t2_taper : float, optional
        Lower bound of the transition time range for right tapering
    t3_taper : float, optional
        Upper bound of the transition time range for right tapering
    right_taper_length : float, optional
        Taper from the right end of the PADDED signal to right_taper_length away
        from the right end. Equivalent to the setting of (t3_taper = time_padded[-1]
        and t2_taper = time_padded[-1]-right_taper_length)

    Detrending Parameters
    ---------------------
    detrend_type : str
        'none' : no detrending
        'zero right end' : bring the right ends of signals to 0; the line to be
            subtracted connects the points (time[0], 0) and (time[-1]+dt, data[:,-1])
        'zero both ends' : bring both ends of signals to 0; the line to be
            subtracted connects the points (time[0], data[:,0]) and
            (time[-1]+dt, data[:,-1])
        where dt := time[1]-time[0] is the time spacing

    Returns
    -------
    data_transitioned : array_like
        Tail-transitioned signals; have the same shape as the data array
    time_padded : array_like
        Padded time array
    data_padded : array_like
        Padded signals; data_padded.shape[-1] matches len(time_padded);
        data_padded.ndim matches data.ndim.
    data_tapered : array_like
        Tapered signals; data_tapered.shape[-1] matches len(time_padded);
        data_tapered.ndim matches data.ndim.
    data_detrended : array_like
        Detrended signals; data_detrended.shape[-1] matches len(time_padded);
        data_detrended.ndim matches data.ndim.
    """
    data = np.array(data)
    assert data.ndim == 1 or data.ndim == 2, "The data array must be one or two dimensional."
    assert len(time) == data.shape[-1], "The data array's last dimension must have the same size as the time array."

    data_transitioned = transition_tail_to_constant(time, data, transition_length)
    time_padded, data_padded = pad(time, data_transitioned, pad_type, left_pad_length, right_pad_length)
    data_tapered = taper(time_padded, data_padded, taper_type, sigmoid_type, **kwargs)
    data_detrended = detrend(time_padded, data_tapered, detrend_type)

    return data_transitioned, time_padded, data_padded, data_tapered, data_detrended


def preprocess(time, data, preprocess_type, transition_length, pad_length, **kwargs):
    """Preprocess time-domain signals in the following steps: tail-transitioning,
    padding, tapering, and detrending. After tail-transitioning and padding, two
    tapering/detrending schemes are supported, specified by the preprocess_type
    argument:

    Scheme 1: tapering both sides of signals, followed by NO detrending
    Scheme 2: tapering the left sides of signals, followed by detrending that
              brings the right ends to 0.

    Signal Parameters
    -----------------
    time : array_like
        Array of time stamps; must have attribute .ndim (e.g., a numpy array)
    data : array_like, 1 or 2 dimensional
        Time-domain signal; data.shape[-1] must match len(time).
    preprocess_type : str
        'taper both sides' : Scheme 1; must specify
            [(both t0_taper and t1_taper) or left_taper_length] and
            [(both t2_taper and t3_taper) or right_taper_length]
        'taper left side and detrend' : Scheme 2; must specify either
            (both t0_taper and t1_taper) or left_taper_length
    transition_length : float
        Length of the signal's tail to be transitioned
    pad_length : float, optional
        Pad the right sides of signals to an additional length of pad_length
    t0_taper : float, optional
        Lower bound of the transition time range for left tapering
    t1_taper : float, optional
        Upper bound of the transition time range for left tapering
    left_taper_length : float, optional
        Taper from the left end of the PADDED signal to left_taper_length away
        from the left end. Equivalent to the setting of (t0_taper = time_padded[0]
        and t1_taper = time_padded[0]+left_taper_length)
    t2_taper : float, optional
        Lower bound of the transition time range for right tapering
    t3_taper : float, optional
        Upper bound of the transition time range for right tapering
    right_taper_length : float, optional
        Taper from the right end of the PADDED signal to right_taper_length away
        from the right end. Equivalent to the setting of (t3_taper = time_padded[-1]
        and t2_taper = time_padded[-1]-right_taper_length)

    Returns
    -------
    data_transitioned : array_like
        Tail-transitioned signals; have the same shape as the data array
    time_padded : array_like
        Padded time array
    data_padded : array_like
        Padded signals; data_padded.shape[-1] matches len(time_padded);
        data_padded.ndim matches data.ndim.
    data_tapered : array_like
        Tapered signals; data_tapered.shape[-1] matches len(time_padded);
        data_tapered.ndim matches data.ndim.
    data_detrended : array_like
        Detrended signals; data_detrended.shape[-1] matches len(time_padded);
        data_detrended.ndim matches data.ndim.
    """
    if preprocess_type == "taper both sides":
        return preprocess_general(
            time,
            data,
            transition_length=transition_length,
            taper_type="both",
            detrend_type="none",
            right_pad_length=pad_length,
            **kwargs,
        )
    elif preprocess_type == "taper left side and detrend":
        return preprocess_general(
            time,
            data,
            transition_length=transition_length,
            taper_type="left",
            detrend_type="zero right end",
            right_pad_length=pad_length,
            **kwargs,
        )
    else:
        raise ValueError("preprocess_type must be 'taper both sides' or 'taper left side and detrend'.")
