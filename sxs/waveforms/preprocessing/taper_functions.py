"""Functions used for tapering time-domain signals
"""

import numpy as np


def smooth_sigmoid(t, t_lower=0., t_upper=1., decreasing=False):
    """A smooth sigmoid function that maps [t_lower, t_upper] to [0, 1]. If the 
    function is increasing, then f(t_lower)=0 and f(t_upper)=1; if decreasing, 
    then f(t_lower)=1 and f(t_upper)=0. The function is smooth: derivatives at 
    all orders are continuous. It is required that t_lower < t_upper and both 
    t_lower and t_upper lie within the time array t. 

    Parameters
    ----------
    t : array_like
        Array of time stamps
    t_lower : float, optional
        Lower bound of the transition time range
    t_upper : float, optional
        Upper bound of the transition time range
    decreasing : boolean, optional
        If True, the sigmoid decreases from 1 to 0; if False, it increases from 
        0 to 1.

    Returns
    -------
    sigmoid_function : array_like
        Values of the sigmoid function

    Notes
    -----
    An increasing smooth sigmoid function has the form 1/[1+exp(Z)], while a 
    decreasing smooth sigmoid function has the form 1/[1+exp(-Z)], where
    Z := 1/T-1/(1-T) and T := (t-t_lower)/(t_upper-t_lower). The implementation 
    below is similar to the function planck_window in the package 
    sxs-collaboration/surrogate_modeling.
    """
    if decreasing:
        return smooth_sigmoid(-t, t_lower=-t_upper, t_upper=-t_lower)

    if t_lower == 0. and t_upper == 1.:
        assert min(t[0], t[-1]) <= 0. and max(t[0], t[-1]) >= 1., \
            "t_lower and t_upper must lie within the time array."
        # To avoid zero division, we calculate the transition in the range
        # [t_lower=tol, t_upper=1-tol] instead of [t_lower=0, t_upper=1].
        # A tol of .005 suffices: for T <= tol, exp(Z) > 1e86; for T <= 1-tol,
        # exp(Z) < 1e-86.
        tol = .005
        transition_zone = (t > tol) * (t < (1.-tol))

        # Replace those t values that lie outside the transition_zone by 0.5 to
        # avoid zero division. This does not affect the final result.
        safe_t = transition_zone * t + (1. - transition_zone) * .5

        safe_exponent = 1./safe_t - 1./(1. - safe_t)
        return transition_zone/(1. + np.exp(safe_exponent)) + (t >= (1.-tol))

    assert t_lower < t_upper, "t_lower must be smaller than t_upper."
    return smooth_sigmoid((t-t_lower)/(t_upper-t_lower))


def cosine_sigmoid(t, t_lower=0., t_upper=1., decreasing=False):
    """A cosine sigmoid function that maps [t_lower, t_upper] to [0, 1]. If the 
    function is increasing, then f(t_lower)=0 and f(t_upper)=1; if decreasing, 
    then f(t_lower)=1 and f(t_upper)=0. The second derivative of this function 
    is NOT continuous. It is required that t_lower < t_upper and both 
    t_lower and t_upper lie within the time array t. 

    Parameters
    ----------
    t : array_like
        Array of time stamps
    t_lower : float, optional
        Lower bound of the transition time range
    t_upper : float, optional
        Upper bound of the transition time range
    decreasing : boolean, optional
        If True, the sigmoid decreases from 1 to 0; if False, it increases from 
        0 to 1.

    Returns
    -------
    sigmoid_function : array_like
        Values of the sigmoid function

    Notes
    -----
    An increasing smooth sigmoid function has the form 0.5*[1-cos(pi*T)], while 
    a decreasing smooth sigmoid function has the form 0.5*[1+cos(pi*T)], where
    T := (t-t_lower)/(t_upper-t_lower).
    """
    if decreasing:
        return cosine_sigmoid(-t, t_lower=-t_upper, t_upper=-t_lower)

    if t_lower == 0. and t_upper == 1.:
        assert min(t[0], t[-1]) <= 0. and max(t[0], t[-1]) >= 1., \
            "t_lower and t_upper must lie within the time array."
        transition_zone = (t > 0.) * (t < 1)
        return transition_zone * 0.5 * (1 - np.cos(np.pi*t)) + (t >= 1.)

    assert t_lower < t_upper, "t_lower must be smaller than t_upper."
    return cosine_sigmoid((t-t_lower)/(t_upper-t_lower))


def taper_function(time, taper_type, sigmoid_type='smooth', **kwargs):
    """A taper function is either a sigmoid function (that can be increasing or
    decreasing) or a window function (that combines an increasing sigmoid and a
    decreasing sigmoid).

    A cosine sigmoid function that maps [t_lower, t_upper] to [0, 1]. If the 
    function is increasing, then f(t_lower)=0 and f(t_upper)=1; if decreasing, 
    then f(t_lower)=1 and f(t_upper)=0. The second derivative of this function 
    is NOT continuous. It is required that t_lower < t_upper and both 
    t_lower and t_upper lie within the time array t. 

    Parameters
    ----------
    time : array_like
        Array of time stamps
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
    output : array_like
        Values of the taper function
    """
    if sigmoid_type == 'smooth':
        sigmoid = smooth_sigmoid
    elif sigmoid_type == 'cosine':
        sigmoid = cosine_sigmoid
    else:
        raise ValueError("sigmoid_type must be 'smooth' or 'cosine'.")

    if taper_type == 'left':
        is_left_tapered = True
        is_right_tapered = False
    elif taper_type == 'right':
        is_left_tapered = False
        is_right_tapered = True
    elif taper_type == 'both':
        is_left_tapered = True
        is_right_tapered = True
    else:
        raise ValueError("taper_type must be 'left', 'right', or 'both'.")

    output = np.ones(len(time))

    if is_left_tapered:
        if 't0_taper' in kwargs and 't1_taper' in kwargs and 'left_taper_length' not in kwargs:
            t0 = kwargs['t0_taper']
            t1 = kwargs['t1_taper']
            assert time[0] <= t0 < t1 <= time[-1], "t0_taper must be smaller than t1_taper, and they both must lie within the time array."
        elif 't0_taper' not in kwargs and 't1_taper' not in kwargs and 'left_taper_length' in kwargs:
            t0 = time[0]
            t1 = time[0]+kwargs['left_taper_length']
            assert t1 <= time[-1], "left_taper_length must be shorter than the time array."
        else:
            raise ValueError("Left tapering is not specified correctly.")

        output *= sigmoid(time, t_lower=t0, t_upper=t1)

    if is_right_tapered:
        if 't2_taper' in kwargs and 't3_taper' in kwargs and 'right_taper_length' not in kwargs:
            t2 = kwargs['t2_taper']
            t3 = kwargs['t3_taper']
            assert time[0] <= t2 < t3 <= time[-1], "t2_taper must be smaller than t3_taper, and they both must lie within the time array."
        elif 't2_taper' not in kwargs and 't3_taper' not in kwargs and 'right_taper_length' in kwargs:
            t2 = time[-1]-kwargs['right_taper_length']
            t3 = time[-1]
            assert time[0] <= t2, "right_taper_length must be shorter than the time array."
        else:
            raise ValueError("Right tapering is not specified correctly.")

        output *= sigmoid(time, t_lower=t2, t_upper=t3, decreasing=True)

    return output
