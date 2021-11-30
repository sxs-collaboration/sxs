import numpy as np
import scipy


def align1d(wa, wb, t1, t2, n_brute_force=None):
    """Align waveforms by shifting in time

    This function determines the optimal time offset to apply to `wb` by minimizing
    the averaged (over time) squared difference in the L² norms (over the sphere)
    of the waveforms:
    
        ∫ [ ‖wa(t)‖ - ‖wb(t+δt)‖ ]² dt

    The integral is taken from time `t1` to `t2`.

    No changes are actually made to the input waveforms, but the result of this
    function can be used to offset the second waveform as in

        δt = align1d(wa, wb, t1, t2)
        wb.t -= δt

    As always, be careful, because `wb.t` is a reference to a numpy array that may
    be shared among copies of the data; use `wb.t = wb.t - δt` if you want to
    create a new copy of that array.

    Note that the input waveforms are assumed to be initially aligned at least well
    enough that:
    
      1) the time span from `t1` to `t2` in the two waveforms will overlap at
         least slightly after the second waveform is shifted in time; and
      2) waveform `wb` contains all the times corresponding to `t1` to `t2` in
         waveform `wa`.

    The first of these can usually be assured by simply aligning the peaks prior to
    calling this function:

        wb.t -= wb.max_norm_time() - wa.max_norm_time()

    The second assumption will be satisfied as long as `t1` is not too close to the
    beginning of `wb` and `t2` is not too close to the end.

    Parameters
    ----------
    wa : WaveformModes
    wb : WaveformModes
        Waveforms to be aligned
    t1 : float
    t2 : float
        Beginning and end of integration interval
    n_brute_force : int, optional
        Number of evenly spaced δt values between (t1-t2) and (t2-t1) to sample
        for the initial guess.  By default, this is just the maximum number of
        time steps in the range (t1, t2) in the input waveforms.  If this is
        too small, an incorrect local minimum may be found.

    Notes
    -----
    Choosing the time interval is usually the most difficult choice to make when
    aligning waveforms.  Assuming you want to align during inspiral, the times
    must span sufficiently long that the waveforms' norm (equivalently, orbital
    frequency changes) significantly from `t1` to `t2`.  This means that you
    cannot always rely on a specific number of orbits, for example.  Also note
    that neither number should be too close to the beginning or end of either
    waveform, to provide some "wiggle room".

    Precession generally causes no problems for this function.  In principle,
    eccentricity, center-of-mass offsets, boosts, or other supertranslations could
    cause problems, but this function begins with a brute-force method of finding
    the optimal time offset that will avoid local minima in all but truly
    outrageous situations.  In particular, as long as `t1` and `t2` are separated
    by enough, there should never be a problem.

    """
    from scipy.interpolate import CubicSpline
    from scipy.integrate import trapezoid
    from scipy.optimize import least_squares

    # Check that (t1, t2) makes sense and is actually contained in both waveforms
    if t2<=t1:
        raise ValueError(f"(t1,t2)=({t1}, {t2}) is out of order")
    if wa.t[0] > t1 or wa.t[-1] < t2:
        raise ValueError(f"(t1,t2)=({t1}, {t2}) not contained in wa.t, which spans ({wa.t[0]}, {wa.t[-1]})")
    if wb.t[0] > t1 or wb.t[-1] < t2:
        raise ValueError(f"(t1,t2)=({t1}, {t2}) not contained in wb.t, which spans ({wb.t[0]}, {wb.t[-1]})")

    # Figure out time offsets to try
    δt_lower = max(t1 - t2, wb.t[0]-t1)
    δt_upper = min(t2 - t1, wb.t[-1]-t2)

    # We'll start by brute forcing, sampling time offsets evenly at as many
    # points as there are time steps in (t1,t2) in the input waveforms
    if n_brute_force is None:
        n_brute_force = max(sum((wa.t>=t1)&(wa.t<=t2)), sum((wb.t>=t1)&(wb.t<=t2)))
    δt_brute_force = np.linspace(δt_lower, δt_upper, num=n_brute_force)

    # Times at which the differences will be evaluated
    t_reference = wa.t[(wa.t>=t1)&(wa.t<=t2)]

    # Define the cost function
    norm_a = wa.norm.ndarray[(wa.t>=t1)&(wa.t<=t2)]
    norm_b = CubicSpline(wb.t, wb.norm.ndarray)
    normalization = trapezoid((norm_a)**2, t_reference)
    def cost(δt):
        # Take the sqrt because least_squares squares the inputs...
        return np.sqrt(trapezoid((norm_a - norm_b(t_reference+δt))**2, t_reference) / normalization)

    # Optimize by brute force
    cost_brute_force = [cost(δt) for δt in δt_brute_force]
    δt = δt_brute_force[np.argmin(cost_brute_force)]

    # Optimize explicitly
    optimum = least_squares(cost, δt, bounds=(δt_lower, δt_upper))

    return optimum.x[0]
