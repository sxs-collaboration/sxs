import numpy as np

import quaternionic

from scipy.integrate import trapezoid

import multiprocessing as mp
from functools import partial


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
    aligning waveforms.  Assuming you want to align during inspiral, the times must
    span sufficiently long that the waveforms' norm (equivalently, orbital
    frequency changes) significantly from `t1` to `t2`.  This means that you cannot
    always rely on a specific number of orbits, for example.  Also note that
    neither number should be too close to the beginning or end of either waveform,
    to provide some "wiggle room".

    Precession generally causes no problems for this function.  In principle,
    eccentricity, center-of-mass offsets, boosts, or other supertranslations could
    cause problems, but this function begins with a brute-force method of finding
    the optimal time offset that will avoid local minima in all but truly
    outrageous situations.  In particular, as long as `t1` and `t2` are separated
    by enough, there should never be a problem.

    """
    from scipy.interpolate import CubicSpline
    from scipy.optimize import least_squares

    # Check that (t1, t2) makes sense and is actually contained in both waveforms
    if t2 <= t1:
        raise ValueError(f"(t1,t2)=({t1}, {t2}) is out of order")
    if wa.t[0] > t1 or wa.t[-1] < t2:
        raise ValueError(f"(t1,t2)=({t1}, {t2}) not contained in wa.t, which spans ({wa.t[0]}, {wa.t[-1]})")
    if wb.t[0] > t1 or wb.t[-1] < t2:
        raise ValueError(f"(t1,t2)=({t1}, {t2}) not contained in wb.t, which spans ({wb.t[0]}, {wb.t[-1]})")

    # Figure out time offsets to try
    δt_lower = max(t1 - t2, t2 - wb.t[-1])
    δt_upper = min(t2 - t1, t1 - wb.t[0])

    # We'll start by brute forcing, sampling time offsets evenly at as many
    # points as there are time steps in (t1,t2) in the input waveforms
    if n_brute_force is None:
        n_brute_force = max(sum((wa.t >= t1) & (wa.t <= t2)), sum((wb.t >= t1) & (wb.t <= t2)))
    δt_brute_force = np.linspace(δt_lower, δt_upper, num=n_brute_force)

    # Times at which the differences will be evaluated
    t_reference = wa.t[(wa.t >= t1) & (wa.t <= t2)]

    # Define the cost function
    norm_a = wa.norm.ndarray[(wa.t >= t1) & (wa.t <= t2)]
    norm_b = CubicSpline(wb.t, wb.norm.ndarray)
    normalization = trapezoid((norm_a) ** 2, t_reference)

    def cost(δt):
        # Take the sqrt because least_squares squares the inputs...
        return np.sqrt(trapezoid((norm_a - norm_b(t_reference + δt)) ** 2, t_reference) / normalization)

    # Optimize by brute force
    cost_brute_force = [cost(δt) for δt in δt_brute_force]
    δt = δt_brute_force[np.argmin(cost_brute_force)]

    # Optimize explicitly
    optimum = least_squares(cost, δt, bounds=(δt_lower, δt_upper))

    return optimum.x[0]


def _cost2d(δt_δϕ, args):
    modes_A, modes_B, t_reference, m, δΨ_factor, normalization = args
    δt, δϕ = δt_δϕ

    # Take the sqrt because least_squares squares the inputs...
    diff = trapezoid(
        np.sum(abs(modes_A(t_reference + δt) * np.exp(1j * m * δϕ) * δΨ_factor - modes_B) ** 2, axis=1),
        t_reference,
    )
    return np.sqrt(diff / normalization)


def align2d(wa, wb, t1, t2, n_brute_force_δt=None, n_brute_force_δϕ=5, max_δt=np.inf, use_δΨ=False, include_modes=None, nprocs=None):
    """Align waveforms by shifting in time and phase

    This function determines the optimal time and phase offset to apply to `wa` by
    minimizing the averaged (over time) L² norm (over the sphere) of the difference
    of the waveforms.

    The integral is taken from time `t1` to `t2`.

    Note that the input waveforms are assumed to be initially aligned at least well
    enough that:

      1) the time span from `t1` to `t2` in the two waveforms will overlap at
         least slightly after the second waveform is shifted in time; and
      2) waveform `wb` contains all the times corresponding to `t1` to `t2` in
         waveform `wa`.

    The first of these can usually be assured by simply aligning the peaks prior to
    calling this function:

        wa.t -= wa.max_norm_time() - wb.max_norm_time()

    The second assumption will be satisfied as long as `t1` is not too close to the
    beginning of `wb` and `t2` is not too close to the end.

    Parameters
    ----------
    wa : WaveformModes
    wb : WaveformModes
        Waveforms to be aligned
    t1 : float
    t2 : float
        Beginning and end of integration interval.
    n_brute_force_δt : int, optional
        Number of evenly spaced δt values between (t1-t2) and (t2-t1) to sample
        for the initial guess.  By default, this is just the maximum number of
        time steps in the range (t1, t2) in the input waveforms.  If this is
        too small, an incorrect local minimum may be found.
    n_brute_force_δϕ : int, optional
        Number of evenly spaced δϕ values between 0 and 2π to sample
        for the initial guess.  By default, this is 5, but if the option
        `None` is provided, then this will be 2 * ell_max + 1. This option
        is not the default because, even though it is formally the right
        thing to do, it takes much longer to produce the same result
        as just using 5, which is much faster.
    max_δt : float, optional
        Max δt to allow for when choosing the initial guess.
    use_δΨ : float, optional
        Whether or not to allow for a h -> -h transformation.
    include_modes: list, optional
        A list containing the (ell, m) modes to be included in the L² norm.
    nprocs: int, optional
        Number of cpus to use. Default is maximum number.
        If -1 is provided, then no multiprocessing is performed.

    Returns
    -------
    error: float
        Cost of scipy.optimize.least_squares
        This is 0.5 ||wa - wb||² / ||wb||²
    wa_prime: WaveformModes
        Resulting waveform after transforming `wa` using `optimum`
    optimum: OptimizeResult
        Result of scipy.optimize.least_squares

    Notes
    -----
    Choosing the time interval is usually the most difficult choice to make when
    aligning waveforms.  Assuming you want to align during inspiral, the times must
    span sufficiently long that the waveforms' norm (equivalently, orbital
    frequency changes) significantly from `t1` to `t2`.  This means that you cannot
    always rely on a specific number of orbits, for example.  Also note that
    neither number should be too close to the beginning or end of either waveform,
    to provide some "wiggle room".

    """
    from scipy.interpolate import CubicSpline
    from scipy.optimize import least_squares
    from .. import WaveformModes

    wa_orig = wa
    wa = wa.copy()
    wb = wb.copy()

    # Check that (t1, t2) makes sense and is actually contained in both waveforms
    if t2 <= t1:
        raise ValueError(f"(t1,t2)=({t1}, {t2}) is out of order")
    if wa.t[0] > t1 or wa.t[-1] < t2:
        raise ValueError(f"(t1,t2)=({t1}, {t2}) not contained in wa.t, which spans ({wa.t[0]}, {wa.t[-1]})")
    if wb.t[0] > t1 or wb.t[-1] < t2:
        raise ValueError(f"(t1,t2)=({t1}, {t2}) not contained in wb.t, which spans ({wb.t[0]}, {wb.t[-1]})")

    # Figure out time offsets to try
    δt_lower = max(-max_δt, max(t1 - t2, t2 - wa.t[-1]))
    δt_upper = min(max_δt, min(t2 - t1, t1 - wa.t[0]))
    
    # We'll start by brute forcing, sampling time offsets evenly at as many
    # points as there are time steps in (t1,t2) in the input waveforms
    if n_brute_force_δt is None:
        n_brute_force_δt = max(sum((wa.t >= t1) & (wa.t <= t2)), sum((wb.t >= t1) & (wb.t <= t2)))
    δt_brute_force = np.array(list(np.linspace(δt_lower, δt_upper, num=n_brute_force_δt)) + [0.0])

    if n_brute_force_δϕ is None:
        n_brute_force_δϕ = 2 * wa.ell_max + 1
    δϕ_brute_force = np.linspace(0, 2 * np.pi, n_brute_force_δϕ, endpoint=False)

    δt_δϕ_brute_force = np.array(np.meshgrid(δt_brute_force, δϕ_brute_force)).T.reshape(-1, 2)

    t_reference = wa.t[np.argmin(abs(wa.t - t1)) : np.argmin(abs(wa.t - t2)) + 1]

    # Remove certain modes, if requested
    ell_max = min(wa.ell_max, wb.ell_max)
    if include_modes is not None:
        for L in range(2, ell_max + 1):
            for M in range(-L, L + 1):
                if not (L, M) in include_modes:
                    wa.data[:, wa.index(L, M)] *= 0
                    wb.data[:, wb.index(L, M)] *= 0

    # Define the cost function
    modes_A = CubicSpline(wa.t, wa[:, wa.index(2, -2) : wa.index(ell_max + 1, -(ell_max + 1))].data)
    modes_B = CubicSpline(wb.t, wb[:, wb.index(2, -2) : wb.index(ell_max + 1, -(ell_max + 1))].data)(t_reference)

    normalization = trapezoid(
        CubicSpline(wb.t, wb[:, wb.index(2, -2) : wb.index(ell_max + 1, -(ell_max + 1))].norm ** 2)(t_reference),
        t_reference,
    )

    m = np.array([M for L in range(2, ell_max + 1) for M in range(-L, L + 1)])

    optimums = []
    wa_primes = []
    δΨ_factors = [1]
    if use_δΨ:
        δΨ_factors = [-1, +1]
    for δΨ_factor in δΨ_factors:
        # Optimize by brute force with multiprocessing
        cost_wrapper = partial(_cost2d, args=[modes_A, modes_B, t_reference, m, δΨ_factor, normalization])

        if nprocs != -1:
            if nprocs is None:
                nprocs = mp.cpu_count()
            pool = mp.Pool(processes=nprocs)
            cost_brute_force = pool.map(cost_wrapper, δt_δϕ_brute_force)
            pool.close()
            pool.join()
        else:
            cost_brute_force = [cost_wrapper(δt_δϕ_brute_force_item) for δt_δϕ_brute_force_item in δt_δϕ_brute_force]

        δt_δϕ = δt_δϕ_brute_force[np.argmin(cost_brute_force)]

        # Optimize explicitly
        optimum = least_squares(cost_wrapper, δt_δϕ, bounds=[(δt_lower, 0), (δt_upper, 2 * np.pi)], max_nfev=50000)
        optimums.append(optimum)
        δt, δϕ = optimum.x

        wa_prime = WaveformModes(
            input_array=(
                wa_orig[:, wa_orig.index(2, -2) : wa_orig.index(ell_max + 1, -(ell_max + 1))].data
                * np.exp(1j * m * δϕ)
                * δΨ_factor
            ),
            time=wa_orig.t - δt,
            time_axis=0,
            modes_axis=1,
            ell_min=2,
            ell_max=ell_max,
        )
        wa_primes.append(wa_prime)

    idx = np.argmin(abs(np.array([optimum.cost for optimum in optimums])))

    return optimums[idx].cost, wa_primes[idx], optimums[idx]


def _cost4d(δt_δso3, args):
    from .. import WaveformModes

    modes_A, modes_B, t_reference, normalization = args
    δt = δt_δso3[0]
    δSO3 = np.exp(quaternionic.array.from_vector_part(δt_δso3[1:]))

    modes_A_at_δt = modes_A(t_reference + δt)
    ell_max = int(np.sqrt(modes_A_at_δt.shape[1] + 4)) - 1

    wa_prime = WaveformModes(
        input_array=(modes_A_at_δt),
        time=t_reference,
        time_axis=0,
        modes_axis=1,
        ell_min=2,
        ell_max=ell_max,
        spin_weight=-2,
    )

    wa_prime = wa_prime.rotate(δSO3)

    # Take the sqrt because least_squares squares the inputs...
    diff = trapezoid(
        np.sum(abs(wa_prime.data - modes_B) ** 2, axis=1),
        t_reference,
    )
    return np.sqrt(diff / normalization)


def align4d(
    wa,
    wb,
    t1,
    t2,
    n_brute_force_δt=1_000,
    n_brute_force_δϕ=None,
    max_δt=np.inf,
    include_modes=None,
    align2d_first=False,
    nprocs=None,
):
    """Align waveforms by optimizing over a time translation and an SO(3) rotation.

    This function determines the optimal transformation to apply to `wa` by 
    minimizing the averaged (over time) L² norm (over the sphere) of the 
    difference of the waveforms.

    The integral is taken from time `t1` to `t2`.

    Note that the input waveforms are assumed to be initially aligned at least 
    well enough that:

      1) the time span from `t1` to `t2` in the two waveforms will overlap at 
         least slightly after the second waveform is shifted in time; and
      2) waveform `wb` contains all the times corresponding to `t1` to `t2` in 
         waveform `wa`.

    The first of these can usually be assured by simply aligning the peaks prior 
    to calling this function:

        wa.t -= wa.max_norm_time() - wb.max_norm_time()

    The second assumption will be satisfied as long as `t1` is not too close to 
    the beginning of `wb` and `t2` is not too close to the end.

    Parameters
    ----------
    wa : WaveformModes
    wb : WaveformModes
        Waveforms to be aligned
    t1 : float
    t2 : float
        Beginning and end of integration interval
    n_brute_force_δt : int, optional
        Number of evenly spaced δt values between (t1-t2) and (t2-t1) to sample 
        for the initial guess.  By default, this is 1,000.  If this is too small, 
        an incorrect local minimum may be found.
    n_brute_force_δϕ : int, optional
        Number of evenly spaced angles about the angular-velocity axis to sample 
        for the initial guess.  By default, this is `2 * (2 * ell_max + 1)`.
    max_δt : float, optional
        Max δt to allow for when choosing the initial guess.
    include_modes: list, optional
        A list containing the (ell, m) modes to be included in the L² norm.
    align2d_first : bool, optional
        Do a 2d align first for an initial guess, with no SO(3) initial guess 
        (besides the phase returned by the 2d solve)
    nprocs: int, optional
        Number of cpus to use.  Default is maximum number.  If -1 is provided, 
        then no multiprocessing is performed.

    Returns
    -------
    error: float
        Cost of scipy.optimize.least_squares.  This is 0.5 ||wa - wb||² / ||wb||²
    wa_prime: WaveformModes
        Resulting waveform after transforming `wa` using `optimum`
    optimum: OptimizeResult
        Result of scipy.optimize.least_squares

    Notes
    -----
    Choosing the time interval is usually the most difficult choice to make when 
    aligning waveforms.  Assuming you want to align during inspiral, the times 
    must span sufficiently long that the waveforms' norm (equivalently, orbital 
    frequency changes) significantly from `t1` to `t2`.  This means that you 
    cannot always rely on a specific number of orbits, for example.  Also note 
    that neither number should be too close to the beginning or end of either 
    waveform, to provide some "wiggle room".

    """
    from scipy.interpolate import CubicSpline
    from scipy.optimize import least_squares
    from .. import WaveformModes

    if wa.spin_weight != wb.spin_weight:
        raise ValueError(f"Waveform spin weights do not match: {wa.spin_weight=} != {wb.spin_weight=}")
    spin_weight = wa.spin_weight
    ell_min = max(wa.ell_min, wb.ell_min)
    ell_max = min(wa.ell_max, wb.ell_max)

    wa_orig = wa
    wa = wa.copy()
    wb = wb.copy()

    # Check that (t1, t2) makes sense and is actually contained in both waveforms
    if t2 <= t1:
        raise ValueError(f"(t1,t2)=({t1}, {t2}) is out of order")
    if wa.t[0] > t1 or wa.t[-1] < t2:
        raise ValueError(f"(t1,t2)=({t1}, {t2}) not contained in wa.t, which spans ({wa.t[0]}, {wa.t[-1]})")
    if wb.t[0] > t1 or wb.t[-1] < t2:
        raise ValueError(f"(t1,t2)=({t1}, {t2}) not contained in wb.t, which spans ({wb.t[0]}, {wb.t[-1]})")

    # Figure out time offsets to try
    δt_lower = max(-max_δt, max(t1 - t2, t2 - wa.t[-1]))
    δt_upper = min(max_δt, min(t2 - t1, t1 - wa.t[0]))

    t_reference = wb.t[np.argmin(abs(wb.t - t1)) : np.argmin(abs(wb.t - t2)) + 1]
    
    if not align2d_first:
        # Get time initial guess
        # Negative sign because align1d aligns wb to wa
        δt_IG = -align1d(wa, wb, t1, t2)
        
        wa_interp = wa.interpolate(t_reference + δt_IG)
        wb_interp = wb.interpolate(t_reference)
        
        # Get rotor initial guess
        omegaa = wa_interp.angular_velocity
        omegab = wb_interp.angular_velocity
        R_IG = quaternionic.align(omegaa, omegab)
    else:
        _, _, res = align2d(
            wa, wb, t1, t2, max_δt=50,
            n_brute_force_δt=n_brute_force_δt, nprocs=nprocs
        )
        δt_IG = res.x[0]

        wa_interp = wa.interpolate(t_reference + δt_IG)
        wb_interp = wb.interpolate(t_reference)

        # Get rotor initial guess
        omegaa = wa_interp.angular_velocity
        omegab = wb_interp.angular_velocity
        R_IG = quaternionic.align(omegaa, omegab)
        
        R_IG = quaternionic.array([0, 0, 0, res.x[1] / 2]) * R_IG

    # Brute force over R_IG * exp(theta * z / 2) with δt_IG
    n_brute_force_δϕ = n_brute_force_δϕ or 2 * (2 * ell_max + 1)
    δt_δso3_brute_force = [
        [δt_IG, *np.log(R_IG * np.exp(quaternionic.array([0, 0, 0, angle / 2]))).vector]
        for angle in np.linspace(0, 2 * np.pi, num=n_brute_force_δϕ, endpoint=False)
    ]

    # Remove certain modes, if requested
    if include_modes is not None:
        for L in range(ell_min, ell_max + 1):
            for M in range(-L, L + 1):
                if not (L, M) in include_modes:
                    wa.data[:, wa.index(L, M)] *= 0
                    wb.data[:, wb.index(L, M)] *= 0

    # Define the cost function
    modes_A = CubicSpline(wa.t, wa[:, wa.index(ell_min, -ell_min) : wa.index(ell_max, ell_max)+1].data)
    modes_B = CubicSpline(wb.t, wb[:, wb.index(ell_min, -ell_min) : wb.index(ell_max, ell_max)+1].data)(t_reference)

    normalization = trapezoid(
        CubicSpline(
            wb.t,
            wb[:, wb.index(ell_min, -ell_min) : wb.index(ell_max, ell_max)+1].norm ** 2
        )(t_reference),
        t_reference,
    )

    # Optimize by brute force with multiprocessing
    cost_wrapper = partial(_cost4d, args=[modes_A, modes_B, t_reference, normalization])

    if nprocs != -1:
        if nprocs is None:
            nprocs = mp.cpu_count()
        pool = mp.Pool(processes=nprocs)
        cost_brute_force = pool.map(cost_wrapper, δt_δso3_brute_force)
        pool.close()
        pool.join()
    else:
        cost_brute_force = [
            cost_wrapper(δt_δso3_brute_force_item)
            for δt_δso3_brute_force_item in δt_δso3_brute_force
        ]

    δt_δso3 = δt_δso3_brute_force[np.argmin(cost_brute_force)]
    
    # Optimize explicitly
    optimum = least_squares(
        cost_wrapper,
        δt_δso3,
        bounds=[(δt_lower, -np.pi/2, -np.pi/2, -np.pi/2), (δt_upper, np.pi/2, np.pi/2, np.pi/2)],
        max_nfev=50000,
    )
    δt = optimum.x[0]
    δSO3 = np.exp(quaternionic.array.from_vector_part(optimum.x[1:]))
    
    wa_prime = WaveformModes(
        input_array=(wa_orig[:, wa_orig.index(ell_min, -ell_min) : wa_orig.index(ell_max, ell_max)+1].data),
        time=wa_orig.t - δt,
        time_axis=0,
        modes_axis=1,
        ell_min=ell_min,
        ell_max=ell_max,
        spin_weight=spin_weight,
    )
    wa_prime = wa_prime.rotate(δSO3)

    return optimum.cost, wa_prime, optimum
