import numpy as np
import multiprocessing as mp
from scipy.integrate import trapezoid
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares

from .waveform_modes import WaveformModes

_MP_STATE = {}


def check_time_constraint(wa, wb, t1, t2):
    """Check that the times t1 and t2 are contained in wa.t and wb.t.

    Parameters
    ----------
    wa : WaveformModes
    wb : WaveformModes
    t1 : float
    t2 : float
    """
    time_intersection = (max(wa.t[0], wb.t[0]), min(wa.t[-1], wb.t[-1]))
    if t1 < time_intersection[0]:
        raise ValueError(f"{t1} is not contained in intersection of waveform time arrays: {time_intersection}.")
    if t2 > time_intersection[1]:
        raise ValueError(f"{t2} is not contained in intersection of waveform time arrays: {time_intersection}.")


def create_unified_waveforms(wa, wb, t1, t2, padding_time_factor=0.2):
    """Convert WaveformModes to common time and modes

    The output waveforms will be interpolated to a set of times
    including — as nearly as possible — the times `t1` and `t2`,
    padded by the given fraction of `(t2-t1)` on other side.

    Parameters
    ----------
    wa : WaveformModes
    wb : WaveformModes
    t1 : float
    t2 : float
    padding_time_factor : float, optional
        Extra time of size (t2 - t1)*padding_time_factor to include
        int the waveforms.  Default is 0.2.

    Returns
    -------
    wa_interp : WaveformModes
    wb_interp : WaveformModes
    """
    check_time_constraint(wa, wb, t1, t2)

    padding_time = (t2 - t1) * padding_time_factor
    t1_padded = max(wa.t[0], wb.t[0], t1 - padding_time)
    t2_padded = min(wa.t[-1], wb.t[-1], t2 + padding_time)

    idx1 = np.argmin(abs(wa.t - t1_padded))
    idx2 = np.argmin(abs(wa.t - t2_padded)) + 1

    wa_interp = wa.interpolate(wa.t[idx1:idx2])
    wb_interp = wb.interpolate(wa_interp.t)

    ell_min = max(wa_interp.ell_min, wb_interp.ell_min)
    ell_max = min(wa_interp.ell_max, wb_interp.ell_max)
    ia1 = wa_interp.index(ell_min, -ell_min)
    ia2 = wa_interp.index(ell_max, ell_max) + 1
    ib1 = wb_interp.index(ell_min, -ell_min)
    ib2 = wb_interp.index(ell_max, ell_max) + 1

    return wa_interp[:, ia1:ia2], wb_interp[:, ib1:ib2]


def L2_difference(wa, wb, t1=-np.inf, t2=np.inf, modes=None, modes_for_norm=None, normalize=True):
    """Compute L² norm of difference between two waveforms

    The norm is integrated over the time window (`t1`, `t2`), and over
    the sphere using the modes specified by `modes`.

    Parameters
    ----------
    wa : WaveformModes
    wb : WaveformModes
        Waveforms to compare.
    t1 : float
        Beginning of L² norm integral.
    t2 : float
        End of L² norm integral.
    modes : list, optional
        Modes (ell, m) to include in numerator of L² norm calculation.
        Default is all modes.
    modes_for_norm : list, optional
        Modes (ell, m) to include in denominator of L² norm
        calculation.  Default is all modes.
    normalize : bool, optional
        Whether or not to divide by sqrt(||wa||²||wb||²) in the sqrt.
        If False, returns the unnormalized L² norm of the residual and
        what would have been the normalization.  Default is True.

    Returns
    -------
    L2_norm: float
        L² norm of the residual between the waveforms.  This is sqrt(
        ||wa - wb||² / sqrt(||wa||² ||wb||²) )
    """
    t1 = max(wa.t[0], wb.t[0], t1)
    t2 = min(wa.t[-1], wb.t[-1], t2)

    already_unified = np.array_equal(wa.t, wb.t) and (wa.ndim == wb.ndim == 1 or np.array_equal(wa.LM, wb.LM))

    if already_unified:
        i1, i2 = np.argmin(abs(wa.t - t1)), np.argmin(abs(wa.t - t2)) + 1
        wa = wa.copy()[i1:i2]
        wb = wb.copy()[i1:i2]
    else:
        wa, wb = create_unified_waveforms(wa, wb, t1, t2, padding_time_factor=0)

    data_diff = wa.data - wb.data

    data_for_norm1 = wa.data
    data_for_norm2 = wb.data

    # Eliminate unwanted modes
    if modes is not None or modes_for_norm is not None:
        for L in range(wa.ell_min, wa.ell_max + 1):
            for M in range(-L, L + 1):
                if modes is not None:
                    if not (L, M) in modes:
                        data_diff[:, wa.index(L, M)] *= 0
                if modes_for_norm is not None:
                    if not (L, M) in modes_for_norm:
                        data_for_norm1[:, wa.index(L, M)] *= 0
                        data_for_norm2[:, wb.index(L, M)] *= 0

    if normalize:
        L2_norm = np.sqrt(
            trapezoid(np.linalg.norm(data_diff, axis=wa.modes_axis) ** 2, wa.t)
            / np.sqrt(
                trapezoid(np.linalg.norm(data_for_norm1, axis=wa.modes_axis) ** 2, wa.t)
                * trapezoid(np.linalg.norm(data_for_norm2, axis=wa.modes_axis) ** 2, wb.t)
            )
        )

        return L2_norm
    else:
        L2_norm_unnormalized = np.sqrt(trapezoid(np.linalg.norm(data_diff, axis=wa.modes_axis) ** 2, wa.t))
        norm = np.sqrt(
            np.sqrt(
                trapezoid(np.linalg.norm(data_for_norm1, axis=wa.modes_axis) ** 2, wa.t)
                * trapezoid(np.linalg.norm(data_for_norm2, axis=wb.modes_axis) ** 2, wb.t)
            )
        )

        return L2_norm_unnormalized, norm


def inner_product(wa, wb, ASD_values=None):
    """Compute the inner product between two waveforms

    No interpolation is performed, including for the `ASD`.  If the
    inputs are `WaveformModes` objects, the modes are assumed to match
    between the two waveforms.

    Note that this can be provided with time-domain or
    frequency-domain waveforms, either over the two-sphere
    (WaveformModes) or at a point on the two-sphere (TimeSeries).

    For frequency-domain waveforms, the .t attribute should be the
    frequency.

    Parameters
    ----------
    wa : WaveformModes or TimeSeries
    wb : WaveformModes or TimeSeries
    ASD_values : ndarray, optioanl
        ASD values for frequency-domain overlaps.
        Default is flat ASD.

    Returns
    -------
    inner_product : float
        Inner product between the two waveforms.
    """
    if ASD_values is None:
        ASD_values = 1

    if not wa.ndim == wb.ndim == 1:
        inner_product = trapezoid(
            np.sum(wa.data * np.conjugate(wb.data), axis=1) / ASD_values**2,
            wa.t,
        )
    else:
        inner_product = trapezoid(
            (wa * np.conjugate(wb)).ndarray / ASD_values**2,
            wa.t,
        )

    return inner_product


def mismatch(wa, wb, x1=-np.inf, x2=np.inf, modes=None, ASD=None):
    """Compute the mismatch between two waveforms

    The mismatch is calculated over the time or frequency window
    (`x1`, `x2`), and — if relevant — over the sphere using the modes
    specified by `modes`.

    Note that this can be provided with time-domain or
    frequency-domain waveforms, either over the two-sphere
    (WaveformModes) or at a point on the two-sphere (TimeSeries).

    For frequency-domain waveforms, the .t attribute should be the
    frequency.

    Parameters
    ----------
    wa : WaveformModes or TimeSeries
    wb : WaveformModes or TimeSeries
    x1 : float, optional
        Beginning of mismatch integral.  Default uses all values.
    x2 : float, optional
        End of mismatch integral.  Default uses all values.
    modes : list, optional
        Modes (ell, m) to include in the mismatch calculation.
        Default is all modes.
    ASD : func, optional
        Function mapping frequencies to the ASD of a detector.
        Default is flat ASD.

    Returns
    -------
    mismatch : float
        Mismatch between the two waveforms.
    """
    x1 = max(wa.t[0], wb.t[0], x1)
    x2 = min(wa.t[-1], wb.t[-1], x2)

    already_unified = np.array_equal(wa.t, wb.t) and (wa.ndim == wb.ndim == 1 or np.array_equal(wa.LM, wb.LM))

    if already_unified:
        i1, i2 = np.argmin(abs(wa.t - x1)), np.argmin(abs(wa.t - x2)) + 1
        wa = wa.copy()[i1:i2]
        wb = wb.copy()[i1:i2]
    else:
        wa, wb = create_unified_waveforms(wa, wb, x1, x2, padding_time_factor=0)

    # Eliminate unwanted modes
    if modes is not None:
        if wa.ndim > 1:
            for L in range(wa.ell_min, wa.ell_max + 1):
                for M in range(-L, L + 1):
                    if not (L, M) in modes:
                        wa.data[:, wa.index(L, M)] *= 0
                        wb.data[:, wb.index(L, M)] *= 0
        else:
            raise ValueError("A `modes` argument was provided, but the input `wa` " "and `wb` only have one dimension.")

    if ASD is not None:
        ASD_values = ASD(wa.t)
    else:
        ASD_values = 1

    wa_wb_overlap = inner_product(wa, wb, ASD_values=ASD_values).real
    wa_norm = inner_product(wa, wa, ASD_values=ASD_values).real
    wb_norm = inner_product(wb, wb, ASD_values=ASD_values).real

    return 1 - wa_wb_overlap / np.sqrt(wa_norm * wb_norm)


def auto_mp_start_method(requested=None):
    if requested is not None:
        return requested
    try:
        methods = mp.get_all_start_methods()
    except Exception:
        methods = []
    return "fork" if "fork" in methods else "spawn"


def limit_worker_threads():
    try:
        from threadpoolctl import threadpool_limits

        return threadpool_limits(limits=1)
    except Exception:
        return None


def init_mismatch_pool(
    wa_t,
    wa_mode_data,
    wb_mode_data,
    t_reference,
    ell_min,
    ell_max,
    spin_weight,
    optimize_polarization,
    dt_lower,
    dt_upper,
    x1,
    x2,
    dt_brute_force,
    dphi_brute_force,
    limit_threads,
):
    import sxs

    _MP_STATE.clear()

    if limit_threads:
        _MP_STATE["_threadpool_limit"] = limit_worker_threads()

    wa_t = np.asarray(wa_t, dtype=float)
    wa_mode_data = np.asarray(wa_mode_data)
    wb_mode_data = np.asarray(wb_mode_data)
    t_reference = np.asarray(t_reference, dtype=float)
    dt_brute_force = np.asarray(dt_brute_force, dtype=float)
    dphi_brute_force = np.asarray(dphi_brute_force, dtype=float)

    m_modes = np.array(
        [m for ell in range(int(ell_min), int(ell_max) + 1) for m in range(-ell, ell + 1)],
        dtype=int,
    )

    _MP_STATE["wa_reference_modes"] = CubicSpline(wa_t, wa_mode_data, axis=0)
    _MP_STATE["wb_reference_modes"] = WaveformModes(
        input_array=wb_mode_data,
        time=t_reference,
        time_axis=0,
        modes_axis=1,
        ell_min=int(ell_min),
        ell_max=int(ell_max),
        spin_weight=spin_weight,
    )
    _MP_STATE["t_reference"] = t_reference
    _MP_STATE["ell_min"] = int(ell_min)
    _MP_STATE["ell_max"] = int(ell_max)
    _MP_STATE["spin_weight"] = spin_weight
    _MP_STATE["m_modes"] = m_modes
    _MP_STATE["optimize_polarization"] = bool(optimize_polarization)
    _MP_STATE["dt_lower"] = float(dt_lower)
    _MP_STATE["dt_upper"] = float(dt_upper)
    _MP_STATE["x1"] = float(x1)
    _MP_STATE["x2"] = float(x2)
    _MP_STATE["dt_brute_force"] = dt_brute_force
    _MP_STATE["dphi_brute_force"] = dphi_brute_force
    _MP_STATE["phase_brute_force"] = np.exp(1j * np.outer(dphi_brute_force, m_modes))
    _MP_STATE["wb_target_cache"] = {}


def get_cached_wb_target(theta, phi):
    key = (float(theta), float(phi))
    cache = _MP_STATE["wb_target_cache"]
    if key not in cache:
        cache[key] = _MP_STATE["wb_reference_modes"].evaluate(theta, phi)
    return cache[key]


def restrict_to_common_modes(wa, wb, include_modes=None):
    ell_min = max(wa.ell_min, wb.ell_min)
    ell_max = min(wa.ell_max, wb.ell_max)
    if ell_max < ell_min:
        raise ValueError("The two waveforms have no overlapping (ell, m) content.")

    ia1 = wa.index(ell_min, -ell_min)
    ia2 = wa.index(ell_max, ell_max) + 1
    ib1 = wb.index(ell_min, -ell_min)
    ib2 = wb.index(ell_max, ell_max) + 1

    wa_sliced = wa[:, ia1:ia2]
    wb_sliced = wb[:, ib1:ib2]

    if include_modes is not None:
        wa_sliced_data = wa_sliced.data
        wb_sliced_data = wb_sliced.data
        for mode in include_modes:
            wa_sliced_data[:, wa_sliced.index(*mode)] *= 0
            wb_sliced_data[:, wb_sliced.index(*mode)] *= 0

            wa_sliced = WaveformModes(
                time=wa_sliced.t,
                input_array=wa_sliced_data,
                time_axis=wa_sliced.time_axis,
                modes_axis=wa_sliced.modes_axis,
                ell_min=wa_sliced.ell_min,
                ell_max=wa_sliced.ell_max,
                spin_weight=wa_sliced.spin_weight,
            )

            wb_sliced = WaveformModes(
                time=wb_sliced.t,
                input_array=wb_sliced_data,
                time_axis=wb_sliced.time_axis,
                modes_axis=wb_sliced.modes_axis,
                ell_min=wb_sliced.ell_min,
                ell_max=wb_sliced.ell_max,
                spin_weight=wb_sliced.spin_weight,
            )

    return wa_sliced, wb_sliced


def choose_points_on_sphere(is_aligned, N_theta=7, N_phi=5):
    """
    N_theta needs to be odd to ensure pi/2 is included.

    If is_aligned == False:
        Distributes (N_theta-2)*N_phi + 2 points uniformly across the sky.

    If is_aligned == True, assumes aligned spins and retains only the upper
    hemisphere plus the equator.
    """
    if type(is_aligned) is not bool:
        raise TypeError("Expected `is_aligned` to be a bool.")
    if N_theta % 2 != 1:
        raise ValueError("N_theta needs to be odd.")

    th_vec = np.arccos(np.linspace(1.0, -1.0, N_theta))
    ph_vec = np.linspace(0.0, 2.0 * np.pi, N_phi, endpoint=False)

    th_phi_list = []
    for th in th_vec:
        if th < np.pi / 2.0 + 1e-5 or not is_aligned:
            if np.isclose(th, 0.0) or np.isclose(th, np.pi):
                th_phi_list.append((float(th), 0.0))
            else:
                for ph in ph_vec:
                    th_phi_list.append((float(th), float(ph)))

    if is_aligned:
        expected = (N_theta // 2 - 1) * N_phi + 1 + N_phi
        if len(th_phi_list) != expected:
            raise RuntimeError("Did not get expected number of aligned-spin points.")
    else:
        expected = (N_theta - 2) * N_phi + 2
        if len(th_phi_list) != expected:
            raise RuntimeError("Did not get expected number of points.")

    return th_phi_list


def build_transformed_point_signal_from_shifted(
    shifted_mode_data,
    theta,
    phi,
    phase_vector,
    wb_target,
    inner_product,
):
    import sxs

    wa_prime = WaveformModes(
        input_array=shifted_mode_data * phase_vector,
        time=_MP_STATE["t_reference"],
        time_axis=0,
        modes_axis=1,
        ell_min=_MP_STATE["ell_min"],
        ell_max=_MP_STATE["ell_max"],
        spin_weight=_MP_STATE["spin_weight"],
    )

    wa_prime_theta_phi = wa_prime.evaluate(theta, phi)

    if _MP_STATE["optimize_polarization"]:
        overlap = inner_product(wa_prime_theta_phi, wb_target)
        if np.abs(overlap) > 0.0:
            wa_prime_theta_phi *= np.exp(-1j * np.angle(overlap))

    return wa_prime_theta_phi


def point_cost_from_shifted(shifted_mode_data, theta, phi, phase_vector, wb_target, t_reference, inner_product):
    wa_trial = build_transformed_point_signal_from_shifted(
        shifted_mode_data=shifted_mode_data,
        theta=theta,
        phi=phi,
        phase_vector=phase_vector,
        wb_target=wb_target,
        inner_product=inner_product,
    )
    data_diff = np.abs(np.asarray(wa_trial - wb_target))
    return float(np.sqrt(trapezoid(data_diff**2, t_reference)))


def solve_single_point(task):
    theta, phi = task
    wb_target = get_cached_wb_target(theta, phi)

    wa_reference_modes = _MP_STATE["wa_reference_modes"]
    t_reference = _MP_STATE["t_reference"]
    dt_brute_force = _MP_STATE["dt_brute_force"]
    dphi_brute_force = _MP_STATE["dphi_brute_force"]
    phase_brute_force = _MP_STATE["phase_brute_force"]

    # brute-force seed, fully local to this sky point
    best_cost = np.inf
    best_x0 = np.array([float(dt_brute_force[0]), float(dphi_brute_force[0])], dtype=float)

    for i, dt in enumerate(dt_brute_force):
        shifted_mode_data = wa_reference_modes(t_reference + dt)
        for j, dphi in enumerate(dphi_brute_force):
            value = point_cost_from_shifted(
                shifted_mode_data=shifted_mode_data,
                theta=theta,
                phi=phi,
                phase_vector=phase_brute_force[j],
                wb_target=wb_target,
                t_reference=t_reference,
                inner_product=inner_product,
            )
            if value < best_cost:
                best_cost = value
                best_x0[:] = (float(dt), float(dphi))

    # least-squares, still local to this sky point
    def cost(x):
        dt, dphi = x
        shifted_mode_data = wa_reference_modes(t_reference + dt)
        phase_vector = np.exp(1j * _MP_STATE["m_modes"] * dphi)
        return point_cost_from_shifted(
            shifted_mode_data=shifted_mode_data,
            theta=theta,
            phi=phi,
            phase_vector=phase_vector,
            wb_target=wb_target,
            t_reference=t_reference,
            inner_product=inner_product,
        )

    optimum = least_squares(
        cost,
        x0=best_x0,
        bounds=(
            [_MP_STATE["dt_lower"], 0.0],
            [_MP_STATE["dt_upper"], 2.0 * np.pi],
        ),
        max_nfev=50000,
    )

    best_dt, best_dphi = map(float, optimum.x)
    shifted_mode_data = wa_reference_modes(t_reference + best_dt)
    phase_vector = np.exp(1j * _MP_STATE["m_modes"] * best_dphi)

    wa_best = build_transformed_point_signal_from_shifted(
        shifted_mode_data=shifted_mode_data,
        theta=theta,
        phi=phi,
        phase_vector=phase_vector,
        wb_target=wb_target,
        inner_product=inner_product,
    )

    return float(
        mismatch(
            wa_best,
            wb_target,
            x1=_MP_STATE["x1"],
            x2=_MP_STATE["x2"],
        )
    )


def compute_optimized_mismatch_at_points(
    wa,
    wb,
    t1,
    t2,
    N_theta=7,
    N_phi=5,
    optimize_polarization=True,
    is_aligned=False,
    max_δt=10,
    include_modes=None,
    nprocs=None,
    mp_start_method=None,
):
    """
    Compute the worst optimized pointwise mismatch over a set of sky locations.

    Parameters
    ----------
    wa : WaveformModes
    wb : WaveformModes
        Waveforms to be aligned
    t1 : float
    t2 : float
        Beginning and end of integration interval.
    N_theta : int, optional
    N_phi : int, optional
        Number of evenly spaced cos(theta) and phi values.
        Default is 7 and 5.
    optimize_polarization : bool, optional
        Whether or not to choose an optimal polarization angle.
        Default is True.
    is_aligned : bool, optional
        Whether or not to only consider points in the northern hemisphere.
    max_δt : float, optional
        Max δt to allow for when choosing the initial guess.
    include_modes: list, optional
        A list containing the (ell, m) modes to be included in the L² norm.
    nprocs: int, optional
        Number of cpus to use. Default is maximum number.
        If 1 is provided, then no multiprocessing is performed.
    mp_start_method : str, optional
        What start method to use ("fork" or "spawn").
        Default is fork.

    Returns
    -------
    max_mismatch : float
        Maximum optimized pointwise mismatch over a set of sky locations.
    """
    if t2 <= t1:
        raise ValueError(f"(t1, t2)=({t1}, {t2}) is out of order.")
    if wa.t[0] > t1 or wa.t[-1] < t2:
        raise ValueError(f"(t1, t2)=({t1}, {t2}) is not contained in wa.t.")
    if wb.t[0] > t1 or wb.t[-1] < t2:
        raise ValueError(f"(t1, t2)=({t1}, {t2}) is not contained in wb.t.")
    if wa.spin_weight != wb.spin_weight:
        raise ValueError(
            f"Spin weights must match, but got wa.spin_weight={wa.spin_weight} " f"and wb.spin_weight={wb.spin_weight}."
        )

    wa, wb = restrict_to_common_modes(wa, wb, include_modes)

    dt_lower = max(-max_δt, max(t1 - t2, t2 - wa.t[-1]))
    dt_upper = min(max_δt, min(t2 - t1, t1 - wa.t[0]))

    wa_peak = wa.max_norm_time()
    wb_peak = wb.max_norm_time()
    count_a = np.count_nonzero((wa.t >= wa_peak + dt_lower) & (wa.t <= wa_peak + dt_upper))
    count_b = np.count_nonzero((wb.t >= wb_peak + dt_lower) & (wb.t <= wb_peak + dt_upper))
    n_brute_force_dt = max(1, int(max(count_a, count_b)))

    dt_brute_force = np.unique(
        np.concatenate(
            [
                np.linspace(dt_lower, dt_upper, num=n_brute_force_dt),
                np.array([0.0], dtype=float),
            ]
        )
    )

    ell_min = wa.ell_min
    ell_max = wa.ell_max
    dphi_brute_force = np.linspace(
        0.0,
        2.0 * np.pi,
        num=2 * ell_max + 1,
        endpoint=False,
    )

    i1 = int(np.argmin(np.abs(wa.t - t1)))
    i2 = int(np.argmin(np.abs(wa.t - t2))) + 1
    t_reference = np.asarray(wa.t[i1:i2], dtype=float)
    x1 = float(t_reference[0])
    x2 = float(t_reference[-1])

    sky_points = [
        (float(theta), float(phi))
        for theta, phi in choose_points_on_sphere(
            is_aligned=is_aligned,
            N_theta=N_theta,
            N_phi=N_phi,
        )
    ]

    wa_mode_data = np.asarray(wa.data)
    wb_mode_data = np.asarray(wb.interpolate(t_reference).data)

    if nprocs is None:
        nprocs = mp.cpu_count()
    nprocs = int(max(1, min(nprocs, len(sky_points))))

    serial_initargs = (
        np.asarray(wa.t, dtype=float),
        wa_mode_data,
        wb_mode_data,
        t_reference,
        ell_min,
        ell_max,
        wa.spin_weight,
        optimize_polarization,
        dt_lower,
        dt_upper,
        x1,
        x2,
        dt_brute_force,
        dphi_brute_force,
        False,  # do not force BLAS threads to 1 in pure serial mode
    )

    parallel_initargs = (
        np.asarray(wa.t, dtype=float),
        wa_mode_data,
        wb_mode_data,
        t_reference,
        ell_min,
        ell_max,
        wa.spin_weight,
        optimize_polarization,
        dt_lower,
        dt_upper,
        x1,
        x2,
        dt_brute_force,
        dphi_brute_force,
        True,  # cap worker BLAS/OpenMP threads
    )

    def run_serial():
        init_mismatch_pool(*serial_initargs)
        mismatches = [solve_single_point(task) for task in sky_points]
        return max(mismatches)

    if nprocs == 1 or len(sky_points) == 1:
        return run_serial()

    ctx = mp.get_context(auto_mp_start_method(mp_start_method))
    with ctx.Pool(
        processes=nprocs,
        initializer=init_mismatch_pool,
        initargs=parallel_initargs,
    ) as pool:
        mismatches = list(pool.imap_unordered(solve_single_point, sky_points, chunksize=1))

    return max(mismatches)
