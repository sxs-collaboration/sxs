"""Functions to enable BMS transformations of waveforms

BMS transformations include the usual Poincaré group — time and space
translations, rotations, and boosts — as well as "supertranslations", which are
a more general form of translations.  Essentially, a supertranslation is a
direction-dependent time translation.

See https://arxiv.org/abs/1509.00862 for a review of BMS transformations and
their computation.

"""


def uprime_generator(u, β):
    """Return u' such that each time step is the smallest in the input time series

    Parameters
    ----------
    u : array_like
        Time steps in the original (stationary) frame
    β : float
        Magnitude of boost velocity as fraction of the speed of light

    Returns
    -------
    uprime : array_like
        Time steps in the boosted frame

    Notes
    -----
    The u' time steps in the boosted frame incorporate data from a range of
    different time slices in the stationary frame.  Because the size of time steps
    can vary as a function of time to represent dynamical data adaptively, the
    timestep sizes of those different slices will differ.  Here, we construct new
    times for the u' series so that the u' steps are no larger than the smallest
    time step on any of those slices from the input data.

    We simplify this somewhat by assuming that the time step sizes in the input `u`
    are fairly smoothly varying, so that the appropriate time step can be
    determined just by looking at the earliest and latest time slices.  This is
    only an approximation, but should be suitable for our data.

    """
    import numpy as np
    from scipy.interpolate import CubicSpline
    # uprime = u / (γ * (1 - v⃗ · n̂))
    # uprime_min = max(min(u) / (γ * (1 - v⃗ · n̂)))
    # uprime_max = min(max(u) / (γ * (1 - v⃗ · n̂)))

    γ = 1 / np.sqrt(1 - β**2)
    uprime_min = max(min(u) / (γ * (1 - β)), min(u) / (γ * (1 + β)))
    uprime_max = min(max(u) / (γ * (1 - β)), max(u) / (γ * (1 + β)))
    if uprime_max < uprime_min:
        raise ValueError(
            f"\n\tThere are no complete slices in the u' coordinate system for u ∈ [{min(u)}, {max(u)}] and β = {β}."
            f"\n\tYou may wish to decrease β or move the origin of the time coordinate closer to (u[0] + u[-1]) / 2."
        )
    uprime = [uprime_min,]
    δuprime_plus = CubicSpline(u[1:], np.diff(u / (γ * (1 + β))), extrapolate=True)
    δuprime_minus = CubicSpline(u[1:], np.diff(u / (γ * (1 - β))), extrapolate=True)
    while uprime[-1] < uprime_max:
        δuprime = min(
            δuprime_plus(uprime[-1] * γ * (1 + β)),
            δuprime_minus(uprime[-1] * γ * (1 - β))
        )
        uprime.append(min(uprime_max, uprime[-1] + δuprime))
    uprime = np.array(uprime)
    return uprime


def Bprime(v⃗, n̂prime):
    """Rotor of aberration spin-weighted fields under boosts

    Implements Equation (2) of arxiv.org/abs/1509.00862

    Parameters
    ----------
    v⃗ : (3,) array_like
        Three-vector representing the velocity of the boosted frame relative to the
        inertial frame, in units where the speed of light is 1
    n̂prime : (..., 3) array_like
        Three-vectors representing the directions in the boosted frame

    Returns
    -------
    Bprm : (..., 4) quaternionic.array
        Quaternions that rotate from the boosted frame to the stationary frame.  In
        addition to rotating n̂prime onto the corresponding n̂ directions, this will
        also rotate tangent vectors appropriately, to account for spin
        transformations.  The shape of this array is n̂prime.shape[:-1]+(4,).

    """
    import numpy as np
    import quaternionic
    ϵ = 1e-15
    β = np.linalg.norm(v⃗)
    if β < ϵ:
        return quaternionic.one
    v̂ = v⃗ / β
    φ = np.arctanh(β)
    Θprime = np.arccos(np.tensordot(v̂, n̂prime, axes=[-1, -1]))
    Θ = 2 * np.arctan(np.exp(-φ) * np.tan(Θprime/2))
    nv = np.cross(n̂prime, v̂)
    nvnorm = np.linalg.norm(nv, axis=-1)
    nv /= nvnorm[..., np.newaxis]
    return np.exp(((Θprime - Θ) / 2) * quaternionic.array.from_vector_part(nv))


def boost(w, v⃗, ell_max):
    """Find modes of waveform boosted by velocity v⃗

    Implements Equation (21) of arxiv.org/abs/1509.00862

    Parameters
    ----------
    w : WaveformModes
        Modes of waveform measured in original frame
    v⃗ : array_like
        Three-vector representing the velocity of the boosted frame relative to the
        inertial frame, in units where the speed of light is 1
    ell_max : int
        Maximum value of `ell` to use while computing the transformation, and to
        provide in the returned object

    Returns
    -------
    wprime : WaveformModes
        Modes of waveform measured in boosted frame or of modes from boosted source
        measured in original frame.  This should have the same properties as the
        input `w`, except with (1) different time data [see Notes, below], (2) a
        minimum `ell` value of 0 even for spin weight other than 0, and (3) a
        maximum `ell` value of `ell_max`.

    Notes
    -----
    Due to the nature of the transformation, some of the information in the input
    waveform must be discarded, because it corresponds to slices of the output
    waveform that are not completely represented in the input.  Thus, the times of
    the output waveform will not just be the Lorentz-transformed times of the input
    waveform.

    Depending on the magnitude β=|v⃗|, a very large value of `ell_max` may be
    needed.  The dominant factor is the translation that builds up over time:
    `β*T`, where `T` is the largest time found in the waveform.  For example, if
    β*T ≈ 1000M, we might need `ell_max=64` to maintain a comparable accuracy as in
    the input data.

    Because of the `β*T` effects, it is usually best to set t=0 at the merger time
    — best approximated as `self.max_norm_time()`.  The largest translation is then
    found early in the waveform, when the waveform is changing slowly.

    """
    import numpy as np
    import quaternionic
    import spherical
    try:
        import spinsfast
    except:
        raise ModuleNotFoundError("You need to install spinsfast manually to use this function")

    if w.data_type.lower() not in ['h', 'psi4']:
        raise NotImplementedError(f"Input waveform `w` has type {w.data_type}, which is not yet implemented")

    ϵ = 1e-15
    β = np.linalg.norm(v⃗)

    if β < ϵ:
        return w.copy()

    γ = 1 / np.sqrt(1 - β**2)
    φ = np.arctanh(β)
    v̂ = v⃗ / β

    nθprime = nϕprime = 2 * ell_max + 1
    θprimeϕprime = spherical.theta_phi(nθprime, nϕprime)
    Rprime = quaternionic.array.from_spherical_coordinates(θprimeφprime)
    n̂prime = (Rprime * quaternionic.z * Rprime.inverse).vector

    R = Bprime(v⃗, n̂prime) * Rprime
    n̂ = (R * quaternionic.z * R.inverse).vector
    doppler_factor = γ * (1 - np.tensordot(v⃗, n̂, axes=[-1, -1]))

    uprime = uprime_generator(w.t, β)
    time_axis, modes_axis = 0, 1
    Hprime = np.zeros((uprime.size, spherical.Ysize(0, ell_max)), dtype=complex)

    # Step through the waveform in segments, so that each segment is small enough to fit
    # comfortably into memory, but large enough to minimize re-computation of SWSHs
    i_step_size = 5_000
    i_padding = 20
    i_outer_1 = 0
    i_inner_1 = 0
    i_inner_2 = min(i_inner_1 + i_step_size, uprime.size)
    i_outer_2 = min(i_inner_2 + i_padding, uprime.size)
    Hprime_grid = np.zeros((i_step_size, nθprime, nϕprime), dtype=complex)
    while True:
        uprime_outer = uprime[i_outer_1:i_outer_2]
        uprime_inner = uprime[i_inner_1:i_inner_2]

        # print(f"Working on ({uprime_inner[[0, -1]]}) of ({uprime[[0, -1]]})")

        # Within each segment, this is the core computation, evaluating the transformed
        # field on a grid, and then converting back to mode weights
        for j in range(Rprime.shape[0]):
            for k in range(Rprime.shape[1]):
                Rjk = R[j, k]
                doppler_factor_jk = doppler_factor[j, k]
                u_outer_1, u_outer_2 = uprime_outer[[0, -1]] * doppler_factor_jk
                i1, i2 = w.index_closest_to(u_outer_1), w.index_closest_to(u_outer_2)
                Hprime_grid[:i_inner_2-i_inner_1, j, k] = (
                    doppler_factor_jk * w[i1:i2].evaluate(Rjk).interpolate(uprime_inner * doppler_factor_jk)
                )
        Hprime[i_inner_1:i_inner_2] = spinsfast.map2salm(Hprime_grid[:i_inner_2-i_inner_1], w.spin_weight, ell_max)

        # Move to the next segment
        if i_inner_2 == uprime.size:
            break
        i_inner_1 = i_inner_2
        i_outer_1 = max(0, i_inner_1 - i_padding)
        i_inner_2 = min(i_inner_2 + i_step_size, uprime.size)
        i_outer_2 = min(i_inner_2 + i_padding, uprime.size)

    return type(w)(
        Hprime, time=uprime, time_axis=time_axis, modes_axis=modes_axis,
        ell_min=0, ell_max=ell_max, data_type=w.data_type, spin_weight=w.spin_weight
    )
