"""Utilities for computing waveform fluxes"""

# This code originally came from scri https://github.com/moble/scri, from the
# file scri/flux.py

import functools
import numpy as np
import numba
from spherical import clebsch_gordan as CG
from spherical import LM_index
from .waveform_modes import WaveformModes
from .. import jit
from .. import TimeSeries


def swsh_indices_to_matrix_indices(matrix_iterator):
    """Convert SWSH-indexed function into caching sparse matrix-indexed function

    This function is designed to decorate functions with signature (ell_min,
    ell_max) that yield tuples of (ellp, mp, ell, m, value), and convert them into
    functions with the same signature but returning tuples of (rows, columns,
    values).  Here, the rows are row indices of a matrix, corresponding to the
    (ellp, mp) in standard `spherical_functions` order.  Similarly the columns are
    column indices corresponding to (ell, m).  The result can be passed into the
    `sparse_expectation_value` function.

    """

    @functools.lru_cache()  # Decorate the decorated function with caching
    def wrapper(ell_min, ell_max, *args, **kwargs):
        rows, columns, values = zip(
            *(  # zip* is the inverse of zip
                (LM_index(ellp, mp, ell_min), LM_index(ell, m, ell_min), value)
                for ellp, mp, ell, m, value in matrix_iterator(ell_min, ell_max, *args, **kwargs)
            )
        )
        return np.array(rows, dtype=int), np.array(columns, dtype=int), np.array(values)

    functools.update_wrapper(wrapper, matrix_iterator)  # Copy over the name, docstring, and such

    return wrapper


@jit
def sparse_expectation_value(abar, rows, columns, values, b):
    """Helper function for the main calculation of `matrix_expectation_value`

    Computes ⟨a|M|b⟩, assuming that ⟨a| and |b⟩ have the same shapes, and that
    (rows, columns, values) are corresponding arrays containing the nonzero
    elements of the matrix M.

    Parameters
    ----------
    abar : ndarray
        The data for ⟨a| (after complex conjugation was applied).

    rows, columns : list
        Lists of same length, containing indices into (respectively)
        abar and b.

    values : ndarray
        Should have shape (n,) where n is the same length as rows and
        columns.  These values are the matrix elements themselves to
        be summed.

    b : ndarray
        The data for |b⟩. Must have the same shape as abar.

    Returns
    -------
    expectation_value : ndarray

    """
    n_times = abar.shape[0]
    n_elements = rows.shape[0]
    expectation_value = np.zeros(n_times, dtype=numba.complex128)
    for i_time in range(n_times):
        for i_element in range(n_elements):
            expectation_value[i_time] += (
                abar[i_time, rows[i_element]] * b[i_time, columns[i_element]] * values[i_element]
            )
    return expectation_value


def matrix_expectation_value(a, M, b, allow_LM_differ=False, allow_times_differ=False):
    """The matrix expectation value ⟨a|M|b⟩(u), where M is a linear operator

    Treat two spin-s waveforms a, b (in the modal representation) as vectors.
    Then we can form the 'expectation value'

      ⟨a|M|b⟩(u) ≡ ∑ⱼₙₗₘ āⱼₙ(u) Mⱼₙₗₘ bₗₘ(u).

    For efficiency with sparse matrices, M should be implemented as a generator
    which yields tuples (j, n, l, m, Mⱼₙₗₘ) for only the non-vanishing matrix
    elements.

    Parameters
    ----------
    a : WaveformModes object
        The waveform |a⟩. This function is "antilinear" in `a`.

    M : callable
        M will be called like
        `for ellp, mp, ell, m, M_el in M(ell_min, ell_max):`
        This is best implemented as a generator yielding tuples
        `(ellp, mp, ell, m, M_{ellp mp ell m})`
        which range over ell, ellp values up to and including ell_max,
        for only the non-vanishing matrix elements of M.

    b : WaveformModes object
        The waveform |b⟩. This function is linear in `b`.

    allow_LM_differ : bool, optional [default: False]
        If True and if the set of (ell,m) modes between a and b
        differ, then the inner product will be computed using the
        intersection of the set of modes.

    allow_times_differ: bool, optional [default: False]
        If True and if the set of times between a and b differ,
        then waveform B will be interpolated to the
        times of waveform A, after intersecting the times.

    Returns
    -------
    TimeSeries
        The timeseries of ⟨a|M|b⟩(u).

    """
    if a.spin_weight != b.spin_weight:
        raise ValueError("Spin weights must match in matrix_expectation_value")

    (ell_min, ell_max) = (a.ell_min, a.ell_max)
    if (a.ell_min != b.ell_min) or (a.ell_max != b.ell_max):
        if allow_LM_differ:
            (ell_min, ell_max) = (max(a.ell_min, b.ell_min), min(a.ell_max, b.ell_max))
            if ell_min > ell_max:
                raise ValueError("Intersection of (ell,m) modes is empty. Assuming this is not desired.")
        else:
            raise ValueError(
                "ell_min and ell_max must match in matrix_expectation_value (use allow_LM_differ=True to override)"
            )

    A = a[:, a.index(ell_min, -ell_min): 1+a.index(ell_max, +ell_max) ]
    B = b[:, b.index(ell_min, -ell_min): 1+b.index(ell_max, +ell_max) ]

    if not np.array_equal(a.t, b.t):

        if not allow_times_differ:
            raise ValueError(
                "Time samples must match in matrix_expectation_value (use allow_times_differ=True to override)"
            )

        # Clip the times
        (t1, t2) = (max(A.t[0], B.t[0]), min(A.t[-1], B.t[-1]))

        idxA1 = np.argmin(abs(A.t - t1))
        idxA2 = np.argmin(abs(A.t - t2)) + 1

        idxB1 = np.argmin(abs(B.t - t1))
        idxB2 = np.argmin(abs(B.t - t2)) + 1

        A = A[idxA1:idxA2, :]
        B = B[idxB1:idxB2, :]

        if not np.array_equal(A.t, B.t):
            # Need to interpolate B onto the times of A
            B = B.interpolate(A.t)

    ##########
    ## Time for actual math!

    rows, columns, values = M(ell_min, ell_max)

    return TimeSeries(
        sparse_expectation_value(np.conj(A), rows, columns, values, B),
        time=A.t)


def energy_flux(h):
    """Compute energy flux from waveform

    This implements Eq. (2.8) from Ruiz et al. (2008) [0707.4654].

    """

    if not isinstance(h, WaveformModes):
        raise ValueError(
            f"Energy flux can only be calculated from a `WaveformModes` object; this object is of type `{type(h)}`."
        )
    if h.data_type == "hdot":
        hdot = h
    elif h.data_type == "h":
        hdot = h.dot
    else:
        raise ValueError(
            f"Input argument is expected to have data_type of 'h' or 'hdot'; this waveform data has {h.data_type=}"
        )

    # No need to use matrix_expectation_value here
    Edot = np.einsum("ij, ij -> i", hdot.conjugate(), hdot).real

    Edot /= 16.0 * np.pi

    return TimeSeries(Edot, time=h.t)


@swsh_indices_to_matrix_indices
def p_z(ell_min, ell_max, s=-2):
    """Generator for pᶻ matrix elements (for use with `matrix_expectation_value`)

    This function is specific to the case where waveforms have s=-2.

      pᶻ = cosθ = 2 √(π/3) Y₁₀

    This is what Ruiz et al. (2008) [0707.4654] call "lᶻ".

    The matrix elements yielded are

      ⟨s,j,n|cosθ|s,l,m⟩ = √[(2l+1)/(2j+1)] ⟨l,m,1,0|j,m⟩ ⟨l,-s,1,0|j,-s⟩

    where the terms on the last line are the ordinary Clebsch-Gordan coefficients.
    Because of the magnetic selection rules, we only have nonzero elements for
    n==m.

    We could have used `_swsh_Y_mat_el` but I am just preemptively combining the
    prefactors.

    """

    for ell in range(ell_min, ell_max + 1):
        ellp_min = max(ell_min, ell - 1)
        ellp_max = min(ell_max, ell + 1)
        for ellp in range(ellp_min, ellp_max + 1):
            for m in range(-ell, ell + 1):
                if (m < -ellp) or (m > ellp):
                    continue
                cg1 = CG(ell, m, 1, 0, ellp, m)
                cg2 = CG(ell, -s, 1, 0, ellp, -s)
                prefac = np.sqrt((2.0 * ell + 1.0) / (2.0 * ellp + 1.0))
                yield ellp, m, ell, m, (prefac * cg1 * cg2)


@swsh_indices_to_matrix_indices
def p_plusminus(ell_min, ell_max, sign, s=-2):
    """Produce the function p_plus or p_minus, based on sign.

      p⁺ = -√(8π/3) Y₁₊₁ = sinθ exp[+iϕ]
      p⁻ = +√(8π/3) Y₁₋₁ = sinθ exp[-iϕ]

    These are what Ruiz et al. (2008) [0707.4654] call "l⁺" and "l⁻".

    We use `swsh_Y_mat_el` to compute the matrix elements.  Notice that since the
    operator has definite m=±1, we only have mp==m±1 nonvanishing in the
    matrix elements.

    """

    if (sign != 1) and (sign != -1):
        raise ValueError("sign must be either 1 or -1 in j_plusminus")

    prefac = -1.0 * sign * np.sqrt(8.0 * np.pi / 3.0)

    def swsh_Y_mat_el(s, l3, m3, l1, m1, l2, m2):
        """Compute a matrix element treating Yₗₘ as a linear operator

        From the rules for the Wigner D matrices, we get the result that

          ⟨s,l₃,m₃|Yₗ₁ₘ₁|s,l₂,m₂⟩ =
            √[(2l₁+1)(2l₂+1)/(4π(2l₃+1))] ⟨l₁,m₁,l₂,m₂|l₃,m₃⟩ ⟨l₁,0,l₂,−s|l₃,−s⟩

        where the terms on the last line are the ordinary Clebsch-Gordan
        coefficients.  See, e.g., Campbell and Morgan (1971).

        """
        cg1 = CG(l1, m1, l2, m2, l3, m3)
        cg2 = CG(l1, 0.0, l2, -s, l3, -s)

        return np.sqrt((2.0 * l1 + 1.0) * (2.0 * l2 + 1.0) / (4.0 * np.pi * (2.0 * l3 + 1))) * cg1 * cg2

    for ell in range(ell_min, ell_max + 1):
        ellp_min = max(ell_min, ell - 1)
        ellp_max = min(ell_max, ell + 1)
        for ellp in range(ellp_min, ellp_max + 1):
            for m in range(-ell, ell + 1):
                mp = round(m + 1 * sign)
                if (mp < -ellp) or (mp > ellp):
                    continue
                yield (ellp, mp, ell, m, (prefac * swsh_Y_mat_el(s, ellp, mp, 1.0, sign, ell, m)))


p_plus = functools.partial(p_plusminus, sign=+1)
p_minus = functools.partial(p_plusminus, sign=-1)
p_plus.__doc__ = p_plusminus.__doc__
p_minus.__doc__ = p_plusminus.__doc__


def momentum_flux(h):
    """Compute momentum flux from waveform

    This implements Eq. (2.11) from Ruiz et al. (2008) [0707.4654] by using
    `matrix_expectation_value` with `p_z`, `p_plus`, and `p_minus`.

    """

    if not isinstance(h, WaveformModes):
        raise ValueError(
            f"Momentum flux can only be calculated from a `WaveformModes` object; this object is of type `{type(h)}`."
        )
    if h.data_type == "hdot":
        hdot = h
    elif h.data_type == "h":
        hdot = h.dot
        hdot._metadata["data_type"] = "hdot"

    else:
        raise ValueError(
            f"Input argument is expected to have data of type `h` or `hdot`; this waveform data has `{h.data_type=}`"
        )

    p_plus_dot = matrix_expectation_value(hdot, functools.partial(p_plus, s=-2), hdot)
    p_minus_dot = matrix_expectation_value(hdot, functools.partial(p_minus, s=-2), hdot)
    p_z_dot = matrix_expectation_value(hdot, functools.partial(p_z, s=-2), hdot)

    # Convert into (x,y,z) basis
    pdot = np.empty((hdot.n_times, 3), dtype=float)
    pdot[:, 0] = 0.5 * (p_plus_dot.real + p_minus_dot.real)
    pdot[:, 1] = 0.5 * (p_plus_dot.imag - p_minus_dot.imag)
    pdot[:, 2] = p_z_dot.real

    pdot /= 16.0 * np.pi

    return TimeSeries(pdot, time=h.t)


@swsh_indices_to_matrix_indices
def j_z(ell_min, ell_max):
    r"""Generator for jᶻ matrix elements (for use with matrix_expectation_value)

    Matrix elements yielded are

      ⟨j,n|jᶻ|l,m⟩ = i m δⱼₗ δₙₘ

    This follows the convention for jᶻ from Ruiz et al. (2008) [0707.4654]

    """
    for ell in range(ell_min, ell_max + 1):
        for m in range(-ell, ell + 1):
            yield ell, m, ell, m, (1.0j * m)


@swsh_indices_to_matrix_indices
def j_plusminus(ell_min, ell_max, sign):
    r"""Produce the function j_plus or j_minus, based on sign.

    The conventions for these matrix elements, to agree with Ruiz et al. (2008)
    [0707.4654], should be:

      ⟨j,n|j⁺|l,m⟩ = i √[(l-m)(l+m+1)] δⱼₗ δₙₘ₊₁
      ⟨j,n|j⁻|l,m⟩ = i √[(l+m)(l-m+1)] δⱼₗ δₙₘ₋₁
    """

    if (sign != 1) and (sign != -1):
        raise ValueError("sign must be either 1 or -1 in j_plusminus")

    for ell in range(ell_min, ell_max + 1):
        for m in range(-ell, ell + 1):
            mp = round(m + 1 * sign)
            if (mp < -ell) or (mp > ell):
                continue
            yield (ell, mp, ell, m, (1.0j * np.sqrt((ell - m * sign) * (ell + m * sign + 1.0))))


j_plus = functools.partial(j_plusminus, sign=+1)
j_minus = functools.partial(j_plusminus, sign=-1)
j_plus.__doc__ = j_plusminus.__doc__
j_minus.__doc__ = j_plusminus.__doc__


def angular_momentum_flux(h, hdot=None):
    """Compute angular momentum flux from waveform

    This implements Eq. (2.24) from Ruiz et al. (2008) [0707.4654] by using
    `matrix_expectation_value` with (internal) helper functions `j_z`, `j_plus`, and
    `j_minus`.

    """

    if not isinstance(h, WaveformModes):
        raise ValueError(
            f"Angular momentum flux can only be calculated from a `WaveformModes` object; `h` is of type `{type(h)}`."
        )
    if (hdot is not None) and (not isinstance(hdot, WaveformModes)):
        raise ValueError(
            f"Angular momentum flux can only be calculated from a `WaveformModes` object; `hdot` is of type `{type(hdot)}`."
        )
    if h.data_type == "h":
        if hdot is None:
            hdot = h.dot
            hdot._metadata["data_type"] = "hdot"
        elif hdot.data_type != "hdot":
            raise ValueError(
                f"Input argument `hdot` is expected to have data_type=='hdot'; this `hdot` waveform data has `{hdot.data_type=}`"
            )
    else:
        raise ValueError(
            f"Input argument `h` is expected to have data of type `h`; this `h` waveform data has `{h.data_type=}`"
        )

    j_plus_dot = matrix_expectation_value(hdot, j_plus, h)
    j_minus_dot = matrix_expectation_value(hdot, j_minus, h)
    j_z_dot = matrix_expectation_value(hdot, j_z, h)

    # Convert into (x,y,z) basis
    jdot = np.empty((hdot.n_times, 3), dtype=float)
    jdot[:, 0] = 0.5 * (j_plus_dot.real + j_minus_dot.real)
    jdot[:, 1] = 0.5 * (j_plus_dot.imag - j_minus_dot.imag)
    jdot[:, 2] = j_z_dot.real

    jdot /= -16.0 * np.pi

    return TimeSeries(jdot, time=h.t)


def poincare_fluxes(h, hdot=None):
    """Compute fluxes of energy, momentum, and angular momentum

    This function will compute the time derivative 1 or 0 times (if an optional
    argument is passed) so is more efficient than separate calls to `energy_flux`,
    `momentum_flux`, and `angular_momentum_flux`.

    Parameters
    ----------
    h : WaveformModes object
        Must have data_type=='h'

    hdot : WaveformModes object, optional [default: None]
        The time derivative of h.  If None, computed from h.  Must
        have data_type=='hdot'.

    Returns
    -------
    edot, pdot, jdot : TimeSeries, TimeSeries, TimeSeries

    """

    if not isinstance(h, WaveformModes):
        raise ValueError(
            f"Poincare fluxes can only be calculated from a `WaveformModes` object; `h` is of type `{type(h)}`."
        )
    if (hdot is not None) and (not isinstance(hdot, WaveformModes)):
        raise ValueError(
            f"Poincare fluxes can only be calculated from a `WaveformModes` object; `hdot` is of type `{type(hdot)}`."
        )
    if h.data_type == "h":
        if hdot is None:
            hdot = h.dot
            hdot._metadata["data_type"] = "hdot"
        elif hdot.data_type != "hdot":
            raise ValueError(
                f"Input argument `hdot` is expected to have data_type='hdot'; this `hdot` waveform data has `{hdot.data_type=}`"
            )
    else:
        raise ValueError(
            f"Input argument `h` is expected to have data of type `h`; this `h` waveform data has `{h.data_type=}`"
        )

    return (energy_flux(hdot), momentum_flux(hdot), angular_momentum_flux(h, hdot))
