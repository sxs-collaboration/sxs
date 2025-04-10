import numpy as np
from scipy.integrate import trapezoid

from .waveform_modes import WaveformModes


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

    return wa_interp[:,ia1:ia2], wb_interp[:,ib1:ib2]


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

    already_unified = (
        np.array_equal(wa.t, wb.t)
        and (wa.ndim == wb.ndim == 1 or np.array_equal(wa.LM, wb.LM))
    )

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

    already_unified = (
        np.array_equal(wa.t, wb.t)
        and (wa.ndim == wb.ndim == 1 or np.array_equal(wa.LM, wb.LM))
    )

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
            raise ValueError(
                "A `modes` argument was provided, but the input `wa` "
                "and `wb` only have one dimension."
            )

    if ASD is not None:
        ASD_values = ASD(wa.t)
    else:
        ASD_values = 1

    wa_wb_overlap = inner_product(wa, wb, ASD_values=ASD_values).real
    wa_norm = inner_product(wa, wa, ASD_values=ASD_values).real
    wb_norm = inner_product(wb, wb, ASD_values=ASD_values).real

    return 1 - wa_wb_overlap / np.sqrt(wa_norm * wb_norm)
