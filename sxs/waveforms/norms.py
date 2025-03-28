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
    """Make two WaveformModes objects share a common time array
    that contains the times t1 and t2.

    Parameters
    ----------
    wa : WaveformModes
    wb : WaveformModes
    t1 : float
    t2 : float
    padding_time_factor : float, optional
        Extra time of size (t2 - t1)*padding_time_factor to include in the waveforms.
        Default is 0.2.

    Returns
    -------
    wa_interp : WaveformModes
    wa_interp : WaveformModes
    """
    check_time_constraint(wa, wb, t1, t2)

    padding_time = (t2 - t1) * padding_time_factor
    t1_padded = max(wa.t[0], wb.t[0], t1 - padding_time)
    t2_padded = min(wa.t[-1], wb.t[-1], t2 + padding_time)

    idx1 = np.argmin(abs(wa.t - t1_padded))
    idx2 = np.argmin(abs(wa.t - t2_padded)) + 1

    wa_interp = wa.interpolate(wa.t[idx1:idx2])
    wb_interp = wb.interpolate(wa_interp.t)

    return wa_interp, wb_interp


def compute_L2_norm(wa, wb, t1=None, t2=None, modes=None, modes_for_norm=None, normalize=True):
    """Compute the L² norm of the residual between two waveforms over the
    time window (t1, t2) using the modes specified by `modes`.

    Parameters
    ----------
    wa : WaveformModes
    wb : WaveformModes
        Interpolated to wa over (t1, t2).
    t1 : float
        Beginning of L² norm integral.
    t2 : float
        End of L² norm integral.
    modes : list, optional
        Modes (ell, m) to include in numerator of L² norm calculation.
        Default is all modes.
    modes_for_norm : list, optional
        Modes (ell, m) to include in denominator of L² norm calculation.
        Default is all modes.
    normalize : bool, optional
        Whether or not to divide by sqrt(||wa||²||wb||²) in the sqrt. If False, returns
        the unnormalized L² norm of the residual and what would have been the normalization.
        Default is True.

    Returns
    -------
    L2_norm: float
        L² norm of the residual between the waveforms
        This is sqrt( ||wa - wb||² / sqrt(||wa||² ||wb||²) )
    """
    wa = wa.copy()
    wb = wb.copy()

    if t1 is None:
        t1 = max(wa.t[0], wb.t[0])

    if t2 is None:
        t2 = min(wa.t[-1], wb.t[-1])

    if not (wa.t.size == wb.t.size and all(wa.t == wb.t)):
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


def compute_inner_product(wa, wb, modes=None, ASD_values=None):
    """Compute the inner product between two waveforms over the
    window (x1, x2) using the modes specified by `modes`.

    Note that this can be provided with time-domain or frequency-domain waveforms,
    either over the two-sphere (WaveformModes) or at a point on the two-sphere (TimeSeries).

    For frequency-domain waveforms, the .t attribute should be the frequency.

    Parameters
    ----------
    wa : WaveformModes or TimeSeries
    wb : WaveformModes or TimeSeries
    modes : list, optional
        Default is all modes.
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
            np.sum(wa.data * wb.bar.data, axis=1) / ASD_values**2,
            wa.t,
        )
    else:
        inner_product = trapezoid(
            (wa * np.conjugate(wb)).ndarray / ASD_values**2,
            wa.t,
        )

    return inner_product


def compute_mismatch(wa, wb, x1=None, x2=None, modes=None, ASD=None):
    """Compute the mismatch between two waveforms over the
    window (x1, x2) using the modes specified by `modes`.

    Note that this can be provided with time-domain or frequency-domain waveforms,
    either over the two-sphere (WaveformModes) or at a point on the two-sphere (TimeSeries).

    For frequency-domain waveforms, the .t attribute should be the frequency.

    Parameters
    ----------
    wa : WaveformModes or TimeSeries
    wb : WaveformModes or TimeSeries
    x1 : float
        Beginning of mismatch integral.
    x2 : float
        End of mismach integral.
    modes : list, optional
        Default is all modes.
    ASD : func, optional
        Function mapping frequencies to the ASD of a detector.
        Default is flat ASD.

    Returns
    -------
    mismatch : float
        Mismatch between the two waveforms.
    """
    wa = wa.copy()
    wb = wb.copy()

    if x1 is None:
        x1 = max(wa.t[0], wb.t[0])

    if x2 is None:
        x2 = min(wa.t[-1], wb.t[-1])

    if not (wa.t.size == wb.t.size and all(wa.t == wb.t)):
        wa, wb = create_unified_waveforms(wa, wb, x1, x2, padding_time_factor=0)

    # Eliminate unwanted modes
    if modes is not None:
        ell_min = min(wa.ell_min, wb.ell_min)
        ell_max = max(wa.ell_max, wb.ell_max)
        for L in range(ell_min, ell_max + 1):
            for M in range(-L, L + 1):
                if not (L, M) in modes:
                    wa.data[:, wa.index(L, M)] *= 0
                    wb.data[:, wb.index(L, M)] *= 0

    wa = wa[np.argmin(abs(wa.t - x1)) : np.argmin(abs(wa.t - x2)) + 1]
    wb = wb[np.argmin(abs(wb.t - x1)) : np.argmin(abs(wb.t - x2)) + 1]

    if ASD is not None:
        ASD_values = ASD(wa.t)
    else:
        ASD_values = 1

    wa_wb_overlap = compute_inner_product(wa, wb, ASD_values=ASD_values).real
    wa_norm = compute_inner_product(wa, wa, ASD_values=ASD_values).real
    wb_norm = compute_inner_product(wb, wb, ASD_values=ASD_values).real

    mismatch = 1 - wa_wb_overlap / np.sqrt(wa_norm * wb_norm)

    return mismatch
