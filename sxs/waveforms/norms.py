import numpy as np
from quaternion.calculus import indefinite_integral as integrate

from .waveform_mts import MTS

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
        raise ValueError(
            f"{t1} is not contained in intersection of waveform time arrays: {time_intersection}."
        )
    if t2 > time_intersection[1]:
        raise ValueError(
            f"{t2} is not contained in intersection of waveform time arrays: {time_intersection}."
        )

def create_unified_waveforms(wa, wb, t1, t2, padding_time_factor=0.2):
    """Make two sxs.WaveformModes objects share a common time array
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
    t1_padded = max(max(wa.t[0], wb.t[0]), t1 - padding_time)
    t2_padded = min(min(wa.t[-1], wb.t[-1]), t2 + padding_time)

    idx1 = np.argmin(abs(wa.t - t1_padded))
    idx2 = np.argmin(abs(wa.t - t2_padded)) + 1

    wa_interp = wa.interpolate(wa.t[idx1:idx2])
    wb_interp = wb.interpolate(wa_interp.t)

    return wa_interp, wb_interp

def compute_L2_norm(
    wa, wb, t1=None, t2=None, modes=None, modes_for_norm=None, normalize=True
):
    """Compute the L2 norm between two waveforms over the
    time window (t1, t2) using the modes specified by `modes`.

    Parameters
    ----------
    wa : WaveformModes
    wb : WaveformModes
        Interpolated to wa over (t1, t2).
    t1 : float
        Beginning of L2 norm integral.
    t2 : float
        End of L2 norm integral.
    modes : list, optional
        Modes (\ell, m) to include in numerator of L2 norm calculation.
        Default is all modes.
    modes_for_norm : list, optional
        Modes (\ell, m) to include in denominator of L2 norm calculation.
        Default is all modes.
    normalize : bool, optional
        Whether or not to divide by ||wa||² in the sqrt. If False, returns
        the unnormalized L2 norm and what would have been the normalization.
        Default is True.

    Returns
    -------
    L2_norm: float
        L2 norm between to the waveforms
        This is sqrt( ||wa - wb||² / ||wa||² )
    """
    wa = wa.copy()
    wb = wb.copy()

    if t1 is None:
        t1 = max(wa.t[0], wb.t[0])
        
    if t2 is None:
        t2 = min(wa.t[-1], wb.t[-1])

    wa, wb = create_unified_waveforms(
        wa, wb, t1, t2, padding_time_factor=0
    )

    data_diff = wa.data - wb.data

    data_for_norm = wa.data
    
    # Eliminate unwanted modes
    if modes is not None or modes_for_norm is not None:
        for L in range(wa.ell_min, wa.ell_max + 1):
            for M in range(-L, L + 1):
                if modes is not None:
                    if not (L, M) in modes:
                        data_diff[:, wa.index(L, M)] *= 0
                if modes_for_norm is not None:
                    if not (L, M) in modes_for_norm:
                        data_for_norm[:, wa.index(L, M)] *= 0
                        
    if normalize:
        L2_norm = np.sqrt(
            integrate(np.linalg.norm(data_diff, axis=wa.modes_axis)**2, wa.t)[-1]
            / integrate(np.linalg.norm(data_for_norm, axis=wa.modes_axis)**2, wa.t)[-1]
        )

        return L2_norm
    else:
        L2_norm_unnormalized = np.sqrt(integrate(np.linalg.norm(data_diff, axis=wa.modes_axis)**2, wa.t)[-1])
        norm = np.sqrt(integrate(np.linalg.norm(data_for_norm, axis=wa.modes_axis)**2, wa.t)[-1])

        return L2_unnormalized, norm

def overlap(wa_mts, wb_mts, t1, t2, modes=None, ASD=None):
    """
    Compute the overlap between two waveforms over the
    time window (t1, t2) using the modes specified by `modes`.

    If an ASD curve is provided, then compute the noise-weighted overlap.

    Parameters
    ----------
    wa_mts : ModesTimeSeries
    wb_mts : ModesTimeSeries
    t1 : float
        Beginning of overlap integral.
    t2 : float
        End of overlap integral.
    modes : list, optional
        Default is all modes.
    ASD : ndarray, optional
        Array containing ASD information for a detector.

    Returns
    -------
    overlap : float
        Overlap between the two waveforms.
    """
    overlap = integrate(
        np.sqrt(4 * np.pi)
        * wa_mts.multiply(wb_mts.bar, truncator=lambda tup: 0).ndarray[:, 0]
        / ASD ** 2,
        wa_mts.t,
    )[-1]

    return overlap


def compute_mismatch(wa, wb, t1=None, t2=None, modes=None, ASD=None):
    """
    Compute the mismatch between two waveforms over the
    time window (t1, t2) using the modes specified by `modes`.

    Parameters
    ----------
    wa : WaveformModes
    wb : WaveformModes
    t1 : float
        Beginning of L2 norm integral.
    t2 : float
        End of L2 norm integral.
    modes : list, optional
        Default is all modes.
    ASD : ndarray, optional
        Array containing requency and ASD information for a detector.
        This should be of shape (..., 2).

    Returns
    -------
    mismatch : float
        Mismatch between the two waveforms.
    """
    wa = wa.copy()
    wb = wb.copy()

    if t1 is None:
        t1 = max(wa.t[0], wb.t[0])
        
    if t2 is None:
        t2 = min(wa.t[-1], wb.t[-1])

    wa, wb = create_unified_waveforms(
        wa, wb, t1, t2, padding_time_factor=0
    )
        
    # Eliminate unwanted modes
    if modes is not None:
        ell_min = min(wa.ell_min, wb.ell_min)
        ell_max = max(wa.ell_max, wb.ell_max)
        for L in range(ell_min, ell_max + 1):
            for M in range(-L, L + 1):
                if not (L, M) in modes:
                    wa.data[:, wa.index(L, M)] *= 0
                    wb.data[:, wb.index(L, M)] *= 0
                    
    wa_mts = MTS(wa)
    wb_mts = MTS(wb)

    if ASD is None:
        ASD = np.ones_like(wa_mts.t)
    
    wa_wb_overlap = overlap(wa_mts, wb_mts, t1, t2, modes, ASD).real
    wa_norm = overlap(wa_mts, wa_mts, t1, t2, modes, ASD).real
    wb_norm = overlap(wb_mts, wb_mts, t1, t2, modes, ASD).real

    mismatch = 1 - wa_wb_overlap / np.sqrt(wa_norm * wb_norm)

    return mismatch
