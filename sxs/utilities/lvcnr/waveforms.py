"""Functions to convert waveform quantities to LVC format"""

import numpy as np
import h5py
import sxs
import spherical_functions as sf

def convert_modes(sxs_format_waveform, metadata, out_filename,
                  modes, extrapolation_order="Extrapolated_N2",
                  log=print, truncation_time=None,
                  tolerance=5e-07, truncation_tol=None):
    """Computes amplitude and phase for an SXS-format waveform, computes
    ROM-spline for each mode, and writes to file.
    """
    from time import perf_counter
    from . import WaveformAmpPhase, Dataset

    extrap = str(extrapolation_order) + ".dir"

    if truncation_time is None:
        truncation_time = metadata['reference_time'] / (metadata["reference_mass1"] + metadata["reference_mass2"])

    h = sxs.load(sxs_format_waveform)
    if not truncation_time is None:
        h = h[np.argmin(abs(h.t - truncation_time)) + 1:]

    start_time = h.t[0]

    peak_time = h.t[np.argmax(h.norm)]

    t = h.t - peak_time

    h_amp_phase = np.empty((h.data.shape[0], 2*h.data.shape[1]))
    for mode in modes:
        L, M = mode
        h_amp_phase[:,2*sf.LM_index(L, M, h.ell_min)] = np.abs(h.data[:,sf.LM_index(L, M, h.ell_min)])
        h_amp_phase[:,2*sf.LM_index(L, M, h.ell_min) + 1] = np.unwrap(np.angle(h.data[:,sf.LM_index(L, M, h.ell_min)]))

    if truncation_tol is True:  # Test for actual identity, not just equality
        truncation_tol = 5e-2 * tolerance
    elif truncation_tol is False:
        truncation_tol = None

    with h5py.File(out_filename, 'w') as out_file:
        out_file.create_dataset('NRtimes', data=t)
        for i, mode in enumerate(modes):
            amp = h_amp_phase[:, 2*i]
            phase = h_amp_phase[:, 2*i+1]

            phase_out = Dataset.from_data(t, phase, tolerance, error_scaling=amp, truncation_tol=truncation_tol)
            amp_out = Dataset.from_data(t, amp, tolerance, truncation_tol=truncation_tol)

            out_group_amp = out_file.create_group('amp_l{0[0]}_m{0[1]}'.format(mode))
            amp_out.write(out_group_amp)
            out_group_phase = out_file.create_group('phase_l{0[0]}_m{0[1]}'.format(mode))
            phase_out.write(out_group_phase)

    return start_time, peak_time, None
