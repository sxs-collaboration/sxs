"""Functions to convert waveform quantities to LVC format"""

import numpy as np
import h5py


def convert_modes(sxs_format_waveform, metadata, out_filename,
                  modes, extrapolation_order="Extrapolated_N2",
                  log=print, truncation_time=None,
                  tolerance=5e-07):
    """Computes amplitude and phase for an SXS-format waveform, computes
    ROM-spline for each mode, and writes to file.
    """
    from time import perf_counter
    from .. import WaveformAmpPhase
    from . import LVCDataset

    extrap = str(extrapolation_order) + ".dir"

    if truncation_time is None:
        truncation_time = metadata['reference_time']

    mode_strings = ["Y_l{0[0]}_m{0[1]}.dat".format(ellm) for ellm in modes]
    h = WaveformAmpPhase.read(sxs_format_waveform, extrap, mode_strings=mode_strings, start_time=truncation_time)

    start_time = h.t[0]
    
    peak_time = h.t[h.peak_index()]

    t = h.t - peak_time

    with h5py.File(out_filename, 'w') as out_file:
        out_file.create_dataset('NRtimes', data=t)
        for i, mode in enumerate(modes):
            amp = h.data[:, 2*i]
            phase = h.data[:, 2*i+1]

            phase_out = LVCDataset.from_data(t, phase, tolerance, error_scaling=amp)
            amp_out = LVCDataset.from_data(t, amp, tolerance)

            out_group_amp = out_file.create_group('amp_l{0[0]}_m{0[1]}'.format(mode))
            amp_out.write(out_group_amp)
            out_group_phase = out_file.create_group('phase_l{0[0]}_m{0[1]}'.format(mode))
            phase_out.write(out_group_phase)

    return start_time, peak_time, h.version_hist
