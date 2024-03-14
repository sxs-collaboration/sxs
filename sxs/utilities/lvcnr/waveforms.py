"""Functions to convert waveform quantities to LVC format"""

import numpy as np
import h5py


def convert_modes(
    sxs_format_waveform,
    metadata,
    out_filename,
    modes,
    extrapolation_order="Extrapolated_N2",
    log=print,
    truncation_time=None,
    tolerance=5e-07,
    truncation_tol=None,
):
    """Computes amplitude and phase for an SXS-format waveform, computes
    ROM-spline for each mode, and writes to file.
    """
    from time import perf_counter
    from . import WaveformAmpPhase, Dataset
    from ... import load
    from ..dicts import KeyPassingDict
    
    extrap = str(extrapolation_order) + ".dir"

    if truncation_time is None:
        truncation_time = metadata["reference_time"] / (metadata["reference_mass1"] + metadata["reference_mass2"])

    h = load(sxs_format_waveform)
    if type(h) == KeyPassingDict:
        h = h[extrap]
    if not truncation_time is None:
        h = h[np.argmin(abs(h.t - truncation_time)) + 1 :]

    start_time = h.t[0]

    peak_time = h.max_norm_time()

    t = h.t - peak_time

    h_amp_phase = np.empty((h.data.shape[0], 2 * h.data.shape[1]))
    h_amp_phase[:, ::2] = h.abs
    h_amp_phase[:, 1::2] = h.arg_unwrapped

    if truncation_tol is True:  # Test for actual identity, not just equality
        truncation_tol = 5e-2 * tolerance
    elif truncation_tol is False:
        truncation_tol = None

    with h5py.File(out_filename, "w") as out_file:
        out_file.create_dataset("NRtimes", data=t)
        for i, mode in enumerate(modes):
            amp = h_amp_phase[:, 2 * i]
            phase = h_amp_phase[:, 2 * i + 1]

            phase_out = Dataset.from_data(t, phase, tolerance, error_scaling=amp, truncation_tol=truncation_tol)
            amp_out = Dataset.from_data(t, amp, tolerance, truncation_tol=truncation_tol)

            out_group_amp = out_file.create_group("amp_l{0[0]}_m{0[1]}".format(mode))
            amp_out.write(out_group_amp)
            out_group_phase = out_file.create_group("phase_l{0[0]}_m{0[1]}".format(mode))
            phase_out.write(out_group_phase)

    return start_time, peak_time, getattr(h, "version_hist", None)
