"""Functions to convert waveform quantities to LVC format"""

import numpy as np
import h5py
import romspline


def first_index_after_time(times, target_time):
    """Returns the index of the first time in a list of times after
    time target_time."""
    return np.abs(times - target_time).argmin() + 1


def first_index_after_relaxation_time(times, metadata, offset=-2):
    """Returns the index of the first time in a list of times after the
    relaxation time, which is given as a key in metadata, i.e.
    metadata['relaxation_time'], except actually return an 
    index offset earlier."""
    relaxation_time = metadata['relaxation_time']
    return first_index_after_time(times, relaxation_time) - offset


def first_index_before_reference_time(times, metadata, offset=2):
    """Returns the index of the first time in a list of times before the
    reference time, which is given as a key in metadata, i.e.
    metadata['reference_time'], except actually return an 
    index offset earlier."""
    reference_time = metadata['reference_time']
    return first_index_after_time(times, reference_time) - offset


def waveform_norm_squared(
        sxs_format_waveform,
        extrapolation_order="Extrapolated_N2"):
    """Takes an SXS-format waveform and returns the sum of the squared
    amplitude of each (l,m) mode of the wave."""
    extrap = str(extrapolation_order) + ".dir"
    times = sxs_format_waveform[extrap]['Y_l2_m2.dat'][:, 0]
    sum_amp_squared = 0.0 * times
    for key in sxs_format_waveform[extrap].keys():
        if "Y_" in key:
            for i in range(1, 3):
                sum_amp_squared += np.square(
                    sxs_format_waveform[extrap][key][:, i])
    return sum_amp_squared


def peak_time_from_sxs(
        sxs_format_waveform,
        metadata,
        extrapolation_order='Extrapolated_N2'):
    """Returns the time when the sum of the squared amplitudes of an
    SXS-format waveform is largest. Note: this is not necessarily the time of
    the peak of the l=m=2 mode."""
    extrap = extrapolation_order + ".dir"
    # All modes have the same time, so just look at the l=m=2 mode to get the
    # times
    times = sxs_format_waveform[extrapolation_order +
                                ".dir"]['Y_l2_m2.dat'][:, 0]
    start = first_index_before_reference_time(times, metadata)
    sum_amp_squared = waveform_norm_squared(
        sxs_format_waveform, extrapolation_order)
    index_peak = start + sum_amp_squared[start:].argmax()
    return sxs_format_waveform[extrap]['Y_l2_m2.dat'][index_peak][0]


def amp_phase_from_sxs(sxs_format_waveform, metadata, modes,
                       extrapolation_order="Extrapolated_N2",
                       log=print, truncation_time=None):
    """Returns amplitude and phase for an SXS-format waveform, for a list of
    Ylm modes. If modes='all', return all modes for l=2 through l=8,
    inclusive."""
    extrap = str(extrapolation_order) + ".dir"

    if modes == "all":
        modes = [[l, m] for l in range(2, 9) for m in range(-l, l+1)]
        
    log("Modes: " + str(modes))
    amps = []
    phases = []
    times_list = []
    l_max = 0
    # All modes have the same time, so just look at the l=m=2 mode to get
    # the times
    times = sxs_format_waveform[extrap]['Y_l2_m2.dat'][:, 0]
    if truncation_time is None:
        start = first_index_before_reference_time(times, metadata)
    else:
        start = first_index_after_time(times, truncation_time)
    peak = peak_time_from_sxs(
        sxs_format_waveform, metadata, extrapolation_order)
    for mode in modes:
        l = mode[0]
        m = mode[1]
        log("Computing mode: l = " + str(l) + ", m = " + str(m))
        mode = "Y_l" + str(l) + "_m" + str(m) + ".dat"
        hlm = sxs_format_waveform[extrap][mode]

        # CHECK ME: is the + sign correct here?
        h = hlm[start:, 1] + 1j * hlm[start:, 2]

        amps.append(np.abs(h))
        phases.append(np.unwrap(np.angle(h)))
        times_list.append(times[start:] - peak)

        if l > l_max:
            l_max = l

    return modes, times_list, amps, phases, times[start], peak, l_max


def spline_amp_phase_from_sxs(sxs_format_waveform, metadata, modes,
                              extrapolation_order="Extrapolated_N2",
                              log=print, truncation_time=None):
    """Returns spline amplitude and phase for an SXS-format waveform, for a
    list of Ylm modes. If modes='all', return all modes for l=2 through l=8,
    inclusive."""
    modes, times, amps, phases, start_time, peak_time, l_max \
        = amp_phase_from_sxs(sxs_format_waveform, metadata, modes,
                             extrapolation_order, log, truncation_time)
    spline_amps = []
    spline_phases = []
    for i, mode in enumerate(modes):
        log("Computing spline for amplitude of mode " + str(mode))
        spline_amps.append(romspline.ReducedOrderSpline(times[i], amps[i]))
        log("Computing spline for phase of mode " + str(mode))
        spline_phases.append(romspline.ReducedOrderSpline(times[i], phases[i]))
    return modes, times, spline_amps, spline_phases, start_time, peak_time, \
        l_max


def write_splines_to_H5(
        out_filename,
        modes,
        spline_amps,
        spline_phases,
        times,
        log=print):
    """Writes spline amplitudes and phases to an HDF5 file
    named out_filename."""
    log("Writing waveform data to " + str(out_filename))
    with h5py.File(out_filename, 'w') as out_file:
        for i, mode in enumerate(modes):
            l = mode[0]
            m = mode[1]
            out_group_amp = out_file.create_group('amp_l%d_m%d' % (l, m))
            out_group_phase = out_file.create_group('phase_l%d_m%d' % (l, m))
            spline_amps[i].write(out_group_amp)
            spline_phases[i].write(out_group_phase)
            if l == 2 and m == 2:
                out_file.create_dataset('NRtimes', data=times[i])
