"""Functions to convert horizon quantities to LVC format"""

import numpy as np
import h5py


def prepare_horizon_quantity(sxs_horizon_quantity, start_time, peak_time):
    """Return times and values of an SXS-format horizon quantity.

    Horizon quantities are things like AhA.dir/ArealMass.dat. This function
    first truncates the horizon data, including only data after the reference
    time. Then, it shifts the time by the same amount as the waveforms (i.e.,
    by the peak time). Then, return the truncated/shifted times and truncated
    values.

    """
    # First, figure out the correct time series
    times_raw_AH = sxs_horizon_quantity[:, 0]
    start_AH = np.argmin(np.abs(times_raw_AH - start_time))
    times_AH = times_raw_AH[start_AH:] - peak_time

    # Loop over remaining components, truncating each one to match times_AH
    quantity_AH_list = []
    for i in range(1, len(sxs_horizon_quantity[0])):
        quantity_AH_list.append(sxs_horizon_quantity[start_AH:, i])
    quantity_AH = np.array(quantity_AH_list)
    return times_AH, quantity_AH


def spline_horizon_quantity(sxs_horizon_quantity, start_time, peak_time, truncation_tol=None):
    """Return spline of a horizon quantity.

    Prepares `sxs_horizon_quantity` by passing it to `prepare_horizon_quantity`
    and then returns a spline of the result.

    """
    from . import Dataset
    tol = 1e-6
    if truncation_tol is True:  # Test for actual identity, not just equality
        truncation_tol = 5e-2 * tol
    elif truncation_tol is False:
        truncation_tol = None
    times_AH, quantity_AH = prepare_horizon_quantity(sxs_horizon_quantity, start_time, peak_time)
    spline_AH_list = [
        Dataset.from_data(times_AH, quantity_AH[i], tol=tol, truncation_tol=truncation_tol)
        for i in range(0, len(sxs_horizon_quantity[0]) - 1)
    ]
    return spline_AH_list


def insert_spline(sxs_horizons, spline_dictionary, spline_keys,
                  horizon_key, quantity_key, start_time, peak_time, log=print, truncation_tol=None):
    """Insert spline of horizon quantity into dictionary of horizon splines.

    Note: spline_keys is a vector of key names (should be length 1 for scalars,
    length 3 for vectors), where each key is the name of a group that will be
    written in the LVC file, such as mass1-vs-time or spin1x-vs-time;
    horizon_key is AhA, AhB, or AhC; and quantity_key is, e.g.,
    ChristodoulouMass, the name of the quantity in Horizons.h5 to be read and
    splined.

    """
    log("Computing " + str(spline_keys))
    ah = str(horizon_key) + ".dir"
    qty = str(quantity_key) + ".dat"
    quantity = sxs_horizons[ah][qty]
    spline = spline_horizon_quantity(quantity, start_time, peak_time, truncation_tol=truncation_tol)
    for i, spline_key in enumerate(spline_keys):
        spline_dictionary[spline_key] = spline[i]


def derived_horizon_quantities_from_sxs(sxs_horizons, start_time, peak_time):
    """Compute nhat, omega_orbit, LNhat, and horizon times from Horizons.h5

    Specifically, returns a tuple containing nhat, a unit vector from the
    secondary black hole to the primary black hole; omega_orbit, the orbital
    frequency; LNhat, a unit vector in the direction of the orbital angular
    momentum; and the horizon times of the primary, secondary, and remnant,
    truncated and shifted.

    """
    t_A, x_A = prepare_horizon_quantity(
        sxs_horizons['AhA.dir']['CoordCenterInertial.dat'], start_time,
        peak_time)
    t_B, x_B = prepare_horizon_quantity(
        sxs_horizons['AhB.dir']['CoordCenterInertial.dat'], start_time,
        peak_time)
    # This is used only for t_C
    t_C, x_C = prepare_horizon_quantity(
        sxs_horizons['AhC.dir']['CoordCenterInertial.dat'], start_time,
        peak_time)

    # n_vec is a unit vector pointing from the secondary (black hole B) to
    # the primary (black hole A). We use np.linalg.norm() to get the Euclidean
    # magnitude of n_vec.
    x_A = x_A.T
    x_B = x_B.T
    n_vec = x_A-x_B
    n_vec_norm = np.linalg.norm(n_vec, axis=-1)
    n_hat = n_vec/n_vec_norm[:, None]

    # We compute dn_vec/dt to get a velocity vector, by computing differences
    # for dn_vec and dt.
    dn_vec = np.diff(n_vec, axis=0)
    dt = np.diff(t_A)
    dn_vec_dt = dn_vec/dt[:, None]

    # The orbital frequency is the magnitude of n_vec x dn_vec/dt
    # This is just from the Newtonian expression |r x v| = r^2 \omega.
    # Because taking a time derivative reduces the length of the array by 1,
    # drop the last time of n_vec so n_vec and dn_vec_dt have the same number
    # of points.
    r_cross_v = np.cross(n_vec[:-1], dn_vec_dt)
    r_cross_v_norm = np.linalg.norm(r_cross_v, axis=-1)
    omega_orbit = r_cross_v_norm / n_vec_norm[:-1]**2

    # Finally, LNhat is a unit vector in the direction of the orbital
    # angular momentum. That is, it is a unit vector in the direction of
    # r x p, which is the same direction as r x v.
    LN_hat = r_cross_v / r_cross_v_norm[:, None]

    # Horizons.h5 stores quantities as functions of time. Append time to the
    # derived quantities.
    n_hat_vs_time = np.c_[t_A, n_hat]
    omega_orbit_vs_time = np.c_[t_A[:-1], omega_orbit]
    LN_hat_vs_time = np.c_[t_A[:-1], LN_hat]
    return n_hat_vs_time, omega_orbit_vs_time, LN_hat_vs_time, t_A, t_B, t_C


def insert_derived_spline(spline_dictionary, spline_keys, derived_quantity, log=print, truncation_tol=None):
    """Inserts a spline into a dictionary of horizon splines

    Derived from Horizons.h5.  Note: spline_keys is a vector of key
    names (should be length 1 for scalars, length 3 for vectors), where
    each key is the name of a group that will be written in the LVC
    file, such as Omega-vs-time or LNhatx-vs-time; and derived_quantity
    is the quantity to be splined.  The derived_quantity should be
    computed using derived_horizon_quantities_from_sxs().

    """
    from . import Dataset
    # N.B. times already shifted, truncated to remove junk by
    # derived_horizon_quantities_from_sxs()
    log("Computing " + str(spline_keys))
    tol = 1e-6
    if truncation_tol is True:  # Test for actual identity, not just equality
        truncation_tol = 5e-2 * tol
    elif truncation_tol is False:
        truncation_tol = None
    times_AH = derived_quantity[:, 0]
    quantity_AH = derived_quantity[:, 1:]
    for i, spline_key in enumerate(spline_keys):
        spline_dictionary[spline_key] = Dataset.from_data(times_AH, quantity_AH[:, i], tol=tol, truncation_tol=truncation_tol)


def horizon_splines_from_sxs(horizons, start_time, peak_time, log=print, truncation_tol=None):
    """Prepare dictionary of horizon-quantity splines

    The LVC format expects such a dictionary.  This function creates
    one, starting with an SXS-format Horizons.h5. The start_time and
    peak_time are determined by waveforms.convert_modes.

    """
    horizon_splines = {}

    # Christodoulou mass
    insert_spline(horizons, horizon_splines, ['mass1-vs-time'],
                  'AhA', 'ChristodoulouMass', start_time, peak_time, log, truncation_tol=truncation_tol)
    insert_spline(horizons, horizon_splines, ['mass2-vs-time'],
                  'AhB', 'ChristodoulouMass', start_time, peak_time, log, truncation_tol=truncation_tol)
    insert_spline(horizons, horizon_splines, ['remnant-mass-vs-time'],
                  'AhC', 'ChristodoulouMass', start_time, peak_time, log, truncation_tol=truncation_tol)

    # Dimensionless spin
    insert_spline(horizons, horizon_splines,
                  ['spin1x-vs-time', 'spin1y-vs-time', 'spin1z-vs-time'],
                  'AhA', 'chiInertial', start_time, peak_time, log, truncation_tol=truncation_tol)
    insert_spline(horizons, horizon_splines,
                  ['spin2x-vs-time', 'spin2y-vs-time', 'spin2z-vs-time'],
                  'AhB', 'chiInertial', start_time, peak_time, log, truncation_tol=truncation_tol)
    insert_spline(horizons, horizon_splines,
                  ['remnant-spinx-vs-time', 'remnant-spiny-vs-time', 'remnant-spinz-vs-time'],
                  'AhC', 'chiInertial', start_time, peak_time, log, truncation_tol=truncation_tol)

    # Position
    insert_spline(
        horizons, horizon_splines,
        ['position1x-vs-time', 'position1y-vs-time', 'position1z-vs-time'],
        'AhA', 'CoordCenterInertial', start_time, peak_time, log, truncation_tol=truncation_tol)
    insert_spline(
        horizons, horizon_splines,
        ['position2x-vs-time', 'position2y-vs-time', 'position2z-vs-time'],
        'AhB', 'CoordCenterInertial', start_time, peak_time, log, truncation_tol=truncation_tol)
    insert_spline(
        horizons, horizon_splines,
        ['remnant-positionx-vs-time', 'remnant-positiony-vs-time', 'remnant-positionz-vs-time'],
        'AhC', 'CoordCenterInertial', start_time, peak_time, log, truncation_tol=truncation_tol)

    # Derived quantities: nhat, omega_orbit, LNhat
    n_hat, omega_orbit, LN_hat, t_A, t_B, t_C = derived_horizon_quantities_from_sxs(horizons, start_time, peak_time)
    insert_derived_spline(horizon_splines, ['nhatx-vs-time', 'nhaty-vs-time', 'nhatz-vs-time'], n_hat, log, truncation_tol=truncation_tol)
    insert_derived_spline(horizon_splines, ['Omega-vs-time'], omega_orbit, log, truncation_tol=truncation_tol)
    insert_derived_spline(horizon_splines, ['LNhatx-vs-time', 'LNhaty-vs-time', 'LNhatz-vs-time'], LN_hat, log, truncation_tol=truncation_tol)

    return horizon_splines, t_A, t_B, t_C


def write_horizon_splines_from_sxs(
        out_filename,
        horizon_splines,
        primary_horizon_times,
        secondary_horizon_times,
        remnant_horizon_times,
        log=print):
    """Takes a dictionary of horizon splines, prepared with
    horizon_splines_from_sxs, and writes each spline into an HDF5 file. Also
    outputs the horizon times for the individual and remnant black holes,
    truncated to remove junk radiation and shifted."""
    log("Writing horizon data to '{0}'".format(out_filename))
    with h5py.File(out_filename, 'a') as out_file:
        for key in horizon_splines.keys():
            out_group = out_file.create_group(key)
            horizon_splines[key].write(out_group)
        out_file.create_dataset('HorizonATimes', data=primary_horizon_times)
        out_file.create_dataset('HorizonBTimes', data=secondary_horizon_times)
        out_file.create_dataset('CommonHorizonTimes', data=remnant_horizon_times)
