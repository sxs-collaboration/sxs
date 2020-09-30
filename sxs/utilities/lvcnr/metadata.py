"""Functions to convert metadata to the LVC format"""

import time
import json
import numpy as np
import h5py

# Compute 1 solar mass * G/c^3
msun_seconds = 4.925491025543576e-06  # LAL value


def sxs_id_from_alt_names(alt_names):
    """Takes an array of alternative names from an SXS metadata.json file
    and returns the SXS ID of the simulation."""
    pattern = 'SXS'
    if not isinstance(alt_names, (list, tuple)):
        alt_names = [alt_names]
    sxs_id = str(next((ss for ss in alt_names if pattern in ss), None))
    return sxs_id


def simulation_type_from_spins(dimensionless_spin_1, dimensionless_spin_2):
    """Classify simulation as non-spinning, aligned-spins, or precessing

    Classifies a simulation with dimensionless_spin_1 on the primary
    (larger) black hole and dimensionless_spin_2 on the secondary black
    hole as nonspinning (if no spin component > 0.01), aligned-spin (if
    only z components are > 0.01), or precessing (otherwise).

    """
    spin_zero_threshold = 0.01  # treat spins smaller than this as zero
    # Types defined in arXiv:1703.01076
    nonspinning_type = "non-spinning"
    aligned_type = "aligned-spins"
    precessing_type = "precessing"

    simulation_type = nonspinning_type

    if (abs(dimensionless_spin_1[2]) > spin_zero_threshold or
            abs(dimensionless_spin_2[2]) > spin_zero_threshold):
        simulation_type = aligned_type

    if (abs(dimensionless_spin_1[0]) > spin_zero_threshold or
            abs(dimensionless_spin_2[0]) > spin_zero_threshold or
            abs(dimensionless_spin_1[1]) > spin_zero_threshold or
            abs(dimensionless_spin_2[1]) > spin_zero_threshold):
        simulation_type = precessing_type

    return simulation_type


def find_comparable_simulations(sxs_id, catalog, catalog_resolutions):
    """Return LVC-format filename containing comparable simulation

    Given an SXD ID sxs_id, a dictionary catalog_data (whose keys are
    SXS ID numbers and whose values are dictionaries including masses
    and spins) and a dictionary catalog_resolutions (whose keys are SXS
    ID numbers and whose values are dictionaries containing lists of
    integers indicating what resolutions are available in the catalog
    for that SXS ID number), return an LVC-format filename that is a
    comparable simulation with multiple resolutions available.

    """
    mass1 = catalog[sxs_id]['initial_mass1']
    mass2 = catalog[sxs_id]['initial_mass2']
    spin1 = catalog[sxs_id]['initial_dimensionless_spin1']
    spin2 = catalog[sxs_id]['initial_dimensionless_spin2']
    mass_ratio = mass1 / mass2
    spin1_magnitude = np.linalg.norm(spin1)
    spin2_magnitude = np.linalg.norm(spin2)

    # Select SXS ID numbers with multiple resolutions available
    has_multiple_resolutions = []
    for key in catalog_resolutions:
        if len(catalog_resolutions[key]) > 1:
            has_multiple_resolutions.append(key)

    # Select SXS ID numbers with multiple resolutions and the same
    # initial data type and spec revision
    same_id_and_spec_revision = []
    for key in has_multiple_resolutions:
        if (catalog[key]['spec_revisions'] == catalog[sxs_id]['spec_revisions']
                and catalog[key]['initial_data_type'] == catalog[sxs_id]['initial_data_type']):
            same_id_and_spec_revision.append(key)

    if len(same_id_and_spec_revision) > 0:
        has_multiple_resolutions = same_id_and_spec_revision

    # Select an SXS ID number with multiple resolutions that is "similar" to
    # sxs_id. First, try to find a case with the same initial data type and
    # SpEC revision. If that is possible, then choose from among the ones with
    # the closest mass ratio and spin magnitudes.  Otherwise, just select a
    # simulation with the closest mass ratios and spin magnitudes.
    #
    # Note: as in the previous script, here I use initial masses and spins, not
    # relaxed masses and spins.
    mass_spin_diff_best = np.inf
    key_best = has_multiple_resolutions[0]
    for key in has_multiple_resolutions:
        current_mass1 = catalog[key]['initial_mass1']
        current_mass2 = catalog[key]['initial_mass2']
        current_spin1 = catalog[key]['initial_dimensionless_spin1']
        current_spin2 = catalog[key]['initial_dimensionless_spin2']
        current_mass_ratio = current_mass1 / current_mass2
        current_spin1_magnitude = np.linalg.norm(current_spin1)
        current_spin2_magnitude = np.linalg.norm(current_spin2)
        mass_spin_diff = (
                np.abs(mass_ratio - current_mass_ratio)
                + np.abs(spin1_magnitude - current_spin1_magnitude)
                + np.abs(spin2_magnitude - current_spin2_magnitude)
        )
        if mass_spin_diff < mass_spin_diff_best:
            mass_spin_diff_best = mass_spin_diff
            key_best = key

    resolution_best = np.max(catalog_resolutions[key_best])
    return key_best.replace(':', '_') + "_Res" + str(resolution_best) + ".h5"


def write_metadata_from_sxs(out_filename, resolution, metadata, catalog,
                            catalog_resolutions, start_time, peak_time, l_max, log=print):
    """Write metadata to an LVC-format file.

    Adds both auxiliary-info/metadata.json and setting attributes
    conforming to the format given by arXiv:1703.01076. Input
    arguments are the output filename out_filename, the resolution of
    the simulation (an int), metadata (read from metadata.json) for
    this simulation/resolution, the start_time, peak time, and l_max
    determined by waveforms.convert_modes.  The argument
    catalog_resolutions is a dictionary (read from a json file) whose
    keys are SXS ID numbers, and whose values are lists of integers
    corresponding to the resolutions available for that SXS ID number
    in the SXS catalog. The argument catalog should point to the SXS
    catalog metadata JSON file, available at
    https://data.black-holes.org/catalog.json or via
    get_sxs_public_metadata.py.

    """
    log("Writing metadata")
    # Get the SXS ID for the current simulation
    names = metadata['alternative_names']
    if not isinstance(names, (list, tuple)):
        names = [names]
    sxs_id = sxs_id_from_alt_names(names)

    with h5py.File(out_filename, 'a') as out_file:

        # Put the metadata for this simulation into the auxiliary-info group
        aux_group = out_file.create_group('auxiliary-info')
        json_this_sim = json.dumps(metadata)
        aux_group.create_dataset('metadata.json', data=json_this_sim)

        # Note: all numerical quantities for LVC attributes
        # are at the reference time unless otherwise indicated
        mass1 = metadata["reference_mass1"]
        mass2 = metadata["reference_mass2"]
        total_mass = mass1 + mass2
        eta = (mass1 * mass2) / np.square(total_mass)

        # Round slightly larger mass ratios to 1
        # CHECK ME: should we do this?
        if 0.25 < eta < 0.2501:
            eta = 0.25

        dimensionless_spin1 = metadata["reference_dimensionless_spin1"]
        dimensionless_spin2 = metadata["reference_dimensionless_spin2"]

        position1 = np.array(metadata["reference_position1"])
        position2 = np.array(metadata["reference_position2"])
        separation = position1 - position2

        n_hat = separation / np.linalg.norm(separation)
        omega_orbit_vec = metadata["reference_orbital_frequency"]
        omega_orbit = np.linalg.norm(omega_orbit_vec)
        omega_grav_22 = 2.0 * omega_orbit

        # Unit vector in direction of orbital angular momentum, which
        # is the same as the direction of the orbital angular velocity
        # at the reference time
        ln_hat = omega_orbit_vec / omega_orbit

        # Eccentricity quantities
        eccentricity = metadata['reference_eccentricity']
        if str(eccentricity)[0] == "<":
            eccentricity = float(str(eccentricity)[1:])
        if str(eccentricity)[0] == ">":
            eccentricity = float(str(eccentricity)[1:])
        elif str(eccentricity) == '[simulation too short]' or str(eccentricity) == '[unknown]':
            # CHECK ME: is this the right thing to do for cases where we can't
            # measure eccentricity?
            log("Warning: eccentricity not measured for this simulation")
            eccentricity = -1.0
        else:
            eccentricity = float(eccentricity)

        mean_anomaly = metadata['reference_mean_anomaly']
        if isinstance(mean_anomaly, str):
            if mean_anomaly == '[unknown]':
                mean_anomaly = -1.0
            else:
                mean_anomaly = float(mean_anomaly)

        # Set metadata attributes of the output file
        out_file.attrs['NR-group'] = "SXS"
        out_file.attrs['name'] = sxs_id
        out_file.attrs['alternative-names'] = ",".join(names)
        out_file.attrs['simulation-type'] \
            = simulation_type_from_spins(dimensionless_spin1, dimensionless_spin2)

        # our metdata uses lowercase for object types (bh, ns), but LVC wants
        # BH or NS, so convert to uppercase with upper()
        out_file.attrs['object1'] = metadata['object1'].upper()
        out_file.attrs['object2'] = metadata['object2'].upper()

        out_file.attrs['PN_approximant'] = "none, NR only"

        out_file.attrs['NR_reference_time'] = float()

        out_file.attrs['mass1'] = mass1
        out_file.attrs['mass2'] = mass2
        out_file.attrs['eta'] = eta

        # Metadata not in arXiv:1703.01076
        out_file.attrs['NR_reference_time'] = metadata['reference_time']
        out_file.attrs['NR_start_time'] = start_time
        out_file.attrs['NR_peak_time'] = peak_time  # not in Patricia's script
        out_file.attrs['NR_frame'] = 'inertial'

        out_file.attrs['f_lower_at_1MSUN'] = omega_grav_22 / (2.0 * np.pi * msun_seconds)
        out_file.attrs['Omega'] = omega_orbit
        out_file.attrs['spin1x'] = dimensionless_spin1[0]
        out_file.attrs['spin1y'] = dimensionless_spin1[1]
        out_file.attrs['spin1z'] = dimensionless_spin1[2]
        out_file.attrs['spin2x'] = dimensionless_spin2[0]
        out_file.attrs['spin2y'] = dimensionless_spin2[1]
        out_file.attrs['spin2z'] = dimensionless_spin2[2]
        out_file.attrs['nhatx'] = n_hat[0]
        out_file.attrs['nhaty'] = n_hat[1]
        out_file.attrs['nhatz'] = n_hat[2]
        out_file.attrs['LNhatx'] = ln_hat[0]
        out_file.attrs['LNhaty'] = ln_hat[1]
        out_file.attrs['LNhatz'] = ln_hat[2]

        out_file.attrs['eccentricity'] = eccentricity
        out_file.attrs['mean_anomaly'] = mean_anomaly

        out_file.attrs['type'] = "NRinjection"
        out_file.attrs['Format'] = 3
        out_file.attrs['NR-code'] = "SpEC"
        out_file.attrs['modification-date'] = time.strftime("%Y-%m-%d")
        out_file.attrs['point-of-contact-email'] = 'questions@black-holes.org'
        # CHECK ME: is this always the right reference to cite?
        out_file.attrs['INSPIRE-bibtex-keys'] = "Boyle:2019kee"

        # CHECK-ME: should this be public?
        out_file.attrs['license'] = "LVC-internal"

        out_file.attrs['NR-techniques'] = (
            "Quasi-Equilibrium-ID, GH, RWZ-h, "
            "Extrapolated-Waveform, ApproxKillingVectorSpin, Christodoulou-Mass"
        )
        out_file.attrs['Lmax'] = l_max

        # If highest resoultion available, 'production-run' = 1, otherwise
        # 'production-run' = 0. If multiple resolutions available, list their
        # file names for 'files-in-error-series'. Otherwise, list filenames of
        # a comparable simulation with multiple resolutions.
        resolutions = catalog_resolutions.get(sxs_id, [])
        if len(resolutions) > 1:
            error_name_base = out_filename.split('/')[-1].split("_Res")[0] + "_Res"
            error_series = [error_name_base + str(i) + ".h5" for i in resolutions]
            out_file.attrs['files-in-error-series'] = ",".join(error_series)
            out_file.attrs['comparable-simulation'] = ""
            if resolution == np.max(resolutions):
                out_file.attrs['production-run'] = 1
            else:
                out_file.attrs['production-run'] = 0
        else:
            out_file.attrs['production-run'] = 1
            out_file.attrs['files-in-error-series'] = ""
            out_file.attrs['comparable-simulation'] = ""
