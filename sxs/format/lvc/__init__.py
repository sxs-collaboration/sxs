# Copyright (c) 2019, Simulating eXtreme Spacetimes Collaboration
# See LICENSE file for details: <https://github.com/sxs-collaboration/sxs/blob/master/LICENSE>

"""Code to convert binary-black-hole simulations from the SXS catalog to
the LVC format described in https://arxiv.org/abs/1703.01076.  Based on
an earlier script by Patricia Schmidt, with input from Geoffrey Lovelace
and Alyssa Garcia.

"""


def bbh_keys_from_simulation_keys(simulation_keys):
    """Extract BBH simulations from a list of all simulations."""
    return [simulation_key for simulation_key in simulation_keys
            if simulation_key.split(':')[-2] == "BBH"]


def convert_simulation(sxs_data_path, resolution, sxs_catalog_metadata_path,
                       sxs_catalog_resolutions_path, out_path, modes=8,
                       truncation_time=None):
    """Convert a simulation from the SXS BBH catalog into the LVC format.
    
    sxs_data_path is a path to a directory that must contain the following:
        i)   rhOverM_Asymptotic_GeometricUnits_CoM.h5
        ii)  Horizons.h5
        iii) metadata.json
    Additionally, the function requires paths to the following:
        iv)  sxs_catalog_metadata_path points to sxs_catalog.json, which
             contains information such as masses and spins. A file in this
             format is available at https://data.black-holes.org/catalog.json.
        v)   sxs_catalog_resolutions_path points to a file whose keys are
             SXS ID numbers and whose values are lists of integers, where each
             integer corresponds to a resolution available in the catalog.
    The option resolution is an integer labeling the resolution of the
    converted waveform. Modes is an array of the format
    [[l1, m1], [l2, m2], ...] listing the l,m modes to convert. This function
    outputs a file in LVC format named SXS_BBH_####_Res#.h5
    in out_path.

    Parameters
    ==========
    sxs_data_path: string
        Path containing rh*h5, Horizons.h5, metadata.json
    resolution: int
        Integer giving the resolution of the data to convert
    sxs_catalog_metadata_path: string
        Path to sxs_catalog.json (via get_sxs_public_metadata.py)
    sxs_catalog_resolutions_path: string
        Path to sxs_catalog_resolutions.json
    out_path: string
        Path where LVC format file is to be output
    modes: '22only' or int for max l [defaults to 8]
        Modes to be placed in the output file.  Passing '22only' results
        in the (2,2) and (2,-2) modes being output.  Otherwise, each
        (l,m) mode up to and including the given l value will be output.
        Note that for backwards compatibility, 'all' is also supported,
        and is equivalent to the default value of `8`.
    truncation_time: None or float [defaults to None]
        If specified, truncate time series at this time instead of at the reference time

    """
    import h5py
    import json
    import os

    from .metadata import sxs_id_from_alt_names, write_metadata_from_sxs
    from .horizons import horizon_splines_from_sxs, write_horizon_splines_from_sxs
    from .waveforms import spline_amp_phase_from_sxs, write_splines_to_H5

    class Log(object):
        """Object to replace `log` function that used global `history`

        Instead of using a global `history` variable, just create an
        instance of this class, and pass it around to any function that
        called the old `log` function.  Just like that function, this
        instance can be called with a string and will print the string while
        storing all the strings passed to it.

        Functions expecting taking an instance of this class can also use
        `print` as a default argument, which will work the same, but not
        store the value.

        """
        def __init__(self):
            self.history = ""

        def __call__(self, string):
            print(string)
            self.history += string + "\n"

        def __str__(self):
            return str(self.history)

        def __repr__(self):
            return repr(self.history)

    if modes == 'all':
        modes = [[l, m] for l in range(2, 9) for m in range(-l, l+1)]
    elif modes == '22only':
        modes = [[2, 2], [2, -2]]
    else:
        l_max = int(modes)
        modes = [[l, m] for l in range(2, l_max+1) for m in range(-l, l+1)]

    log = Log()

    with open(sxs_data_path + "/metadata.json", 'r') as f:
        metadata = json.load(f)
    with open(sxs_catalog_metadata_path, 'r') as f:
        sxs_catalog_metadata = json.load(f)
    with open(sxs_catalog_resolutions_path, 'r') as f:
        sxs_catalog_resolutions = json.load(f)

    sxs_id = sxs_id_from_alt_names(metadata['alternative_names'])
    log("Converting " + sxs_id)
    log("Running script " + os.path.realpath(__file__) + " in " + os.getcwd())
    log("convert_sxs_to_lvc.py called with the following parameters:")
    log("  sxs_data_path: " + sxs_data_path)
    log("  resolution: " + str(resolution))
    log("  sxs_catalog_metadata_path: " + sxs_catalog_metadata_path)
    log("  sxs_catalog_resolutions_path: " + sxs_catalog_resolutions_path)
    log("  modes: " + str(modes))
    log("  out_path: " + str(out_path))

    extrapolation_order = "Extrapolated_N2"
    log("Extrapolation order: " + extrapolation_order)

    out_name = out_path + "/" + sxs_id.replace(':', '_') + "_Res" \
                        + str(resolution) + ".h5"
    log("Output filename is " + out_name)

    with h5py.File(sxs_data_path + "/rhOverM_Asymptotic_GeometricUnits_CoM.h5", 'r') as rhOverM:
        version_hist = rhOverM.get('VersionHist.ver', None)
        modes, times, spline_amps, spline_phases, start_time, peak_time, l_max \
            = spline_amp_phase_from_sxs(rhOverM, metadata, modes, 
                                        extrapolation_order, log, truncation_time)
    write_splines_to_H5(out_name, modes, spline_amps, spline_phases, times, log)

    with h5py.File(sxs_data_path + "/Horizons.h5", 'r') as horizons:
        horizon_splines_to_write, t_A, t_B, t_C = horizon_splines_from_sxs(horizons, start_time, peak_time, log)
    write_horizon_splines_from_sxs(out_name, horizon_splines_to_write, t_A, t_B, t_C, log)

    write_metadata_from_sxs(out_name, resolution, metadata,
                            sxs_catalog_metadata, sxs_catalog_resolutions,
                            start_time, peak_time, l_max, log)

    with h5py.File(out_name, 'a') as out_file:
        # Save a copy of this script in auxiliary-info
        with open(os.path.realpath(__file__), 'r') as f:
            out_file["auxiliary-info"].create_dataset('convert_sxs_to_lvc.py', data=f.read())

        if version_hist is not None:
            log("Writing VersionHist.ver")
            out_file["auxiliary-info"].create_dataset('VersionHist.ver', data=version_hist)
        else:
            log("No VersionHist.ver found. Data being converted is version 0.")
        log("Writing log")
        out_file["auxiliary-info"].create_dataset('ConversionLog.txt', data=log.history)
