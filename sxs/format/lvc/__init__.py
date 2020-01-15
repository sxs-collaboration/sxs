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


def convert_simulation(sxs_data_path, out_path,
                       sxs_catalog_path='~/.sxs/catalog', resolution=None, modes=8,
                       truncation_time=None):
    """Convert a simulation from the SXS BBH catalog into the LVC format.
    
    This function outputs a file in LVC format named SXS_BBH_####_Res#.h5 in out_path.

    Parameters
    ==========
    sxs_data_path: string
        Path to directory containing rhOverM_Asymptotic_GeometricUnits_CoM.h5, Horizons.h5,
        and metadata.json files.
    out_path: string
        Path where LVC format file is to be output
    sxs_catalog_path: string [defaults to '~/.sxs/catalog']
        Path passed to `sxs.zenodo.catalog.read`.  See that function's docstring for details.
    resolution: int or None [defaults to None]
        Integer giving the resolution (Lev) of the data to convert.  If this is not given,
        the resolution is determined automatically from sxs_data_path.
    modes: '22only' or int for max l [defaults to 8]
        Modes to be placed in the output file.  Passing '22only' results
        in the (2,2) and (2,-2) modes being output.  Otherwise, each
        (l,m) mode up to and including the given l value will be output.
        Note that for backwards compatibility, 'all' is also supported,
        and is equivalent to the default value of `8`.
    truncation_time: None or float [defaults to None]
        If specified, truncate time series at this time instead of at the reference time

    """
    import os
    import time
    import json
    import platform
    import textwrap
    import numpy
    import scipy
    import h5py
    import romspline
    import sxs

    from ... import lev_number, zenodo
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

    code_versions = textwrap.dedent("""\
        python=={python}
        numpy=={numpy}
        scipy=={scipy}
        h5py=={h5py}
        # h5py.version.api_version: {h5py_api}
        # h5py.version.hdf5_version: {h5py_hdf5}
        romspline=={romspline}
        sxs=={sxs}""".format(
            python=platform.python_version(),
            numpy=numpy.version.version,
            scipy=scipy.version.full_version,
            h5py=h5py.version.version,
            h5py_api=h5py.version.api_version,
            h5py_hdf5=h5py.version.hdf5_version,
            romspline=romspline.__version__,
            sxs=sxs.__version__
        ))

    log = Log()

    # Load metadata.json from this simulation
    with open(os.path.join(sxs_data_path, "metadata.json"), 'r') as f:
        metadata = json.load(f)

    # Load catalog metadata
    if sxs_catalog_path is None:
        sxs_catalog = zenodo.catalog.read()
    else:
        sxs_catalog = zenodo.catalog.read(sxs_catalog_path)

    # Determine the resolution of the input simulation, if needed
    if resolution is None:
        resolution = lev_number(sxs_data_path)
    if resolution is None:
        raise ValueError('No `resolution` value found in input arguments or data path.')

    sxs_catalog_resolutions = zenodo.catalog.resolutions_for_simulations(sxs_catalog)

    sxs_id = sxs_id_from_alt_names(metadata['alternative_names'])
    log("Converting " + sxs_id)
    log("sxs.format.lvc.convert_simulation called with the following parameters:")
    log("  sxs_data_path: " + sxs_data_path)
    log("  out_path: " + str(out_path))
    log("  sxs_catalog_path: " + sxs_catalog_path)
    log("  resolution: " + str(resolution))
    log("  modes: " + str(modes))
    log("  truncation_time: " + str(truncation_time))

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
                            sxs_catalog, sxs_catalog_resolutions,
                            start_time, peak_time, l_max, log)

    with h5py.File(out_name, 'a') as out_file:
        # Save information about versions of code used in this function
        out_file["auxiliary-info"].create_dataset('CodeVersions.txt', data=code_versions)

        # Copy VersionHist.ver into the new file, if available
        if version_hist is not None:
            log("Writing VersionHist.ver")
            out_file["auxiliary-info"].create_dataset('VersionHist.ver', data=version_hist)
        else:
            log("No VersionHist.ver found. Data being converted is version 0.")

        # Store the log output by this script as a dataset
        log("Finishing at "+time.strftime('%l:%M%p %Z on %b %d, %Y'))
        log("Writing log")
        out_file["auxiliary-info"].create_dataset('ConversionLog.txt', data=log.history)
