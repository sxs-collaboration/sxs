"""Class and function to convert SXS data to LVC-NR format"""


class SimulationConverter(object):
    class Log(object):
        """Object to replace `log` function that used global `history`

        Instead of using a global `history` variable, just create an instance of this
        class, and pass it around to any function that called the old `log` function.
        Just like that function, this instance can be called with a string and will
        print the string while storing all the strings passed to it.

        Functions expecting an instance of this class can also use `print` as a default
        argument, which will work the same, but not store the value.

        """

        def __init__(self, quiet):
            self.history = ""
            self.quiet = quiet

        def __call__(self, string):
            if not self.quiet:
                print(string)
            self.history += string + "\n"

        def __str__(self):
            return str(self.history)

        def __repr__(self):
            return repr(self.history)

    def __init__(self, modes=8, tolerance=1e-06, quiet=False):
        """Create an object to be used for converting many waveforms to LVC format

        Parameters
        ----------
        modes : {int, '22only'}, optional
            Modes to be placed in the output file.  Passing '22only' results in the
            (2,2) and (2,-2) modes being output.  Otherwise, each (l,m) mode up to and
            including the given integer value will be output.  Note that for backwards
            compatibility, 'all' is also supported, and is equivalent to the default
            value of `8`.
        tolerance : float, optional
            Target tolerance used in `sxs.utilities.greedy_spline.minimal_indices`.
        quiet : bool, optional
            If False (the default), echo each line of the log as it is created;
            otherwise just store the final log in the output file.

        """
        import platform
        import numpy
        import scipy
        import h5py
        from ... import __version__, load, zenodo

        self.modes = modes
        self.tolerance = tolerance
        self.quiet = quiet

        self.code_versions = (
            f"python=={platform.python_version()}\n"
            f"numpy=={numpy.version.version}\n"
            f"scipy=={scipy.version.full_version}\n"
            f"h5py=={h5py.version.version}\n"
            f"# h5py_api=={h5py.version.api_version}\n"
            f"# h5py_hdf5=={h5py.version.hdf5_version}\n"
            f"sxs=={__version__}\n"
        )

        self.command = (
            f"sxs.utilities.lvcnr.convert_simulation(\n"
            f"    sxs_data_path={{sxs_data_path!r}},\n"
            f"    out_path={{out_path!r}},\n"
            f"    truncation_time={{truncation_time!r}},\n"
            f"    resolution={{resolution!r}},\n"
            f"    modes={modes!r},\n"
            f"    tolerance={tolerance!r},\n"
            f"    quiet={quiet!r}\n"
            f")"
        )

        # Make sense of the `modes` parameter
        if modes == "all":
            self.modes = [[l, m] for l in range(2, 9) for m in range(-l, l + 1)]
        elif modes == "22only":
            self.modes = [[2, 2], [2, -2]]
        else:
            l_max = int(modes)
            self.modes = [[l, m] for l in range(2, l_max + 1) for m in range(-l, l + 1)]
        self.ell_max = max(lm[0] for lm in self.modes)

        # Load catalog metadata
        catalog = load("catalog")
        self.sxs_catalog = {
            "simulations": catalog.simulations,
            "records": catalog.records,
        }

        self.sxs_catalog_resolutions = zenodo.catalog.resolutions_for_simulations(self.sxs_catalog)

    def convert(
        self,
        sxs_data_path,
        out_path,
        waveform_name=None,
        out_waveform_name=None,
        truncation_time=None,
        resolution=None,
        truncation_tol=None,
    ):
        """Convert a simulation from the SXS BBH catalog into the LVC format.

        This function outputs a file in LVC format named SXS_BBH_####_Res#.h5 in
        out_path.

        Parameters
        ----------
        sxs_data_path : string
            Path to directory containing waveform file,
            Horizons.h5, and metadata.txt or metadata.json files.
        out_path : string
            Path where LVC-format file is to be output
        waveform_name : string
            Name of waveform to load in sxs_data_path. If not specified,
            will try to load Strain_N2.h5 or rhOverM_Asymptotic_GeometricUnits_CoM.h5
        out_waveform_name : string
            Name of LVC-format file to output in out_path. If not specified,
            will use SXS_BBH_####_Res#.h5
        truncation_time : {None, float}
            If specified, truncate time series at this time instead of at the reference
            time
        resolution : {None, int}
            Integer giving the resolution (Lev) of the data to convert.  If this is not
            given, the resolution is determined automatically from sxs_data_path.
        truncation_tol : {None, bool, callable, float, array_like}, optional
            If None (the default) or False, nothing happens.  If True, the waveform
            data (amplitude and phase) are "truncated" so that bits with significance
            lower than `5e-2 * self.tolerance` are set to zero, for improved
            compression.  Any other input is passed to `sxs.TimeSeries.truncate`.  Note
            that this is not typically a very effective setting — perhaps providing
            another 10% compression; the output file sizes are dominated by fairly
            redundant time data unaffected by this parameter.

        """
        import os
        import time
        import h5py
        from ... import lev_number, Metadata

        from .metadata import sxs_id_from_alt_names, write_metadata_from_sxs
        from .horizons import horizon_splines_from_sxs, write_horizon_splines_from_sxs
        from .waveforms import convert_modes

        log = self.Log(self.quiet)
        log(
            self.command.format(
                sxs_data_path=sxs_data_path, out_path=out_path, truncation_time=truncation_time, resolution=resolution
            )
        )
        log("Starting at " + time.strftime("%H:%M%p %Z on %b %d, %Y"))

        # Load metadata from this simulation
        try:
            metadata = Metadata.from_file(os.path.join(sxs_data_path, "metadata"))
        except:
            raise ValueError("Cannot load metadata.")

        # Determine the resolution of the input simulation, if needed
        if resolution is None:
            resolution = lev_number(sxs_data_path)
        if resolution is None:
            raise ValueError("No `resolution` value found in input arguments or data path.")

        sxs_id = sxs_id_from_alt_names(metadata["alternative_names"])
        log("Converting " + sxs_id)

        extrapolation_order = "Extrapolated_N2"
        log("Extrapolation order: " + extrapolation_order)

        if not out_waveform_name is None:
            out_name = out_path + "/" + out_waveform_name
        else:
            out_name = out_path + "/" + sxs_id.replace(":", "_") + "_Res" + str(resolution) + ".h5"
        log("Output filename is '{0}'".format(out_name))

        if waveform_name is None:
            if os.path.exists(os.path.join(sxs_data_path, "Strain_N2.h5")):
                waveform_path = sxs_data_path + "/Strain_N2"
            elif os.path.exists(os.path.join(sxs_data_path, "rhOverM_Asymptotic_GeometricUnits_CoM.h5")):
                waveform_path = sxs_data_path + "/rhOverM_Asymptotic_GeometricUnits_CoM.h5"
            else:
                raise ValueError("Cannot find a default waveform file.")
        else:
            waveform_path = sxs_data_path + "/" + waveform_name

        start_time, peak_time, version_hist = convert_modes(
            waveform_path,
            metadata,
            out_name,
            self.modes,
            extrapolation_order,
            log,
            truncation_time,
            tolerance=self.tolerance / 2.0,
            truncation_tol=truncation_tol,
        )

        with h5py.File(sxs_data_path + "/Horizons.h5", "r") as horizons:
            horizon_splines_to_write, t_A, t_B, t_C = horizon_splines_from_sxs(
                horizons, start_time, peak_time, log, truncation_tol=truncation_tol
            )
        write_horizon_splines_from_sxs(out_name, horizon_splines_to_write, t_A, t_B, t_C, log)

        write_metadata_from_sxs(
            out_name,
            resolution,
            metadata,
            self.sxs_catalog,
            self.sxs_catalog_resolutions,
            start_time,
            peak_time,
            self.ell_max,
            log,
        )

        with h5py.File(out_name, "a") as out_file:
            # Save information about versions of code used in this function
            out_file["auxiliary-info"].create_dataset("CodeVersions.txt", data=self.code_versions)

            # Copy VersionHist.ver into the new file, if available
            if version_hist is not None:
                log("Writing VersionHist.ver")
                out_file["auxiliary-info"].create_dataset("VersionHist.ver", data=version_hist)
            else:
                log("No VersionHist.ver found. Data being converted is version 0.")

            # Store the log output by this script as a dataset
            log("Finishing at " + time.strftime("%H:%M%p %Z on %b %d, %Y"))
            log("Writing log")
            out_file["auxiliary-info"].create_dataset("ConversionLog.txt", data=log.history)


def convert_simulation(
    sxs_data_path,
    out_path,
    waveform_name=None,
    out_waveform_name=None,
    truncation_time=None,
    resolution=None,
    modes=8,
    tolerance=1e-06,
    quiet=False,
):
    """Convert a simulation from the SXS BBH catalog into the LVC format.

    This function outputs a file in LVC format named SXS_BBH_####_Res#.h5 in
    out_path.

    Note that this function is essentially a wrapper for
    `SimulationConverter.convert`.  If you have very many systems to convert, it is
    significantly faster to create the SimulationConverter object once, and then
    call the `convert` method for each system.

    Parameters
    ----------
    sxs_data_path : string
        Path to directory containing waveform file,
        Horizons.h5, and metadata.txt or metadata.json files.
    out_path : string
        Path where LVC-format file is to be output
    waveform_name : string
        Name of waveform to load in sxs_data_path. If not specified,
        will try to load Strain_N2.h5 or rhOverM_Asymptotic_GeometricUnits_CoM.h5
    out_waveform_name : string
        Name of LVC-format file to output in out_path. If not specified,
        will use SXS_BBH_####_Res#.h5
    truncation_time : {None, float}, optional
        If specified, truncate time series at this time instead of at the reference time
    resolution : {None, int}, optional
        Integer giving the resolution (Lev) of the data to convert.  If this is not given,
        the resolution is determined automatically from sxs_data_path.
    modes : {int, '22only'}, optional
        Modes to be placed in the output file.  Passing '22only' results in the (2,2)
        and (2,-2) modes being output.  Otherwise, each (l,m) mode up to and including
        the given l value will be output.  Note that for backwards compatibility, 'all'
        is also supported, and is equivalent to the default value of `8`.
    tolerance : float, optional
        Target tolerance used in `sxs.utilities.greedy_spline.minimal_indices`.
    quiet : bool, optional
        If False (the default), echo each line of the log as it is created; otherwise
        just store the final log in the output file.

    """
    lvc_converter = SimulationConverter(modes, tolerance, quiet)
    return lvc_converter.convert(sxs_data_path, out_path, waveform_name, out_waveform_name, truncation_time, resolution)
