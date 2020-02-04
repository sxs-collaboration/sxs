# Copyright (c) 2019, Simulating eXtreme Spacetimes Collaboration
# See LICENSE file for details: <https://github.com/sxs-collaboration/sxs/blob/master/LICENSE>

"""Code to convert binary-black-hole simulations from the SXS catalog to
the LVC format described in https://arxiv.org/abs/1703.01076.  Based on
an earlier script by Patricia Schmidt, with input from Geoffrey Lovelace
and Alyssa Garcia.

"""

from . import waveforms, horizons, metadata, comparisons
from .comparisons import compare


class LVCDataset(object):
    """Represents a single dataset of a group in an LVC-format HDF5 file

    Note that the de facto LVC format definition is in the code that
    checks it, which seems to be the lvcnr module.  More
    specifically, the classes in lvcnrpy.format.specs at
    https://git.ligo.org/waveforms/lvcnrpy.  In particular, fields
    for format 3 are required to have subfields 'X', 'Y', 'deg',
    'tol', and 'errors'.

    The first four of these make perfect sense.  The last is not so
    relevant for LVC purposes; it is inherited from romspline, and is
    actually the L1 convergence measure of romspline.  In particular,
    it is nearly -- but distinctly not -- the same size as 'X' and
    'Y'.  This is naturally useful for investigating the algorithm
    itself, but is irrelevant to using its results.  Moreover, since
    this class uses the "peak-greed" variant of the algorithm, that
    dataset no longer means the same thing.  But since it is
    required, we include here an `errors` member, containing a single
    element, being the largest error (scaled, if relevant) in the
    final result.

    """
    def __init__(self, *args, **kwargs):
        if args or kwargs:
            raise ValueError("This is an empty constructor; use `from_data` or `read` to add data.")

    @classmethod
    def from_data(cls, x, y, tol, rel=False, error_scaling=None):
        """Construct reduced-order dataset from (x, y) data"""
        import numpy as np
        from ...utilities.decimation.peak_greed import minimal_grid
        lvc_dataset = LVCDataset()
        lvc_dataset.tol = tol
        if error_scaling is None:
            indices = minimal_grid(x, y, tol=tol)
        else:
            indices = minimal_grid(x, y, tol=tol, error_scale=error_scaling)
        lvc_dataset.deg = 3
        lvc_dataset.X = x[indices].copy()
        lvc_dataset.Y = y[indices].copy()
        if error_scaling is None:
            lvc_dataset.errors = np.array([np.max(np.abs(y - lvc_dataset.spline(x)))])
        else:
            lvc_dataset.errors = np.array([np.max(np.abs(error_scaling * (y - lvc_dataset.spline(x))))])
        # lvc_dataset.compression_ratio = x.size/lvc_dataset.X.size
        return lvc_dataset

    def write(self, output_group):
        import h5py
        if not isinstance(output_group, h5py.Group):
            raise Exception("Parameter `output_group` must be an h5py.Group (or File) object.")
        output_group.create_dataset('deg', data=self.deg, dtype='int')
        output_group.create_dataset('tol', data=self.tol, dtype='double')
        output_group.create_dataset('errors', data=self.errors, dtype='double', compression='gzip', shuffle=True)
        output_group.create_dataset('X', data=self.X, dtype='double', compression='gzip', shuffle=True)
        output_group.create_dataset('Y', data=self.Y, dtype='double', compression='gzip', shuffle=True)

    @classmethod
    def read(cls, input_group):
        import h5py
        if not isinstance(input_group, h5py.Group):
            raise Exception("Parameter `input_group` must be an h5py.Group (or File) object.")
        lvc_dataset = LVCDataset()
        lvc_dataset.deg = input_group['deg'][()]
        lvc_dataset.tol = input_group['tol'][()]
        lvc_dataset.errors = input_group['errors'][:]
        lvc_dataset.X = input_group['X'][:]
        lvc_dataset.Y = input_group['Y'][:]
        return lvc_dataset

    def spline(self, xprime):
        from scipy.interpolate import InterpolatedUnivariateSpline as spline
        return spline(self.X, self.Y, k=self.deg)(xprime)


def bbh_keys_from_simulation_keys(simulation_keys):
    """Extract BBH simulations from a list of all simulations."""
    return [simulation_key for simulation_key in simulation_keys
            if simulation_key.split(':')[-2] == "BBH"]


class SimulationConverter(object):
    class Log(object):
        """Object to replace `log` function that used global `history`

        Instead of using a global `history` variable, just create an
        instance of this class, and pass it around to any function that
        called the old `log` function.  Just like that function, this
        instance can be called with a string and will print the string while
        storing all the strings passed to it.

        Functions expecting an instance of this class can also use `print`
        as a default argument, which will work the same, but not store the
        value.

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

    def __init__(self, sxs_catalog_path='~/.sxs/catalog', modes=8, tolerance=1e-06, quiet=False):
        """Create an object to be used for converting many waveforms to LVC format

        Parameters
        ==========
        sxs_catalog_path: string [defaults to '~/.sxs/catalog']
            Path passed to `sxs.zenodo.catalog.read`.  See that function's docstring for details.
        modes: '22only' or int for max l [defaults to 8]
            Modes to be placed in the output file.  Passing '22only' results
            in the (2,2) and (2,-2) modes being output.  Otherwise, each
            (l,m) mode up to and including the given l value will be output.
            Note that for backwards compatibility, 'all' is also supported,
            and is equivalent to the default value of `8`.
        tolerance: float [defaults to 1e-6]
            Target tolerance used in `sxs.utilities.greedy_spline.minimal_indices`.
        quiet: bool [defaults to False]
            If False, echo each line of the log as it is created; otherwise just store the final
            log in the output file.

        """
        import os
        import time
        import json
        import platform
        import textwrap
        import numpy
        import scipy
        import h5py
        import sxs

        self.sxs_catalog_path=sxs_catalog_path
        self.modes = modes
        self.tolerance = tolerance
        self.quiet = quiet

        self.code_versions = textwrap.dedent("""\
            python=={python}
            numpy=={numpy}
            scipy=={scipy}
            h5py=={h5py}
            # h5py.version.api_version: {h5py_api}
            # h5py.version.hdf5_version: {h5py_hdf5}
            sxs=={sxs}""".format(
                python=platform.python_version(),
                numpy=numpy.version.version,
                scipy=scipy.version.full_version,
                h5py=h5py.version.version,
                h5py_api=h5py.version.api_version,
                h5py_hdf5=h5py.version.hdf5_version,
                sxs=sxs.__version__
            ))

        self.command = textwrap.dedent("""\
            sxs.format.lvc.convert_simulation(
                sxs_data_path={{sxs_data_path!r}},
                out_path={{out_path!r}},
                truncation_time={{truncation_time!r}},
                resolution={{resolution!r}},
                sxs_catalog_path={sxs_catalog_path!r},
                modes={modes!r},
                tolerance={tolerance!r},
                quiet={quiet!r}
            )""".format(
                sxs_catalog_path=sxs_catalog_path,
                modes=modes,
                tolerance=tolerance,
                quiet=quiet
            ))

        # Make sense of the `modes` parameter
        if modes == 'all':
            self.modes = [[l, m] for l in range(2, 9) for m in range(-l, l+1)]
        elif modes == '22only':
            self.modes = [[2, 2], [2, -2]]
        else:
            l_max = int(modes)
            self.modes = [[l, m] for l in range(2, l_max+1) for m in range(-l, l+1)]
        self.ell_max = max(lm[0] for lm in self.modes)

        # Load catalog metadata
        if sxs_catalog_path is None:
            self.sxs_catalog = sxs.zenodo.catalog.read()
        else:
            self.sxs_catalog = sxs.zenodo.catalog.read(sxs_catalog_path)

        self.sxs_catalog_resolutions = sxs.zenodo.catalog.resolutions_for_simulations(self.sxs_catalog)


    def convert(self, sxs_data_path, out_path, truncation_time=None, resolution=None):
        """Convert a simulation from the SXS BBH catalog into the LVC format.

        This function outputs a file in LVC format named SXS_BBH_####_Res#.h5 in out_path.

        Parameters
        ==========
        sxs_data_path: string
            Path to directory containing rhOverM_Asymptotic_GeometricUnits_CoM.h5, Horizons.h5,
            and metadata.json files.
        out_path: string
            Path where LVC format file is to be output
        truncation_time: None or float [defaults to None]
            If specified, truncate time series at this time instead of at the reference time
        resolution: int or None [defaults to None]
            Integer giving the resolution (Lev) of the data to convert.  If this is not given,
            the resolution is determined automatically from sxs_data_path.

        """
        import os
        import time
        import json
        import h5py
        import sxs

        from .metadata import sxs_id_from_alt_names, write_metadata_from_sxs
        from .horizons import horizon_splines_from_sxs, write_horizon_splines_from_sxs
        from .waveforms import convert_modes

        log = self.Log(self.quiet)
        log(self.command.format(sxs_data_path=sxs_data_path, out_path=out_path,
                                truncation_time=truncation_time, resolution=resolution))
        log("Starting at "+time.strftime('%H:%M%p %Z on %b %d, %Y'))

        # Load metadata.json from this simulation
        with open(os.path.join(sxs_data_path, "metadata.json"), 'r') as f:
            metadata = json.load(f)

        # Determine the resolution of the input simulation, if needed
        if resolution is None:
            resolution = sxs.lev_number(sxs_data_path)
        if resolution is None:
            raise ValueError('No `resolution` value found in input arguments or data path.')

        sxs_id = sxs_id_from_alt_names(metadata['alternative_names'])
        log("Converting " + sxs_id)

        extrapolation_order = "Extrapolated_N2"
        log("Extrapolation order: " + extrapolation_order)

        out_name = out_path + "/" + sxs_id.replace(':', '_') + "_Res" + str(resolution) + ".h5"
        log("Output filename is '{0}'".format(out_name))

        start_time, peak_time, version_hist = convert_modes(sxs_data_path + "/rhOverM_Asymptotic_GeometricUnits_CoM.h5",
                                                            metadata, out_name, self.modes, extrapolation_order, log,
                                                            truncation_time, tolerance=self.tolerance/2.0)

        with h5py.File(sxs_data_path + "/Horizons.h5", 'r') as horizons:
            horizon_splines_to_write, t_A, t_B, t_C = horizon_splines_from_sxs(horizons, start_time, peak_time, log)
        write_horizon_splines_from_sxs(out_name, horizon_splines_to_write, t_A, t_B, t_C, log)

        write_metadata_from_sxs(out_name, resolution, metadata,
                                self.sxs_catalog, self.sxs_catalog_resolutions,
                                start_time, peak_time, self.ell_max, log)

        with h5py.File(out_name, 'a') as out_file:
            # Save information about versions of code used in this function
            out_file["auxiliary-info"].create_dataset('CodeVersions.txt', data=self.code_versions)

            # Copy VersionHist.ver into the new file, if available
            if version_hist is not None:
                log("Writing VersionHist.ver")
                out_file["auxiliary-info"].create_dataset('VersionHist.ver', data=version_hist)
            else:
                log("No VersionHist.ver found. Data being converted is version 0.")

            # Store the log output by this script as a dataset
            log("Finishing at "+time.strftime('%H:%M%p %Z on %b %d, %Y'))
            log("Writing log")
            out_file["auxiliary-info"].create_dataset('ConversionLog.txt', data=log.history)


def convert_simulation(sxs_data_path, out_path, truncation_time=None, resolution=None,
                       sxs_catalog_path='~/.sxs/catalog', modes=8, tolerance=1e-06, quiet=False):
    """Convert a simulation from the SXS BBH catalog into the LVC format.
    
    This function outputs a file in LVC format named SXS_BBH_####_Res#.h5 in out_path.

    Parameters
    ==========
    sxs_data_path: string
        Path to directory containing rhOverM_Asymptotic_GeometricUnits_CoM.h5, Horizons.h5,
        and metadata.json files.
    out_path: string
        Path where LVC format file is to be output
    truncation_time: None or float [defaults to None]
        If specified, truncate time series at this time instead of at the reference time
    resolution: int or None [defaults to None]
        Integer giving the resolution (Lev) of the data to convert.  If this is not given,
        the resolution is determined automatically from sxs_data_path.
    sxs_catalog_path: string [defaults to '~/.sxs/catalog']
        Path passed to `sxs.zenodo.catalog.read`.  See that function's docstring for details.
    modes: '22only' or int for max l [defaults to 8]
        Modes to be placed in the output file.  Passing '22only' results
        in the (2,2) and (2,-2) modes being output.  Otherwise, each
        (l,m) mode up to and including the given l value will be output.
        Note that for backwards compatibility, 'all' is also supported,
        and is equivalent to the default value of `8`.
    tolerance: float [defaults to 1e-6]
        Target tolerance used in `sxs.utilities.greedy_spline.minimal_indices`.
    quiet: bool [defaults to False]
        If False, echo each line of the log as it is created; otherwise just store the final
        log in the output file.

    """
    lvc_converter = SimulationConverter(sxs_catalog_path, modes, tolerance, quiet)
    return lvc_converter.convert(sxs_data_path, out_path, truncation_time, resolution)
