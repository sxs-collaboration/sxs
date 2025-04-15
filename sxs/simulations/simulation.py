from pathlib import Path
from warnings import warn
from .. import doi_url, Metadata
from ..utilities import (
    sxs_id_and_version, lev_number, sxs_path_to_system_path,
    download_file, sxs_directory, read_config
)


def search_prefixes(file, lev, files, ending=""):
    """Find the actual filename present in the list of files

    Different versions of Zenodo and CaltechDATA place different
    restrictions on filenames — and specifically, things I would
    consider to be directories.  If there are multiple Levs for a
    simulation, we have to somehow specify the Lev in the filename.
    The permitted strings for doing so have changed over the years, so
    the directory separator was allowed to be a "/" for old systems,
    but must be something else now; we settled on ":".  However, if
    there is just one Lev in a simulation, and we try to save the
    files in a consistent way — including the Lev with the appropriate
    separator, different versions of Zenodo and CaltechDATA will
    either allow or silently remove that prefix.  The simplest
    approach is to just search over all possibilities, for which
    filename is actually present in the data.  That's what this
    function does.
    
    """
    for prefix in [f"{lev}:", f"{lev}/", ""]:
        fn = f"{prefix}{file}"
        if f"{fn}{ending}" in files:
            return fn
    raise ValueError(f"{file}{ending} not found in any form in files")


def Simulation(location, *args, **kwargs):
    """Construct a Simulation object from a location string

    The location string should be an SXS ID with an optional version
    number, and possibly a Lev specification, as in
    "SXS:BBH:1234v2.0/Lev5".  The default version number is "v2.0",
    and the default Lev number is the highest available.  Note that
    the version number must be in the form "v2.0", not "2.0", and it
    must exist, or no files will be found to load the data from.

    The returned object will be either a `Simulation_v1` or
    `Simulation_v2` object, depending on the version number.
    Hopefully, most details of the two versions will be hidden from
    the user, so that the interface is identical.

    Note that some simulations are deprecated and/or superseded by
    other simulations.  By default, this function will raise an error
    if you try to load a deprecated or superseded simulation.  There
    are several ways to override this behavior:

    1. Pass `ignore_deprecation=True` to completely bypass even
        checking for deprecation or supersession.  No warnings or
        errors will be issued.
    2. Include an explicit version number in the `location`
        string, as in "SXS:BBH:0001v2.0".  A warning will be issued
        that the simulation is deprecated, but it will be loaded
        anyway.
    3. Pass `auto_supersede=True` to automatically load the
        superseding simulation, if there is only one.  Because no
        superseding simulation can be *precisely* the same as the
        deprecated one, there may be multiple superseding simulations
        that have very similar parameters, in which case an error will
        be raised and you must explicitly choose one.  If there is
        only one, a warning will be issued, but the superseding
        simulation will be loaded.

    Otherwise, a `ValueError` will be raised, with an explanation and
    suggestions on what you might want to do.

    Parameters
    ----------
    location : str
        The location string for the simulation.  This must include the
        SXS ID, but can also include the version and Lev number, as
        described above.

    Keyword Arguments
    -----------------
    ignore_deprecation : bool, optional
        If `True`, completely bypass checking for deprecation or
        supersession.  No warnings or errors will be issued.  Default
        is `False`.
    auto_supersede : bool or float, optional
        If present, automatically load the closest undeprecated
        simulation.  If this is a float and the distance to the
        closest undeprecated simulation (see the following argument)
        is larger, a ValueError will be raised.  Otherwise, a warning
        will be issued with the ID of the new simulation and the
        distance.  Note that this can also be set in the configuration
        file by calling, e.g.,
        `sxs.write_config(auto_supersede=True)`.
    metadata_metric : MetadataMetric, optional
        Metric to use for comparing simulations when automatically
        superseding deprecated a simulation.  If not provided, the
        default metric will be used.
    extrapolation : str, optional
        The extrapolation order to use for the strain and Psi4 data.
        This is only relevant for versions 1 and 2 of the data format,
        both of which default to "N2".  Other options include "N3",
        "N4", and "Outer".  "Nx" refers to extrapolation by
        polynomials in 1/r with degree `x`, while "Outer" refers to
        data extracted at the outermost extraction radius but
        corrected for time-dilation and areal-radius effects.
    download_file_info : bool, optional
        If `True`, download the information about the files from the
        CaltechDATA record.  If `False`, only use the file information
        that is already available (which will raise an error if the
        file information has not previously been downloaded).  If not
        present, use the value from the configuration file, defaulting
        to `True` if it is not configured.

    Returns
    -------
    simulation : SimulationBase
        A `Simulation_v1` or `Simulation_v2` object, depending on the
        version of the simulation data.

    Note that all remaining arguments (including keyword arguments)
    are passed on to the `SimulationBase`, `Simulation_v1`, and/or
    `Simulation_v2` constructors, though none currently recognize any
    arguments other than those listed above.
    
    """
    import numpy as np
    from packaging.version import Version
    from .. import load, sxs_directory
    from ..metadata.metric import MetadataMetric

    # Load the simulation catalog
    simulations = load("simulations")
    v = Version(simulations.tag)
    latest_version = f"v{v.major}.{v.minor}"

    # Extract the simulation ID, version, and Lev from the location string
    simulation_id, input_version = sxs_id_and_version(location)
    if not simulation_id:
        if (sim_id := location.split("/Lev")[0]) in simulations:
            simulation_id = sim_id
            input_version = latest_version
        else:
            raise ValueError(f"Invalid SXS ID in '{simulation_id}'")
    input_lev_number = lev_number(location)  # Will be `None` if not present

    # Check if simulation ID exists in the catalog
    if simulation_id not in simulations:
        raise ValueError(f"Simulation '{simulation_id}' not found in simulation catalog")

    # Attach metadata to this object
    metadata = Metadata(simulations[simulation_id])
    series = simulations.dataframe.loc[simulation_id]

    # If input_version is not the default, remove "files" from metadata
    version_is_not_default = (
        input_version
        and input_version != max(metadata.get("DOI_versions", []), default="")
    )
    if version_is_not_default:
        metadata = type(metadata)({
            key: value for key, value in metadata.items() if key != "files"
        })

    # Check if the specified version exists in the simulation catalog
    if not hasattr(metadata, "DOI_versions"):
        input_version = latest_version
    if input_version != latest_version and input_version not in metadata.DOI_versions:
        raise ValueError(f"Version '{input_version}' not found in simulation catalog for '{simulation_id}'")

    # Set various pieces of information about the simulation
    version = input_version or max(metadata.DOI_versions)
    if not version.startswith("v"):
        raise ValueError(f"Invalid version string '{version}'")
    sxs_id_stem = simulation_id
    sxs_id = f"{sxs_id_stem}{version}"
    url = f"{doi_url}{sxs_id}"

    # Deal with deprecations
    deprecated = "deprecated" in metadata.get("keywords", [])
    if deprecated and not kwargs.get("ignore_deprecation", False):
        auto_supersede = kwargs.get("auto_supersede", read_config("auto_supersede", False))
        if not bool(auto_supersede):
            if not input_version:
                raise ValueError(
                    f"Simulation '{location}' is deprecated.  You could\n"
                    +  "  1. pass `ignore_deprecation=True` to load the latest available version,\n"
                    +  "  2. manually choose a different simulation from the catalog,\n"
                    +  "  3. pass `auto_supersede=0.01` to load a match closer than 0.01 in the catalog if one exists,\n"
                    +  "  4. pass `auto_supersede=True` to load the closest match in the catalog, or\n"
                    + f"  5. include the version number, as in '{sxs_id_stem}v2.0', to load a specific version.\n"
                )
            else:
                message = ("\n"
                    + f"Simulation '{sxs_id_stem}' is deprecated, but you explicitly\n"
                    + f"requested version '{input_version}', so it is being used.\n"
                    + f"Pass `ignore_deprecation=True` to quiet this warning.\n"
                )
                warn(message)
        else:
            if input_version:
                message = ("\n"
                    + f"\nSimulation '{sxs_id}' is deprecated.  You explicitly requested.\n"
                    + f"version '{input_version}', but you also passed the `auto_supersede` option.\n"
                    + f"Using the specified version, as that takes precedence.\n"
                )
                warn(message)
            else:
                original_kwargs = kwargs.copy()
                original_kwargs["ignore_deprecation"] = True
                original = Simulation(location, *args, **original_kwargs)
                metadata_metric = kwargs.pop("metadata_metric", MetadataMetric())
                superseding, distance = original.closest_simulation(
                    dataframe=simulations.dataframe,
                    metadata_metric=metadata_metric
                )
                if isinstance(auto_supersede, str):
                    try:
                        auto_supersede = float(auto_supersede)
                    except:
                        pass
                if isinstance(auto_supersede, float) and distance > auto_supersede:
                    raise ValueError(
                        f"Simulation '{sxs_id}' is deprecated, but the closest undeprecated\n"
                        f"simulation is more distant than allowed by `auto_supersede` argument:\n"
                        f"{distance:.3g} > {auto_supersede:.3g}.\n"
                    )
                message = (
                    f"\nSimulation '{sxs_id}' is being automatically superseded by '{superseding}'."
                    f"\nThe distance between them in the given metadata metric is {distance:.3g}."
                )
                warn(message)
                new_location = f"{superseding}{input_version}"
                if input_lev_number:
                    new_location += f"/Lev{input_lev_number}"
                return Simulation(new_location, *args, **kwargs)

    # Note the deprecation status in the kwargs, even if ignoring deprecation
    kwargs["deprecated"] = deprecated

    # We want to do this *after* deprecation checking, to avoid possibly unnecessary web requests
    if version_is_not_default:
        # The default metadata points to files with a different version, so delete this info
        if "files" in metadata:
            del metadata["files"]
    files = get_file_info(metadata, sxs_id, download=kwargs.get("download_file_info", None))

    # If Lev is given as part of `location`, use it; otherwise, use the highest available
    lev_numbers = metadata.get(
        "lev_numbers",
        sorted({lev_num for f in metadata.get("files", []) if (lev_num:=lev_number(f))})
    )
    if not lev_numbers:
        raise ValueError(f"Could not find Levs for {location}")
    if input_lev_number is not None and input_lev_number not in lev_numbers:
        raise ValueError(
            f"Lev number '{input_lev_number}' not found in simulation files for {sxs_id}"
        )
    max_lev_number = max(lev_numbers)
    output_lev_number = input_lev_number or max_lev_number
    if output_lev_number is None:
        raise ValueError(
            f"No Lev number found for {location}"
        )
    location = f"{sxs_id_stem}{version}/Lev{output_lev_number}"

    # Keep the metadata around unless we're asking for an old version
    # or a less-than-maximal Lev
    if (
        version_is_not_default
        or (lev_numbers and output_lev_number != max_lev_number)
    ):
        metadata = None

    # Finally, figure out which version of the simulation to load and dispatch
    version_number = float(version[1:])
    if 1 <= version_number < 2.0:
        sim = Simulation_v1(
            metadata, series, version, sxs_id_stem, sxs_id, url, files, lev_numbers, output_lev_number, location, *args, **kwargs
        )
    elif 2 <= version_number < 3.0:
        sim = Simulation_v2(
            metadata, series, version, sxs_id_stem, sxs_id, url, files, lev_numbers, output_lev_number, location, *args, **kwargs
        )
    elif 3 <= version_number < 4.0 or version == latest_version:
        sim = Simulation_v3(
            metadata, series, version, sxs_id_stem, sxs_id, url, files, lev_numbers, output_lev_number, location, *args, **kwargs
        )
    else:
        raise ValueError(f"Version '{version}' not yet supported")
    sim.__file__ = str(sxs_directory("cache") / sxs_path_to_system_path(sim.sxs_id))
    return sim


class SimulationBase:
    """Base class for Simulation objects

    Note that users almost certainly never need to call this function;
    see the `Simulation` function or `sxs.load` function instead.
    
    Attributes
    ----------
    metadata : Metadata or None
        Metadata object for the simulation.  If `None`, the metadata
        will be loaded automatically.
    series : pandas.Series
        The metadata, as extracted from the `simulations.dataframe`,
        meaning that it has columns consistent with other simulations,
        even when the underlying Metadata objects do not.  Note that
        `metadata` is an alias for this attribute, just based on the
        use of that name for `simulations`, but technically `pandas`
        distinguishes a single row like this as a `Series` object.
    version : str
        Version number of the simulation
    sxs_id_stem : str
        SXS ID without the version number or Lev
    sxs_id : str
        SXS ID with the version number
    location : str
        Location string for the simulation, including the SXS ID,
        version number, and Lev number.
    url : str
        URL for the DOI of the simulation
    files : dict
        Dictionary of file information for the simulation.  The keys
        are of the form "Lev5:Horizons.h5", and the values are
        dictionaries with keys "checksum", "size", and "link".
    lev_numbers : list
        List of available Lev numbers for the simulation.
    lev_number : int
        Chosen Lev number for the simulation.
    horizons : Horizons
        Horizons object for the simulation
    strain : Waveform
        Strain Waveform object for the simulation.  Note that `h` is
        an alias for this attribute, both of which are case
        insensitive: `Strain` and `H` are also acceptable.
    psi4 : Waveform
        Psi4 Waveform object for the simulation.  Note that this
        attribute is also case insensitive: `Psi4` is also acceptable.
    psi3 : Waveform
        Psi3 Waveform object for the simulation.  Note that this
        attribute is also case insensitive: `Psi3` is also acceptable.
        In versions 1 and 2, this attribute will raise an error
        because the data is not available.
    psi2 : Waveform
        Psi2 Waveform object for the simulation.  Note that this
        attribute is also case insensitive: `Psi2` is also acceptable.
        In versions 1 and 2, this attribute will raise an error
        because the data is not available.
    psi1 : Waveform
        Psi1 Waveform object for the simulation.  Note that this
        attribute is also case insensitive: `Psi1` is also acceptable.
        In versions 1 and 2, this attribute will raise an error
        because the data is not available.
    psi0 : Waveform
        Psi0 Waveform object for the simulation.  Note that this
        attribute is also case insensitive: `Psi0` is also acceptable.
        In versions 1 and 2, this attribute will raise an error
        because the data is not available.
    """
    def __init__(self,
        metadata, series, version, sxs_id_stem, sxs_id, url, files, lev_numbers, lev_number, location,
        *args, **kwargs
    ):
        self.series = series
        self.version = version
        self.sxs_id_stem = sxs_id_stem
        self.sxs_id = sxs_id
        self.url = url
        self.files = files
        self.lev_numbers = lev_numbers
        self.lev_number = lev_number
        self.location = location
        self.deprecated = kwargs.get("deprecated", False)
        self.metadata = metadata or self.load_metadata()

    def __repr__(self):
        chi1 = self.series["reference_dimensionless_spin1"]
        chi2 = self.series["reference_dimensionless_spin2"]
        e = self.metadata.reference_eccentricity
        construction = f"""{type(self).__qualname__}("{self.location}")\n# """
        if self.deprecated:
            construction += "DEPRECATED "
        construction += f"n_orbits={self.metadata.number_of_orbits:.3g} "
        construction += f"q={self.metadata.reference_mass_ratio:.3g} "
        construction += f"""chi1=[{", ".join(f"{c:.3g}" for c in chi1)}] """
        construction += f"""chi2=[{", ".join(f"{c:.3g}" for c in chi2)}] """
        construction += f"e={e:.3g} simulation" if type(e) is float else f"{e=} simulation"
        return construction

    def __str__(self):
        return repr(self)
    
    def distances(self, dataframe=None, metadata_metric=None, drop_deprecated=False):
        """Measure the distance from this simulation to others

        Parameters
        ----------
        dataframe : pandas.DataFrame, optional
            DataFrame of simulations to compare to.  If not provided,
            the full catalog of simulations will be loaded as
            `sxs.load("simulations").dataframe`.
        metadata_metric : MetadataMetric, optional
            Metric to use for comparing simulations.  If not provided,
            the default metric will be used.
        drop_deprecated : bool, optional
            If `True`, remove deprecated simulations from the
            `dataframe` before measuring distances.

        Returns
        -------
        distances : pandas.Series
            Distance from this simulation to each element of the
            `dataframe`.  This series will be indexed by the index of
            the `dataframe`.  If this simulation is in the
            `dataframe`, it will have a distance of 0.

        See Also
        --------
        simulations_sorted_by_distance : Sort dataframe of simulations
            by "distance" to this one
        sxs.metadata.metric.MetadataMetric : Metric for comparing
            metadata

        """
        from numpy import sqrt
        from ..metadata.metric import MetadataMetric
        from .. import load
        if dataframe is None:
            dataframe = load("simulations").dataframe
        metadata_metric = metadata_metric or MetadataMetric()
        if drop_deprecated:
            dataframe = dataframe[~dataframe["deprecated"]]
        return dataframe.apply(
            lambda m: sqrt(metadata_metric(self.metadata, m)),
            axis=1
        )

    def closest_simulation(self, dataframe=None, metadata_metric=None, warning_threshold=1e-2):
        """Return the closest undeprecated simulation to this one

        Note that any simulation in `dataframe` with zero distance
        from this one will be ignored; the returned index will not
        refer to this simulation, even if it is undeprecated.

        Parameters
        ----------
        dataframe : pandas.DataFrame, optional
            DataFrame of simulations to compare to.  If not provided,
            the full catalog of simulations will be loaded as
            `sxs.load("simulations").dataframe`.
        metadata_metric : MetadataMetric, optional
            Metric to use for comparing simulations.  If not provided,
            the default metric will be used.
        warning_threshold : float, optional
            Threshold distance above which a warning will be issued
            that the closest simulation is fairly distant.  Default is
            1e-2.

        Returns
        -------
        closest_index : str
            Index of the closest undeprecated simulation in the
            `dataframe`.
        distance : float

        """
        d = self.distances(
            dataframe=dataframe,
            metadata_metric=metadata_metric,
            drop_deprecated=True
        )
        d = d[d > 0].sort_values()
        if d.iloc[0] > warning_threshold:
            warn(f"\nClosest simulation ({d.index[0]}) is fairly distant: {d.iloc[0]:.3g}")
        return d.index[0], d.iloc[0]
    
    @property
    def dataframe(self):
        return self.series

    @property
    def versions(self):
        return self.metadata.DOI_versions

    @property
    def lev(self):
        if self.lev_number is None:
            return ""
        else:
            return f"Lev{self.lev_number}"

    @property
    def Lev(self):
        return self.lev

    @property
    def metadata_path(self):
        for separator in [":", "/"]:
            for ending in [".json", ".txt"]:
                for prefix in ["", f"{self.lev}{separator}" if self.lev else ""]:
                    if (fn := f"{prefix}metadata{ending}") in self.files:
                        return fn
        raise ValueError(
            f"Metadata file not found in simulation files for {self.location}"
        )

    def load_metadata(self):
        from .. import load
        metadata_path = self.metadata_path
        metadata_location = self.files.get(metadata_path)["link"]
        sxs_id_path = Path(self.sxs_id)
        metadata_truepath = Path(sxs_path_to_system_path(sxs_id_path / metadata_path))
        return Metadata(load(metadata_location, truepath=metadata_truepath))

    def load_horizons(self):
        from .. import load
        sxs_id_path = Path(self.sxs_id)
        horizons_path = self.horizons_path
        if horizons_path in self.files:
            horizons_location = self.files.get(horizons_path)["link"]
        else:
            # Some simulations used the SXS ID as a prefix in file paths
            # within the Zenodo upload in version 1.x of the catalog.
            if (extended_location := f"{self.sxs_id_stem}/{horizons_path}") in self.files:
                horizons_location = self.files.get(extended_location)["link"]
            else:
                raise ValueError(
                    f"File '{horizons_path}' not found in simulation files for {self.location}"
                )
        horizons_truepath = Path(sxs_path_to_system_path(sxs_id_path / horizons_path))
        return load(horizons_location, truepath=horizons_truepath)

    @property
    def horizons(self):
        if not hasattr(self, "_horizons"):
            self._horizons = self.load_horizons()
        return self._horizons
    Horizons = horizons

    @property
    def strain(self):
        if not hasattr(self, "_strain"):
            self._strain = self.load_waveform(*self.strain_path)
        return self._strain
    Strain = strain
    h = strain
    H = strain

    # I'm not entirely sure about the conjugations and factors of 2 in
    # shear and news in our conventions.  These will have to wait for
    # later.
    #
    # @property
    # def shear(self):
    #     if not hasattr(self, "_shear"):
    #         self._shear = self.strain.bar / 2
    #     return self._shear
    # sigma = shear
    # σ = shear
    # Shear = shear
    # Sigma = shear
    # Σ = shear
    #
    # @property
    # def news(self):
    #     if not hasattr(self, "_news"):
    #         self._news = self.strain.dot
    #     return self._news
    # News = news

    @property
    def psi4(self):
        if not hasattr(self, "_psi4"):
            self._psi4 = self.load_waveform(*self.psi4_path)
        return self._psi4
    Psi4 = psi4

    @property
    def psi3(self):
        raise AttributeError(f"Psi3 is not available for version {self.version} of the data")
    Psi3 = psi3

    @property
    def psi2(self):
        raise AttributeError(f"Psi2 is not available for version {self.version} of the data")
    Psi2 = psi2

    @property
    def psi1(self):
        raise AttributeError(f"Psi1 is not available for version {self.version} of the data")
    Psi1 = psi1

    @property
    def psi0(self):
        raise AttributeError(f"Psi0 is not available for version {self.version} of the data")
    Psi0 = psi0

    def to_lvk(self, **kwargs):
        r"""Convert an SXS simulation to LVK convention.

        Returns an SXS waveform (modes or polarizations) and dynamics
        (including angular velocities, frame quaternions, and spins) in
        the inertial frame that coincides with the waveform-defined frame
        defined at a reference time `t_ref` or reference frequency
        `f_ref`.

        Parameters
        ----------
        t_ref : float, optional
            The reference time at which the waveform frame is specified.
            This is measured in units of M, and defined relative to the
            epoch time (see below).  Either `t_ref` or `f_ref` must be
            specified.  If `t_ref` is given, it is used to compute
            `f_ref`.
        f_ref : float, optional
            The reference frequency, in units of cycles/M, at which the
            waveform frame is specified.  Either `t_ref` or `f_ref` must
            be specified.  If `f_ref` is given, it is used to compute
            `t_ref`.
        dt : float, optional
            The time step, in units of M, to which to interpolate the
            waveform.
        f_low : float, optional
            The lower frequency bound, in units of cycles/M, for the
            waveform.
        ell_max : int, optional
            The maximum ell to include in the waveform.
        phi_ref : float, optional
            The binary's phase in the coprecessing frame, measured at
            `t_ref`. Should be between 0 and $2\pi$.
        inclination : float, optional
            Angle between the binary's angular momentum and the line of
            sight of the observer, measured at `t_ref`.  Should be between
            0 and $\pi$.
        ell_max_epoch : int, optional
            The maximum ell to include in the epoch time calculation,
            which sets t=0 at the maximum of the L^2 norm, calculated by
            including all modes up to and including this ell value.

        Returns
        -------
        times : float array
            Uniformly spaced 1D array of times, in units of M, at which
            the waveform and dynamics quantities are returned.  Aligned
            such that peak of waveform modes with ell=2 is at t=0.
        hlm_dict : dict [optional]
            Dictionary of waveform modes in the inertial frame that
            coincides with the waveform-defined coprecessing frame at
            `f_ref`.  Each mode in the dictionary is a 1D array of
            complex-valued floats with values corresponding to each time
            node.  Keys:[(ell,m)] for all ell<=ell_max and -ell<=m<=+ell.
            This is returned only if the input values of `phi_ref` and
            `inclination` are both None.
        hp, hc : float arrays [optional]
            1D-arrays of real-valued GW polarizations evaluated in the
            frame of the observer at each time node.  Polarizations are
            computed using all modes up to ell=ell_max, with one value at
            each time node.  These are returned only if either of the
            input values of `phi_ref` and `inclination` is not None.
        dynamics_dict : dict
            Dictionary of real-valued arrays of dynamics quantities:
                * "t_ref": The reference time at which the waveform frame
                is specified.
                * "f_ref": The waveform's frequency at `t_ref`.
                * "t_low": The earliest time in the waveform.
                * "f_low": The waveform's frequency at `t_low`.
                * "chi1_ref": Cartesian spin components for the more
                massive object, evaluated at `t_ref` in the
                waveform-defined inertial frame.
                * "chi2_ref": Cartesian spin components for the less
                massive object, evaluated at `t_ref` in the
                waveform-defined inertial frame.
                * "frame_quat": Quaternions describing the instantaneous
                transformation from the inertial frame to the corotating
                frame at each time node as described in arXiv:1905.09300
                and using conventions in Appendix B of arXiv:1110.2965.
                Array of shape (len(times),4) with four quaternion
                components at each time node.  The first element at each
                time node is the scalar part.
                * "frame_omega": Angular velocity vector of the corotating
                frame at each time node.  Array of shape (len(times),3)
                with three components at each time node.
                * "times_spins": The times, in units of M, at which `chi`
                and `chi2` are returned.  Note that these are coordinate
                times deep within the dynamic region of the simulation,
                and so cannot be precisely related to the times in the
                asymptotic waveform.
                * "chi1": Cartesian spin components for the more massive
                object, evaluated in the waveform-defined inertial
                frame.  Array of shape (len(times_spins),3) which
                contains the Cartesian spin components
                \{chi_{1x},chi_{1y},chi_{1z}\} at each time node.
                * "chi2": Cartesian spin components for the less massive
                object, evaluated in the waveform-defined inertial
                frame.  Array of shape (len(times_spins),3) which
                contains the Cartesian spin components
                \{chi_{2x},chi_{2y},chi_{2z}\} at each time node.
                * "times_remnant": The times, in units of M, at which
                `chi_remnant` and `mass_remnant` are returned.  Note
                that these are coordinate times deep within the dynamic
                region of the simulation, and so cannot be precisely
                related to the times in the asymptotic waveform.
                * "chi_remnant": Cartesian spin components for the remnant
                black hole, evaluated in the waveform-defined inertial
                frame.  Array of shape (len(times_remnant),3) which
                contains the Cartesian spin components
                \{chi_{rx},chi_{ry},chi_{rz}\}.
                * "mass_remnant": The Christodoulou mass of the remnant
                black hole as a function of time.  Array of shape
                (len(times_remnant),).

        See also
        ========
        sxs.load : General-purpose function to load SXS data in native
            format
        sxs.waveforms.to_lvc_conventions : Inner function that does all
            the work for this function

        Conventions
        ===========
        We assume geometric units for time (units of M) and frequency
        (units of cycles/M), with total mass M equal to 1.
        
        Epoch time is defined by the peak of the $L^2$ norm of the modes:

            $$t_e = \argmax_t \sum_\ell \sum_m |h_{\ell,m}(t)|^2,$$

        where the sum over $\ell$ ranges from 2 to `ell_max_epoch`.  All
        time axes are then adjusted so that $t_e = 0$.

        Frequencies are measured from the waveform, rather than orbital
        trajectories, in terms of the angular-velocity vector given by
        equation (7) of arXiv:1302.2919.  The frequency is defined as the
        magnitude of this vector divided by $2\pi$.

        Waveforms, spins, dynamics, times, inclination angle and GW
        reference phase are all defined in the inertial frame that
        coincides with the "waveform-defined frame" at the reference time.
        This frame is chosen so that the $z$ axis is aligned with the
        dominant principal axis of the matrix given by equation (2a) of
        arXiv:1109.5224, except that the strain is used in place of
        $\psi_4$.  That axis is only determined up to a sign, so we choose
        the positive $z$ axis to be more parallel than antiparallel to the
        angular velocity.  Rotation about the $z$ axis is chosen to
        approximate the condition that the more massive black hole is
        located on the positive $x$ axis at the reference time, but can be
        written solely in terms of the waveform modes:

            * $\Im{h_{2,2} + \bar{h}_{2,-2}} = 0$
            * $\Re{h_{2,2} + \bar{h}_{2,-2}} < 0$
            * $\Im{h_{2,1} + \bar{h}_{2,-1}} < 0$

        The first two conditions are necessary for cases of symmetric
        systems in which $h_{2,\pm 1}=0$; the last condition breaks the
        degeneracy of the first two under rotation about the $z$ axis by
        $\pi$.  For configurations that are symmetric under that rotation,
        $h_{2,1}$ will be zero, so this condition will be impossible to
        apply, but the symmetry of the system will mean that there is no
        difference in the result.

        Quaternion conventions are described in Appendix B of
        arXiv:1110.2965.  In particular, we use the convention that

            Q = q_0 + vec(q) -> (q_0, q_1, q_2, q_3)

        where q_0 is the scalar part.

        """
        from ..waveforms.format_handlers.lvc import to_lvc_conventions
        strain = self.load_waveform(
            *self.strain_path,
            transform_to_inertial=False,
        )
        return to_lvc_conventions(strain, self.horizons, **kwargs)


class Simulation_v1(SimulationBase):
    """Simulation object for version 1 of the data format
    
    Note that users almost certainly never need to call this function;
    see the `Simulation` function or `sxs.load` function instead.  See
    also `SimulationBase` for the base class that this class inherits
    from.
    """
    # We have to deal with the fact that some early file paths on
    # Zenodo included the SXS ID as a prefix, while others did not.
    # This means that we have to check for both possibilities in
    # `load_horizons` and `load_waveform`.

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.extrapolation = kwargs.get("extrapolation", "N2")

    @property
    def horizons_path(self):
        return search_prefixes("Horizons.h5", self.lev, self.files)

    @property
    def strain_path(self):
        extrapolation = (
            f"Extrapolated_{self.extrapolation}.dir"
            if self.extrapolation != "Outer"
            else "OutermostExtraction.dir"
        )
        return (
            search_prefixes("rhOverM_Asymptotic_GeometricUnits_CoM.h5", self.lev, self.files),
            extrapolation
        )

    @property
    def psi4_path(self):
        extrapolation = (
            f"Extrapolated_{self.extrapolation}.dir"
            if self.extrapolation != "Outer"
            else "OutermostExtraction.dir"
        )
        return (
            search_prefixes("rMPsi4_Asymptotic_GeometricUnits_CoM.h5", self.lev, self.files),
            extrapolation
        )

    def load_waveform(self, file_name, group, transform_to_inertial=True):
        # Note: `transform_to_inertial` is a slightly unnatural argument in version 1,
        # since the data are actually stored in the inertial frame.  If the value is
        # `False`, we will transform the data to the corotating frame, which is what
        # the corresponding argument in version 2 achieves naturally (since the data
        # are stored in the corotating frame in that version).
        from .. import load
        if file_name in self.files:
            location = self.files.get(file_name)["link"]
        else:
            # Some simulations used the SXS ID as a prefix in file paths
            # within the Zenodo upload in version 1.x of the catalog.
            if (extended_file_name := f"{self.sxs_id_stem}/{file_name}") in self.files:
                location = self.files.get(extended_file_name)["link"]
            else:
                raise ValueError(f"File '{file_name}' not found in simulation files")
        sxs_id_path = Path(self.sxs_id)
        truepath = Path(sxs_path_to_system_path(sxs_id_path / file_name))
        w = load(
            location, truepath=truepath, extrapolation_order=group,
            transform_to_inertial=transform_to_inertial
        )
        w.metadata = self.metadata
        return w


class Simulation_v2(SimulationBase):
    """Simulation object for version 2 of the data format
        
    Note that users almost certainly never need to call this function;
    see the `Simulation` function or `sxs.load` function instead.  See
    also `SimulationBase` for the base class that this class inherits
    from.
    """
    # Default extrapolation order for this simulation version
    default_extrapolation = "N2"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.extrapolation = kwargs.get("extrapolation", self.default_extrapolation)

    @property
    def horizons_path(self):
        return search_prefixes("Horizons.h5", self.lev, self.files)

    @property
    def strain_path(self):
        return (
            search_prefixes(f"Strain_{self.extrapolation}", self.lev, self.files, ".h5"),
            "/"
        )

    @property
    def psi4_path(self):
        extrapolation = (
            f"Extrapolated_{self.extrapolation}.dir"
            if self.extrapolation != "Outer"
            else "OutermostExtraction.dir"
        )
        prefix = f"{self.lev}:" if len(self.lev_numbers)>1 else ""
        return (
            search_prefixes(f"ExtraWaveforms", self.lev, self.files, ".h5"),
            f"/rMPsi4_Asymptotic_GeometricUnits_CoM_Mem/{extrapolation}"
        )

    def load_waveform(self, file_name, group, transform_to_inertial=True):
        from .. import load
        # Note that `name` should not have the file ending on input,
        # but we will replace it regardless with `.with_suffix`.
        file_name = Path(file_name)
        sxs_id_path = Path(self.sxs_id)
        h5_path = str(file_name.with_suffix(".h5"))
        json_path = str(file_name.with_suffix(".json"))
        h5_location = self.files.get(h5_path)["link"]
        json_location = self.files.get(json_path)["link"]
        h5_truepath = Path(sxs_path_to_system_path(sxs_id_path / h5_path))
        json_truepath = Path(sxs_path_to_system_path(sxs_id_path / json_path))
        json_truepath = sxs_directory("cache") / json_truepath
        if not Path(json_location).exists() and not json_truepath.exists():
            if not read_config("download", True):
                raise ValueError(f"{json_truepath} not found and download is disabled")
            download_file(json_location, json_truepath)
        return load(
            h5_location, truepath=h5_truepath, group=group, metadata=self.metadata,
            transform_to_inertial=transform_to_inertial
        )


class Simulation_v3(Simulation_v2):
    # Default extrapolation order for this simulation version
    default_extrapolation = "N2"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.extrapolation = kwargs.get("extrapolation", self.default_extrapolation)

    @property
    def strain_path(self):
        return (
            search_prefixes(f"Strain_{self.extrapolation}", self.lev, self.files, ".h5"),
            "/"
        ) if self.extrapolation == self.default_extrapolation else (
            search_prefixes("ExtraWaveforms", self.lev, self.files, ".h5"),
            f"/Strain_{self.extrapolation}.dir"
        )

    @property
    def psi4_path(self):
        return (
            search_prefixes("ExtraWaveforms", self.lev, self.files, ".h5"),
            f"/Psi4_{self.extrapolation}.dir"
        )


def get_file_info(metadata, sxs_id, download=None):
    from .. import load, load_via_sxs_id
    if "files" in metadata:
        return metadata["files"]
    truepath = Path(sxs_path_to_system_path(sxs_id)) / "zenodo_metadata.json"
    record = load_via_sxs_id(sxs_id, "export/json", truepath=truepath, download=download)
    if not record.get("files", {}).get("entries", []):
        # CaltechDATA files should generally be already stored in the `simulations`,
        # but just in case, we get to this point
        truepath = Path(sxs_path_to_system_path(sxs_id)) / "caltechdata_files.json"
        url = record["links"]["files"]
        record = load(url, truepath=truepath, download=download)
        entries = record["entries"]
        return {
            entry["key"]: {
                "checksum": entry["checksum"],
                "size": entry["size"],
                "link": entry["links"]["content"],
            }
            for entry in entries
        }
    else:
        entries = record["files"]["entries"]
        return {
            str(filename): {
                "checksum": entry["checksum"],
                "size": entry["size"],
                "link": entry["links"]["content"],
            }
            for filename, entry in entries.items()
        }
