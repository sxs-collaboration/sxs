from pathlib import Path
from warnings import warn
from .. import doi_url, Metadata
from ..utilities import (
    sxs_id_and_version, lev_number, sxs_path_to_system_path,
    download_file, sxs_directory, read_config
)


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

    Other parameters
    ----------------
    ignore_deprecation : bool
        If `True`, completely bypass checking for deprecation or
        supersession.  No warnings or errors will be issued.
    auto_supersede : bool
        If `True`, automatically load the superseding simulation, if
        there is only one.  If there are multiple superseding
        simulations, an error will be raised, and you must explicitly
        choose one.  If there is only one, a warning will still be
        issued, but the superseding simulation will be loaded.  Note
        that this can also be set in the configuration file with
        `sxs.write_config(auto_supersede=True)`.
    extrapolation : str
        The extrapolation order to use for the strain and Psi4 data.
        This is only relevant for versions 1 and 2 of the data format,
        both of which default to "N2".  Other options include "N3",
        "N4", and "Outer".  "Nx" refers to extrapolation by
        polynomials in 1/r with degree `x`, while "Outer" refers to
        data extracted at the outermost extraction radius but
        corrected for time-dilation and areal-radius effects.

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
    from .. import load

    # Extract the simulation ID, version, and Lev from the location string
    simulation_id, input_version = sxs_id_and_version(location)
    if not simulation_id:
        raise ValueError(f"Invalid SXS ID in '{simulation_id}'")
    input_lev_number = lev_number(location)  # Will be `None` if not present

    # Load the simulation catalog and check if simulation ID exists in the catalog
    simulations = load("simulations")
    if simulation_id not in simulations:
        raise ValueError(f"Simulation '{simulation_id}' not found in simulation catalog")

    # Attach metadata to this object
    metadata = Metadata(simulations[simulation_id])
    series = simulations.dataframe.loc[simulation_id]

    # Check if the specified version exists in the simulation catalog
    if input_version not in metadata.DOI_versions:
        raise ValueError(f"Version '{input_version}' not found in simulation catalog for '{simulation_id}'")

    # Set various pieces of information about the simulation
    version = input_version or max(metadata.DOI_versions)
    if not version.startswith("v"):
        raise ValueError(f"Invalid version string '{version}'")
    sxs_id_stem = simulation_id
    sxs_id = f"{sxs_id_stem}{version}"
    url = f"{doi_url}{sxs_id}"

    # Deal with "superseded_by" field, or "deprecated" keyword in the metadata
    if not kwargs.get("ignore_deprecation", False):
        auto_supersede = kwargs.get("auto_supersede", read_config("auto_supersede", False))
        if (
            input_version
            and not auto_supersede
            and ("deprecated" in metadata.get("keywords", []) or metadata.get("superseded_by", False))
        ):
            message = ("\n"
                + f"Simulation '{sxs_id_stem}' is deprecated and/or superseded.\n"
                + "Normally, this simulation should no longer be used, but you\n"
                + f"explicitly requested version '{input_version}', so it is being used.\n"
            )
            warn(message)
        else:
            if "superseded_by" in metadata:
                superseded_by = metadata["superseded_by"]
                if auto_supersede and isinstance(superseded_by, list):
                    raise ValueError(
                        f"`auto_supersede` is enabled, but simulation '{sxs_id}' is\n"
                        + "superseded by multiple simulations.  You must choose one\n"
                        + "explicitly from the list:\n"
                        + "\n".join(f"  {s}" for s in superseded_by)
                        + "\nAlternatively, you could pass `ignore_deprecation=True` or\n"
                        + "specify a version to load this waveform anyway."
                    )
                elif auto_supersede and isinstance(superseded_by, str):
                    message = f"\nSimulation '{sxs_id}' is being automatically superseded by '{superseded_by}'."
                    warn(message)
                    new_location = f"{superseded_by}{input_version}"
                    if input_lev_number:
                        new_location += f"/Lev{input_lev_number}"
                    return Simulation(new_location, *args, **kwargs)
                elif isinstance(superseded_by, list):
                    raise ValueError(
                        f"Simulation '{sxs_id}' is superseded by multiple simulations.\n"
                        + "Even if you enable `auto_supersede`, with multiple options, you\n"
                        + "must choose one explicitly from the list:\n"
                        + "\n".join(f"  {s}" for s in superseded_by)
                        + "\nAlternatively, you could pass `ignore_deprecation=True` or\n"
                        + "specify a version to load this waveform anyway."
                    )
                elif isinstance(superseded_by, str):
                    raise ValueError(
                        f"Simulation '{sxs_id}' is superseded by '{superseded_by}'.\n"
                        + "Note that you could enable `auto_supersede` to automatically\n"
                        + "load the superseding simulation.  Alternatively, you could\n"
                        + "pass `ignore_deprecation=True` or specify a version to load\n"
                        + "this waveform anyway."
                    )
                else:
                    raise ValueError(
                        f"Simulation '{sxs_id}' is superseded by '{superseded_by}'.\n"
                        + "Note that you could pass `ignore_deprecation=True` or\n"
                        + "specify a version to load this waveform anyway."
                    )
            if "deprecated" in metadata.get("keywords", []):
                raise ValueError(
                    f"Simulation '{sxs_id}' is deprecated but has no superseding simulation.\n"
                    + "Note that you could pass `ignore_deprecation=True` or specify a version\n"
                    + "to  to load this waveform anyway."
                )

    # We want to do this *after* deprecation checking, to avoid possibly unnecessary web requests
    files = get_file_info(metadata, sxs_id)

    # If Lev is given as part of `location`, use it; otherwise, use the highest available
    lev_numbers = sorted({lev for f in files if (lev:=lev_number(f))})
    output_lev_number = input_lev_number or max(lev_numbers)
    location = f"{sxs_id_stem}{version}/Lev{output_lev_number}"

    # Finally, figure out which version of the simulation to load and dispatch
    version_number = float(version[1:])
    if 1 <= version_number < 2.0:
        return Simulation_v1(
            metadata, series, version, sxs_id_stem, sxs_id, url, files, lev_numbers, output_lev_number, location, *args, **kwargs
        )
    elif 2 <= version_number < 3.0:
        return Simulation_v2(
            metadata, series, version, sxs_id_stem, sxs_id, url, files, lev_numbers, output_lev_number, location, *args, **kwargs
        )
    else:
        raise ValueError(f"Version '{version}' not yet supported")


class SimulationBase:
    """Base class for Simulation objects

    Note that users almost certainly never need to call this function;
    see the `Simulation` function or `sxs.load` function instead.
    
    Attributes
    ----------
    metadata : Metadata
        Metadata object for the simulation
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
        self.metadata = metadata
        self.series = series
        self.version = version
        self.sxs_id_stem = sxs_id_stem
        self.sxs_id = sxs_id
        self.url = url
        self.files = files
        self.lev_numbers = lev_numbers
        self.lev_number = lev_number
        self.location = location

    def __repr__(self):
        return f"""{type(self).__qualname__}("{self.sxs_id}")"""

    def __str__(self):
        return repr(self)
    
    @property
    def dataframe(self):
        return self.series

    @property
    def versions(self):
        return self.metadata.DOI_versions

    @property
    def lev(self):
        return f"Lev{self.lev_number}"

    @property
    def Lev(self):
        return self.lev

    def load_horizons(self):
        from .. import load
        sxs_id_path = Path(sxs_path_to_system_path(self.sxs_id))
        horizons_path = self.horizons_path
        horizons_location = self.files.get(horizons_path)["link"]
        horizons_truepath = sxs_id_path / sxs_path_to_system_path(horizons_path)
        return load(horizons_location, truepath=horizons_truepath)

    @property
    def horizons(self):
        if not hasattr(self, "_horizons"):
            self._horizons = self.load_horizons()
        return self._horizons

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
        return f"{self.lev}/Horizons.h5"
    
    def load_horizons(self):
        from .. import load
        sxs_id_path = Path(sxs_path_to_system_path(self.sxs_id))
        horizons_path = self.horizons_path
        if horizons_path in self.files:
            horizons_location = self.files.get(horizons_path)["link"]
        else:
            if (extended_horizons_path := f"{self.sxs_id_stem}/{horizons_path}") in self.files:
                horizons_location = self.files.get(extended_horizons_path)["link"]
            else:
                raise ValueError(f"File '{horizons_path}' not found in simulation files")
        horizons_truepath = sxs_id_path / sxs_path_to_system_path(horizons_path)
        return load(horizons_location, truepath=horizons_truepath)

    @property
    def strain_path(self):
        extrapolation = (
            f"Extrapolated_{self.extrapolation}.dir"
            if self.extrapolation != "Outer"
            else "OutermostExtraction.dir"
        )
        return (
            f"{self.lev}/rhOverM_Asymptotic_GeometricUnits_CoM.h5",
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
            f"{self.lev}/rMPsi4_Asymptotic_GeometricUnits_CoM.h5",
            extrapolation
        )

    def load_waveform(self, file_name, group):
        from .. import load
        if file_name in self.files:
            location = self.files.get(file_name)["link"]
        else:
            if (extended_file_name := f"{self.sxs_id_stem}/{file_name}") in self.files:
                location = self.files.get(extended_file_name)["link"]
            else:
                raise ValueError(f"File '{file_name}' not found in simulation files")
        sxs_id_path = Path(sxs_path_to_system_path(self.sxs_id))
        truepath = sxs_id_path / sxs_path_to_system_path(file_name)
        w = load(location, truepath=truepath, extrapolation_order=group)
        w.metadata = self.metadata
        return w


class Simulation_v2(SimulationBase):
    """Simulation object for version 2 of the data format
        
    Note that users almost certainly never need to call this function;
    see the `Simulation` function or `sxs.load` function instead.  See
    also `SimulationBase` for the base class that this class inherits
    from.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.extrapolation = kwargs.get("extrapolation", "N2")

    @property
    def horizons_path(self):
        return f"{self.lev}:Horizons.h5"

    @property
    def strain_path(self):
        return (
            f"{self.lev}:Strain_{self.extrapolation}",
            "/"
        )

    @property
    def psi4_path(self):
        extrapolation = (
            f"Extrapolated_{self.extrapolation}.dir"
            if self.extrapolation != "Outer"
            else "OutermostExtraction.dir"
        )
        return (
            f"{self.lev}:ExtraWaveforms",
            f"/rMPsi4_Asymptotic_GeometricUnits_CoM_Mem/{extrapolation}"
        )

    def load_waveform(self, file_name, group):
        from .. import load
        # Note that `name` should not have the file ending on input,
        # but we will strip it regardless with `.stem`.
        file_name = Path(file_name).stem
        sxs_id_path = Path(sxs_path_to_system_path(self.sxs_id))
        h5_path = f"{file_name}.h5"
        json_path = f"{file_name}.json"
        h5_location = self.files.get(h5_path)["link"]
        json_location = self.files.get(json_path)["link"]
        h5_truepath = sxs_id_path / sxs_path_to_system_path(h5_path)
        json_truepath = sxs_id_path / sxs_path_to_system_path(json_path)
        if not json_truepath.exists():
            download_file(json_location, sxs_directory("cache") / json_truepath)
        return load(h5_location, truepath=h5_truepath, group=group, metadata=self.metadata)


def get_file_info(metadata, sxs_id):
    from .. import load_via_sxs_id
    if "files" in metadata:
        return metadata["files"]
    truepath = Path(sxs_path_to_system_path(sxs_id)) / "zenodo_metadata.json"
    record = load_via_sxs_id(sxs_id, "export/json", truepath=truepath)
    entries = record["files"]["entries"]
    return {
        str(filename): {
            "checksum": entry["checksum"],
            "size": entry["size"],
            "link": entry["links"]["content"],
        }
        for filename, entry in entries.items()
    }
