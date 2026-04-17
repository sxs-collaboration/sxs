from pathlib import Path
from .simulation import SimulationBase
from .. import Metadata
from ..utilities import download_file, sxs_directory


def _not_defined_property(name):
    @property
    def prop(self):
        raise AttributeError(f"'{name}' is not supported for version {self.version} of the catalog.")

    return prop


def _not_implemented_function(name):
    def function(self, *args, **kwargs):
        raise NotImplementedError(f"'{name}' cannot be implemented for version {self.version} of the data.")

    return function


def RITSimulation(location, *args, **kwargs):
    """Construct a RITSimulation object from a location string.

    The location string should be an RIT ID, as in "RIT:BBH:xxxx" or
    "RIT:eBBH:xxxx". It must exist, or no files will be found to load the data
    from.

    The returned object will be a `RITSimulation_v4` object.

    Parameters
    ----------
    location : str
        The location string for the simulation.  This must include the
        RIT ID as described above.

    Returns
    -------
    simulation : SimulationBase
        A `RITSimulation_v4` object, depending on the version of the simulation
        data.

    Note that all remaining arguments (including keyword arguments)
    are passed on to the `SimulationBase` constructors.

    """
    from .. import load

    # Load the simulation catalog
    simulations = load("RITsimulations")

    # Check if simulation ID exists in the catalog
    if location not in simulations:
        raise ValueError(f"Simulation '{location}' not found in RIT simulations catalog.")

    # Attach metadata to this object
    series = simulations.dataframe.loc[location]
    metadata = Metadata(dict(series))
    files = {
        "strain": {"link": metadata.extrap_strain_url},
        "psi4": {"link": metadata.extrap_psi4_url},
        "metadata": {"link": metadata.metadata_url},
    }
    version = "v4.0"

    sim = RITSimulation_v4(
        metadata=metadata,
        series=series,
        version=version,
        files=files,
        location=location,
    )
    sim.__file__ = str(sxs_directory("cache") / location)
    return sim


class RITSimulation_v4(SimulationBase):
    """Simulation object for version 4 of the RIT data format

    Note that users almost certainly never need to call this function;
    see the `RITSimulation` function or `sxs.load` function instead.  See
    also `SimulationBase` for the base class that this class inherits
    from.
    """

    def __init__(self, metadata, series, version, files, location, *args, **kwargs):
    # Note that we want just a subset of the constructor of SimulationBase here.
    # Therefore, we don't use the super of SimulationsBase.
    # This keeps the constructor clean and avoids setting positional arguments -
    # sxs_id_stem, sxs_id, url, lev_numbers and lev_number to None.
        self.metadata = metadata
        self.series = series
        self.version = version
        self.location = location
        self.files = files

    def __repr__(self):
        chi1 = self.series["relaxed_chi1"]
        chi2 = self.series["relaxed_chi2"]
        e = self.series.eccentricity
        construction = f"""{type(self).__qualname__}("{self.location}")\n# """
        construction += f"n_orbits={self.metadata.number_of_orbits:.3g} "
        construction += f"q={(1/self.metadata.relaxed_mass_ratio_1_over_2):.3g} "
        construction += f"""chi1=[{", ".join(f"{c:.3g}" for c in chi1)}] """
        construction += f"""chi2=[{", ".join(f"{c:.3g}" for c in chi2)}] """
        construction += f"e={e:.3g} simulation" if type(e) is float else f"{e=} simulation"
        return construction

    def load_waveform(self, location):

        from ..waveforms import lvcnr

        strain_url = self.metadata.extrap_strain_url
        filename = Path(strain_url).name

        path = sxs_directory("cache") / location / filename

        if not path.exists():
            download_file(strain_url, path)

        w = lvcnr.load(path)
        w.metadata = self.metadata

        return w

    distances = _not_implemented_function("distances")
    closest_simulation = _not_implemented_function("closest_simulation")
    load_horizons = _not_implemented_function("load_horizons")

    @property
    def strain(self):
        if not hasattr(self, "_strain"):
            self._strain = self.load_waveform(self.location)
        return self._strain

    Strain = strain
    h = strain
    H = strain

    versions = _not_defined_property("versions")
    lev = _not_defined_property("lev")
    Lev = _not_defined_property("Lev")
    metadata_path = _not_defined_property("metadata_path")
    horizons = _not_defined_property("horizons")
