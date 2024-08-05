from pathlib import Path
from .. import doi_url, Metadata
from ..utilities import (
    sxs_id_and_version, lev_number, sxs_path_to_system_path,
    download_file, sxs_directory,
)
from . import default_version

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
    
    """
    _, version = sxs_id_and_version(location)
    if not version:
        version = default_version
    if not version.startswith("v"):
        raise ValueError(f"Invalid version '{version}'")
    version_number = float(version[1:])
    if 1 <= version_number < 2.0:
        return Simulation_v1(location)
    elif 2 <= version_number < 3.0:
        return Simulation_v2(location)
    else:
        raise ValueError(f"Version '{version}' not yet supported")


class SimulationBase:
    def __init__(self, location, *args, **kwargs):
        from .. import load
        simulation_id, version = sxs_id_and_version(location)
        if not simulation_id:
            raise ValueError(f"Invalid SXS ID in '{simulation_id}'")
        simulations = load("simulations")
        if simulation_id not in simulations:
            raise ValueError(f"Simulation '{simulation_id}' not found in simulation catalog")
        self.metadata = Metadata(simulations[simulation_id])
        if version not in self.metadata.DOI_versions:
            raise ValueError(f"Version '{version}' not found in simulation catalog for '{simulation_id}'")
        self.version = version or default_version
        self.id = simulation_id
        self.sxs_id = f"{self.id}{self.version}"
        self.url = f"{doi_url}{self.sxs_id}"
        self.files = get_file_info(self.metadata, self.sxs_id)
        self.highest_lev_number = max(
            lev for f in self.files
            if (lev:=lev_number(f))
        )
        self.lev_number = lev_number(location) or self.highest_lev_number

    def __repr__(self):
        return f"""{type(self).__qualname__}("{self.sxs_id}")"""

    def __str__(self):
        return repr(self)

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

    ## I'm not entirely sure about the conjugations and factors of 2 in shear and news
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
    # We have to deal with the fact that some early file paths on
    # Zenodo included the SXS ID as a prefix, while others did not.
    # This means that we have to check for both possibilities in
    # `load_horizons` and `load_waveform`.

    def __init__(self, location, *args, **kwargs):
        super().__init__(location, *args, **kwargs)
        self.extrapolation = kwargs.get("extrapolation", "N4")

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
            if (extended_horizons_path := f"{self.id}/{horizons_path}") in self.files:
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
            if (extended_file_name := f"{self.id}/{file_name}") in self.files:
                location = self.files.get(extended_file_name)["link"]
            else:
                raise ValueError(f"File '{file_name}' not found in simulation files")
        sxs_id_path = Path(sxs_path_to_system_path(self.sxs_id))
        truepath = sxs_id_path / sxs_path_to_system_path(file_name)
        w = load(location, truepath=truepath, extrapolation_order=group)
        w.metadata = self.metadata
        return w


class Simulation_v2(SimulationBase):
    def __init__(self, location, *args, **kwargs):
        super().__init__(location, *args, **kwargs)
        self.extrapolation = kwargs.get("extrapolation", "N4")

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
