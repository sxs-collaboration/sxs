"""Container interface to the catalog of SXS simulations"""

import collections
import numpy as np
import pandas as pd

class SimulationsDataFrame(pd.DataFrame):
    @property
    def BBH(self):
        """Restrict dataframe to just binary black hole systems"""
        return type(self)(self[self["object_types"] == "BHBH"])
    BHBH = BBH
    
    @property
    def BHNS(self):
        """Restrict dataframe to just black hole-neutron star systems"""
        return type(self)(self[self["object_types"] == "BHNS"])
    NSBH = BHNS
    
    @property
    def NSNS(self):
        """Restrict dataframe to just binary neutron star systems"""
        return type(self)(self[self["object_types"] == "NSNS"])
    BNS = NSNS

    @property
    def noneccentric(self):
        """Restrict dataframe to just non-eccentric systems (e<1e-3)"""
        return type(self)(self[self["reference_eccentricity_bound"] < 1e-3])

    @property
    def eccentric(self):
        """Restrict dataframe to just eccentric systems (e>=1e-3)"""
        return type(self)(self[self["reference_eccentricity_bound"] >= 1e-3])
    
    @property
    def nonprecessing(self):
        """Restrict dataframe to just nonprecessing systems
        
        The criterion used here is that the sum of the x-y components
        of the spins is less than 1e-3 at the reference time.
        """
        return type(self)(self[
            (self["reference_chi1_perp"] + self["reference_chi2_perp"]) < 1e-3
        ])
    
    @property
    def precessing(self):
        """Restrict dataframe to just precessing systems
        
        The criterion used here is that the sum of the x-y components
        of the spins is at least 1e-3 at the reference time.
        """
        return type(self)(self[
            (self["reference_chi1_perp"] + self["reference_chi2_perp"]) >= 1e-3
        ])
    
    @property
    def IMR(self):
        """Restrict dataframe to just BBH IMR systems
        
        "IMR" stands for inspiral, merger, and ringdown.  Systems that
        will *not* be in this group include simulations that
        correspond to physical IMR systems, but were not continued
        through the merger.

        The criteria used here are just that the reference
        eccentricity and remnant mass are actual (finite) numbers.
        Currently, at least, the existence of a measured eccentricity
        means that the system is not hyperbolic or head-on.
        """
        df = self.BBH
        return type(df)(df[
            np.isfinite(df["reference_eccentricity"])
            & np.isfinite(df["remnant_mass"])
        ])
    
    @property
    def hyperbolic(self):
        """Restrict dataframe to just hyperbolic systems
        
        The criterion used here is that the (normalized) ADM mass is
        greater than 1.
        """
        total_mass = self["initial_mass1"] + self["initial_mass2"]
        normalized_ADM = self["initial_ADM_energy"] / total_mass
        return type(self)(self[
            np.isfinite(total_mass) & (total_mass > 0) & (normalized_ADM > 1)
        ])


class Simulations(collections.OrderedDict):
    """Interface to the catalog of SXS simulations
    
    Creation
    --------
    You probably don't need to create this object yourself.  The
    easiest way to create this object is just to use the `sxs.load`
    function:

    ```python
    import sxs

    simulations = sxs.load("simulations")
    ```

    Note SXS members may also wish to read a local copy of the
    simulation annex, which can be done with
    ```python
    simulations = sxs.load("simulations", annex_dir="/path/to/SimulationAnnex.git")
    ```
    which will re-read the annex (which may take about a minute), or
    ```python
    simulations = sxs.load("simulations", local=True)
    ```
    if the annex has not been updated since the last time you
    used the `annex_dir` argument.  Once you have done this,
    calls to `sxs.load` will automatically use this local copy
    of the simulations.
    """
    last_modified_url = "https://api.github.com/repos/sxs-collaboration/sxs/contents/simulations.json?ref=simulations"
    url = "https://github.com/sxs-collaboration/sxs/raw/simulations/simulations.json"

    def __init__(self, sims):
        """Initialize the Simulations dictionary

        Note that the constructor is not generally useful from outside
        this class.  See `Simulations.load` for a more useful
        initialization function, or simply call
        `sxs.load("simulations")`.

        """
        from .. import Metadata
        super(Simulations, self).__init__(
            (k, Metadata(sims[k])) for k in sorted(sims)
        )

    @classmethod
    def remote_timestamp(cls, download):
        import requests
        from datetime import datetime, timezone
        if not download:
            return datetime.min.replace(tzinfo=timezone.utc)
        failed = False
        try:
            response = requests.head(
                Simulations.last_modified_url,
                headers={"X-GitHub-Api-Version": "2022-11-28"},
            )
            if response.status_code != 200 or "Last-Modified" not in response.headers:
                failed = True
            else:
                remote_timestamp = datetime.strptime(
                    response.headers["Last-Modified"], "%a, %d %b %Y %H:%M:%S GMT"
                ).replace(tzinfo=timezone.utc)
        except Exception as e:
            print("Got exception while trying to get the remote timestamp:", e)
            failed = True
        if failed:
            print(
                f"Failed to get the remote timestamp from <{Simulations.last_modified_url}>.\n"
                + "Assuming it is old."
            )
            return datetime.min.replace(tzinfo=timezone.utc)
        return remote_timestamp
    
    @classmethod
    def local(cls, directory=None, *, download=None):
        """Load the local catalog of SXS simulations
        
        This function loads the standard public catalog, but also
        includes any local simulations found in the given directory.
        If no directory is provided, it will look for the local
        simulations file in the sxs cache directory.

        Parameters
        ----------
        directory : {None, str, Path}, optional
            A directory containing subdirectories of SXS simulations.
            See `sxs.local_simulations` for details about what is
            expected in this directory.  If None (the default), it
            will look for the local simulations file in the sxs cache
            directory.
        download : {None, bool}, optional
            Passed to `Simulations.load` when loading the public set
            of simulations.

        See Also
        --------
        sxs.local_simulations : Search for local simulations
        sxs.write_local_simulations : Write local simulations to a file
        
        """
        import json
        from .local import write_local_simulations
        from .. import sxs_directory

        local_path = sxs_directory("cache") / "local_simulations.json"
        if directory is not None:
            write_local_simulations(directory)
        if not local_path.exists():
            if directory is not None:
                raise ValueError(f"Writing local simulations for {directory=} failed")
            else:
                raise ValueError(
                    f"Local simulations file not found, but no `directory` was provided.\n"
                    + "If called from `sxs.load`, just pass the name of the directory."
                )
        with local_path.open("r") as f:
            local_simulations = json.load(f)
        simulations = cls.load(download)
        doi_versions = {
            k: v["DOI_versions"]
            for k,v in simulations.items()
            if "DOI_versions" in v
        }
        simulations.update(local_simulations)
        for k,v in doi_versions.items():
            simulations[k]["DOI_versions"] = v
        simulations.__file__ = str(local_path)
        return simulations

    @classmethod
    def load(cls, download=None, *, local=False, annex_dir=None):
        """Load the catalog of SXS simulations

        Note that — unlike most SXS data files — the simulations file
        is updated frequently.  As a result, this function — unlike
        the loading functions for most SXS data files — will download
        the simulations by default each time it is called.  However,
        also note that this function is itself cached, meaning that
        the same dict will be returned on each call in a given python
        session.  If you want to avoid that behavior, use
        `Simulations.reload`.

        Parameters
        ----------
        download : {None, bool}, optional
            If False, this function will look for the simulations in
            the sxs cache and raise an error if it is not found.  If
            True, this function will download the simulations and
            raise an error if the download fails.  If None (the
            default), it will try to download the file, warn but fall
            back to the cache if that fails, and only raise an error
            if the simulations is not found in the cache.  Note that
            this ignores the sxs configuration file entirely.

        Keyword-only Parameters
        -----------------------
        local : {None, bool}, optional
            If True, this function will load local simulations from
            the sxs cache.  To prepare the cache, you may wish to call
            `sxs.write_local_simulations`.
        annex_dir : {None, str, Path}, optional
            If provided and `local=True`, this function will load
            local simulations from the given directory.  This is
            equivalent to calling `Simulations.local(directory)`.

        See Also
        --------
        sxs.sxs_directory : Locate cache directory Simulations.reload
        : Avoid caching the result of this function

        """
        from datetime import datetime, timezone
        import json
        import zipfile
        from .. import sxs_directory, read_config
        from ..utilities import download_file

        if hasattr(cls, "_simulations"):
            return cls._simulations

        if local or annex_dir is not None:
            cls._simulations = cls.local(annex_dir, download=download)
            return cls._simulations

        progress = read_config("download_progress", True)

        remote_timestamp = cls.remote_timestamp(download is not False)  # Test for literal `False`

        cache_path = sxs_directory("cache") / "simulations.zip"

        if cache_path.exists():
            local_timestamp = datetime.fromtimestamp(cache_path.stat().st_mtime, timezone.utc)
        elif download is False:
            raise ValueError(f"Simulations not found in '{cache_path}' and downloading was turned off")
        else:
            local_timestamp = datetime.min.replace(tzinfo=timezone.utc)

        download_failed = False
        if (download or download is None) and remote_timestamp > local_timestamp:
            # 1. Download the full json file (zipped in flight, but auto-decompressed on arrival)
            # 2. Zip to a temporary file (using bzip2, which is better than the in-flight compression)
            # 3. Replace the original simulations.zip with the temporary zip file
            # 4. Remove the full json file
            # 5. Make sure the temporary zip file is gone too
            temp_json = cache_path.with_suffix(".temp.json")
            temp_zip = cache_path.with_suffix(".temp.zip")
            try:
                try:
                    download_file(cls.url, temp_json, progress=progress, if_newer=False)
                except Exception as e:
                    if download:
                        raise RuntimeError(f"Failed to download '{cls.url}'; try setting `download=False`") from e
                    download_failed = e  # We'll try the cache
                else:
                    if temp_json.exists():
                        with zipfile.ZipFile(temp_zip, "w", compression=zipfile.ZIP_BZIP2) as simulations_zip:
                            simulations_zip.write(temp_json, arcname="simulations.json")
                        temp_zip.replace(cache_path)
            finally:
                temp_json.unlink(missing_ok=True)
                temp_zip.unlink(missing_ok=True)

        if not cache_path.exists():
            if download is False:  # Test if it literally *is* False, rather than just casts to False
                raise ValueError(f"The simulations file was not found in '{cache_path}', and downloading was turned off")
            elif download_failed:
                raise ValueError(f"Simulations not found in '{cache_path}' and download failed") from download_failed
            else:
                raise ValueError(f"Simulations not found in '{cache_path}' for unknown reasons")

        try:
            with zipfile.ZipFile(cache_path, "r") as simulations_zip:
                try:
                    with simulations_zip.open("simulations.json") as simulations_json:
                        try:
                            simulations = json.load(simulations_json)
                        except Exception as e:
                            raise ValueError(f"Failed to parse 'simulations.json' in '{cache_path}'") from e
                except Exception as e:
                    raise ValueError(f"Failed to open 'simulations.json' in '{cache_path}'") from e
        except Exception as e:
            raise ValueError(f"Failed to open '{cache_path}' as a ZIP file") from e

        sims = cls(simulations)
        sims.__file__ = str(cache_path)

        cls._simulations = sims
        return sims

    @classmethod
    def reload(cls, download=True):
        """Reload the catalog of SXS simulations, without caching

        Clears the cache of `Simulations.load` and returns the result of calling it again.
        Note that in this function, the default value of `download` is `True`, rather
        than `None` as in `Simulations.load` — though both behaviors are available.

        Parameters
        ----------
        download : {None, bool}, optional
            If False, this function will look for the simulations in the sxs cache and
            raise an error if it is not found.  If True (the default), this function
            will download the simulations and raise an error if the download fails.  If
            None (the default), it will try to download the file, warn but fall back to
            the cache if that fails, and only raise an error if the simulations is not
            found in the cache.  Note that this ignores the sxs configuration file
            entirely.

        See Also
        --------
        sxs.sxs_directory : Locate cache directory
        Simulations.load : Caching version of this function

        """
        cls.load.cache_clear()
        return cls.load(download=download)

    @property
    def dataframe(self):
        """Create pandas.DataFrame containing metadata for all
        simulations

        Note that `pandas` is the standard Python interface for
        heterogeneous data tables, like the one we have here.  This
        interface allows for more convenient slicing and querying of
        data than the list of `dict`s provided by the `Simulations`
        object.
        
        This can also be a more convenient way to access the metadata
        because the raw metadata has missing keys and mixed formats.
        Iif a key is missing from the metadata for a particular key,
        the dataframe will just have a `NaN` in that entry, rather
        than raising an exception.  Other keys may have unexpected
        entries — such as the `"reference_eccentricity"` field, which
        is *usually* a float but may be a string like "<0.0001" if the
        eccentricity is not known precisely, but is only bounded.  The
        dataframe introduces a new column called
        `"reference_eccentricity_bound"` that is always a float giving
        an upper bound on the eccentricity.

        See the `pandas` documentation for more information on how to
        use the resulting dataframe, or the `Simulations` tutorial for
        examples.
        
        """
        import numpy as np
        import pandas as pd

        if hasattr(self, "_dataframe"):
            return self._dataframe

        simulations = pd.DataFrame.from_dict(self, orient="index")

        def floater(x):
            try:
                f = float(x)
            except:
                f = np.nan
            return f

        def floaterbound(x):
            try:
                f = float(x)
            except:
                try:
                    f = float(x.replace("<", ""))
                except:
                    f = np.nan
            return f

        def norm(x):
            try:
                n = np.linalg.norm(x)
            except:
                n = np.nan
            return n

        def three_vec(x):
            try:
                a = np.array(x, dtype=float)
                if a.shape != (3,):
                    raise ValueError("Don't understand input as a three-vector")
            except:
                a = np.array([np.nan, np.nan, np.nan])
            return a

        def datetime_from_string(x):
            try:
                dt = pd.to_datetime(x).tz_convert("UTC")
            except:
                dt = pd.to_datetime("1970-1-1").tz_localize("UTC")
            return dt

        sims_df = SimulationsDataFrame(pd.concat((
            simulations["reference_time"].map(floater),
            simulations["reference_mass_ratio"].map(floater),
            simulations["reference_dimensionless_spin1"].map(three_vec),
            simulations["reference_dimensionless_spin1"].map(norm).rename("reference_dimensionless_spin1_mag"),
            simulations["reference_dimensionless_spin2"].map(three_vec),
            simulations["reference_dimensionless_spin2"].map(norm).rename("reference_dimensionless_spin2_mag"),
            simulations["reference_chi_eff"].map(floater),
            simulations["reference_chi1_perp"].map(floater),
            simulations["reference_chi2_perp"].map(floater),
            simulations["reference_eccentricity"].map(floater),
            simulations["reference_eccentricity"].map(floaterbound).rename("reference_eccentricity_bound"),
            simulations["reference_mean_anomaly"].map(floater),
            simulations["reference_orbital_frequency"].map(three_vec),
            simulations["reference_orbital_frequency"].map(norm).rename("reference_orbital_frequency_mag"),
            (
                simulations["reference_position1"].map(three_vec)
                -simulations["reference_position2"].map(three_vec)
            ).map(norm).rename("reference_separation"),
            simulations["reference_position1"].map(three_vec),
            simulations["reference_position2"].map(three_vec),
            simulations["reference_mass1"].map(floater),
            simulations["reference_mass2"].map(floater),
            simulations["reference_dimensionless_spin1"].map(norm).rename("reference_chi1_mag"),
            simulations["reference_dimensionless_spin2"].map(norm).rename("reference_chi2_mag"),
            simulations["relaxation_time"].map(floater),
            #simulations["merger_time"].map(floater),
            simulations["common_horizon_time"].map(floater),
            simulations["remnant_mass"].map(floater),
            simulations["remnant_dimensionless_spin"].map(three_vec),
            simulations["remnant_dimensionless_spin"].map(norm).rename("remnant_dimensionless_spin_mag"),
            simulations["remnant_velocity"].map(three_vec),
            simulations["remnant_velocity"].map(norm).rename("remnant_velocity_mag"),
            #simulations["final_time"].map(floater),
            simulations["EOS"].fillna(simulations["eos"]),
            simulations["disk_mass"].map(floater),
            simulations["ejecta_mass"].map(floater),
            simulations["object_types"].astype("category"),
            simulations["initial_data_type"].astype("category"),
            simulations["initial_separation"].map(floater),
            simulations["initial_orbital_frequency"].map(floater),
            simulations["initial_adot"].map(floater),
            simulations["initial_ADM_energy"].map(floater),
            simulations["initial_ADM_linear_momentum"].map(three_vec),
            simulations["initial_ADM_linear_momentum"].map(norm).rename("initial_ADM_linear_momentum_mag"),
            simulations["initial_ADM_angular_momentum"].map(three_vec),
            simulations["initial_ADM_angular_momentum"].map(norm).rename("initial_ADM_angular_momentum_mag"),
            simulations["initial_mass1"].map(floater),
            simulations["initial_mass2"].map(floater),
            simulations["initial_mass_ratio"].map(floater),
            simulations["initial_dimensionless_spin1"].map(three_vec),
            simulations["initial_dimensionless_spin1"].map(norm).rename("initial_dimensionless_spin1_mag"),
            simulations["initial_dimensionless_spin2"].map(three_vec),
            simulations["initial_dimensionless_spin2"].map(norm).rename("initial_dimensionless_spin2_mag"),
            simulations["initial_position1"].map(three_vec),
            simulations["initial_position2"].map(three_vec),
            #simulations["object1"].astype("category"),
            #simulations["object2"].astype("category"),
            # simulations["url"],
            #simulations["simulation_name"],
            #simulations["alternative_names"],
            # simulations["metadata_path"],
            # simulations["end_of_trajectory_time"].map(floater),
            # simulations["merger_time"].map(floater),
            simulations["number_of_orbits"].map(floater),
            simulations["superseded_by"],
            simulations["DOI_versions"],
            simulations["keywords"],
            simulations["date_link_earliest"].map(datetime_from_string),
            simulations["date_run_earliest"].map(datetime_from_string),
            simulations["date_run_latest"].map(datetime_from_string),
            simulations["date_postprocessing"].map(datetime_from_string),
        ), axis=1))

        sims_df.insert(0, "deprecated", (
            ~sims_df.superseded_by.isna()
            | sims_df["keywords"].map(lambda ks: "deprecated" in ks)
        ))

        # We have ignored the following fields present in the
        # simulations.json file (as of 2024-08-04), listed here with
        # the number of non-null entries:
        #
        # alternative_names                2778
        # point_of_contact_email           2778
        # authors_emails                   2776
        # simulation_bibtex_keys           2778
        # code_bibtex_keys                 2778
        # initial_data_bibtex_keys         2778
        # quasicircular_bibtex_keys        2778
        # metadata_version                 2778
        # spec_revisions                   2778
        # spells_revision                  2778
        # merger_time                         9
        # final_time                         12
        # reference_spin1                     2
        # reference_spin2                     1
        # nitial_spin1                        2
        # initial_spin2                       2
        # remnant_spin                        2
        # initial_mass_withspin2              2
        # end_of_trajectory_time              3

        self._dataframe = sims_df
        return sims_df

    table = dataframe
