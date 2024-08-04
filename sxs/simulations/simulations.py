"""Container interface to the catalog of SXS simulations"""

import functools
import collections



class Simulations(collections.OrderedDict):
    """Interface to the catalog of SXS simulations"""
    url = "https://gist.githubusercontent.com/moble/d3fa38db9f1257c76e61006f5e182884/raw/f2b0939c2520ec3377c9f6fb5ee5236d938430f5/simulations.json"

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
    @functools.lru_cache()
    def load(cls, download=None):
        """Load the catalog of SXS simulations

        Note that — unlike most SXS data files — the simulations file is updated
        frequently.  As a result, this function — unlike the loading functions for most
        SXS data files — will download the simulations by default each time it is called.
        However, also note that this function is itself cached, meaning that the same
        dict will be returned on each call in a given python session.  If you want to
        avoid that behavior, use `Simulations.reload`.

        Parameters
        ----------
        download : {None, bool}, optional
            If False, this function will look for the simulations in the sxs cache and
            raise an error if it is not found.  If True, this function will download
            the simulations and raise an error if the download fails.  If None (the
            default), it will try to download the file, warn but fall back to the cache
            if that fails, and only raise an error if the simulations is not found in the
            cache.  Note that this ignores the sxs configuration file entirely.

        See Also
        --------
        sxs.sxs_directory : Locate cache directory
        Simulations.reload : Avoid caching the result of this function

        """
        import json
        import zipfile
        from .. import sxs_directory, read_config
        from ..utilities import download_file

        progress = read_config("download_progress", True)

        cache_path = sxs_directory("cache") / "simulations.zip"

        if cache_path.exists():
            if_newer = cache_path
        else:
            if_newer = False

        if download or download is None:
            # 1. Download the full json file (zipped in flight, but auto-decompressed on arrival)
            # 2. Zip to a temporary file (using bzip2, which is better than the in-flight compression)
            # 3. Replace the original simulations.zip with the temporary zip file
            # 4. Remove the full json file
            # 5. Make sure the temporary zip file is gone too
            temp_json = cache_path.with_suffix(".temp.json")
            temp_zip = cache_path.with_suffix(".temp.zip")
            try:
                try:
                    download_file(cls.url, temp_json, progress=progress, if_newer=if_newer)
                except Exception as e:
                    if download:
                        raise RuntimeError(f"Failed to download '{cls.url}'; try setting `download=False`") from e
                    download_failed = e  # We'll try the cache
                else:
                    download_failed = False
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

        return cls(simulations)

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
    @functools.lru_cache()
    def dataframe(self):
        """Return a pandas.DataFrame containing the metadata for all simulations"""
        import numpy as np
        import pandas as pd
        simulations = pd.DataFrame.from_dict(self._dict, orient="index")

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

        def space_translation(x):
            try:
                a = np.array(x["space_translation"])
            except:
                a = np.array([np.nan, np.nan, np.nan])
            return a

        def space_translation_norm(x):
            try:
                n = np.linalg.norm(np.array(x["space_translation"]))
            except:
                n = np.nan
            return n

        def boost_velocity(x):
            try:
                a = np.array(x["boost_velocity"])
            except:
                a = np.array([np.nan, np.nan, np.nan])
            return a

        def boost_velocity_norm(x):
            try:
                n = np.linalg.norm(np.array(x["boost_velocity"]))
            except:
                n = np.nan
            return n

        sims_df = pd.concat((
            simulations["object_types"].astype("category"),
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
            simulations["com_parameters"].map(space_translation).rename("com_correction_space_translation"),
            simulations["com_parameters"].map(space_translation_norm).rename("com_correction_space_translation_mag"),
            simulations["com_parameters"].map(boost_velocity).rename("com_correction_boost_velocity"),
            simulations["com_parameters"].map(boost_velocity_norm).map(norm).rename("com_correction_boost_velocity_mag"),
            simulations["reference_time"].map(floater),
            (
                simulations["reference_position1"].map(three_vec)
                -simulations["reference_position2"].map(three_vec)
            ).map(norm).rename("reference_separation"),
            simulations["reference_orbital_frequency"].map(norm).rename("reference_orbital_frequency_mag"),
            simulations["reference_mass_ratio"].map(floater),
            simulations["reference_dimensionless_spin1"].map(norm).rename("reference_chi1_mag"),
            simulations["reference_dimensionless_spin2"].map(norm).rename("reference_chi2_mag"),
            simulations["reference_chi_eff"].map(floater),
            simulations["reference_chi1_perp"].map(floater),
            simulations["reference_chi2_perp"].map(floater),
            simulations["reference_eccentricity"].map(floater),
            simulations["reference_eccentricity"].map(floaterbound).rename("reference_eccentricity_bound"),
            simulations["reference_mean_anomaly"].map(floater),
            simulations["reference_mass1"].map(floater),
            simulations["reference_mass2"].map(floater),
            simulations["reference_dimensionless_spin1"].map(three_vec),
            simulations["reference_dimensionless_spin1"].map(norm).rename("reference_dimensionless_spin1_mag"),
            simulations["reference_dimensionless_spin2"].map(three_vec),
            simulations["reference_dimensionless_spin2"].map(norm).rename("reference_dimensionless_spin2_mag"),
            simulations["reference_orbital_frequency"].map(three_vec),
            simulations["reference_position1"].map(three_vec),
            simulations["reference_position2"].map(three_vec),
            simulations["relaxation_time"].map(floater),
            #simulations["merger_time"].map(floater),
            simulations["common_horizon_time"].map(floater),
            simulations["remnant_mass"].map(floater),
            simulations["remnant_dimensionless_spin"].map(three_vec),
            simulations["remnant_dimensionless_spin"].map(norm).rename("remnant_dimensionless_spin_mag"),
            simulations["remnant_velocity"].map(three_vec),
            simulations["remnant_velocity"].map(norm).rename("remnant_velocity_mag"),
            #simulations["final_time"].map(floater),
            simulations["eos"],
            simulations["initial_data_type"].astype("category"),
            #simulations["object1"].astype("category"),
            #simulations["object2"].astype("category"),
            simulations["disk_mass"].map(floater),
            simulations["ejecta_mass"].map(floater),
            simulations["url"],
            #simulations["simulation_name"],
            #simulations["alternative_names"],
            simulations["metadata_path"],
        ), axis=1)

        #  57  reference_spin1                2 non-null      object
        #  58  reference_spin2                1 non-null      object
        #  59  nitial_spin1                   2 non-null      object
        #  60  initial_spin2                  2 non-null      object
        #  61  remnant_spin                   2 non-null      object
        #  62  initial_mass_withspin2         2 non-null      float64

        return sims_df

    table = dataframe
