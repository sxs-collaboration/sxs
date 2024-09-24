"""Functions to facilitate generic handling of SXS-format data files"""

import contextlib
from . import waveforms, doi_url


class JSONHandler:
    """Utility for loading and saving JSON files"""
    @classmethod
    def load(cls, file, **kwargs):
        import json
        if hasattr(file, "read"):
            return json.load(file, **kwargs)
        with open(file, "r") as f:
            return json.load(f, **kwargs)
    @classmethod
    def save(cls, obj, file, **kwargs):
        import json
        if hasattr(file, "write"):
            return json.dump(obj, file, **kwargs)
        with open(file, "w") as f:
            return json.dump(obj, f, **kwargs)


def sxs_handler(format_string):
    """Find an object to load from or save to files in the given format.

    Parameters
    ----------
    format_string : str

    Returns
    -------
    handler : object
        This object will have (at least) two attributes: `load` and `save`, which
        can be called as `handler.load(file, **kwargs)` and `handler.save(obj,
        file, **kwargs)`.

    See Also
    --------
    sxs_loader : Returns the function that will load a given file
    sxs.utilities.file_format : Returns just the string found in the file

    """
    import itertools
    import re
    from . import catalog, metadata, horizons, waveforms

    if not format_string:
        raise ValueError("Empty string cannot be associated with a handler")
    elif format_string.lower().startswith("catalog"):
        format_string = re.sub(r"^catalog\.?", "", format_string, count=1, flags=re.IGNORECASE)
        return catalog.formats.get(format_string, catalog.formats[None])
    elif format_string.lower().startswith("metadata"):
        format_string = re.sub(r"^metadata\.?", "", format_string, count=1, flags=re.IGNORECASE)
        return metadata.formats.get(format_string, metadata.formats[None])
    elif format_string.lower().startswith("horizons"):
        format_string = re.sub(r"^horizons\.?", "", format_string, count=1, flags=re.IGNORECASE)
        return horizons.formats.get(format_string, horizons.formats[None])
    elif format_string.lower().startswith("waveforms"):
        format_string = re.sub(r"^waveforms\.?", "", format_string, count=1, flags=re.IGNORECASE)
        return waveforms.formats.get(format_string, waveforms.formats[None])
    elif format_string.lower() == "json":
        return JSONHandler
    else:
        format_list = [
            catalog.formats,
            metadata.formats,
            horizons.formats,
            waveforms.formats,
        ]
        format_cycler = itertools.cycle(format_list)
        for _ in range(len(format_list)):
            format_dict = next(format_cycler)
            if format_string in format_dict:
                if any(format_string in next(format_cycler) for _ in range(len(format_list)-1)):
                    raise ValueError(f"Format string '{format_string}' found in multiple sxs format groups")
                return format_dict[format_string]
    raise ValueError(f"Format '{format_string}' is unknown to the `sxs` package; maybe you need to update `sxs`")


def sxs_loader(file, group=None):
    """Find the function that will load the given file

    If a format is specified, we assume that it is a top-level attribute or member
    named "sxs_format" or just "format".  If neither exists, return None; the
    calling function should check for this possibility.

    Parameters
    ----------
    file : file-like object, string, or pathlib.Path
    group : str, optional
        If the file is an HDF5 or JSON file, this is the group to inspect.

    Returns
    -------
    load : callable
        This function can be called as `load(file, **kwargs)`.

    See Also
    --------
    sxs_handler : Returns an object to load from and save to a given format
    sxs.utilities.file_format : Returns just the string found in the file

    Notes
    -----
    Note that this function has a very different signature than the related
    functions `sxs_handler`, which takes only the format string, not an actual
    file.

    """
    import re
    import pathlib
    from .utilities import file_format
    format_string = file_format(file, group)
    if format_string is None:
        file_string = str(pathlib.Path(file).name).lower()
        if "catalog" in file_string:
            format_string = "catalog"
        elif file_string.startswith("metadata"):
            format_string = "metadata"
        elif "horizons" in file_string:
            format_string = "horizons"
        elif re.match("(rh_|rhoverm_|rpsi4_|rmpsi4_)", file_string):
            format_string = "waveforms"
        elif file_string.endswith("json"):
            format_string = "json"
        else:
            raise ValueError(f"File '{file}' contains no recognized format information")
    handler = sxs_handler(format_string)
    # noinspection PyUnresolvedReferences
    return handler.load


def _safe_resolve_exists(path):
    """Evaluate `path.resolve().exists()`  without throwing exception
    
    This is just here to work around a bug that turned up in Windows
    on python 3.8.  It's not clear if it turns up in other versions.
    """
    try:
        return path.resolve().exists()
    except:
        return False

def load(location, download=None, cache=None, progress=None, truepath=None, **kwargs):
    """Load an SXS-format dataset, optionally downloading and caching

    The dataset can be the full catalog of all SXS simulations, or
    metadata, horizon data, or a waveform from an individual
    simulation.

    Parameters
    ----------
    location : {str, pathlib.Path}
        A local file path, URL, SXS path, or SXS path pattern.  See
        Notes below.
    download : {None, bool}, optional
        If this is True and the data is recognized as starting with an
        SXS ID but cannot be found in the cache, the data will be
        downloaded automatically.  If this is None (the default) and
        an SXS configuration file is found with a `download` key, that
        value will be used.  If this is False, any configuration will
        be ignored, and no files will be downloaded.  Note that if
        this is True but `cache` is None, `cache` will automatically
        be switched to True.
    cache : {None, bool}, optional
        The cache directory is determined by `sxs.sxs_directory`, and
        any downloads will be stored in that directory.  If this is
        None (the default) and `download` is True it will be set to
        True.  If this is False, any configuration will be ignored and
        any files will be downloaded to a temporary directory that
        will be deleted when python exits.
    progress : {None, bool}, optional
        If True, full file names will be shown and, if a nonzero
        Content-Length header is returned, a progress bar will be
        shown during any downloads.  Default is None, which just reads
        the configuration value with `read_config("download_progress",
        True)`, defaulting to True.
    truepath : {None, str}, optional
        If the file is downloaded, this allows the output path to be
        overridden, rather than selected automatically.  The output
        path will be stored in `truepath` relative to the cache
        directory.

    Keyword Parameters
    ------------------
    All remaining parameters are passed to the `load` function
    responsible for the requested data.

    See Also
    --------
    sxs.sxs_directory : Locate configuration and cache files
    sxs.write_config : Set defaults for `download` and `cache`
        parameters

    Notes
    -----
    This function can load data in various ways.

      1) Given an absolute or relative path to a local file, it just
         loads the data directly.

      2) If `truepath` is set, and points to a file that exists —
         whether absolute, relative to the current working directory,
         or relative to the cache directory — that file will be
         loaded.

      3) If `location` is a valid URL including the scheme (https://,
         or http://), it will be downloaded regardless of the
         `download` parameter and optionally cached.

      4) Given an SXS simulation specification — like "SXS:BBH:1234",
         "SXS:BBH:1234v2.0", "SXS:BBH:1234/Lev5", or
         "SXS:BBH:1234v2.0/Lev5" — the simulation is loaded as an
         `sxs.Simulation` object.

      5) Given an SXS path — like
         "SXS:BBH:1234/Lev5/h_Extrapolated_N2.h5" — the file is
         located in the catalog for details.  This function then looks
         in the local cache directory and loads it if present.

      6) If the SXS path is not found in the cache directory and
         `download` is set to `True` (when this function is called, or
         in the sxs config file) this function attempts to download
         the data.  Note that `download` must be explicitly set in
         this case, or a ValueError will be raised.

    If the file is downloaded, it will be stored in the cache
    according to the `location`, unless `truepath` is set as noted
    above, in which case it is stored there.  Note that downloading is
    switched off by default, but if it is switched on (set to True),
    the cache is also switched on by default.

    """
    import pathlib
    import urllib.request
    from . import Simulations, Simulation, read_config, sxs_directory, Catalog
    from .utilities import url, download_file, sxs_path_to_system_path, sxs_id_version_lev_exact_re

    # Note: `download` and/or `cache` may still be `None` after this
    if download is None:
        download = read_config("download", True)
    if cache is None:
        cache = read_config("cache")
    if progress is None:
        progress = read_config("download_progress", True)

    # We set the cache path to be persistent if `cache` is `True` or `None`.  Thus,
    # we test for whether or not `cache` literally *is* `False`, rather than just
    # if it casts to `False`.
    cache_path = sxs_directory("cache", persistent=(cache is not False))

    path = pathlib.Path(sxs_path_to_system_path(location)).expanduser()  # .resolve()
    h5_path = path.with_suffix('.h5')
    json_path = path.with_suffix('.json')

    if not path.exists():
        if truepath and (testpath := pathlib.Path(sxs_path_to_system_path(truepath)).expanduser()).exists():
            path = testpath

        elif truepath and (testpath := cache_path / sxs_path_to_system_path(truepath)).exists():
            path = testpath

        elif _safe_resolve_exists(h5_path):
            path = h5_path

        elif _safe_resolve_exists(json_path):
            path = json_path

        elif "scheme" in url.parse(location):
            m = url.parse(location)
            truepath = truepath or urllib.request.url2pathname(f"{m['host']}/{m['port']}/{m['resource']}")
            path = cache_path / sxs_path_to_system_path(truepath)
            if not path.resolve().exists():
                if download is False:  # Again, we want literal False, not casting to False
                    raise ValueError(f"File '{truepath}' not found in cache, but downloading turned off")
                download_file(location, path, progress=progress)

        elif location == "catalog":
            return Catalog.load(download=download)

        elif location == "simulations":
            return Simulations.load(download=download)

        elif sxs_id_version_lev_exact_re.match(location):
            return Simulation(location, download=download, cache=cache, progress=progress, **kwargs)

        else:
            # Try to find an appropriate SXS file in the catalog
            catalog = Catalog.load(download=download)
            selections = catalog.select_files(location)
            if not selections:
                raise ValueError(f"Nothing found matching '{location}'")
            if progress:
                print("Found the following files to load from the SXS catalog:")
                print("    " + "\n    ".join(selections))
            paths = []
            for sxs_path, file_info in selections.items():
                truepath = truepath or sxs_path_to_system_path(file_info.get("truepath", sxs_path))
                path = cache_path / sxs_path_to_system_path(truepath)
                if not path.resolve().exists():
                    download_url = file_info["download"]
                    download_file(download_url, path, progress=progress)
                paths.append(path)
            loaded = [load(path, download=False, progress=progress, **kwargs) for path in paths]
            if len(loaded) == 1:
                return loaded[0]
            else:
                return loaded

    loader = sxs_loader(path, kwargs.get("group", None))

    loaded = loader(path, **kwargs)
    try:
        loaded.__file__ = str(path)
    except:
        pass
    return loaded


def load_via_sxs_id(sxsid, location, *, download=None, cache=None, progress=None, truepath=None, **kwargs):
    """Load a path via a (possibly versioned) SXS ID
    
    Given some SXS ID like "SXS:BBH:1234" or a versioned ID like
    "SXS:BBH:1234v2.0", we may wish to first resolve the DOI to a
    specific Zenodo record, and then load a specific path under that
    record — for example, we may want to export the Zenodo record as
    JSON by appending "export/json" to the Zenodo URL, which we
    *could* get as

        load("https://zenodo.org/records/13152488/export/json")

    because that's the record "SXS:BBH:1234v2.0" resolves to.
    However, we don't want to keep track of the Zenodo URL, so we
    just use this function instead, as

        load_via_sxs_id("SXS:BBH:1234v2.0", "export/json")
    
    """
    from pathlib import Path
    import requests
    from .utilities import sxs_path_to_system_path
    url = f"{doi_url}{sxsid}"
    response = requests.head(url, allow_redirects=True)
    if response.status_code != 200:
        raise ValueError(f"Could not load via DOI {url=}")
    final_url = f"{response.url}/{location}"
    truepath = truepath or Path(sxs_path_to_system_path(sxsid)) / location
    return load(final_url, download, cache, progress, truepath, **kwargs)


@contextlib.contextmanager
def loadcontext(*args, **kwargs):
    """Context manager for backwards compatibility

    This context manager takes precisely the same arguments as `sxs.load` and
    yields precisely the same results; essentially it is a trivial wrapper around
    the `load` function.  The benefit of this approach is that it can be used in
    precisely the same way as `h5py.File` would have been used previously.  For
    example, the old approach would be to open an HDF5 file like this:

        with h5py.File("Horizons.h5", "r") as horizons:
            areal_mass = horizons["AhA.dir/ArealMass.dat"]

    With this function, the same effect can be achieved as

        with sxs.loadcontext("Horizons.h5") as horizons:
            areal_mass = horizons["AhA.dir/ArealMass.dat"]

    Each of the datasets found in Horizons.h5 as well as the old NRAR-style files
    will be available through this interface, even when using newer files in
    different formats.  Thus, only one line of code would need to change to use the
    new interface.

    However, be aware that this may not an be efficient use of memory, and is
    almost certainly slower than the newer interfaces.  Wherever possible, you
    should update your code to use newer interfaces.  Failing to do so will leave
    you open to ridicule from your peers and loved ones.

    See Also
    --------
    load

    """
    yield load(*args, **kwargs)


def load_lvc(
        sxs_id, *,
        t_ref=None, f_ref=None,
        dt=None,
        f_low=None,
        ell_max=None,
        phi_ref=None, inclination=None,
        ell_max_epoch=None,
        **kwargs
):
    r"""Load an SXS waveform in LVC convention.

    This is a deprecated function that is a thin wrapper around the
    method `sxs.WaveformModes.to_lvk`.  It is recommended that you use
    that method directly.
    """
    lev = kwargs.pop("lev", "Lev")  # If number is unspecified, load chooses highest
    waveform_name = kwargs.pop("waveform_name", "Strain_N2.h5")
    h = load(f"{sxs_id}/{lev}/{waveform_name}", transform_to_inertial=False)
    horizons = load(f"{sxs_id}/{lev}/Horizons.h5")
    return waveforms.to_lvc_conventions(
        h, horizons,
        t_ref=t_ref, f_ref=f_ref,
        dt=dt,
        f_low=f_low,
        ell_max=ell_max,
        phi_ref=phi_ref, inclination=inclination,
        ell_max_epoch=ell_max_epoch,
        **kwargs
    )
