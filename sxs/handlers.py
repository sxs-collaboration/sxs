"""Functions to facilitate generic handling of SXS-format data files"""

import contextlib


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


def sxs_loader(file):
    """Find the function that will load the given file

    If a format is specified, we assume that it is a top-level attribute or member
    named "sxs_format" or just "format".  If neither exists, return None; the
    calling function should check for this possibility.

    Parameters
    ----------
    file : file-like object, string, or pathlib.Path

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
    format_string = file_format(file)
    if format_string is None:
        file_string = str(pathlib.Path(file).name).lower()
        if "catalog" in file_string:
            format_string = "catalog"
        elif "metadata" in file_string:
            format_string = "metadata"
        elif "horizons" in file_string:
            format_string = "horizons"
        elif re.match("(rh_|rhoverm_|rpsi4_|rmpsi4_)", file_string):
            format_string = "waveforms"
        else:
            raise ValueError(f"File '{file}' contains no recognized format information")
    handler = sxs_handler(format_string)
    # noinspection PyUnresolvedReferences
    return handler.load


def load(location, download=None, cache=None, progress=None, **kwargs):
    """Load an SXS-format dataset, optionally downloading and caching

    The dataset can be the full catalog of all SXS simulations, or metadata,
    horizon data, or a waveform from an individual simulation.

    Parameters
    ----------
    location : {str, pathlib.Path}
        A local file path, URL, SXS path, or SXS path pattern.  See Notes below.
    download : {None, bool}, optional
        If this is True and the data is recognized as starting with an SXS ID but
        cannot be found in the cache, the data will be downloaded automatically.
        If this is None (the default) and an SXS configuration file is found with a
        `download` key, that value will be used.  If this is False, any
        configuration will be ignored, and no files will be downloaded.  Note that
        if this is True but `cache` is None, `cache` will automatically be switched
        to True.
    cache : {None, bool}, optional
        The cache directory is determined by `sxs.sxs_directory`, and any downloads
        will be stored in that directory.  If this is None (the default) and
        `download` is True it will be set to True.  If this is False, any
        configuration will be ignored and any files will be downloaded to a
        temporary directory that will be deleted when python exits.
    progress : {None, bool}, optional
        If True, full file names will be shown and, if a nonzero Content-Length
        header is returned, a progress bar will be shown during any downloads.
        Default is None, which just reads the configuration value with
        `read_config("download_progress", True)`, defaulting to True.

    Keyword Parameters
    ------------------
    All remaining parameters are passed to the `load` function responsible for the
    requested data.

    See Also
    --------
    sxs.sxs_directory : Locate configuration and cache files
    sxs.write_config : Set defaults for `download` and `cache` parameters

    Notes
    -----
    This function can load data in various ways.

      1) Given an absolute or relative path to a local file, it just loads the data
         directly.

      2) If `location` is a valid URL including the scheme (https://, or http://),
         it will be downloaded regardless of the `download` parameter and
         optionally cached.

      3) Given an SXS path — like 'SXS:BBH:1234/Lev5/h_Extrapolated_N2.h5' — the
         file is located in the catalog for details.  This function then looks in
         the local cache directory and loads it if present.

      4) If the SXS path is not found in the cache directory and `download` is set
         to `True` (when this function is called, or in the sxs config file) this
         function attempts to download the data.  Note that `download` must be
         explicitly set in this case, or a ValueError will be raised.

    Note that downloading is switched off by default, but if it is switched on (set
    to True), the cache is also switched on by default.

    """
    import pathlib
    import urllib.request
    from . import Catalog, read_config, sxs_directory
    from .utilities import url, download_file, sxs_path_to_system_path

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

    path = pathlib.Path(sxs_path_to_system_path(location)).expanduser().resolve()
    h5_path = path.with_suffix('.h5')
    json_path = path.with_suffix('.json')

    if not path.exists():
        if h5_path.exists():
            path = h5_path

        elif json_path.exists():
            path = json_path

        elif "scheme" in url.parse(location):
            m = url.parse(location)
            path_name = urllib.request.url2pathname(f"{m['host']}/{m['port']}/{m['resource']}")
            path = cache_path / path_name
            if not path.exists():
                if download is False:  # Again, we want literal False, not casting to False
                    raise ValueError(f"File '{path_name}' not found in cache, but downloading turned off")
                download_file(location, path, progress=progress)

        elif location == "catalog":
            return Catalog.load(download=download)

        else:
            # Try to find an appropriate SXS file
            catalog = Catalog.load(download=download)
            selections = catalog.select_files(location)
            if not selections:
                raise ValueError(f"Nothing found matching '{location}'")
            if progress:
                print("Found the following files to load from the SXS catalog:")
                print("    " + "\n    ".join(selections))
            paths = []
            for sxs_path, file_info in selections.items():
                truepath = sxs_path_to_system_path(file_info.get("truepath", sxs_path))
                path = cache_path / truepath
                if not path.exists():
                    download_url = file_info["download"]
                    download_file(download_url, path, progress=progress)
                paths.append(path)
            loaded = [load(path, download=False, progress=progress, **kwargs) for path in paths]
            if len(loaded) == 1:
                return loaded[0]
            else:
                return loaded

    loader = sxs_loader(path)

    return loader(path, **kwargs)


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
