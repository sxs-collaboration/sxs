"""Functions to find user-specific config and cache directories"""

import re
import platform
import functools
from .sxs_identifiers import sxs_path_re

_platform_system = platform.system()


def read_config(key=None, default=None):
    """Read variable from SXS configuration file

    The configuration file is named `config.json` and is stored in the directory
    returned by `sxs_directory("config")`.

    Parameters
    ----------
    key : {None, str}
        If None (the default), the entire configuration file is returned as a
        dictionary.  Otherwise, this is used as a key for this dictionary.
    default : Any
        An arbitrary object that is returned if the `key` is not found in the
        dictionary.

    Returns
    -------
    config : {dict, Any}
        If `key` is None (the default), the entire configuration dictionary is
        returned.  Otherwise, the value indexed by that key in the configuration
        file is returned.

    """
    import json
    config_path = sxs_directory("config") / "config.json"
    if config_path.exists():
        config = json.load(config_path.open("r"))
    else:
        config = {}
    if key is None:
        return config
    else:
        return config.get(key, default)


def write_config(**kwargs):
    """Write variables to SXS configuration file

    The configuration file is named `config.json` and is stored in the directory
    returned by `sxs_directory("config")`.

    All parameters must be passed as keyword arguments, as in

        write_config(key=value)

    which are then inserted as `key:value` pairs into the config dictionary and
    written into the `config.json` file.

    Useful settings include

      * `write_config(download=True)`, to download data whenever necessary
      * `write_config(cache=True)`, to ensure that caching is used for downloads
      * `write_config(cache_directory="/some/path/you/choose")`, to manually choose
        the directory to use for caching

    """
    import json
    config_path = sxs_directory("config") / "config.json"
    if config_path.exists():
        config = json.load(config_path.open("r"))
    else:
        config = {}
    config.update(**kwargs)
    with config_path.open("w") as c:
        json.dump(config, c, indent=4, separators=(',', ': '))


@functools.lru_cache()
def sxs_directory(directory_type, persistent=True):
    # noinspection SpellCheckingInspection
    """Return the SXS directory location, creating it if necessary

    Parameters
    ----------
    directory_type : {"cache", "config"}
        The type of user directory to be found.
    persistent : bool
        If True (the default) try to return a persistent directory; otherwise,
        return a temporary directory that will be deleted when the calling python
        process exits (via the standard-library function `atexit.register`).

    Returns
    -------
    directory : pathlib.Path
        The full path to the config or cache directory.  The directory is created
        if it does not exist, and it is checked to make sure it can be written to.

    Notes
    -----
    This function's return value is cached after the first call in a given python
    session.  To clear that cache, execute `sxs_directory.cache_clear()`.

    In order of priority, if `persistent` is True, this function will choose one of
    these directories:

      1) The value found by `read_config("cache_directory")` if `directory_type` is
         "cache"
      2) Environment variable 'SXSCACHEDIR' or 'SXSCONFIGDIR'
      3) On 'linux' or 'freebsd' platforms
         a) Environment variable 'XDG_CACHE_HOME' or 'XDG_CONFIG_HOME'
         b) The '.cache/sxs' or '.config/sxs' directory in the user's home directory
      4) The '.sxs' directory in the user's home directory

    Whichever directory is chosen, this function will try to create the directory
    if necessary, and then check that it can be accessed and written to.  If that
    check fails, or if `persistent` is False, the function will create and return

      5) A temporary directory that will be removed when python exits

    This is based on matplotlib._get_config_or_cache_dir.

    """
    import warnings
    import sys
    import os
    import atexit
    import shutil
    import tempfile
    from pathlib import Path

    if directory_type not in ["cache", "config"]:
        raise ValueError(f"Can only find 'cache' or 'config' directories, not '{directory_type}'")

    suffix = Path("cache") if directory_type == "cache" else Path()

    if persistent:
        # Try to read config file first
        if directory_type == "cache":
            sxs_dir = read_config("cache_directory")
            if sxs_dir is not None:
                sxs_dir = Path(sxs_dir).expanduser().resolve()
                try:
                    sxs_dir.mkdir(parents=True, exist_ok=True)
                except OSError:
                    pass
                else:
                    if os.access(str(sxs_dir), os.W_OK) and sxs_dir.is_dir():
                        return sxs_dir
                message = (
                    f"\nThe `sxs` module failed to find or create a writable directory at {sxs_dir},\n"
                    f"even though {sxs_directory('config')} specified that directory for the cache."
                )
                warnings.warn(message)

        # Use SXSCONFIGDIR or SXSCACHEDIR if they are set
        sxs_dir = os.getenv(f'SXS{directory_type.upper()}DIR', default=False)
        if sxs_dir:
            sxs_dir = Path(sxs_dir).expanduser().resolve()

        # On linux/freebsd
        #     a) Use ${XDG_CONFIG_HOME}/sxs or ${XDG_CACHE_HOME}/sxs if they are set
        #     b) Default to ~/.config/sxs or ~/.cache/sxs
        elif sys.platform.startswith(('linux', 'freebsd')):
            xdg_base = os.environ.get(f'XDG_{directory_type.upper()}_HOME')
            if xdg_base is None:
                xdg_base = Path.home() / f".{directory_type}"
            sxs_dir = Path(xdg_base).expanduser().resolve() / "sxs"

        # Elsewhere, default to ~/.sxs ~/.sxs/cache
        else:
            sxs_dir = Path.home() / ".sxs" / suffix

        # Ensure that we have a writable directory
        try:
            sxs_dir.mkdir(parents=True, exist_ok=True)
        except OSError:
            pass
        else:
            if os.access(str(sxs_dir), os.W_OK) and sxs_dir.is_dir():
                return sxs_dir

        unwritable_dir = sxs_dir

    # We've fallen through to creating a temporary directory.  We want the cache directory to be a
    # subdirectory of the main directory, so we fetch it through the lru_cache.  We can't do this
    # above, because it should be possible to get different results for certain settings.
    if directory_type == "cache":
        sxs_dir = sxs_directory("config", persistent=persistent) / "cache"
        try:
            sxs_dir.mkdir(exist_ok=True)
        except OSError:
            pass
        else:
            if os.access(str(sxs_dir), os.W_OK) and sxs_dir.is_dir():
                return sxs_dir

    # If the config or cache directory cannot be created or is not a writable
    # directory, create a temporary one.
    tmpdir = os.environ[f'SXS{directory_type.upper()}DIR'] = tempfile.mkdtemp(prefix="sxs-")
    atexit.register(shutil.rmtree, tmpdir)
    sxs_dir = Path(tmpdir) / suffix
    if persistent:
        # noinspection PyUnboundLocalVariable
        message = (
            f"\nThe `sxs` module created a temporary {directory_type} directory at {sxs_dir}\n"
            f"because the default path ({unwritable_dir}) is not a writable directory;\n"
            f"it is highly recommended to set the SXS{directory_type.upper()}DIR environment\n"
            f"variable to a writable directory to enable caching of downloaded waveforms."
        )
        warnings.warn(message)
    if directory_type == "cache":
        sxs_dir.mkdir(exist_ok=True)
    return sxs_dir


def sxs_path_to_system_path(path):
    r"""Translate SXS path to a system-compatible path

    Parameters
    ----------
    path : str
        SXS-style path to a file — for example, r"SXS:BBH:0123\Lev4:Horizons.h5"
        becomes r"SXS_BBH_0123\Lev4_Horizons.h5" on Windows.  Other systems can
        handle the original path, so are not changed.

    Notes
    -----
    SXS-style paths begin with SXS IDs, which contain colon characters.  These
    colons are incompatible with Windows file systems, so we simply replace the
    colons with underscores.

    """
    if _platform_system == "Windows":
        return sxs_path_re.sub(lambda s: s.group(0).replace(":", "_"), str(path))
    else:
        return path


def cached_path(path_pattern):
    """Return path to file in local cache

    Parameters
    ----------
    path_pattern : str
        A pattern to search for among the catalog files.  See the docstring of
        `sxs.Catalog.select` for details about how this pattern is used.

    Returns
    -------
    path : pathlib.Path
        Full path to file on this system.  The file may not exist.

    Raises
    ------
    StopIteration
        When no result is found for the input `path_pattern`.

    Notes
    -----
    This function returns just one path, corresponding to the first result returned
    by `catalog.select_files(path_pattern)`.

    """
    import contextlib
    from .. import load
    with contextlib.redirect_stdout(None):
        catalog = load("catalog")
    file_infos = catalog.select_files(path_pattern)
    sxs_path, file_info = next(iter(file_infos.items()))
    cache_path = sxs_directory("cache")
    truepath = sxs_path_to_system_path(file_info.get("truepath", sxs_path))
    path = cache_path / truepath
    return path
