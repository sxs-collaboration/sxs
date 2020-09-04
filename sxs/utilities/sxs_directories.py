import functools


@functools.lru_cache()
def get_sxs_directory(directory_type):
    """Return the SXS directory location

    Parameters
    ----------
    directory_type : {"cache", "config"}
        The type of user directory to be found.

    Returns
    -------
    directory : str
        The full path to the directory.  This is created if it does not exist, and
        it is checked to make sure it can be written to.

    Notes
    -----
    This function determines the location of the SXS directory of the given type —
    either 'cache' or 'config' — and ensures that it exists and is writeable.

    1) Environment variable 'SXSCACHEDIR' or 'SXSCONFIGDIR'
    2) On 'linux' or 'freebsd' platforms
       a) Environment variable 'XDG_CACHE_HOME' or 'XDG_CONFIG_HOME'
       b) The '.cache/sxs' or '.config/sxs' directory in the user's home directory
    3) The '.sxs' directory in the user's home directory
    4) A temporary directory that will be removed when python exits

    This is based on matplotlib's _get_config_or_cache_dir function.

    """
    import warnings
    import sys
    import os
    import atexit
    import shutil
    import tempfile
    from pathlib import Path

    suffix = Path("cache") if directory_type == "cache" else Path()
    if directory_type == "cache":
        suffix = ""

    configdir = os.environ.get(f'SXS{directory_type.upper()}DIR')
    if configdir:
        configdir = Path(configdir).expanduser().resolve()
    elif sys.platform.startswith(('linux', 'freebsd')):
        xdg_base = os.environ.get(f'XDG_{directory_type.upper()}_HOME') or Path.home() / f".{directory_type}"
        configdir = Path(xdg_base).expanduser().resolve() / "sxs"
    else:
        configdir = Path.home() / ".sxs"

    try:
        configdir.mkdir(parents=True, exist_ok=True)
    except OSError:
        pass
    else:
        if os.access(str(configdir), os.W_OK) and configdir.is_dir():
            return str(configdir)

    # If the config or cache directory cannot be created or is not a writable
    # directory, create a temporary one.
    tmpdir = os.environ[f'SXS{directory_type.upper()}DIR'] = tempfile.mkdtemp(prefix="sxs-")
    atexit.register(shutil.rmtree, tmpdir)
    warnings.warn(
        f"SXS created a temporary config/cache directory at {configdir} because "
        f"the default path ({tmpdir}) is not a writable directory; it is highly "
        f"recommended to set the SXS{directory_type.upper()}DIR environment variable to a "
        f"writable directory to enable caching of downloaded waveforms."
    )
    return tmpdir
