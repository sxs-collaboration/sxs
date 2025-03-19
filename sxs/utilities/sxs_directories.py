"""Functions to find user-specific config and cache directories"""

import re
import platform
import functools
from sxscatalog.utilities import read_config, write_config, sxs_directory
from .sxs_identifiers import sxs_path_re

_platform_system = platform.system()


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
