"""Various utilities used the sxs package"""

import functools
import numba

jit = functools.partial(numba.njit, cache=True)
vectorize = functools.partial(numba.vectorize, nopython=True, cache=True)
guvectorize = functools.partial(numba.guvectorize, nopython=True, cache=True)

# This is a generically reliable set of widths to feed to multishuffle when using 64-bit floats
default_shuffle_widths = (8, 8, 4, 4, 4, 2,) + (1,) * 34

from . import url, inspire, monotonicity, decimation, lvcnr, references
from .downloads import download_file
from .bitwise import xor, multishuffle
from .sxs_identifiers import (
    sxs_identifier_regex, sxs_identifier_re, lev_regex, lev_re, sxs_id, lev_number, simulation_title
)
from .sxs_directories import sxs_directory, read_config, write_config, sxs_path_to_system_path, cached_path
from .select import select_by_path_component
from .formats import file_format
from .pretty_print import fit_to_console
from .files import md5checksum, lock_file_manager, find_simulation_directories, find_files
from .dicts import KeyPassingDict


def version_info():
    """Find all relevant package versions

    This function attempts to import each of the various packages relevant to this
    package — such as numpy, quaternionic, etc. — and returns a dictionary mapping
    the package names to their version strings.  It attempts to import some
    packages — like spinsfast and scri — which are not required; if they are not
    found, they simply are not included in the result.

    """
    import importlib
    import sys
    info = {
        "python": sys.version,
    }
    modules = [
        "numpy", "scipy", "h5py", "numba",
        "quaternion", "spherical_functions", "spinsfast", "scri",
        "quaternionic", "spherical", "sxs"
    ]
    for module in modules:
        try:
            m = importlib.import_module(module)
            info[module] = m.__version__
        except Exception:
            pass
    return info
