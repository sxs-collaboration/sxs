"""Various utilities used the sxs package"""

import functools
import numba

jit = functools.partial(numba.njit, cache=True)
vectorize = functools.partial(numba.vectorize, nopython=True, cache=True)
guvectorize = functools.partial(numba.guvectorize, nopython=True, cache=True)

# This is a generically reliable set of widths to feed to multishuffle when using 64-bit floats
default_shuffle_widths = (16, 4, 2, 2, 2) + (1,) * 22 + (4,) * 4
default_shuffle_widths_old = (8, 8, 4, 4, 4, 2,) + (1,) * 34

from . import url, inspire, monotonicity, decimation, lvcnr, references
from .downloads import download_file
from .bitwise import diff, xor, multishuffle
from .sxs_identifiers import (
    sxs_identifier_regex, sxs_identifier_re,
    lev_regex, lev_re,
    sxs_id_version_lev_regex, sxs_id_version_lev_re,
    sxs_id_version_lev_exact_regex, sxs_id_version_lev_exact_re,
    sxs_path_regex, sxs_path_re,
    sxs_id, sxs_id_and_version,
    lev_number, simulation_title, sxs_id_to_url,
)
from .sxs_directories import (
    sxs_directory, read_config, write_config, sxs_path_to_system_path, cached_path
)
from .select import select_by_path_component
from .formats import file_format
from .pretty_print import fit_to_console
from .files import (
    md5checksum, lock_file_manager, find_simulation_directories, find_files
)
from .dicts import KeyPassingDict
from .smooth_functions import (
    transition_function, transition_function_inplace,
    transition_function_derivative, transition_function_derivative_inplace,
    transition_to_constant, transition_to_constant_inplace,
    bump_function, bump_function_inplace
)
from .inspire import inspire2doi


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


class SimpleVersion:
    """Basic versioning class

    This class uses a simple two-element version string, enabling printing,
    incrementing, and comparison.

    Note that this is just intended for use SXS data versions; for more general
    use, you probably want something like
    [`packaging.version`](https://packaging.pypa.io/en/latest/version.html).

    """
    def __init__(self, v):
        if "." not in v:
            v = f"1.{v}"
        self.major, self.minor = map(int, v.split("."))
    def __str__(self):
        return f"{self.major}.{self.minor}"
    def __repr__(self):
        return f"{self.major}.{self.minor}"
    def increment(self, level="minor"):
        if level=="major":
            self.major += 1
        elif level=="minor":
            self.minor += 1
        else:
            raise ValueError("""Can only increment major or minor levels, not "{level}".""")
        return self

    def __hash__(self):
        return hash(str(self))

    def __lt__(self, other):
        if not isinstance(other, SimpleVersion):
            return NotImplemented
        return self.major < other.major or (self.major == other.major and self.minor < other.minor)

    def __le__(self, other):
        if not isinstance(other, SimpleVersion):
            return NotImplemented
        return self.major <= other.major or (self.major == other.major and self.minor <= other.minor)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SimpleVersion):
            return NotImplemented
        return self.major == other.major and self.minor == other.minor

    def __ge__(self, other):
        if not isinstance(other, SimpleVersion):
            return NotImplemented
        return self.major >= other.major or (self.major == other.major and self.minor >= other.minor)

    def __gt__(self, other):
        if not isinstance(other, SimpleVersion):
            return NotImplemented
        return self.major > other.major or (self.major == other.major and self.minor > other.minor)

    def __ne__(self, other: object) -> bool:
        if not isinstance(other, SimpleVersion):
            return NotImplemented
        return self.major != other.major or self.minor != other.minor
