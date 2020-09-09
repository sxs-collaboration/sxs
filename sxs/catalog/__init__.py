"""Interface to the catalog of SXS data

"""

import functools
from .catalog import Catalog


@functools.lru_cache()
def load():
    """Load SXS catalog from file, optionally downloading

    """
    raise NotImplementedError()


formats = {
    None: Catalog,
    "": Catalog,
}
