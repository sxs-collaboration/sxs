"""Interface to the catalog of SXS data

"""

import functools
from .catalog import Catalog


formats = {
    None: Catalog,
    "": Catalog,
}
