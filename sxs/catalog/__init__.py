"""Interface to the catalog of SXS data

"""

from .catalog import Catalog


formats = {
    None: Catalog,
    "": Catalog,
}
