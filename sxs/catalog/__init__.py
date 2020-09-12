"""Interface to the catalog of SXS data

"""

from .catalog import Catalog
from .create import _create
from .description import catalog_file_description


formats = {
    None: Catalog,
    "": Catalog,
}
