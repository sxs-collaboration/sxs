"""Interface to SXS metadata files"""

from .metadata import Metadata

formats = {
    None: Metadata,
    "": Metadata,
    "metadata": Metadata,
}
