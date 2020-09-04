"""Interface to SXS metadata files"""

from .container import Metadata

formats = {
    None: Metadata,
    "metadata": Metadata,
}
