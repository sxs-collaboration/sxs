"""Interface to SXS metadata files"""

from .metadata import Metadata
from .metric import MetadataMetric

formats = {
    None: Metadata,
    "": Metadata,
    "metadata": Metadata,
}
