from __future__ import division, print_function, absolute_import

__doc_title__ = "SXS python code"
__doc__ = "A collection of python code used by the SXS collaboration"

from ._version import __version__

import sxs.doxygen
import sxs.metadata
import sxs.references
import sxs.zenodo

sxs_identifier_regex = r'(?P<sxs_identifier>SXS:(?P<simulation_type>BBH|BHNS|NSNS):(?P<sxs_number>[0-9]*))'
