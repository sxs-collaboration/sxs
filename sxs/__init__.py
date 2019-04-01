from __future__ import division, print_function, absolute_import

__doc_title__ = "SXS python code"
__doc__ = "A collection of python code used by the SXS collaboration"

from ._version import __version__

import sxs.doxygen
import sxs.metadata
import sxs.references
import sxs.validate
import sxs.zenodo

sxs_identifier_regex = r'(?P<sxs_identifier>SXS:(?P<simulation_type>BBH|BHNS|NSNS):(?P<sxs_number>[0-9]*))'


def sxs_id(s):
    """Return the SXS ID contained in the input string

    An SXS ID is anything that matches the following regular expression:

        SXS:(BBH|BHNS|NSNS):[0-9]*

    If no match is found, the empty string is returned.  If multiple matches are found, only the
    first is returned.

    """
    import re
    m = re.search(sxs_identifier_regex, s)
    if m:
        return m['sxs_identifier']
    else:
        return ''


