from __future__ import division, print_function, absolute_import

__doc_title__ = "SXS python code"
__doc__ = "A collection of python code used by the SXS collaboration"

from ._version import __version__

import sxs.doxygen
import sxs.format
import sxs.metadata
import sxs.references
import sxs.utilities
import sxs.validate
import sxs.zenodo

sxs_identifier_regex = r'(?P<sxs_identifier>SXS:(?P<simulation_type>BBH|BHNS|NSNS):(?P<sxs_number>[0-9]*))'
lev_regex = r'Lev(?P<lev>[-0-9]*)'


def sxs_id(s):
    """Return the SXS ID contained in the input string

    An SXS ID is anything that matches the following regular expression:

        SXS:(BBH|BHNS|NSNS):[0-9]*

    If no match is found, the empty string is returned.  If multiple matches are found, only the
    first is returned.

    The input parameter may be
    1) A string believed to contain the SXS ID
    2) The path to a file (like common-metadata.txt) containing such a string
    3) An object with a 'title' attribute
    4) An object with a 'title' item

    """
    import os.path
    import re
    try:
        s = s['title']
    except (TypeError, KeyError):
        pass
    if hasattr(s, 'title') and not isinstance(s.title, type('a'.title)):
        s = s.title
    try:
        if os.path.isfile(s):
            with open(s, 'r') as f:
                s = [l.strip() for l in f.splitlines()]
            for line in s:
                sxs_id_line = sxs_id(line)
                if sxs_id_line:
                    return sxs_id_line
            return ''
    except TypeError:
        pass
    m = re.search(sxs_identifier_regex, s)
    if m:
        return m['sxs_identifier']
    else:
        return ''


def lev_number(s):
    """Return the integer Lev number contained in the input string

    The Lev number is anything that matches the following regular expression:

        Lev(?P<lev>[-0-9]*)

    If no match is found, `None` is returned.  If multiple matches are found,
    only the first is returned.

    """
    import os.path
    import re
    m = re.search(lev_regex, s)
    if m:
        return int(m['lev'])
    else:
        return None
