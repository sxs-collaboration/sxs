"""Simple regexes to understand SXS IDs"""

import re
import os

sep_regex = r"(:|/)" if os.sep == "/" else r"(:|/|\\)"

sxs_identifier_regex = (
    r"(?P<sxs_identifier>SXS:(?P<simulation_type>BBH|BHNS|NSNS)(?:_ExtCCE)?:"
    r"(?P<sxs_number>[0-9]+))(?:(v|V)(?P<version>[0-9.]+))?"
)
lev_regex = r"Lev(?P<lev>-?[0-9]+)"
sxs_id_version_lev_regex = sxs_identifier_regex + rf"(?:{sep_regex}{lev_regex})?"
sxs_id_version_lev_exact_regex = f"^{sxs_id_version_lev_regex}$"

file_regex = r"(?P<file>[a-zA-Z0-9_]+\.[a-zA-Z0-9]+)"
sxs_path_regex = sxs_id_version_lev_regex + rf"(?:{sep_regex}{file_regex})?"

sxs_identifier_re = re.compile(sxs_identifier_regex)
lev_re = re.compile(lev_regex)
sxs_id_version_lev_re = re.compile(sxs_id_version_lev_regex)
sxs_id_version_lev_exact_re = re.compile(sxs_id_version_lev_exact_regex)
sxs_path_re = re.compile(sxs_path_regex)

def sxs_id(s, default="", include_version=False):
    """Return the SXS ID contained in the input string

    An SXS ID is anything that matches the following regular expression:

        SXS:(BBH|BHNS|NSNS)(_ExtCCE)?:[0-9]+(v[0-9.]+)?

    If no match is found, the empty string is returned.  If multiple matches are
    found, only the first is returned.

    The input parameter may be:
    1) A string believed to contain the SXS ID
    2) The path to a file (like common-metadata.txt) containing such a string
    3) An object with a 'title' attribute
    4) An object with a 'title' item

    """
    id, version = sxs_id_and_version(s, default)
    if include_version:
        return f"{id}{version}"
    else:
        return id


def sxs_id_and_version(s, default=""):
    """Return the SXS ID and version contained in the input string"""
    import os.path
    import re
    try:
        s = s["title"]
    except (TypeError, KeyError):
        pass
    if hasattr(s, "title") and not isinstance(s.title, type("a".title)):
        s = s.title
    try:
        if os.path.isfile(s):
            with open(s, "r") as f:
                s = [l.strip() for l in f.splitlines()]
            for line in s:
                sxs_id_line = sxs_id_and_version(line)
                if sxs_id_line[0]:
                    return sxs_id_line
            return default, ""
    except TypeError:
        pass
    m = re.search(sxs_identifier_regex, s)
    if m:
        return m["sxs_identifier"], (f"v{m['version']}" if m["version"] else "")
    else:
        return default, ""


def simulation_title(sxs_id):
    """Create simulation title given just the SXS ID

    Currently supported types return output that looks like

        Binary black-hole simulation SXS:BBH:0001
        Black-hole neutron-star binary simulation SXS:BHNS:0001
        Binary neutron-star simulation SXS:NSNS:0001

    """
    import re
    sxs_system_re = re.compile(sxs_identifier_regex)
    m = sxs_system_re.search(sxs_id)
    if not m:
        raise ValueError(f"No SXS identifier found in '{sxs_id}'")
    sxs_identifier = m["sxs_identifier"]
    simulation_type = m["simulation_type"]
    if simulation_type == "BBH":
        title = f"Binary black-hole simulation {sxs_identifier}"
    elif simulation_type == "BHNS":
        title = f"Black-hole neutron-star binary simulation {sxs_identifier}"
    elif simulation_type == "NSNS":
        title = f"Binary neutron-star simulation {sxs_identifier}"
    else:
        raise ValueError(f"Did not recognize SXS system type '{simulation_type}'; should be BBH, BHNS, or NSNS.")
    return title


def lev_number(s):
    """Return the integer Lev number contained in the input string

    The Lev number is anything that matches the following regular expression:

        Lev(?P<lev>-?[0-9]+)

    If no match is found, `None` is returned.  If multiple matches are found,
    only the first is returned.

    """
    import re
    m = re.search(lev_regex, s)
    if m:
        return int(m["lev"])
    else:
        return None


def sxs_id_to_url(sxs_id):
    """Return the DOI URL for the SXS ID"""
    return f"https://doi.org/10.26138/{sxs_id}"
