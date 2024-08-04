"""Regular expressions for validating and parsing URLs"""

import re

# Much of this regex is quite similar to the django URL validator
# <https://github.com/django/django/blob/master/django/core/validators.py>
scheme = r"(?P<scheme>https?)://"
user_pass = r"(?P<user_pass>(?:[^\s:@/]+(?::[^\s:@/]*)?(?=@))?)@?"
characters = "a-z" + "\u00a1-\uffff" + "0-9"
ipv4 = r"(?:25[0-5]|2[0-4]\d|[0-1]?\d?\d)(?:\.(?:25[0-5]|2[0-4]\d|[0-1]?\d?\d)){3}"
ipv6 = r"\[[0-9a-f:.]+\]"
hostname = rf"[{characters}](?:[{characters}-]{{0,61}}[{characters}])?"
domain = rf"(?:\.(?!-)[{characters}-]{{1,63}}(?<!-))*"
tld = r"\.(?!-)(?:[a-z" + "\u00a1-\uffff" + r"-]{2,63}|xn--[a-z0-9]{1,59})(?<!-)\.?"
host = f"(?P<host>{hostname}{domain}{tld}|localhost)"
port = r":?(?P<port>(?:(?<=:)\d{2,5})?)"
resource = r"(?P<resource>(?:[/?#][^\s]*)?)"
url = scheme + user_pass + host + port + resource + r"\Z"

url_regex = re.compile(url, re.IGNORECASE)


def parse(url):
    """Parse the parts of a URL

    Parameters
    ----------
    url : str

    Returns
    -------
    match : dict
        If the URL is valid, this returns a dictionary with keys "scheme",
        "user_pass", "host", "port", and "resource" describing the various parts of
        the URL.  Note that if any of these parts are missing, the values of those
        keys will be empty strings.  If the URL is not valid, the dictionary will
        be empty.

    See Also
    --------
    validate : just return True or False depending on whether the URL is valid

    """
    match = url_regex.match(url)
    if hasattr(match, "groupdict"):
        return match.groupdict()
    else:
        return {}


def validate(url):
    """Validate a URL

    Parameters
    ----------
    url : str

    Returns
    -------
    valid : bool

    See Also
    --------
    parse : parse the parts of a URL into a dictionary

    """
    return bool(url_regex.match(url))
