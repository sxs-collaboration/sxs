"""Low-level interface to Invenio-backed APIs

"""

from .login import Login
from .deposit import Deposit
from .records import Records

url_standard = 'https://zenodo.org/'
url_sandbox = 'https://sandbox.zenodo.org/'
