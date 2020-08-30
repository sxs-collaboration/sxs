# Copyright (c) 2020, Michael Boyle
# See LICENSE file for details:
# <https://github.com/sxs-collaboration/sxs/blob/master/LICENSE>

try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:  # pragma: no cover
    import importlib_metadata

__version__ = importlib_metadata.version(__name__)

from . import catalog, data, utilities


def load(*args, **kwargs):
    """Load SXS data

    

    """
    if len(args) == 0:
        raise ValueError("No arguments passed")
    if args[0] == 'catalog':
        return catalog.load(args[1:], **kwargs)
    raise NotImplementedError()
