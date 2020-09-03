"""Collection of utilities for the sxs package"""

import functools
import numba

jit = functools.partial(numba.njit, cache=True)
vectorize = functools.partial(numba.vectorize, nopython=True, cache=True)
guvectorize = functools.partial(numba.guvectorize, nopython=True, cache=True)

# This is a generically reliable set of widths to feed to multishuffle when using 64-bit floats
default_shuffle_widths = (8, 8, 4, 4, 4, 2,) + (1,) * 34

from .bitwise import xor, multishuffle
