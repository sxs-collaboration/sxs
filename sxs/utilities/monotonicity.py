"""Functions to check and enforce monotonicity of array values"""

import numpy as np
from . import jit


@jit
def index_is_monotonic(y):
    length = y.size
    monotonic = np.ones_like(y, dtype=np.bool_)
    direction = y[-1] - y[0]
    if direction > 0.0:
        max_value = y[0]
        for i in range(1, length):
            if y[i] <= max_value:
                monotonic[i] = False
            else:
                max_value = y[i]
    else:
        min_value = y[0]
        for i in range(1, length):
            if y[i] >= min_value:
                monotonic[i] = False
            else:
                min_value = y[i]
    return monotonic


def monotonic_indices(y):
    indices = np.arange(y.size)
    return indices[index_is_monotonic(y)]


def monotonize(y):
    return y[monotonic_indices(y)]
