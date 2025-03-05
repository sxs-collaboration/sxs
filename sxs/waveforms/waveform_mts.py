import spherical
from . import WaveformModes

class _ModesTimeSeries(WaveformModes):

    @property
    def s(self):
        return self.spin_weight

    @property
    def u(self):
        return self.time

    from spherical.modes.algebra import (
        conjugate, bar, _real_func, real, _imag_func, imag, norm,
        add, subtract, multiply, divide
    )

    conj = conjugate

    from spherical.modes.utilities import (
        truncate_ell, _check_broadcasting
    )

    from spherical.modes.ufuncs import __array_ufunc__


def MTS(*args, **kwargs):
    kwargs.setdefault('multiplication_truncator', max)
    return _ModesTimeSeries(*args, **kwargs)
