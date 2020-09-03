"""Interface to SXS data sets"""

from .time_series import TimeSeries
from .metadata import Metadata
from .horizons import Horizons, HorizonQuantities
from .waveforms import WaveformModes, WaveformGrid, WaveformSignal
from . import time_series, metadata, horizons, waveforms


# This is a generically reliable set of widths to feed to multishuffle when using 64-bit floats
default_shuffle_widths = (8, 8, 4, 4, 4, 2,) + (1,) * 34


def load(sxs_id, version='latest', data='h', lev='+', extrapolation='N4', inertial=True, **kwargs):
    """

    """
    raise NotImplementedError()
