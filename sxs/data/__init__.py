"""Interface to SXS data sets"""

from .time_series import TimeSeries
from .metadata import Metadata
from .horizons import Horizons, HorizonQuantities
from .waveforms import WaveformModes, WaveformGrid, WaveformSignal
from . import time_series, metadata, horizons, waveforms


def load(sxs_id, version='latest', data='h', lev='+', extrapolation='N4', inertial=True, **kwargs):
    """

    """
    raise NotImplementedError()
