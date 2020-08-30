"""Interface to SXS data sets


"""



from .metadata import Metadata
from .time_series import TimeSeries
from .waveforms import WaveformModes, WaveformGrid, WaveformSignal
from . import formats, metadata, time_series, horizons, waveforms


def load(sxs_id, version='latest', data='h', lev='+', extrapolation='N4', inertial=True, *args, **kwargs):
    """

    """
    raise NotImplementedError()
