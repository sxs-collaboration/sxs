"""Interface to SXS data sets


"""



from . import formats, metadata, time_series, horizon_quantities, waveforms
from .metadata import Metadata
from .time_series import TimeSeries
from .waveforms import WaveformModes, WaveformGrid, WaveformSignal


def load(sxs_id, data='h', lev='+', extrapolation='N4', inertial=True, *args, **kwargs):
    """

    """
    raise NotImplementedError()
