from .. import TimeSeries
from . import WaveformMixin


class WaveformSignal(TimeSeries, WaveformMixin):
    def __init__(self, *args, **kwargs):
        raise NotImplementedError()
