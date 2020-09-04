from .. import TimeSeries
from . import WaveformMixin


class WaveformGrid(TimeSeries, WaveformMixin):
    def __init__(self, *args, **kwargs):
        raise NotImplementedError()
