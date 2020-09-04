from .. import TimeSeries
from . import WaveformMixin


class WaveformModes(TimeSeries, WaveformMixin):
    def __init__(self, *args, **kwargs):
        raise NotImplementedError()
