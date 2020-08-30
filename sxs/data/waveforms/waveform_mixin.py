import abc


class WaveformMixin(abc.ABC):
    def __init__(self, *args, **kwargs):
        raise NotImplementedError()
