"""Base class for all waveform-like objects"""

import abc
import numpy as np
import quaternionic


class WaveformMixin(abc.ABC):
    @property
    @abc.abstractmethod
    def data(self):  # Handy alias for backwards compatibility
        return self.ndarray

    @property
    @abc.abstractmethod
    def data_type(self):
        return self._metadata.get("data_type", "unknown")

    @property
    @abc.abstractmethod
    def spin_weight(self):
        return self._metadata.get("spin_weight", None)

    @property
    @abc.abstractmethod
    def frame(self):
        return self._metadata.get("frame", np.atleast_2d(quaternionic.one))

    @property
    @abc.abstractmethod
    def frame_type(self):
        if "frame" in self._metadata:
            return self._metadata.get("frame_type", "unknown")
        else:
            return self._metadata.get("frame_type", "inertial")

    # metadata
    # boost_weight
    # type
    #     name
    #     radius_scaling
    #     mass_scaling
    #     conformal_scaling
    #     frame ∈ {'inertial', 'corotating', ...}
    # copy()
    # evaluate(units, total_mass, theta, phi, orientation, time, ...)  # Return WaveformSignal
