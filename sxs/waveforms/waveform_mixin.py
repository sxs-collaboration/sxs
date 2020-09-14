import abc


class WaveformMixin(abc.ABC):
    # Note: This is currently a pretty trivial class, but as the code develops
    # this will be increasingly useful.

    @property
    @abc.abstractmethod
    def t(self):  # Handy alias for backwards compatibility
        return self.time

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
    def frame(self):
        return self._metadata.get("frame", np.atleast_2d(quaternionic.one))

    @property
    @abc.abstractmethod
    def frame_type(self):
        return self._metadata.get("frame_type", "unknown")

    # metadata
    # spin_weight
    # boost_weight
    # type
    #     name
    #     radius_scaling
    #     mass_scaling
    #     conformal_scaling
    #     frame ∈ {'inertial', 'corotating', ...}
    # copy()
    # evaluate(units, total_mass, theta, phi, orientation, time, ...)  # Return WaveformSignal
