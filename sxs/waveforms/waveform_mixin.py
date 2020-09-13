import abc


class WaveformMixin(abc.ABC):
    @property
    @abc.abstractmethod
    def data(self):
        return self.ndarray

    # metadata
    # spin_weight
    # boost_weight
    # type
    #     name
    #     radius_scaling
    #     mass_scaling
    #     conformal_scaling
    #     frame ['inertial', 'corotating', ...]
    # copy()
    # evaluate(units, total_mass, theta, phi, orientation, time, ...)  # Return WaveformSignal
