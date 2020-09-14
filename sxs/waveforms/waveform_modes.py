import numpy as np
import spherical
from .. import TimeSeries
from . import WaveformMixin


class WaveformModes(TimeSeries, WaveformMixin):
    import functools

    def __new__(cls, input_array, *args, **kwargs):
        for requirement in ["mode_axis", "ell_min", "ell_max"]:
            if requirement not in kwargs and requirement not in getattr(input_array, "_metadata", {}):
                raise ValueError(f"{cls} could not find required argument '{requirement}'")
        self = super().__new__(cls, input_array, *args, **kwargs)
        return self

    @property
    def t(self):  # A handy alias for backwards-compatibility
        return self.time

    @property
    def data(self):  # A handy alias for backwards-compatibility
        return self.ndarray

    @property
    def modes_axis(self):
        return self._metadata["modes_axis"]

    @property
    def ell_min(self):
        return self._metadata["ell_min"]

    @property
    def ell_max(self):
        return self._metadata["ell_max"]

    @property
    def n_modes(self):
        return self.shape[self.modes_axis]

    @property
    @functools.lru_cache()
    def LM(self):
        """Array of (ell, m) values in the data

        This array is just a flat array of `[ell, m]` pairs.  It is automatically
        recomputed each time `ell_min` or `ell_max` is changed.  Specifically, it is

            np.array([
                [ell, m]
                for ell in range(self.ell_min, self.ell_max+1)
                for m in range(-ell, ell+1)
            ])

        """
        return np.array([[ell, m] for ell in range(self.ell_min, self.ell_max+1) for m in range(-ell, ell+1)])

    def index(self, ell, m):
        """Index of given (ell,m) mode in the data

        Parameters
        ----------
        ell : int
        m : int

        Returns
        -------
        idx : int
            Index such that self.LM[idx] == [ell, m]

        """
        return spherical.LM_index(ell, m, self.ell_min)

    @property
    def data_type(self):
        return self._metadata["data_type"]

    @property
    def abs(self):
        """Absolute value of the data"""
        return np.abs(self)

    @property
    def arg(self):
        """Complex phase angle of the data

        Note that the result is not "unwrapped", meaning that there may be
        discontinuities as the phase approaches ±π.

        Returns
        -------
        phase : TimeSeries
            Values are in the interval (-π, π].

        See Also
        --------
        numpy.angle
        WaveformModes.arg_unwrapped

        """
        return np.angle(self).view(TimeSeries)

    @property
    def arg_unwrapped(self):
        """Complex phase angle of the data, unwrapped along the time axis

        The result is "unwrapped", meaning that discontinuities as the phase approaches
        ±π are removed

        Returns
        -------
        phase : TimeSeries
            Values at the initial time are in the interval (-π, π], but may evolve to
            arbitrary real values.

        See Also
        --------
        numpy.angle
        numpy.unwrap
        WaveformModes.arg

        """
        return np.unwrap(self.arg, axis=self.time_axis)

    def norm(self, take_sqrt=False, indices=slice(None, None, None)):
        """Compute a norm of the waveform

        Parameters
        ----------
        take_sqrt : bool, optional
            If True, return the standard L² norm, which involves taking the square
            root.  If False (the default), the square root is not taken, meaning that
            the returned value is the square of the usual L² norm.
        indices : array_like, optional
            If present, the data is sliced along the time axis with this quantity —
            specifically, using `numpy.take`.  See that function's docstring for more
            explanation.

        Returns
        -------
        n : TimeSeries

        See Also
        --------
        numpy.linalg.norm
        numpy.take

        Notes
        -----
        The integral of the (squared) magnitude of the data equals the sum of the
        (squared) magnitude of the modes for orthonormal basis functions, meaning that
        the L² norm of the function equals the basic Euclidean norm of the modes.  We
        assume that these modes are expanded in a band-limited but otherwise complete
        orthonormal basis.

        """
        # Because of MKL and similar fanciness, np.linalg.norm is faster than most
        # other built-in things, even with its sqrt.  Numba can be faster; taking this
        # shortcut out of laziness for now.
        if indices == slice(None, None, None):
            n = np.linalg.norm(self.ndarray, axis=self.modes_axis)
        else:
            n = np.linalg.norm(np.take(self.ndarray, indices, axis=self.time_axis), axis=self.modes_axis)
        if take_sqrt:
            return n
        else:
            return n**2

    def max_norm_index(self, skip_fraction_of_data=4):
        """Index of time step with largest norm

        The optional argument skips a fraction of the data.  The default is
        4, which means that it only searches the last three-fourths of the
        data for the max.  If 0 or 1 is input, this is ignored, and all the
        data is searched.

        """
        if skip_fraction_of_data == 0 or skip_fraction_of_data == 1:
            indices = slice(None, None, None)
            return np.argmax(self.norm(indices=indices))
        else:
            indices = slice(self.n_times // skip_fraction_of_data, None, None)
            return np.argmax(self.norm(indices=indices)) + (self.n_times // skip_fraction_of_data)

    def max_norm_time(self, skip_fraction_of_data=4):
        """Return time at which largest norm occurs in data

        See `help(max_norm_index)` for explanation of the optional argument.

        """
        return self.t[self.max_norm_index(skip_fraction_of_data=skip_fraction_of_data)]

    def truncate(self, tol=1e-10):
        """Truncate the precision of this object's data in place

        This function sets bits in the data to 0 when they have lower significance than
        will alter the norm of the Waveform by a fraction `tol` at that instant in
        time.

        """
        if tol != 0.0:
            tol_per_mode = tol / np.sqrt(self.n_modes)
            absolute_tolerance = np.linalg.norm(self.ndarray, axis=self.modes_axis) * tol_per_mode
            power_of_2 = (2.0 ** np.floor(-np.log2(absolute_tolerance)))[:, np.newaxis]
            self.ndarray *= power_of_2
            np.round(self.ndarray, out=self.ndarray)
            self.ndarray /= power_of_2

    def convert_to_conjugate_pairs(self):
        """Convert modes to conjugate-pair format in place

        This function alters this object's modes to store the sum and difference of
        pairs with opposite `m` values.  If we denote the modes `f[l, m]`, then we
        define

            s[l, m] = (f[l, m] + f̄[l, -m]) / √2
            d[l, m] = (f[l, m] - f̄[l, -m]) / √2

        For m<0 we replace the mode data with `d[l, -m]`, for m=0 we do nothing, and
        for m>0 we replace the mode data with `s[l, m]`.  That is, the mode data on
        output look like this:

            [d[2, 2], d[2, 1], f[2, 0], s[2, 1], s[2, 2], d[3, 3], d[3, 2], ...]

        The factor of √2 is chosen so that the norm (sum of the magnitudes squared) at
        each time for this data is the same as it is for the original data.

        """
        for ell in range(self.ell_min, self.ell_max + 1):
            for m in range(1, ell + 1):
                i_plus = self.index(ell, m)
                i_minus = self.index(ell, -m)
                mode_plus = self.ndarray[..., i_plus].copy()
                mode_minus = self.ndarray[..., i_minus].copy()
                self.ndarray[..., i_plus] = (mode_plus + np.conjugate(mode_minus)) / np.sqrt(2)
                self.ndarray[..., i_minus] = (mode_plus - np.conjugate(mode_minus)) / np.sqrt(2)

    def convert_from_conjugate_pairs(self):
        """Convert modes from conjugate-pair format in place

        This function reverses the effects of `convert_to_conjugate_pairs`.  See that
        function's docstring for details.

        """
        for ell in range(self.ell_min, self.ell_max + 1):
            for m in range(1, ell + 1):
                i_plus = self.index(ell, m)
                i_minus = self.index(ell, -m)
                mode_plus = self.ndarray[..., i_plus].copy()
                mode_minus = self.ndarray[..., i_minus].copy()
                self.ndarray[..., i_plus] = (mode_plus + mode_minus) / np.sqrt(2)
                self.ndarray[..., i_minus] = np.conjugate(mode_plus - mode_minus) / np.sqrt(2)
