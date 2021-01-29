"""The main container for waveform objects with mode weights"""

import re
import numpy as np
import quaternionic
import spherical
from .. import TimeSeries
from . import WaveformMixin
from .mode_utilities import _expectation_value_LL, _expectation_value_Ldt

NRAR_mode_regex = re.compile(r"""Y_l(?P<L>[0-9]+)_m(?P<M>[-+0-9]+)\.dat""")


class WaveformModes(WaveformMixin, TimeSeries):
    """Array-like object representing time-series data of SWSH modes

    We generally assume that the data represent mode weights in expansions of
    functions in terms of spin-weighted spherical harmonics (where the standard
    spherical harmonics just happen to have spin weight 0).

    This object is based on the TimeSeries object, but has many additional methods
    for manipulating modes.

    Parameters
    ----------
    input_array : (..., N, ..., M, ...) array_like
        Input data representing the dependent variable, in any form that can be
        converted to a numpy array.  This includes scalars, lists, lists of tuples,
        tuples, tuples of tuples, tuples of lists, and numpy ndarrays.  It can have
        an arbitrary number of dimensions, but the length `N` along `time_axis`
        must match the length of `time`, and the length `M` along `modes_axis` (see
        below) must be consistent with `ell_min` and `ell_max`.  Values must be
        finite.
    time : (N,) array_like
        1-D array containing values of the independent variable.  Values must be
        real, finite, and in strictly increasing order.
    time_axis : int, optional
        Axis along which `input_array` is assumed to be varying in time, meaning
        that for `time[i]` the corresponding values are `np.take(input_array, i,
        axis=time_axis)`.  If this is not given, the first axis of `input_array`
        that has the same length as `time` is chosen as the time axis — which may
        be prone to errors.
    modes_axis : int
        Axis of the array along which the modes are stored.  See Notes below.
    ell_min : int
        Smallest value of ℓ stored in the data
    ell_max : int
        Largest value of ℓ stored in the data

    Notes
    -----
    We assume that the modes at a given time, say, are stored in a flat array, in
    order of increasing `m` and `ell`, with `m` varying most rapidly — something
    like the following:

        [f(ell, m) for ell in range(ell_min, ell_max+1) for m in range(-ell,ell+1)]

    The total size is implicitly `ell_max * (ell_max + 2) - ell_min ** 2 + 1`.

    For backwards compatibility, it is possible to retrieve individual modes in the
    same way as the old NRAR-format HDF5 files would be read, as in

        h_22 = waveform["Y_l2_m2.dat"]

    Note that "History.txt" may not contain anything but an empty string, because
    history is not retained in more recent data formats.  Also note that — while
    not strictly a part of this class — the loaders that open waveform files will
    return a dict-like object when the extrapolation order is not specified.  That
    object can also be accessed in a backwards-compatible way much like the root
    directory of the NRAR-format HDF5 files.  For example:

        with sxs.loadcontext("rhOverM_Asymptotic_GeometricUnits_CoM.h5") as f:
            h_22 = f["Extrapolated_N2.dir/Y_l2_m2.dat"]

    This code is identical to the equivalent code using `h5py` except that the call
    to `h5py.File` is replaced with the call to `sxs.loadcontext`.  The `.dat`
    datasets are reconstructed on the fly, but should be bitwise-identical to the
    output from the HDF5 file whenever the underlying format is NRAR.

    """
    import functools

    def __new__(cls, input_array, *args, **kwargs):
        for requirement in ["modes_axis", "ell_min", "ell_max"]:
            if requirement not in kwargs and requirement not in getattr(input_array, "_metadata", {}):
                raise ValueError(f"{cls} could not find required argument '{requirement}'")
        self = super().__new__(cls, input_array, *args, **kwargs)
        return self

    def __getitem__(self, key):
        if isinstance(key, str):
            if key == "History.txt":
                return self._metadata.get("history", "")
            # Assume we're asking for Y_l2_m2.dat or something
            if self.ndarray.shape != (self.n_times, self.n_modes):
                raise ValueError(f"Data has shape {self.ndarray.shape}, which is incompatible with NRAR format")
            match = NRAR_mode_regex.match(key)
            if not match:
                raise ValueError(f"Key '{key}' did not match mode format")
            ell, m = int(match["L"]), int(match["M"])
            if ell < self.ell_min or ell > self.ell_max:
                raise ValueError(f"Key '{key}' requested ell={ell} value not found in this data")
            if abs(m) > ell:
                raise ValueError(f"Key '{key}' requested (ell, m)=({ell}, {m}) value does not make sense")
            index = self.index(ell, m)
            data = np.take(self.ndarray, index, axis=self.modes_axis).view(float).reshape((-1, 2))
            return np.hstack((self.time[:, np.newaxis], data))
        obj, time_key = self._slice(key)
        if time_key is None:
            raise ValueError(f"Fancy indexing (with {key}) is not supported")
        obj = obj.view(type(self))
        if "frame" in obj._metadata and obj.frame.shape == (self.n_times, 4):
            obj._metadata["frame"] = obj.frame[time_key, :]
        return obj

    @property
    def modes_axis(self):
        """Axis of the array storing the various modes

        See the documentation of this class for an explanation of what this means.

        See Also
        --------
        time_axis : Axis of the array along which time varies

        """
        return self._metadata["modes_axis"]

    @property
    def ell_min(self):
        """Smallest value of ℓ stored in the data"""
        return self._metadata["ell_min"]

    @property
    def ell_max(self):
        """Largest value of ℓ stored in the data"""
        return self._metadata["ell_max"]

    @property
    def n_modes(self):
        """Total number of mode weights stored in the data"""
        return self.shape[self.modes_axis]

    @property
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
        """Mode index of given (ell,m) mode in the data

        Parameters
        ----------
        ell : int
        m : int

        Returns
        -------
        idx : int
            Index such that self.LM[idx] == [ell, m]

        """
        return spherical.Yindex(ell, m, self.ell_min)

    @property
    def abs(self):
        """Absolute value of the data

        Returns
        -------
        absolute : TimeSeries
            Because the absolute values make no sense as mode weights, this is just a
            plain TimeSeries object.

        See Also
        --------
        arg

        """
        return np.abs(self.view(TimeSeries))

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
        arg_unwrapped

        """
        return np.angle(self.view(TimeSeries))

    @property
    def arg_unwrapped(self):
        """Complex phase angle of the data, unwrapped along the time axis

        The result is "unwrapped", meaning that discontinuities as the phase approaches
        ±π are removed by adding an appropriate amount to all following data points.

        Returns
        -------
        phase : TimeSeries
            Values at the initial time are in the interval (-π, π], but may evolve to
            arbitrary real values.

        See Also
        --------
        numpy.angle
        numpy.unwrap
        arg

        """
        return TimeSeries(np.unwrap(self.arg, axis=self.time_axis), self.time)

    @property
    def norm(self):
        """Compute the L² norm of the waveform

        Returns
        -------
        n : TimeSeries

        See Also
        --------
        numpy.linalg.norm
        numpy.take
        norm2 : squared version of this

        Notes
        -----
        The integral of the (squared) magnitude of the data equals the sum of the
        (squared) magnitude of the modes for orthonormal basis functions, meaning that
        the L² norm of the function equals the basic Euclidean norm of the modes.  We
        assume that these modes are expanded in a band-limited but otherwise complete
        orthonormal basis.

        """
        return TimeSeries(np.linalg.norm(self, axis=self.modes_axis), self.time)

    @property
    def bar(self):
        """Return waveform modes of function representing conjugate of this function

        N.B.: This property is different from the `.conjugate` method; see below.

        See Also
        --------
        re : Return modes of function representing the real part of this function
        im : Return modes of function representing the imaginary part of this function

        Notes
        -----
        This property is different from the `.conjugate` (or `.conj`) method, in that
        `.conjugate` returns the conjugate of the mode weights of the function, whereas
        this property returns the mode weights of the conjugate of the function.  That
        is, `.conjugate` treats the data as a generic numpy array, and simply returns
        the complex conjugate of the raw data without considering what the data
        actually represents.  This property treats the data as a function represented
        by its mode weights.

        The resulting function has the negative spin weight of the input function.

        We have

            conjugate(f){s, l, m} = (-1)**(s+m) * conjugate(f{-s, l, -m})

        """
        return spherical.modes.algebra.bar(self)

    @property
    def re(self):
        """Return waveform modes of function representing real part of this function

        N.B.: This property is different from the `.real` method; see below.

        See Also
        --------
        im : Equivalent method for the imaginary part
        bar : Return modes of function representing the conjugate of this function

        Notes
        -----
        This property is different from the `.real` method, in that `.real` returns the
        real part of the mode weights of the function, whereas this property returns
        the mode weights of the real part of the function.  That is, `.real` treats the
        data as a generic numpy array, and simply returns the real part of the raw data
        without considering what the data actually represents.  This property treats
        the data as a function represented by its mode weights.

        Note that this only makes sense for functions of spin weight zero; taking the
        real part of functions with nonzero spin weight will depend too sensitively on
        the orientation of the coordinate system to make sense.  Therefore, this
        property raises a ValueError for other spins.

        The condition that a function `f` be real is that its modes satisfy

            f{l, m} = conjugate(f){l, m} = (-1)**(m) * conjugate(f{l, -m})

        [Note that conjugate(f){l, m} != conjugate(f{l, m}).]  As usual, we enforce
        that condition by essentially averaging the two modes:

            f{l, m} = (f{l, m} + (-1)**m * conjugate(f{l, -m})) / 2

        """
        return spherical.modes.algebra._real_func(self, False)

    @property
    def im(self):
        """Return waveform modes of function representing imaginary part of this function

        N.B.: This property is different from the `.imag` method; see below.

        See Also
        --------
        re : Equivalent method for the real part
        bar : Return modes of function representing the conjugate of this function

        Notes
        -----
        This property is different from the `.imag` method, in that `.imag` returns the
        imaginary part of the mode weights of the function, whereas this property
        returns the mode weights of the imaginary part of the function.  That is,
        `.imag` treats the data as a generic numpy array, and simply returns the
        imaginary part of the raw data without considering what the data actually
        represents.  This property treats the data as a function represented by its
        mode weights.

        Note that this only makes sense for functions of spin weight zero; taking the
        imaginary part of functions with nonzero spin weight will depend too
        sensitively on the orientation of the coordinate system to make sense.
        Therefore, this property raises a ValueError for other spins.

        The condition that a function `f` be imaginary is that its modes satisfy

            f{l, m} = -conjugate(f){l, m} = (-1)**(m+1) * conjugate(f{l, -m})

        [Note that conjugate(f){l, m} != conjugate(f{l, m}).]  As usual, we enforce
        that condition by essentially averaging the two modes:

            f{l, m} = (f{l, m} + (-1)**(m+1) * conjugate(f{l, -m})) / 2

        """
        return spherical.modes.algebra._imag_func(self, False)

    from spherical.modes.derivatives import (
        Lsquared, Lz, Lplus, Lminus,
        Rsquared, Rz, Rplus, Rminus,
        eth, ethbar
    )

    @property
    def eth_GHP(self):
        """Spin-raising derivative operator defined by Geroch-Held-Penrose

        The operator ð is defined in https://dx.doi.org/10.1063/1.1666410

        See Also
        --------
        eth : Related operator in the Newman-Penrose convention
        ethbar : Similar operator in the Newman-Penrose convention
        ethbar_GHP : Conjugate of this operator

        Notes
        -----
        We assume that the Ricci rotation coefficients satisfy β=β'=0, meaning that
        this operator equals the Newman-Penrose operator ð multiplied by 1/√2.

        """
        return self.eth / np.sqrt(2)

    @property
    def ethbar_GHP(self):
        """Spin-lowering derivative operator defined by Geroch-Held-Penrose

        The operator ð̄ is defined in https://dx.doi.org/10.1063/1.1666410

        See Also
        --------
        eth : Related operator in the Newman-Penrose convention
        ethbar : Similar operator in the Newman-Penrose convention
        eth_GHP : Conjugate of this operator

        Notes
        -----
        We assume that the Ricci rotation coefficients satisfy β=β'=0, meaning that
        this operator equals the Newman-Penrose operator ð̄ multiplied by 1/√2.

        """
        return self.ethbar / np.sqrt(2)

    def max_norm_index(self, skip_fraction_of_data=4):
        """Index of time step with largest norm

        The optional argument skips a fraction of the data.  The default is 4, which
        means that it skips the first 1/4 of the data, and only searches the last 3/4
        of the data for the max.  This must be strictly greater than 1, or the entire
        data is searched for the maximum of the norm.

        """
        if skip_fraction_of_data <= 1:
            return np.argmax(self.norm)
        else:
            i = int(self.n_times // skip_fraction_of_data)
            return np.argmax(self[i:].norm) + i

    def max_norm_time(self, skip_fraction_of_data=4):
        """Return time at which largest norm occurs in data

        See `help(max_norm_index)` for explanation of the optional argument.

        """
        return self.t[self.max_norm_index(skip_fraction_of_data=skip_fraction_of_data)]

    def interpolate(self, new_time, derivative_order=0, out=None):
        """Interpolate this object to a new set of times

        Note that if this object has "frame" data and the derivative order is nonzero,
        it is not entirely clear what is desired.  In those cases, the frame is just
        interpolated to the new times, but no derivative or antiderivative is taken.

        Parameters
        ----------
        new_time : array_like
            Points to evaluate the interpolant at
        derivative_order : int, optional
            Order of derivative to evaluate.  If negative, the antiderivative is
            returned.  Default value of 0 returns the interpolated data without
            derivatives or antiderivatives.  Must be between -3 and 3, inclusive.

        See Also
        --------
        scipy.interpolate.CubicSpline :
            The function that this function is based on.
        antiderivative :
            Calls this funtion with `new_time=self.time` and
            `derivative_order=-antiderivative_order` (defaulting to a single
            antiderivative).
        derivative :
            Calls this function `new_time=self.time` and
            `derivative_order=derivative_order` (defaulting to a single derivative).
        dot :
            Property calling `self.derivative(1)`.
        ddot :
            Property calling `self.derivative(2)`.
        int :
            Property calling `self.antiderivative(1)`.
        iint :
            Property calling `self.antiderivative(2)`.

        """
        result = TimeSeries.interpolate(self, new_time, derivative_order=derivative_order, out=out)
        if self.frame.shape == (self.n_times, 4) and not np.array_equal(self.time, result.time):
            self._metadata["frame"] = quaternionic.squad(self.frame, self.time, result.time)
        return result

    def truncate(self, tol=1e-10):
        """Truncate the precision of this object's data in place

        This function sets bits in the data to 0 when they have lower significance than
        will alter the norm of the Waveform by a fraction `tol` at that instant in
        time.

        """
        if tol != 0.0:
            tol_per_mode = tol / np.sqrt(self.n_modes)
            abs_tolerance = np.linalg.norm(self.ndarray, axis=self.modes_axis, keepdims=True) * tol_per_mode
            super().truncate(abs_tolerance)

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
        mode_plus = np.empty_like(self.ndarray[..., 0])
        mode_minus = np.empty_like(mode_plus)
        for ell in range(self.ell_min, self.ell_max + 1):
            for m in range(1, ell + 1):
                i_plus = self.index(ell, m)
                i_minus = self.index(ell, -m)
                mode_plus[:] = self.ndarray[..., i_plus]
                mode_minus[:] = self.ndarray[..., i_minus]
                self.ndarray[..., i_plus] = (mode_plus + np.conjugate(mode_minus)) / np.sqrt(2)
                self.ndarray[..., i_minus] = (mode_plus - np.conjugate(mode_minus)) / np.sqrt(2)

    def convert_from_conjugate_pairs(self):
        """Convert modes from conjugate-pair format in place

        This function reverses the effects of `convert_to_conjugate_pairs`.  See that
        function's docstring for details.

        """
        mode_plus = np.empty_like(self.ndarray[..., 0])
        mode_minus = np.empty_like(mode_plus)
        for ell in range(self.ell_min, self.ell_max + 1):
            for m in range(1, ell + 1):
                i_plus = self.index(ell, m)
                i_minus = self.index(ell, -m)
                mode_plus[:] = self.ndarray[..., i_plus]
                mode_minus[:] = self.ndarray[..., i_minus]
                self.ndarray[..., i_plus] = (mode_plus + mode_minus) / np.sqrt(2)
                self.ndarray[..., i_minus] = np.conjugate(mode_plus - mode_minus) / np.sqrt(2)

    def evaluate(self, *directions):
        """Evaluate waveform in a particular direction or set of directions

        Parameters
        ----------
        directions : array_like
            Directions of the observer relative to the source may be specified using
            the usual spherical coordinates, and an optional polarization angle (see
            Notes below).  These can be expressed as 2 or 3 floats (where the third is
            the polarization angle), or as an array with final dimension of size 2 or
            3.  Alternatively, the input may be a `quaternionic.array` (see Notes and
            arxiv.org/abs/1604.08140).  Input arrays can have multiple leading
            dimensions; the final dimension is always considered to hold the
            directions, and the other dimensions are retained in the output.

        Returns
        -------
        signal : array_like
            Note that this is complex-valued, meaning that it represents both
            polarizations.  To get the signal measured by a single detector, just take
            the real part.

        Notes
        -----
        To evaluate mode weights and obtain values, we need to evaluate spin-weighted
        spherical harmonics (SWSHs).  Though usually treated as functions of just the
        angles (θ, ϕ), a mathematically correct treatment (arxiv.org/abs/1604.08140)
        defines SWSHs as functions of the rotation needed to rotate the basis (x̂, ŷ, ẑ)
        onto the usual spherical-coordinate basis (θ̂, ϕ̂, n̂).  This function can take
        quaternionic arrays representing such a rotation directly, or the (θ, ϕ)
        coordinates themselves, and optionally the polarization angle ψ, which are
        automatically converted to quaternionic arrays.

        We define the spherical coordinates (θ, ϕ) such that θ is the polar angle
        (angle between the z axis and the point) and ϕ is the azimuthal angle (angle
        between x axis and orthogonal projection of the point into the x-y plane).
        This gives rise to the standard unit tangent vectors (θ̂, ϕ̂).

        We also define the polarization angle ψ as the angle through which we must
        rotate the vector θ̂ in a positive sense about n̂ to line up with the vector
        defining the legs of the detector.  If not given, this angle is treated as 0.

        Examples
        --------
        We can evaluate the signal in a single direction:

        >>> θ, ϕ, ψ = 0.1, 0.2, 0.3
        >>> w.evaluate(θ, ϕ)  # Default polarization angle
        >>> w.evaluate(θ, ϕ, ψ)  # Specified polarization angle

        Or we can evaluate in a set of directions:

        >>> w.evaluate([[θ, ϕ], [θ+0.4, ϕ], [θ+0.8, ϕ]])

        We can also evaluate on a more extensive set of directions.  Here, we construct
        an equi-angular grid to evaluate the waveform on (though irregular grids are
        also acceptable as long as you can pack them into a numpy array).

        >>> n_theta = n_phi = 2 * w.ell_max + 1
        >>> equiangular = np.array([
            [
                [theta, phi]
                for phi in np.linspace(0.0, 2*np.pi, num=n_phi, endpoint=False)
            ]
            for theta in np.linspace(0.0, np.pi, num=n_theta, endpoint=True)
        ])
        >>> w.evaluate(equiangular)

        """
        if len(directions) == 1:
            directions = directions[0]
        if isinstance(directions, quaternionic.array):
            R = directions
        else:
            directions = np.asarray(directions, dtype=float)
            if directions.shape[-1] == 2:
                R = quaternionic.array.from_spherical_coordinates(directions)
            elif directions.shape[-1] == 3:
                R = quaternionic.array.from_euler_angles(directions[..., 1], directions[..., 0], directions[..., 2])
            else:
                raise ValueError(
                    f"Input `directions` array must be quaternionic, or "
                    f"final dimension of must have size 2 or 3, not {directions.shape[-1]}"
                )

        # Compute the shape of the output
        modes_axis = self.modes_axis
        out_shape = (
            self.shape[:modes_axis]
            + R.shape[:-1]
            + tuple() if (modes_axis % self.ndim) == (-1 % self.ndim) else self.shape[modes_axis+1:]
        )

        # For now, we'll keep the new dimensions flat
        Rflat = R.ndarray.reshape(-1, 4)
        signal_shape = list(self.shape)
        signal_shape[modes_axis] = Rflat.shape[0]
        signal = np.zeros(tuple(signal_shape), dtype=complex)

        # Now, loop through, evaluating for each input R value
        wigner = spherical.Wigner(self.ell_max, ell_min=self.ell_min, mp_max=abs(self.spin_weight))
        sYlm = np.empty(wigner.Ysize, dtype=complex)
        slices = [slice(None) for _ in range(signal.ndim)]
        if np.array_equal(self.frame, np.atleast_2d(quaternionic.one)):  # frame is time-independent
            for i_R in range(Rflat.shape[0]):
                slices[modes_axis] = i_R
                wigner.sYlm(self.spin_weight, Rflat[i_R], out=sYlm)
                signal[tuple(slices)] = np.dot(self.ndarray, sYlm)
        else:  # Need to account for time-dependent frame
            data_slices = [slice(None) for _ in range(self.ndarray.ndim)]
            time_axis = self.time_axis
            for i_t in range(self.n_times):
                slices[time_axis] = i_t
                data_slices[time_axis] = i_t
                R_t = self.frame[i_t].inverse * Rflat
                for i_R, R_i in enumerate(Rflat):
                    slices[modes_axis] = i_R
                    wigner.sYlm(self.spin_weight, R_i, out=sYlm)
                    signal[tuple(slices)] = np.dot(self.ndarray[tuple(data_slices)], sYlm)

        return TimeSeries(signal.reshape(out_shape), self.time)

    # TODO:
    # expectation_value_L
    # # Don't bother with inner_product_LL, as it doesn't appear to be used; maybe a more general version?
    # inner_product
    # mode_frame : Minimally rotating O'Shaughnessy et al. frame
    # to_mode_frame
    # corotating_frame

    @property
    def expectation_value_LL(self):
        """Compute the matrix expectation value ⟨w|LᵃLᵇ|w⟩

        Here, Lᵃ is the usual angular-momentum operator familiar from quantum physics,
        and

            ⟨w|LᵃLᵇ|w⟩ = ℜ{Σₗₘₙ w̄ˡᵐ ⟨l,m|LᵃLᵇ|l,n⟩ wˡⁿ}

        This quantity is important for computing the angular velocity of a waveform.

        See Also
        --------
        expectation_value_Ldt
        angular_velocity

        """
        mode_weights = np.moveaxis(self.ndarray, self.modes_axis, -1)
        output_shape = mode_weights.shape[:-1] + (3, 3)
        mode_weights = mode_weights.reshape(-1, mode_weights.shape[-1])
        LL = np.zeros((mode_weights.shape[0], 3, 3), dtype=float)
        _expectation_value_LL(mode_weights, self.LM, LL)
        return LL.reshape(output_shape)

    @property
    def expectation_value_Ldt(self):
        """Compute the matrix expectation value ⟨w|Lᵃ∂ₜ|w⟩

        Here, Lᵃ is the usual angular-momentum operator familiar from quantum physics,
        ∂ₜ is the partial derivative with respect to time, and

            ⟨w|Lᵃ∂ₜ|w⟩ = ℑ{Σₗₘₙ w̄ˡᵐ ⟨l,m|Lᵃ|l,n⟩ ∂ₜwˡⁿ}

        This quantity is important for computing the angular velocity of a waveform.

        See Also
        --------
        expectation_value_LL
        angular_velocity

        """
        mode_weights = np.moveaxis(self.ndarray, self.modes_axis, -1)
        output_shape = mode_weights.shape[:-1] + (3,)
        mode_weights = mode_weights.reshape(-1, mode_weights.shape[-1])
        mode_weights_dot = np.moveaxis(self.dot.ndarray, self.modes_axis, -1)
        mode_weights_dot = mode_weights_dot.reshape(-1, mode_weights_dot.shape[-1])
        Ldt = np.zeros((mode_weights.shape[0], 3), dtype=float)
        _expectation_value_Ldt(mode_weights, mode_weights_dot, self.LM, Ldt)
        return Ldt.reshape(output_shape)

    @property
    def angular_velocity(self):
        """Angular velocity of waveform

        This function calculates the angular velocity of a WaveformModes object from
        its modes — essentially, the angular velocity of the rotating frame in which
        the time dependence of the modes is minimized.  This was introduced in Sec. II
        of "Angular velocity of gravitational radiation and the corotating frame"
        <http://arxiv.org/abs/1302.2919>.

        It can be calculated in terms of the expectation values ⟨w|Lᵃ∂ₜ|w⟩ and
        ⟨w|LᵃLᵇ|w⟩ according to the relation

            ⟨w|LᵇLᵃ|w⟩ ωₐ = -⟨w|Lᵇ∂ₜ|w⟩

        For each set of modes (e.g., at each instant of time), this is a simple linear
        equation in 3 dimensions to be solved for ω.

        See Also
        --------
        expectation_value_LL
        expectation_value_Ldt

        """
        # Calculate the <L∂ₜ> vector and <LL> matrix at each instant
        ldt = self.expectation_value_Ldt
        ll = self.expectation_value_LL

        # Solve ⟨w|LᵇLᵃ|w⟩ ωₐ = -⟨w|Lᵇ∂ₜ|w⟩ for ω
        ω = -np.linalg.solve(ll, ldt)

        return ω

    def boost(self, v⃗, ell_max):
        """Find modes of waveform boosted by velocity v⃗

        Implements Equation (21) of arxiv.org/abs/1509.00862

        Parameters
        ----------
        v⃗ : array_like
            Three-vector representing the velocity of the boosted frame relative to the
            inertial frame, in units where the speed of light is 1
        ell_max : int
            Maximum value of `ell` to use while computing the transformation, and to
            provide in the returned object.  See Notes, below.

        Returns
        -------
        wprime : WaveformModes
            Modes of waveform measured in boosted frame or of modes from boosted source
            measured in original frame.  This should have the same properties as the
            input waveform, except with (1) different time data [see Notes, below], (2)
            a minimum `ell` value of 0 even for spin weight other than 0, and (3) a
            maximum `ell` value of `ell_max`.

        Notes
        -----
        Due to the nature of the transformation, some of the information in the input
        waveform must be discarded, because it corresponds to slices of the output
        waveform that are not completely represented in the input.  Thus, the times of
        the output waveform will not just be the Lorentz-transformed times of the input
        waveform.

        Depending on the magnitude β=|v⃗|, a very large value of `ell_max` may be
        needed.  The dominant factor is the translation that builds up over time:
        `β*T`, where `T` is the largest time found in the waveform.  For example, if
        β*T ≈ 1000M, we might need `ell_max=64` to maintain a comparable accuracy as in
        the input data.

        Because of the `β*T` effects, it is usually best to set t=0 at the merger time
        — best approximated as `self.max_norm_time()`.  The largest translation is then
        found early in the waveform, when the waveform is changing slowly.

        """
        from .transformations import boost
        return boost(self, v⃗, ell_max)

    def rotate(self, quat):
        """Rotate decomposition basis of modes represented by this waveform

        This returns a new waveform object, with missing "frame" data.

        Parameters
        ----------
        quat : quaternionic.array
            This must have one quaternion or the same number of quaternions as the
            number of times in the waveform.

        """
        from spherical.wigner import _rotate

        R = quaternionic.array(quat)
        wigner = spherical.Wigner(self.ell_max, ell_min=self.ell_min)  #, mp_max=abs(self.spin_weight))
        D = np.zeros(wigner.Dsize, dtype=complex)
        mode_weights = self.ndarray
        rotated_mode_weights = np.zeros_like(mode_weights)
        mode_weights = np.moveaxis(mode_weights, self.modes_axis, -1)
        rotated_mode_weights = np.moveaxis(rotated_mode_weights, self.modes_axis, -1)
        shape = rotated_mode_weights.shape
        if quat.shape == (4,) or quat.shape == (1, 4):
            wigner.D(R, out=D)
            mode_weights = mode_weights.reshape(-1, mode_weights.shape[-1])
            rotated_mode_weights = rotated_mode_weights.reshape(-1, mode_weights.shape[-1])
            _rotate(
                mode_weights, rotated_mode_weights,
                wigner.ell_min, wigner.ell_max, wigner.mp_max,
                self.ell_min, self.ell_max, self.spin_weight,
                D
            )
        elif quat.shape == (self.n_times, 4):
            slices = [slice(None) for _ in range(self.ndim)]
            time_axis = self.time_axis if self.time_axis < self.modes_axis else self.time_axis - 1
            for i_t in range(self.n_times):
                wigner.D(R[i_t], out=D)
                slices[time_axis] = i_t
                s = tuple(slices)
                m = mode_weights[s]
                r = rotated_mode_weights[s]
                m = m.reshape(-1, m.shape[-1])
                r = r.reshape(-1, r.shape[-1])
                _rotate(
                    m, r,
                    wigner.ell_min, wigner.ell_max, wigner.mp_max,
                    self.ell_min, self.ell_max, self.spin_weight,
                    D
                )
        else:
            raise ValueError(
                f"Quaternionic array shape {R.shape} not understood; expected {(4,)}, {(1, 4)}, or {(self.n_times, 4)}"
            )
        rotated_mode_weights = rotated_mode_weights.reshape(shape)
        rotated_mode_weights = np.moveaxis(rotated_mode_weights, -1, self.modes_axis)
        new_metadata = self._metadata.copy()
        new_metadata.pop("frame", None)
        return type(self)(rotated_mode_weights, **new_metadata)

    def to_inertial_frame(self):
        """Return a copy of this waveform in the inertial frame"""
        if "frame" not in self._metadata:
            raise ValueError("This waveform has no frame information")
        if self.frame.shape[0] == 1:
            raise ValueError("This waveform appears to already be in an inertial frame")
        if self.frame.shape != (self.n_times, 4):
            raise ValueError(f"Frame shape {frame.shape} not understood; expected {(self.n_times, 4)}")
        w = self.rotate(~self.frame)
        w._metadata["frame_type"] = "inertial"
        return w

    def to_corotating_frame(
        self, R0=None, tolerance=1e-12, z_alignment_region=None, return_omega=False, truncate_log_frame=False
    ):
        """Return a copy of this waveform in the corotating frame

        The corotating frame is defined to be a rotating frame for which the (L² norm
        of the) time-dependence of the modes expressed in that frame is minimized.
        This leaves the frame determined only up to an overall rotation.  In this 

        Parameters
        ----------
        R0 : quaternionic, optional
            Initial value of frame when integrating angular velocity.  Defaults to the
            identity.
        tolerance : float, optional
            Absolute tolerance used in integration of angular velocity
        z_alignment_region : {None, 2-tuple of floats}, optional
            If not None, the dominant eigenvector of the <LL> matrix is aligned with
            the z axis, averaging over this portion of the data.  The first and second
            elements of the input are considered fractions of the inspiral at which to
            begin and end the average.  For example, (0.1, 0.9) would lead to starting
            10% of the time from the first time step to the max norm time, and ending
            at 90% of that time.
        return_omega : bool, optional
            If True, return a 2-tuple consisting of the waveform in the corotating
            frame (the usual returned object) and the angular-velocity data.  That is
            frequently also needed, so this is just a more efficient way of getting the
            data.
        truncate_log_frame : bool, optional
            If True, set bits of log(frame) with lower significance than `tolerance` to
            zero, and use exp(truncated(log(frame))) to rotate the waveform.  Also
            returns `log_frame` along with the waveform (and optionally `omega`)

        """
        raise NotImplementedError()
        frame, omega = corotating_frame(
            self, R0=R0, tolerance=tolerance, z_alignment_region=z_alignment_region, return_omega=True
        )
        if truncate_log_frame:
            log_frame = np.log(frame).ndarray
            power_of_2 = 2 ** int(-np.floor(np.log2(2 * tolerance)))
            log_frame = np.round(log_frame * power_of_2) / power_of_2
            frame = np.exp(quaternionic.array(log_frame))
        w = self.rotate_decomposition_basis(frame)
        w._metadata["frame_type"] = "corotating"
        w._metadata["frame"] = frame
        if return_omega:
            if truncate_log_frame:
                return (w, omega, log_frame)
            else:
                return (w, omega)
        else:
            if truncate_log_frame:
                return (w, log_frame)
            else:
                return w
