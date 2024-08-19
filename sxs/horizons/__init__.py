"""Interfaces for horizon data

This module provides two basic container classes that provide a uniform
interface to data describing horizons in numerical-relativity simulations, as
well as functions for saving and loading horizon data in standardized formats.

"""

import numpy as np
import inflection
from . import spec_horizons_h5, xor_multishuffle_bzip2
from .. import TimeSeries

xmb = xor_multishuffle_bzip2

formats = {
    None: spec_horizons_h5,
    "": spec_horizons_h5,
    "nrar": spec_horizons_h5,
    "spec_horizons_h5": spec_horizons_h5,
    "xmb": xmb,
    "xor_multishuffle_bzip2": xor_multishuffle_bzip2,
}


class HorizonQuantities(object):
    """Container object for various TimeSeries related to an individual horizon

    Parameters
    ----------
    time : (N,) array_like
        Times at which the horizon quantities are measured.
    areal_mass : (N,) array_like
        The areal (or irreducible) mass of the horizon, defined as the square-root
        of its surface area divided by 16π, where area is measured as a function of
        time.
    christodoulou_mass : (N,) array_like
        The Christodoulou mass `mᵪ` of the horizon is related to the areal (or
        irreducible) mass `mₐ` and the dimensionful spin magnitude `s` by the
        expression `mᵪ² = mₐ² + s²/4mₐ²`.
    coord_center_inertial : (N, 3) array_like
        Cartesian coordinates of the center of the apparent horizon, in the
        "inertial frame," the asymptotically inertial frame in which the
        gravitational waves are measured.
    dimensionful_inertial_spin : (N, 3) array_like
        Cartesian vector components of the spin angular momentum measured on the
        apparent horizon in the "inertial frame".
    chi_inertial : (N, 3) array_like
        Cartesian components of the spin angular momentum measured in the "inertial
        frame", made dimensionless by dividing by the square of the Christodoulou
        mass.

    Attributes
    ----------
    All of the above parameters are converted to TimeSeries objects and accessible
    as attributes, along with the following:

    dimensionful_inertial_spin_mag : (N,) TimeSeries
        Euclidean norm of the `dimensionful_inertial_spin` quantity
    chi_inertial_mag : (N,) TimeSeries
        Euclidean norm of the `chi_inertial` quantity
    chi_mag_inertial : (N,) TimeSeries
        Euclidean norm of the `chi_inertial` quantity (alias for the above,
        maintained for backwards compatibility)

    Methods
    -------
    __getitem__(key)
        Another interface for accessing the attributes, with more flexibility.  See
        below for explanation.

    Notes
    -----
    In addition to the standard attribute access, as in

        horizon.coord_center_inertial

    it is also possible to access that attribute equivalently via indexing as

        horizon["coord_center_inertial"]
        horizon["CoordCenterInertial"]

    Here, the latter is simply an alias for the former.  For backwards
    compatibility, it is also possible to access the attributes "horizontally
    stacked" (via `np.hstack`) with the time data.  That is, rather than being an
    Nx3 vector-valued function of time as returned above, when an attribute ends
    with ".dat" it returns an Nx4 array:

       horizon["CoordCenterInertial.dat"]

    The result of that call can be sliced as `[:, 0]` to access the time data and
    `[:, 1:]` to access the 3-d vector as a function of time.  Together with
    related behavior in the `Horizons` class, this provides full backward
    compatibility with SpEC-format Horizons.h5 files, in the sense that a
    `Horizons` object can be indexed in exactly the same way as a Horizons.h5 file.

    """

    def __init__(self, **kwargs):
        self.time = kwargs["time"]
        self.areal_mass = TimeSeries(kwargs["areal_mass"], time=self.time, time_axis=0)
        self.christodoulou_mass = TimeSeries(kwargs["christodoulou_mass"], time=self.time, time_axis=0)
        self.coord_center_inertial = TimeSeries(kwargs["coord_center_inertial"], time=self.time, time_axis=0)
        self.dimensionful_inertial_spin = TimeSeries(kwargs["dimensionful_inertial_spin"], time=self.time, time_axis=0)
        self.chi_inertial = TimeSeries(kwargs["chi_inertial"], time=self.time, time_axis=0)

    @property
    def dimensionful_inertial_spin_mag(self):
        return np.linalg.norm(self.dimensionful_inertial_spin, axis=1)

    @property
    def chi_inertial_mag(self):
        return np.linalg.norm(self.chi_inertial, axis=1)

    chi_mag_inertial = chi_inertial_mag  # backwards-compatible alias because the name is inconsistent

    def __getitem__(self, key):
        dat = key.endswith(".dat")
        standardized_key = inflection.underscore(key.replace(".dat", ""))
        attribute = getattr(self, standardized_key)
        if dat:
            if getattr(attribute, "ndim", 0) == 1:
                attribute = attribute[:, np.newaxis]
            if hasattr(attribute, "ndarray"):
                attribute = attribute.ndarray
            return np.hstack((self.time[:, np.newaxis], attribute))
        return attribute


class Horizons(object):
    """Container object for several HorizonQuantities objects

    Parameters
    ----------
    A, B, C : HorizonQuantities, optional
        If these are not given, they will be `None`.

    Attributes
    ----------
    A, B, C, a, b, c : {HorizonQuantities, None}
        The lowercase versions are simply aliases for the uppercase ones.

    Methods
    -------
    __getitem__(key)
        Indexes the individual HorizonQuantities objects, and optionally passes
        indexes through to the underlying object.  See below for explanation.

    See also
    --------
    HorizonQuantities : Containers for data pertaining to each of the horizons

    Notes
    -----
    This is a small container to provide an interface for several HorizonQuantities
    objects, which can be accessed in several ways.  Up to three horizons are
    supported, and are named A, B, and C.  Typically A and B will represent objects
    in a merging binary and C will represent the remnant, though this convention is
    not enforced.  It is expected that components that do not have horizons (e.g.,
    neutron stars) will be represented as `None` rather than HorizonQuantities
    objects.  If this object is named `horizons`, each individual horizon can be
    accessed in any of these ways:

        horizons.A
        horizons.a
        horizons["A"]
        horizons["a"]
        horizons["AhA.dir"]

    and similarly for B and C.  These different ways of accessing `A` are
    essentially aliases; they return precisely the same object.  In addition, the
    attributes of those horizon objects are passed through — for example, as

        horizons.A.time
        horizons["A/time"]

    to access the time data for horizon A, or

        horizons.A.coord_center_inertial
        horizons["A/coord_center_inertial"]
        horizons["AhA.dir/CoordCenterInertial"]
        horizons["AhA.dir/CoordCenterInertial.dat"]

    These are equivalent and return precisely the same thing, except for the last
    one.  If an attribute ends with ".dat", the returned quantity will be the
    quantity that appears in a SpEC-format Horizons.h5 file, which is "horizontally
    stacked" (via `np.hstack`) with the time data.  That is, rather than being an
    Nx3 vector-valued function of time, when this attribute ends with ".dat" it
    returns an Nx4 array.  For scalar-valued functions of time, the returned object
    has shape Nx2, rather than just N.  This provides full backward compatibility
    with SpEC-format Horizons.h5 files, in the sense that a `Horizons` object can
    be indexed in exactly the same way as a Horizons.h5 file.  Also note that the
    function `sxs.loadcontext` provides a context manager just like `h5py`:

        with sxs.loadcontext("Horizons.h5") as f:
            time = f["AhA.dir/ArealMass.dat"][:, 0]
            areal_mass = f["AhA.dir/ArealMass.dat"][:, 1]

    This code is identical to the equivalent code using `h5py` except that the call
    to `h5py.File` is replaced with the call to `sxs.loadcontext`.  The `.dat`
    datasets are re-computed on the fly.

    """
    def __init__(self, **kwargs):
        self.A = kwargs.pop("A", None)
        self.B = kwargs.pop("B", None)
        self.C = kwargs.pop("C", None)

    def __getitem__(self, key):
        shorter_key = key.replace("AhA.dir", "A").replace("AhB.dir", "B").replace("AhC.dir", "C")
        if shorter_key.upper() in "ABC":
            return getattr(self, shorter_key.upper())
        elif len(shorter_key.split("/", maxsplit=1)) == 2:
            horizon, sub_key = shorter_key.split("/", maxsplit=1)
            return getattr(self, horizon)[sub_key]
        else:
            raise ValueError(f"Cannot find key '{key}' in this {type(self).__name__} object")

    @property
    def a(self):
        return self.A

    @property
    def b(self):
        return self.B

    @property
    def c(self):
        return self.C

    @property
    def newtonian_com(self):
        """Newtonian center of mass as function of time

        This returns only the center of mass of the binary components; the center of
        mass of the common horizon is just `horizons.C.coord_center_inertial`.

        Returns
        -------
        com : ndarray
            This has shape (self.A.n_times, 3), representing the components of the
            vector as a function of time.

        See Also
        --------
        average_com_motion : fit uniform motion to this result

        Notes
        -----
        This just evaluates the simple formula

            (m_A * x_A + m_B * x_B) / (m_A + m_B)

        where the masses are the respective Christodoulou masses, and the positions are
        taken from the `coord_center_inertial` properties of the respective horizons.
        This is highly susceptible to the vagaries of gauge, and must always be taken
        with plentiful grains of salt.

        """
        m_A = self.A.christodoulou_mass[:, np.newaxis]
        x_A = self.A.coord_center_inertial
        m_B = self.B.christodoulou_mass[:, np.newaxis]
        x_B = self.B.coord_center_inertial
        m = m_A + m_B
        com = ((m_A * x_A) + (m_B * x_B)) / m
        return com

    def average_com_motion(self, skip_beginning_fraction=0.01, skip_ending_fraction=0.10):
        """Fit uniform motion to measured Newtonian center of mass

        Parameters
        ----------
        skip_beginning_fraction : float, optional
            Exclude this portion from the beginning of the data.  Note that this is
            a fraction, rather than a percentage.  The default value is 0.01,
            meaning the first 1% of the data will be ignored.
        skip_ending_fraction : float, optional
            Exclude this portion from the end of the data.  Note that this is a
            fraction, rather than a percentage.  The default value is 0.10, meaning
            the last 10% of the data will be ignored.

        Returns
        -------
        x_i : length-3 array of floats
            Best-fit initial position of the center of mass
        v_i : length-3 array of floats
            Best-fit initial velocity of the center of mass
        t_i : float
            Initial time used.  This is determined by the `skip_beginning_fraction`
            input parameter.
        t_f : float
            Final time used.  This is determined by the `skip_ending_fraction` input
            parameter.

        See Also
        --------
        newtonian_com : measured quantity as function of time

        Notes
        -----
        See the docstring of `newtonian_com` for some relevant caveats.  The
        translation to be applied to the data should be calculated given the values
        returned by this function as

            com_average = sxs.TimeSeries(
                x_i[np.newaxis] + v_i[np.newaxis] * horizons.A.time[:, np.newaxis],
                horizons.A.time
            )

        """
        from scipy.integrate import simpson

        t = self.A.time
        com = self.newtonian_com

        # We will be skipping the beginning and end of the data;
        # this gives us the initial and final indices
        t_i, t_f = t[0] + (t[-1] - t[0]) * skip_beginning_fraction, t[-1] - (t[-1] - t[0]) * skip_ending_fraction
        i_i, i_f = np.argmin(np.abs(t - t_i)), np.argmin(np.abs(t - t_f))

        # Find the optimum analytically
        com_0 = simpson(com[i_i : i_f + 1], x=t[i_i : i_f + 1], axis=0)
        com_1 = simpson((t[:, np.newaxis] * com)[i_i : i_f + 1], x=t[i_i : i_f + 1], axis=0)
        x_i = 2 * (com_0 * (2 * t_f ** 3 - 2 * t_i ** 3) + com_1 * (-3 * t_f ** 2 + 3 * t_i ** 2)) / (t_f - t_i) ** 4
        v_i = 6 * (com_0 * (-t_f - t_i) + 2 * com_1) / (t_f - t_i) ** 3

        return x_i, v_i, t_i, t_f

    @property
    def n⃗(self):
        """Vector pointing from horizon A to horizon B

        This function can be spelled `n⃗`, `nvec`, or `separation`, interchangeably.

        Returns
        -------
        n⃗ : ndarray
            This has shape (self.A.n_times, 3), representing the components of the
            vector as a function of time.

        See Also
        --------
        n̂, nhat : Normalized version of this vector
        λ̂, lambdahat : Normalized time-derivative of n̂
        ℓ̂, ellhat : Normalized angular-velocity vector of n̂

        """
        return self.B.coord_center_inertial - self.A.coord_center_inertial

    nvec = n⃗
    separation = n⃗

    @property
    def n̂(self):
        """Unit vector pointing from horizon A to horizon B

        This function can be spelled `n̂` or `nhat`, interchangeably.

        Returns
        -------
        n̂ : ndarray
            This has shape (self.A.n_times, 3), representing the components of the
            vector as a function of time.

        See Also
        --------
        n⃗, nvec, separation : Non-normalized version of this vector
        λ̂, lambdahat : Normalized time-derivative of n̂
        ℓ̂, ellhat : Normalized angular-velocity vector

        Notes
        -----
        Note that (n̂, λ̂, ℓ̂) forms a right-handed frame, which is commonly used in
        post-Newtonian theory and similar treatments.

        """
        n⃗ = self.separation
        return n⃗ / np.linalg.norm(n⃗, axis=1)[:, np.newaxis]

    nhat = n̂

    @property
    def λ⃗(self):
        """Time-derivative of separation vector

        This function can be spelled `λ⃗` or `lambdavec`, interchangeably.

        Returns
        -------
        λ⃗ : ndarray
            This has shape (self.A.n_times, 3), representing the components of the
            vector as a function of time.

        See Also
        --------
        n⃗, nvec, separation : (Non-normalized) separation vector between two horizons
        n̂, nhat : Normalized separation vector
        λ̂, lambdahat : Normalized version of this vector
        ℓ̂, ellhat : Normalized angular-velocity vector

        """
        return self.n⃗.dot
    
    lambdavec = λ⃗

    @property
    def λ̂(self):
        """Time-derivative of normalized separation vector

        This function can be spelled `λ̂` or `lambdahat`, interchangeably.

        Returns
        -------
        λ̂ : ndarray
            This has shape (self.A.n_times, 3), representing the components of the
            vector as a function of time.

        See Also
        --------
        n⃗, nvec, separation : (Non-normalized) separation vector between two horizons
        n̂, nhat : Normalized separation vector
        ℓ̂, ellhat : Normalized angular-velocity vector

        Notes
        -----
        Note that (n̂, λ̂, ℓ̂) forms a right-handed frame, which is commonly used in
        post-Newtonian theory and similar treatments.

        """
        λ⃗ = self.λ⃗
        return λ⃗ / np.linalg.norm(λ⃗, axis=1)[:, np.newaxis]

    lambdahat = λ̂

    @property
    def ℓ̂(self):
        """Normalized angular-velocity vector

        This function can be spelled `ℓ̂` or `ellhat`, interchangeably.

        Returns
        -------
        ℓ̂ : ndarray
            This has shape (self.A.n_times, 3), representing the components of the
            vector as a function of time.

        See Also
        --------
        n⃗, nvec, separation : (Non-normalized) separation vector between two horizons
        n̂, nhat : Normalized separation vector
        λ̂, lambdahat : Normalized time-derivative of n̂

        Notes
        -----
        Note that (n̂, λ̂, ℓ̂) forms a right-handed frame, which is commonly used in
        post-Newtonian theory and similar treatments.

        """
        return TimeSeries(np.cross(self.n̂, self.λ̂), time=self.n̂.time)

    ellhat = ℓ̂

    @property
    def R(self):
        """Frame Rotor taking (x̂, ŷ, ẑ) onto (n̂, λ̂, ℓ̂) at each instant
        
        For example, if λ̂ᵢ is the value of λ̂ at time tᵢ, and Rᵢ the
        corresponding output from this function, then we have
        
            λ̂ᵢ = Rᵢ * quaternionic.y / Rᵢ

        Returns
        -------
        R : quaternionic.array
            This has shape (self.A.n_times, 4), representing the rotor
            at each time.  Note that this is *not* a TimeSeries
            object, as is returned by several other functions in this
            class.  However, the corresponding times are available as
            `R.time`.
        """
        import numpy as np
        import quaternionic
        reference_frame = np.array([quaternionic.x.vector, quaternionic.y.vector, quaternionic.z.vector])
        target_frame = np.array([self.n̂, self.λ̂, self.ℓ̂])
        R = quaternionic.array([
            quaternionic.align(
                target_frame[:, i, :],
                reference_frame
            )
            for i in range(target_frame.shape[1])
        ])
        quaternionic.unflip_rotors(R, inplace=True)
        R.time = self.A.time
        return R

    @property
    def Ω⃗(self):
        """Angular velocity vector of the binary

        This function can be spelled `Ω⃗`, `ω⃗`, `OmegaVec`, or
        `omegaVec`, interchangeably.

        Returns
        -------
        Ω⃗ : TimeSeries
            This represents the angular velocity as a function of time.
        """
        R = self.R
        return TimeSeries(R.to_angular_velocity(R.time), time=R.time)
    
    ω⃗ = Ω⃗
    OmegaVec = Ω⃗
    omegaVec = Ω⃗

    @property
    def Ω(self):
        """Magnitude of the angular velocity of the binary

        This function can be spelled `Ω`, `ω`, `Omega`, or `omega`,
        interchangeably.

        Returns
        -------
        Ω : TimeSeries
            This represents the magnitude of the angular velocity as a
            function of time.
        """
        return TimeSeries(np.linalg.norm(self.Ω⃗, axis=1), time=self.A.time)
    
    ω = Ω
    Omega = Ω
    omega = Ω
