"""Interfaces for horizon data

This module provides two basic container classes that provide a uniform
interface to data describing horizons in numerical-relativity simulations:

    HorizonQuantities
        Container for various TimeSeries related to an individual horizon
    Horizons
        Container for several HorizonQuantities related to a merger

There are also functions for saving and loading horizon data in standardized
formats:

    spec_horizons_h5.save
    spec_horizons_h5.load
    xor_multishuffle_bzip2.save
    xor_multishuffle_bzip2.load

"""

from . import spec_horizons_h5, xor_multishuffle_bzip2
xmb = xor_multishuffle_bzip2

formats = {
    None: spec_horizons_h5,
    "nrar": spec_horizons_h5,
    "spec_horizons_h5": spec_horizons_h5,
    "xmb": xmb,
    "xor_multishuffle_bzip2": xor_multishuffle_bzip2,
}


def load(file, **kwargs):
    """Load horizon data from an SXS-format file

    This is the highest-level loader for horizon data, and should generally be
    preferred over more specific loaders like those in `spec_horizons_h5` or
    `xor_multishuffle_bzip2`; this function automatically detects the format
    and dispatches as needed.

    Parameters
    ----------
    file : file-like object, string, or pathlib.Path
        Path to the file on disk or a file-like object (such as an open file
        handle) to be opened by h5py.File.

    Returns
    -------
    horizons : sxs.Horizons
        This is a container for the horizon objects.  See Notes below.

    See also
    --------
    sxs.Horizons : Container object for all of the horizons
    sxs.HorizonQuantities : Container objects for each of the horizons

    Notes
    -----
    The returned object can be indexed just like the original SpEC-format
    Horizons.h5 file:

        horizons["AhA.dir/CoordCenterInertial.dat"]

    However, the `horizons` object also has a more natural and general interface
    that should be preferred for compatibility with other formats, in which the
    same vector-valued function of time can be accessed as

        horizons.A.coord_center_inertial

    See the documentation of `Horizons` and `HorizonQuantities` for more details.

    """
    raise NotImplementedError()
    horizons.full_path = full_path
    return horizons


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
        The Christodoulou mass `mᵪ` of the horizon is related to the irreducible
        mass `mᵢ` and the spin magnitude `s` by `mᵪ² = mᵢ² + s²/4mᵢ²`.
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
    as attributes, along with the following

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
        from .. import TimeSeries
        self.time = kwargs["time"]
        self.areal_mass = TimeSeries(kwargs["areal_mass"], time=self.time, time_axis=0)
        self.christodoulou_mass = TimeSeries(kwargs["christodoulou_mass"], time=self.time, time_axis=0)
        self.coord_center_inertial = TimeSeries(kwargs["coord_center_inertial"], time=self.time, time_axis=0)
        self.dimensionful_inertial_spin = TimeSeries(kwargs["dimensionful_inertial_spin"], time=self.time, time_axis=0)
        self.chi_inertial = TimeSeries(kwargs["chi_inertial"], time=self.time, time_axis=0)

    @property
    def dimensionful_inertial_spin_mag(self):
        import numpy as np
        return np.linalg.norm(self.dimensionful_inertial_spin, axis=1)

    @property
    def chi_inertial_mag(self):
        import numpy as np
        return np.linalg.norm(self.chi_inertial, axis=1)

    chi_mag_inertial = chi_inertial_mag  # backwards-compatible alias because the name is inconsistent

    def __getitem__(self, key):
        import numpy as np
        import inflection
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
    be indexed in exactly the same way as a Horizons.h5 file.

    """
    def __init__(self, **kwargs):
        self.A = kwargs.pop("A", None)
        self.B = kwargs.pop("B", None)
        self.C = kwargs.pop("C", None)

    @property
    def a(self):
        return self.A

    @property
    def b(self):
        return self.B

    @property
    def c(self):
        return self.C

    def __getitem__(self, key):
        shorter_key = key.replace("AhA.dir", "A").replace("AhB.dir", "B").replace("AhC.dir", "C")
        if shorter_key.upper() in "ABC":
            return getattr(self, shorter_key.upper())
        elif len(shorter_key.split("/", maxsplit=1)) == 2:
            horizon, subkey = shorter_key.split("/", maxsplit=1)
            return getattr(self, horizon)[subkey]
        else:
            raise ValueError(f"Cannot find key '{key}' in this {type(self).__name__} object")
