


class HorizonQuantities(object):
    import numpy as np

    def __init__(self, **kwargs):
        from . import TimeSeries
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

    This is a fairly small container to provide an interface for several
    HorizonQuantities objects, which can be accessed in several ways.  Up to
    three horizons are supported, and are named A, B, and C.  Typically A and B
    will represent objects in a merging binary and C will represent the
    remnant, though this convention is not enforced.  If this object is named
    `horizons`, each individual horizon can be accessed in any of these ways:

        horizons.A
        horizons.a
        horizons["A"]
        horizons["a"]
        horizons["AhA.dir"]

    and similarly for B and C.  These different ways of accessing `A` are
    essentially aliases; they return precisely the same object.  In addition,
    the attributes of those horizon objects are passed through — for example,
    as

        horizons.A.time
        horizons["A/time"]

    to axis the time data for horizon A, or

        horizons.A.coord_center_inertial
        horizons["A/coord_center_inertial"]
        horizons["AhA.dir/CoordCenterInertial"]
        horizons["AhA.dir/CoordCenterInertial.dat"]

    Note that these are equivalent and return precisely the same thing, except
    for the last one.  If an attribute ends with ".dat", the returned quantity
    will be the quantity that appears in a SpEC-format Horizons.h5 file, which
    is "stacked" with the time data.  That is, rather than being an Nx3
    vector-valued function of time, when this attribute ends with ".dat" it
    returns an Nx4 array.  For scalar-valued functions of time, the returned
    object has shape Nx2, rather than just N.

    """
    import numpy as np

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
