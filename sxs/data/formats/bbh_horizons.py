

def save():
    raise NotImplementedError()


def load(file_name):
    """Load SpEC Horizons.h5 file as `Horizons` object

    Parameters
    ----------
    file_name : str or file-like
        Name of the file on disk, or file-like object (such as an open file
        handle) to be opened by h5py.File.

    Returns
    -------
    horizons : sxs.data.Horizons object
        This is a container for the horizon objects.  See Notes below.

    Notes
    -----
    The returned object can be indexed just like the original h5 file — as in

        horizons["AhA.dir/CoordCenterInertial.dat"]

    This returns an array of shape Nx4, where N corresponds to the number of
    time steps, and each time step contains the time itself and the three
    components of the coordinate vector.  Alternatively, the object can be
    accessed more naturally, so that `horizons.A` returns an object
    encapsulating all the data for horizon A

        horizons.A  # returns an object encapsulating data for horizon A
        horizons["A"]  # equivalent to the above
        horizons.A.time  # array of all the time values
        horizons.A["time"]  # equivalent to the above
        horizons["A/time"]  # equivalent to the above
        horizons.A.coord_center_inertial  # Nx3 array of coordinate values
        horizons.A["coord_center_inertial"]  # equivalent to the above

    Note that, unless the requested quantity ends with ".dir", the returned
    values do not come with the redundant time information; you are expected to
    access that separately.

    """
    import h5py
    from .. import Horizons, HorizonQuantities

    hqs = {}
    with h5py.File(file_name, "r") as f:
        for horizon_name in "ABC":
            dir_name = f"Ah{horizon_name}.dir"
            if dir_name in f:
                g = f[dir_name]
                horizon = HorizonQuantities(
                    time=g["ArealMass.dat"][:, 0],
                    areal_mass=g["ArealMass.dat"][:, 1],
                    christodoulou_mass=g["ChristodoulouMass.dat"][:, 1],
                    coord_center_inertial=g["CoordCenterInertial.dat"][:, 1:],
                    dimensionful_inertial_spin=g["DimensionfulInertialSpin.dat"][:, 1:],
                    chi_inertial=g["chiInertial.dat"][:, 1:],
                )
                hqs[horizon_name] = horizon
    return Horizons(**hqs)
