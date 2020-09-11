"""File I/O for SpEC-format Horizons.h5 files"""


def save(horizons, file):
    """Save `Horizons` object as SpEC Horizons.h5 file

    Parameters
    ----------
    horizons : sxs.Horizons
        Horizons object to be written to file.
    file : file-like object, string, or pathlib.Path
        Path to the file on disk or a file-like object (such as an open file
        handle) to be written by h5py.File.

    See Also
    --------
    sxs.horizons.spec_horizons_h5.load : load the output file format

    """
    import h5py

    with h5py.File(file, "w") as f:
        f.attrs["sxs_format"] = "horizons.spec_horizons_h5"
        for horizon_name in "ABC":
            horizon = horizons[horizon_name]
            if horizon is not None:
                g = f.create_group(f"Ah{horizon_name}.dir")
                for dat in [
                    "ArealMass.dat", "ChristodoulouMass.dat", "DimensionfulInertialSpinMag.dat", "chiMagInertial.dat",
                    "CoordCenterInertial.dat", "DimensionfulInertialSpin.dat", "chiInertial.dat"
                ]:
                    data = horizon[dat]
                    chunks = (data.shape[0], 1)
                    g.create_dataset(dat, data=data, shuffle=True, compression="gzip", chunks=chunks)


def load(file):
    """Load SpEC-format Horizons.h5 file as `Horizons` object

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
    sxs.horizons.spec_horizons_h5.save : save to this file format

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
    import h5py
    from .. import Horizons, HorizonQuantities

    hqs = {}
    with h5py.File(file, "r") as f:
        sxs_format = f.attrs.get("sxs_format", "horizons.spec_horizons_h5")
        if sxs_format != "horizons.spec_horizons_h5":
            raise ValueError(
                f"\nAttribute 'sxs_format' found in '{file}' is '{sxs_format}'.\n"
                f"This function only accepts 'horizons.spec_horizons_h5' formats.\n"
                f"Use a higher-level `load` function to auto-detect the format."
            )
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
