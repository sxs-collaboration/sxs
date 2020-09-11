"""File I/O for XMB-format horizons.h5 files"""

from ..utilities import default_shuffle_widths


def save(horizons, file, truncate=lambda x: x, shuffle_widths=default_shuffle_widths):
    """Save `Horizons` object as XMB-format horizons.h5 file

    This function saves the input horizons object in a more compact form than the
    older SpEC-format Horizons.h5 file.  Unlike the older format, only one copy of
    the time data is stored and magnitude datasets are not stored.

    Parameters
    ----------
    horizons : sxs.Horizons
        Horizons object to be written to file.
    file : file-like object, string, or pathlib.Path
        Path to the file on disk or a file-like object (such as an open file
        handle) to be written by h5py.File.
    truncate : callable, optional
        Function that truncates the data to a desired accuracy.  This should set
        bits beyond a desired precision to 0 so that they will compress
        effectively.  By default, no truncation is performed.
    shuffle_widths : array_like, optional
        This must be a series of integers totalling 64 (or the bit width of the
        data in the input `horizons`).  This argument is passed to the multishuffle
        routine.

    See Also
    --------
    sxs.horizons.xor_multishuffle_bzip2.load : load the output file format
    sxs.utilities.xor_multishuffle_bzip2 : compresses data

    """
    import bz2
    import numpy as np
    import h5py
    from ..utilities import xor, multishuffle

    shuffle = multishuffle(tuple(shuffle_widths))

    with h5py.File(file, "w") as f:
        f.attrs["sxs_format"] = "horizons.xor_multishuffle_bzip2"
        f.attrs["shuffle_widths"] = np.array(shuffle_widths, dtype=np.uint8)
        for horizon_name in "ABC":
            horizon = horizons[horizon_name]
            if horizon is not None:
                compressor = bz2.BZ2Compressor()
                # Time
                data = compressor.compress(shuffle(xor(horizon.time)).tobytes())
                # 1-d data sets
                for dat in ["areal_mass", "christodoulou_mass"]:
                    data = data + compressor.compress(
                        shuffle(xor(truncate(horizon[dat].ndarray).view(np.uint64))).tobytes()
                    )
                # 2-d data sets
                for dat in ["coord_center_inertial", "dimensionful_inertial_spin", "chi_inertial"]:
                    data = data + compressor.compress(
                        shuffle(xor(truncate(horizon[dat].ndarray).view(np.uint64).flatten("F"))).tobytes()
                    )
                data = data + compressor.flush()
                d = f.create_dataset(horizon_name, data=np.void(data))
                d.attrs["n_times"] = horizon.time.size


def load(file, ignore_format=False):
    """Load SpEC-format Horizons.h5 file as `Horizons` object

    Parameters
    ----------
    file : file-like object, string, or pathlib.Path
        Path to the file on disk or a file-like object (such as an open file
        handle) to be opened by h5py.File.
    ignore_format : bool, optional
        If True, the format attribute in the input file will be ignored, and this
        function will attempt to load the data regardless.  Note that this will
        return incorrect data or raise an exception if the format is incompatible.
        The default is False, meaning that a ValueError will be raised if the
        format is incorrect.

    Returns
    -------
    horizons : sxs.Horizons
        This is a container for the horizon objects.  See Notes below.

    See also
    --------
    sxs.Horizons : Container object for all of the horizons
    sxs.HorizonQuantities : Container objects for each of the horizons
    sxs.horizons.xor_multishuffle_bzip2.save : save to this file format

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
    import bz2
    import numpy as np
    import h5py
    from .. import Horizons, HorizonQuantities
    from ..utilities import xor, multishuffle

    hqs = {}
    with h5py.File(file, "r") as f:
        if not ignore_format:
            if "sxs_format" not in f.attrs:
                raise ValueError(f"Attribute 'sxs_format' not found in '{file}'")
            sxs_format = f.attrs["sxs_format"]
            if sxs_format != "horizons.xor_multishuffle_bzip2":
                raise ValueError(
                    f"\nAttribute 'sxs_format' found in '{file}' is '{sxs_format}'.\n"
                    f"This function only accepts 'horizons.spec_horizons_h5' formats.\n"
                    f"Use a higher-level `load` function to auto-detect the format."
                )

        shuffle_widths = tuple(f.attrs["shuffle_widths"])
        shuffle = multishuffle(shuffle_widths, forward=False)

        for horizon_name in "ABC":
            if horizon_name in f:
                d = f[horizon_name]
                uncompressed_data = bz2.decompress(d[...])
                n_times = d.attrs["n_times"]
                sizeof_float = 8  # in bytes
                bytes_per_series = n_times * sizeof_float
                hq_kwargs = {}
                i = 0

                # scalar data sets
                for dat in ["time", "areal_mass", "christodoulou_mass"]:
                    hq_kwargs[dat] = xor(shuffle(
                        np.frombuffer(uncompressed_data[i: i + bytes_per_series], dtype=np.uint64)
                    ), reverse=True).view(np.float64)
                    i += bytes_per_series

                # 3-vector data sets
                for dat in ["coord_center_inertial", "dimensionful_inertial_spin", "chi_inertial"]:
                    hq_kwargs[dat] = np.asarray(xor(shuffle(
                        np.frombuffer(uncompressed_data[i: i + 3 * bytes_per_series], dtype=np.uint64)
                    ), reverse=True, preserve_dtype=True).reshape((3, -1)).T.view(np.float64), order="C")
                    i += 3 * bytes_per_series

                horizon = HorizonQuantities(**hq_kwargs)

                hqs[horizon_name] = horizon

    return Horizons(**hqs)
