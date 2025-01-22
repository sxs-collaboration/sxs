from ... import sxs_directory
from . import nrar


def load(file, **kwargs):
    """Load a waveform from a GR-Athena++ `tar` file.

    Parameters
    ==========
    file : str
        The path to the `tar` file containing the waveform data.

    Other Parameters
    ================
    subfile : str
        The name of the subfile within the `tar` file to read.
        Default is "rh_CCE_GeometricUnits.h5".  Other potentially
        useful values replace "rh" with "rPsi4" and/or "CCE" with
        "Asymptotic" or "FiniteRadii".
    radius : str
        The extraction radius to use.  Default is "50.00".
    
    Notes
    =====
    Waveforms from GR-Athena++ are distributed as `tar` files
    containing a set of HDF5 files, containing CCE data,
    "extrapolated" data (where "extrapolation" refers to the
    single-radius PN-based correction method that Nakano introduced),
    or finite-radius data, for either Psi4 or h.  Each of those files
    is NRAR-formatted exactly like the old SXS data.  We could extract
    the `tar` file, and then use the `nrar.load` function to read the
    data.  Alternatively — as is done here — we could just use the
    `tar` file as a file-like object, and pass that to the `nrar.load`
    function.  This is a bit more efficient, and is preferred if we
    intend to directly use the `tar` files.

    The `tar` file is generally named with the resolution of the
    simulation, which — for the initial catalog at least — may be any
    of 128, 192, 256, 320, or 384.

    Within the `tar` file, we have files like the following — though
    "384" may be replaced by any resolution:

    * "384/rPsi4_Asymptotic_GeometricUnits.h5"
    * "384/rPsi4_CCE_GeometricUnits.h5"
    * "384/rPsi4_FiniteRadii_GeometricUnits.h5"
    * "384/rh_Asymptotic_GeometricUnits.h5"
    * "384/rh_CCE_GeometricUnits.h5"
    * "384/rh_FiniteRadii_GeometricUnits.h5"

    And finally, within each of those h5 files, we have waveforms
    corresponding to a series of extraction radii, which may include

    * "50.00"
    * "60.00"
    * "70.00"
    * "80.00"
    * "90.00"
    * "100.00"
    * "120.00"
    * "140.00"

    (Note the two 0s after the decimal point.)  The `FiniteRadii` and
    `Asymptotic` data generally contain all of the above radii, while
    the `CCE` data generally contains only the "50.00" and "100.00"
    radius.  The catalog paper uses "50.00" as the radius when
    discussing properties of the waveforms, possibly because CCE only
    "exhibits convergence behavior" for data from that radius.  They
    also use "100.00" for "extrapolated" waveforms some times.  The
    default used here is "50.00".

    """
    from pathlib import Path
    import tarfile

    resolution = Path(file).stem
    subfile = kwargs.pop("subfile", "rh_CCE_GeometricUnits.h5")
    radius = kwargs.pop("radius", "50.00")

    if subfile.startswith("rPsi4"):
        data_type = nrar.psi4
    elif subfile.startswith("rh"):
        data_type = nrar.h

    with tarfile.open(file, "r") as tf:
        tf_names = [tfi.name for tfi in tf]
        index = tf_names.index(f"{resolution}/{subfile}")
        h5file = tf.extractfile(list(tf)[index])
        w = nrar.load(
            h5file,
            h5_group=radius,
            frame_type=nrar.Inertial,
            data_type=data_type,
            m_is_scaled_out=True,
            r_is_scaled_out=True,
        )

    w._metadata["frame_type"] = "inertial"
    if subfile.startswith("rPsi4"):
        w._metadata["data_type"] = "psi4"
    elif subfile.startswith("rh"):
        w._metadata["data_type"] = "h"

    return w
