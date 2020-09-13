import warnings
import os
import re
import ast
import numpy as np
from . import WaveformModes
from ..metadata import Metadata
from ..utilities.monotonicity import index_is_monotonic


def load(file_name, **kwargs):
    """Read data from an H5 file in NRAR format

    Parameters
    ----------
    file_name : str
        Path to H5 file containing the data, optionally including the path within
        the file itself to the directory containing the data.  For example, the
        standard SXS data with N=2 might be obtained with the file name
        `'rhOverM_Asymptotic_GeometricUnits.h5/Extrapolated_N2.dir'`.

    Keyword parameters
    ------------------
    frameType : int, optional
    dataType : int, optional

    r_is_scaled_out : bool, optional
        True if the radius is scaled out, meaning that the asymptotic form of this
        data could be finite and nonzero.
    m_is_scaled_out : bool, optional
        True if the mass is scaled out of the data, so that the data does not
        change when the total mass changes (presumably for vacuum systems only).

    """
    import os.path
    import re
    import h5py
    import quaternionic
    import spherical

    # This unfortunate concoction is needed to determine the (ell,m) values of the various mode data sets
    pattern_Ylm = re.compile(r"""Y_l(?P<L>[0-9]+)_m(?P<M>[-+0-9]+)\.dat""")

    # Initialize an empty object to be filled with goodies
    w = WaveformModes(constructor_statement=f'scri.SpEC.read_from_h5("{file_name}", **{kwargs})')

    # Get an h5py handle to the desired part of the h5 file
    try:
        file_name, root_group = file_name.rsplit(".h5", 1)
        file_name += ".h5"
    except ValueError:
        root_group = ""  # FileName is just a file, not a group in a file

    with h5py.File(file_name, "r") as f_h5:
        if root_group:
            f = f_h5[root_group]
        else:
            f = f_h5

        # If it exists, add the metadata file to `w` as an object.  So, for example, the initial spin on
        # object 1 can be accessed as `w.metadata.initial_spin1`.  See the documentation of
        # `sxs.metadata.Metadata` for more details.  And in IPython, tab completion works on the
        # `w.metadata` object.
        try:
            w.metadata = Metadata.from_file(
                os.path.join(os.path.dirname(file_name), "metadata"), ignore_invalid_lines=True, cache_json=False
            )
        except:
            pass  # Probably couldn't find the metadata.json/metadata.txt file

        try:  # Make sure the h5py.File gets closed, even in the event of an exception

            # Get the frame data, converting to quaternion objects
            try:
                w.frame = quaternionic.array(f["Frame"])
            except KeyError:
                pass  # There was no frame data

            # Get the descriptive items
            try:
                w.frameType = int(f.attrs["FrameType"])
            except KeyError:
                if "frameType" in kwargs:
                    w.frameType = int(kwargs.pop("frameType"))
                else:
                    warning = (
                        f"\n`frameType` was not found in '{file_name}' or the keyword arguments.\n"
                        + "Using default value `{}`.  You may want to set it manually.\n\n".format(FrameNames[1])
                    )
                    warnings.warn(warning)
                    w.frameType = 1
            try:
                w.dataType = translate_data_types_GWFrames_to_waveforms(int(f.attrs["DataType"]))
            except KeyError:
                if "dataType" in kwargs:
                    w.dataType = int(kwargs.pop("dataType"))
                else:
                    found = False
                    for type_int, type_name in zip(reversed(DataType), reversed(DataNames)):
                        if type_name.lower() in file_name.lower():
                            found = True
                            w.dataType = type_int
                            warning = (
                                f"\n`dataType` was not found in '{file_name}' or the keyword arguments.\n"
                                + f"Using default value `{type_name}`.  You may want to set it manually."
                            )
                            warnings.warn(warning)
                            break
                    if not found:
                        warning = (
                            f"\n`dataType` was not found in '{file_name}' or the keyword arguments.\n"
                            + "You may want to set it manually."
                        )
                        warnings.warn(warning)
            try:
                w.r_is_scaled_out = bool(f.attrs["RIsScaledOut"])
            except KeyError:
                if "r_is_scaled_out" in kwargs:
                    w.r_is_scaled_out = bool(kwargs.pop("r_is_scaled_out"))
                else:
                    warning = (
                        f"\n`r_is_scaled_out` was not found in '{file_name}' or the keyword arguments.\n"
                        + "Using default value `True`.  You may want to set it manually.\n\n"
                    )
                    warnings.warn(warning)
                    w.r_is_scaled_out = True
            try:
                w.m_is_scaled_out = bool(f.attrs["MIsScaledOut"])
            except KeyError:
                if "m_is_scaled_out" in kwargs:
                    w.m_is_scaled_out = bool(kwargs.pop("m_is_scaled_out"))
                else:
                    warning = (
                        f"\n`m_is_scaled_out` was not found in '{file_name}' or the keyword arguments.\n"
                        + "Using default value `True`.  You may want to set it manually.\n\n"
                    )
                    warnings.warn(warning)
                    w.m_is_scaled_out = True

            # Get the names of all the data sets in the h5 file, and check for matches
            YLMdata = [data_set for data_set in list(f) for m in [pattern_Ylm.search(data_set)] if m]
            if len(YLMdata) == 0:
                raise ValueError(
                    f"Couldn't understand data set names in '{file_name}'.\n"
                    + "Maybe you need to add the directory within the h5 file.\n"
                    + f"E.g.: '{file_name}/Extrapolated_N2.dir'."
                )

            # Sort the data set names by increasing ell, then increasing m
            YLMdata = sorted(
                YLMdata,
                key=lambda data_set: [
                    int(pattern_Ylm.search(data_set).group("L")),
                    int(pattern_Ylm.search(data_set).group("M")),
                ],
            )
            LM = np.array(
                sorted(
                    [
                        [int(m.group("L")), int(m.group("M"))]
                        for data_set in YLMdata
                        for m in [pattern_Ylm.search(data_set)]
                        if m
                    ]
                )
            )
            ell_min, ell_max = min(LM[:, 0]), max(LM[:, 0])
            if not np.array_equal(LM, spherical.LM_range(ell_min, ell_max)):
                raise ValueError(f"Input [ell,m] modes are not complete.  Found modes:\n{LM}\n")
            n_modes = len(LM)

            # Get the time data (assuming all are equal)
            T = f[YLMdata[0]][:, 0]
            monotonic = index_is_monotonic(T)
            w.t = T[monotonic]
            n_times = len(w.t)

            # Loop through, setting data in each mode
            w.data = np.empty((n_times, n_modes), dtype=complex)
            for m, DataSet in enumerate(YLMdata):
                if f[DataSet].shape[0] != n_times:
                    raise ValueError(
                        f"The number of time steps in this dataset should be {n_times}; "
                        + "it is {} in '{}'.".format(f[DataSet].shape[0], DataSet)
                    )
                w.data[:, m] = f[DataSet][:, 1:3].view(dtype=np.complex)[monotonic, 0]

            # Now that the data is set, we can set these
            w.ells = ell_min, ell_max

            # If possible, retrieve the CoM-correction parameters
            try:
                if hasattr(f, "attrs") and "space_translation" in f.attrs:
                    w.space_translation = np.array(list(f.attrs["space_translation"]))
                elif old_history:
                    pattern = r"'{}': array\((.*?)\)".format("space_translation")
                    matches = re.search(pattern, old_history)
                    if matches:
                        w.space_translation = np.array(ast.literal_eval(matches.group(1)))
            except:
                pass
            try:
                if hasattr(f, "attrs") and "boost_velocity" in f.attrs:
                    w.boost_velocity = np.array(list(f.attrs["boost_velocity"]))
                elif old_history:
                    pattern = r"'{}': array\((.*?)\)".format("boost_velocity")
                    matches = re.search(pattern, old_history)
                    if matches:
                        w.boost_velocity = np.array(ast.literal_eval(matches.group(1)))
            except:
                pass

            # If possible, retrieve the CoM-correction parameters
            try:
                if "VersionHist.ver" in f_h5:
                    w.version_hist = [
                        [git_hash.decode("ascii"), message.decode("ascii")]
                        for git_hash, message in f_h5["VersionHist.ver"][()].tolist()
                    ]
            except:
                pass

            # Check up on the validity of the waveform
            if not w.ensure_validity(alter=True, assertions=False):
                raise ValueError("The data resulting from this input is invalid")

        except KeyError:
            print("\nThis H5 file appears to have not stored all the required information.\n\n")
            raise  # Re-raise the exception after adding our information

    if kwargs:
        import pprint
        warnings.warn("\nUnused kwargs passed to this function:\n{}".format(pprint.pformat(kwargs, width=1)))

    return w


def save(w, file_name, file_write_mode="w", attributes={}, use_NRAR_format=True):
    """Output to an HDF5 file

    Note that the file_name is prepended with some descriptive information
    involving the data type and the frame type, such as 'rhOverM_Corotating_' or
    'rMpsi4_Aligned_'.

    """

    import os.path
    import h5py
    import warnings

    group = None
    if ".h5" in file_name and not file_name.endswith(".h5"):
        file_name, group = file_name.split(".h5")
        file_name += ".h5"
    # Add descriptive prefix to file_name
    base_name = w.descriptor_string + "_" + os.path.basename(file_name)
    if not os.path.dirname(file_name):
        file_name = base_name
    else:
        file_name = os.path.join(os.path.dirname(file_name), base_name)
    # Open the file for output
    try:
        f = h5py.File(file_name, file_write_mode)
    except OSError:  # If that did not work...
        print(f"write_to_h5 was unable to open the file '{file_name}'.\n\n")
        raise  # re-raise the exception after the informative message above
    try:
        # If we are writing to a group within the file, create it
        if group:
            g = f.create_group(group)
        else:
            g = f
        # Now write all the data to various groups in the file
        g.attrs["OutputFormatVersion"] = "scri.SpEC"
        g.create_dataset("History.txt", data="\n".join(w.history) + f"\n\nwrite_to_h5({w}, {file_name})\n")
        g.attrs["FrameType"] = w.frameType
        g.attrs["DataType"] = translate_data_types_waveforms_to_GWFrames(w.dataType)
        g.attrs["RIsScaledOut"] = int(w.r_is_scaled_out)
        g.attrs["MIsScaledOut"] = int(w.m_is_scaled_out)
        if len(w.version_hist) > 0:
            try:
                version_hist = [[e.encode("ascii", "ignore") for e in hm] for hm in w.version_hist]
            except AttributeError:
                version_hist = w.version_hist
            g.create_dataset("VersionHist.ver", data=version_hist)
        for attr in attributes:
            try:
                g.attrs[attr] = attributes[attr]
            except:
                warning = f"scri.SpEC.write_to_h5 unable to output attribute {attr}={attributes[attr]}"
                warnings.warn(warning)
        if use_NRAR_format:
            for i_m in range(w.n_modes):
                ell, m = w.LM[i_m]
                Data_m = g.create_dataset(
                    f"Y_l{ell}_m{m}.dat",
                    data=[[t, d.real, d.imag] for t, d in zip(w.t, w.data[:, i_m])],
                    compression="gzip",
                    shuffle=True,
                )
                Data_m.attrs["ell"] = ell
                Data_m.attrs["m"] = m
        else:
            g.create_dataset("Time", data=w.t.tolist(), compression="gzip", shuffle=True)
            if len(w.frame) > 0:
                g.create_dataset("Frame", data=[[r.w, r.x, r.y, r.z] for r in w.frame])
            else:
                g.create_dataset("Frame", shape=())
            Data = g.create_group("Data")
            for i_m in range(w.n_modes):
                ell, m = w.LM[i_m]
                Data_m = Data.create_dataset(
                    "l{}_m{:+}".format(int(ell), int(m)), data=w.data[:, i_m], compression="gzip", shuffle=True
                )
                Data_m.attrs["ell"] = ell
                Data_m.attrs["m"] = m
    finally:  # Use `finally` to make sure this happens:
        f.close()
