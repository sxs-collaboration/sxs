"""Handlers for waveforms in the old NRAR data format"""

import sys
import warnings
import ast
import numpy as np
from .. import WaveformModes
from ...metadata import Metadata
from ...utilities.monotonicity import index_is_monotonic


FrameType = [UnknownFrameType, Inertial, Coprecessing, Coorbital, Corotating] = range(5)
FrameNames = ["UnknownFrameType", "Inertial", "Coprecessing", "Coorbital", "Corotating"]

DataType = [UnknownDataType, psi0, psi1, psi2, psi3, psi4, sigma, h, hdot, news, psin] = range(11)
DataNames = ["UnknownDataType", "Psi0", "Psi1", "Psi2", "Psi3", "Psi4", "sigma", "h", "hdot", "news", "psin"]
SpinWeights = [sys.maxsize, 2, 1, 0, -1, -2, 2, -2, -2, -2, sys.maxsize]
ConformalWeights = [sys.maxsize, 2, 1, 0, -1, -2, 1, 0, -1, -1, -3]
RScaling = [sys.maxsize, 5, 4, 3, 2, 1, 2, 1, 1, 1, 0]
MScaling = [sys.maxsize, 2, 2, 2, 2, 2, 0, 0, 1, 1, 2]
DataNamesLaTeX = [
    r"\mathrm{unknown data type}",
    r"\psi_0",
    r"\psi_1",
    r"\psi_2",
    r"\psi_3",
    r"\psi_4",
    r"\sigma",
    r"h",
    r"\dot{h}",
    r"\mathrm{n}",
    r"\psi_n",
]


def translate_frame_type_to_sxs_string(f):
    return ["unknown", "inertial", "coprecessing", "coorbital", "corotating"][f]

def translate_data_type_to_spin_weight(d):
    return [sys.maxsize, 2, 1, 0, -1, -2, 2, -2, -2, -2, sys.maxsize][d]

def translate_data_type_to_sxs_string(d):
    return ["unknown", "psi0", "psi1", "psi2", "psi3", "psi4", "sigma", "h", "hdot", "news", "psin"][d]

def translate_data_types_GWFrames_to_waveforms(d):
    if d < 8:
        return {0: UnknownDataType, 1: h, 2: hdot, 3: psi4, 4: psi3, 5: psi2, 6: psi1, 7: psi0}[d]
    else:
        return DataType[d - 8]


def translate_data_types_waveforms_to_GWFrames(d):
    if d in [UnknownDataType, h, hdot, psi4, psi3, psi2, psi1, psi0]:
        return {UnknownDataType: 0, h: 1, hdot: 2, psi4: 3, psi3: 4, psi2: 5, psi1: 6, psi0: 7}[d]
    else:
        return d + 8


def frame_type_string(w):
    return FrameNames[w._metadata.get("frame_type", UnknownFrameType)]


def data_type_string(w):
    return DataNames[w._metadata.get("data_type", UnknownDataType)]


def data_type_latex(w):
    return DataNamesLaTeX[w._metadata.get("data_type", UnknownDataType)]


def descriptor_string(w):
    """Create a simple string describing the content of the waveform

    This string will be suitable for file names.  For example, 'rMpsi4' or
    'rhOverM'.  It uses the waveform's knowledge of itself, so if this is
    incorrect, the result will be incorrect.

    """
    dataType = w._metadata.get("data_type", UnknownDataType)
    r_is_scaled_out = w._metadata.get("r_is_scaled_out", 1)
    m_is_scaled_out = w._metadata.get("m_is_scaled_out", 1)
    data_type_string = DataNames[dataType]
    if dataType == UnknownDataType:
        return data_type_string
    descriptor = ""
    if r_is_scaled_out:
        if RScaling[dataType] == 1:
            descriptor = "r"
        elif RScaling[dataType] > 1:
            descriptor = "r" + str(RScaling[dataType])
    if m_is_scaled_out:
        Mexponent = MScaling[dataType] - (RScaling[dataType] if r_is_scaled_out else 0)
        if Mexponent < -1:
            descriptor = descriptor + data_type_string + "OverM" + str(-Mexponent)
        elif Mexponent == -1:
            descriptor = descriptor + data_type_string + "OverM"
        elif Mexponent == 0:
            descriptor = descriptor + data_type_string
        elif Mexponent == 1:
            descriptor = descriptor + "M" + data_type_string
        elif Mexponent > 1:
            descriptor = descriptor + "M" + str(Mexponent) + data_type_string
    else:
        descriptor = descriptor + data_type_string
    return descriptor


def load(file, **kwargs):
    """Read data from an H5 file in NRAR format

    Note that either `h5_group` or `extrapolation_order` must be given.

    Parameters
    ----------
    file : {str, pathlib.Path}
        Path to H5 file containing the data, optionally including the path within
        the file itself to the directory containing the data.  For example, the
        standard SXS data with N=2 might be obtained with the file name
        `'rhOverM_Asymptotic_GeometricUnits.h5/Extrapolated_N2.dir'`.

    Keyword parameters
    ------------------
    h5_group : str, optional
        HDF5 group in which to find the data within the H5 file
    extrapolation_order : {None, str, int, Ellipsis}, optional
        Extrapolation order to load from the file.  This can be a string naming the
        HDF5 group to find the data in — like "OutermostExtraction.dir" or
        "Extrapolated_N2.dir" — or it can be an integer that is translated into
        such a string, where -1 is translated into "OutermostExtraction.dir".  The
        default value of `None` will return a dict-like object mapping names like
        the above to their corresponding waveforms, but will also issue a warning
        that it is doing so.  The final option of `Ellipsis` (which can also be
        input as `...`) will also return this dict-like object, but will not issue
        a warning.
    frame_type : int, optional
    data_type : int, optional
    r_is_scaled_out : bool, optional
        True if the radius is scaled out, meaning that the asymptotic form of this
        data could be finite and nonzero.
    m_is_scaled_out : bool, optional
        True if the mass is scaled out of the data, so that the data does not
        change when the total mass changes (presumably for vacuum systems only).
    transform_to_inertial : bool, optional
        This is a slightly unnatural argument for this format, because the data is
        already in the inertial frame.  It is included for compatibility with the
        RPDMB loader.  Here, the default value is `True`, meaning that nothing is
        done.  If set to `False`, the data is transformed to the corotating frame.

    """
    import pathlib
    import re
    import h5py
    import quaternionic
    import spherical
    from ...utilities import KeyPassingDict

    # This unfortunate concoction is needed to determine the (ell,m) values of the various mode data sets
    pattern_Ylm = re.compile(r"""Y_l(?P<L>[0-9]+)_m(?P<M>[-+0-9]+)\.dat""")

    # Initialize an empty object to be filled with goodies
    w_attributes = {}

    # Get an h5py handle to the desired part of the h5 file
    file_str = str(file)
    split = file_str.rsplit(".h5", 1)  # Raises ValueError if ".h5" is not in the string
    if len(split) != 2:
        split = file_str.rsplit(".hdf5", 1)
        file_ending = ".hdf5"
    else:
        file_ending = ".h5"
    if len(split) != 2:
        raise ValueError(f"Could not find a valid HDF5 filename ending in '{file_str}'")
    file_str, root_group = split
    file_str = file_str + file_ending
    if not root_group:
        root_group = kwargs.pop("h5_group", "")
    if not root_group:
        extrapolation_order = kwargs.pop("extrapolation_order", None)
        if extrapolation_order is None:
            warning = "\nCould not find root group as `h5_group` or as `extrapolation_order`; returning all groups"
            warnings.warn(warning)
        elif extrapolation_order is Ellipsis:
            pass
        elif isinstance(extrapolation_order, str):
            root_group = extrapolation_order
        elif extrapolation_order == -1:
            root_group = "OutermostExtraction.dir"
        else:
            root_group = f"Extrapolated_N{extrapolation_order}.dir"

    file_path = pathlib.Path(file_str).expanduser().resolve()
    file_name = file_path.name
    file_dir = file_path.parent

    with h5py.File(file_path, "r") as f_h5:
        if root_group:
            if root_group not in f_h5:
                raise ValueError(f"Input root group '{root_group}' was not found in '{file_path}'")
            f = f_h5[root_group]
        else:
            return KeyPassingDict(**{
                new_h5_group: load(file, **dict(h5_group=new_h5_group, **kwargs))
                for new_h5_group in f_h5 if "VersionHist.ver" not in new_h5_group
            }, **{
                new_h5_group: f_h5[new_h5_group][:]
                for new_h5_group in f_h5 if "VersionHist.ver" in new_h5_group
            })

        # If it exists, add the metadata file to `w` as an object.  So, for example, the initial spin on
        # object 1 can be accessed as `w.metadata.initial_spin1`.  See the documentation of
        # `sxs.metadata.Metadata` for more details.  And in IPython, tab completion works on the
        # `w.metadata` object.
        try:
            w_attributes["metadata"] = Metadata.from_file(
                file_dir / "metadata", ignore_invalid_lines=True, cache_json=False
            )
        except:
            pass  # Probably couldn't find the metadata.json/metadata.txt file

        try:  # Print a more explanatory error message if this fails
            # Extract the old History to search for COM-correction parameters
            if "History.txt" in f:
                try:
                    old_history = f["History.txt"][()].decode()
                except AttributeError:
                    old_history = f["History.txt"][()]
            else:
                old_history = ""
            w_attributes["history"] = old_history

            # Get the frame data, converting to quaternion objects
            if "Frame" in f:
                w_attributes["frame"] = quaternionic.array(f["Frame"])

            # Get the descriptive items
            if "FrameType" in f.attrs:
                w_attributes["frame_type"] = translate_frame_type_to_sxs_string(int(f.attrs["FrameType"]))
            elif "frame_type" in kwargs:
                w_attributes["frame_type"] = kwargs.pop("frame_type")
            else:
                warning = (
                    f"\n`frameType` was not found in '{file_str}' or the keyword arguments.\n"
                    + "Using default value `inertial`.  You may want to set it manually.\n\n"
                )
                warnings.warn(warning)
                w_attributes["frame_type"] = "inertial"

            if "DataType" in f.attrs:
                w_attributes["data_type"] = translate_data_types_GWFrames_to_waveforms(int(f.attrs["DataType"]))
            elif "data_type" in kwargs:
                w_attributes["data_type"] = int(kwargs.pop("data_type"))
            elif any(file_name.startswith(s) for s in {"rhOverM_", "rh_"}):
                w_attributes["data_type"] = h
            elif any(file_name.startswith(s) for s in {"rMPsi4_", "rPsi4_"}):
                w_attributes["data_type"] = psi4
            else:
                found = False
                for type_int, type_name in zip(reversed(DataType), reversed(DataNames)):
                    if type_name.lower() in file_name.lower():
                        found = True
                        w_attributes["data_type"] = type_int
                        break
                if not found:
                    warning = (
                        f"\n`dataType` was not found in '{file_str}' or the keyword arguments.\n"
                        + f"Using default value `{type_name}`.  You may want to set it manually."
                    )
                    warnings.warn(warning)
                    w_attributes["data_type"] = 0
            w_attributes["spin_weight"] = translate_data_type_to_spin_weight(w_attributes["data_type"])
            w_attributes["data_type"] = translate_data_type_to_sxs_string(w_attributes["data_type"])

            if "RIsScaledOut" in f.attrs:
                w_attributes["r_is_scaled_out"] = bool(f.attrs["RIsScaledOut"])
            elif "r_is_scaled_out" in kwargs:
                w_attributes["r_is_scaled_out"] = bool(kwargs.pop("r_is_scaled_out"))
            elif any(file_name.startswith(s) for s in {"rhOverM_", "rh_", "rMPsi4_", "rPsi4_"}):
                w_attributes["r_is_scaled_out"] = True
            else:
                warning = (
                    f"\n`r_is_scaled_out` was not found in '{file_str}' or the keyword arguments.\n"
                    + "Using default value `True`.  You may want to set it manually.\n\n"
                )
                warnings.warn(warning)
                w_attributes["r_is_scaled_out"] = True

            if "MIsScaledOut" in f.attrs:
                w_attributes["m_is_scaled_out"] = bool(f.attrs["MIsScaledOut"])
            elif "m_is_scaled_out" in kwargs:
                w_attributes["m_is_scaled_out"] = bool(kwargs.pop("m_is_scaled_out"))
            elif any(file_name.startswith(s) for s in {"rhOverM_", "rMPsi4_"}):
                w_attributes["m_is_scaled_out"] = True
            else:
                warning = (
                    f"\n`m_is_scaled_out` was not found in '{file_str}' or the keyword arguments.\n"
                    + "Using default value `True`.  You may want to set it manually.\n\n"
                )
                warnings.warn(warning)
                w_attributes["m_is_scaled_out"] = True

            # Get the names of all the data sets in the h5 file, and check for matches
            YLMdata = [data_set for data_set in list(f) for m in [pattern_Ylm.search(data_set)] if m]
            if len(YLMdata) == 0:
                raise ValueError(
                    f"Couldn't understand data set names in '{file_str}'.\n"
                    + "Maybe you need to add the directory within the h5 file.\n"
                    + f"E.g.: '{file_str}/Extrapolated_N2.dir'."
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
            if not np.array_equal(LM, spherical.Yrange(ell_min, ell_max)):
                raise ValueError(f"Input [ell,m] modes are not complete.  Found modes:\n{LM}\n")
            n_modes = len(LM)

            # Get the time data (assuming all are equal)
            T = f[YLMdata[0]][:, 0]
            monotonic = index_is_monotonic(T)
            w_attributes["t"] = T[monotonic]
            n_times = len(w_attributes["t"])

            # Loop through, setting data in each mode
            w_attributes["data"] = np.empty((n_times, n_modes), dtype=complex)
            for m, DataSet in enumerate(YLMdata):
                if f[DataSet].shape[0] != n_times:
                    raise ValueError(
                        f"The number of time steps in this dataset should be {n_times}; "
                        + "it is {} in '{}'.".format(f[DataSet].shape[0], DataSet)
                    )
                w_attributes["data"][:, m] = f[DataSet][:, 1:3].view(dtype=complex)[monotonic, 0]

            # Now that the data is set, we can set these
            w_attributes["ells"] = ell_min, ell_max

            # If possible, retrieve the CoM-correction parameters
            try:
                if hasattr(f, "attrs") and "space_translation" in f.attrs:
                    w_attributes["space_translation"] = np.array(list(f.attrs["space_translation"]))
                elif old_history:
                    pattern = r"'{}': array\((.*?)\)".format("space_translation")
                    matches = re.search(pattern, old_history)
                    if matches:
                        w_attributes["space_translation"] = np.array(ast.literal_eval(matches.group(1)))
            except:
                pass
            try:
                if hasattr(f, "attrs") and "boost_velocity" in f.attrs:
                    w_attributes["boost_velocity"] = np.array(list(f.attrs["boost_velocity"]))
                elif old_history:
                    pattern = r"'{}': array\((.*?)\)".format("boost_velocity")
                    matches = re.search(pattern, old_history)
                    if matches:
                        w_attributes["boost_velocity"] = np.array(ast.literal_eval(matches.group(1)))
            except:
                pass

            # If possible, retrieve the CoM-correction parameters
            try:
                if "VersionHist.ver" in f_h5:
                    w_attributes["version_hist"] = [
                        [git_hash.decode("ascii"), message.decode("ascii")]
                        for git_hash, message in f_h5["VersionHist.ver"][()].tolist()
                    ]
            except:
                pass

        except KeyError as e:
            raise ValueError("\nThis H5 file appears to have not stored all the required information.\n\n") from e

    transform_to_inertial = kwargs.pop("transform_to_inertial", True)

    if kwargs:
        import pprint
        warnings.warn("\nUnused kwargs passed to this function:\n{}".format(pprint.pformat(kwargs, width=1)))

    data = w_attributes.pop("data")
    time = w_attributes.pop("t")
    ell_min, ell_max = w_attributes.pop("ells")
    metadata = w_attributes.pop("metadata", {})
    w = WaveformModes(
        data,
        time=time,
        time_axis=0,
        modes_axis=1,
        ell_min=ell_min,
        ell_max=ell_max,
        **w_attributes
    )
    if metadata:
        w.metadata = metadata

    if not transform_to_inertial:
        w = w.to_corotating_frame()

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
        file_name, group = file_name.rsplit(".h5", 1)
        file_name += ".h5"

    # Add descriptive prefix to file_name
    base_name = descriptor_string(w) + "_" + os.path.basename(file_name)
    if not os.path.dirname(file_name):
        file_name = base_name
    else:
        file_name = os.path.join(os.path.dirname(file_name), base_name)

    # Open the file for output
    with h5py.File(file_name, file_write_mode) as f:
        # If we are writing to a group within the file, create it
        if group:
            g = f.create_group(group)
        else:
            g = f

        # Now write all the data to various groups in the file
        g.attrs["OutputFormatVersion"] = "scri.SpEC"
        g.attrs["FrameType"] = w._metadata["frame_type"]
        g.attrs["DataType"] = translate_data_types_waveforms_to_GWFrames(w._metadata["data_type"])
        g.attrs["RIsScaledOut"] = int(w._metadata.get("r_is_scaled_out", 1))
        g.attrs["MIsScaledOut"] = int(w._metadata.get("m_is_scaled_out", 1))
        if "version_hist" in w._metadata:
            try:
                version_hist = [[e.encode("ascii", "ignore") for e in hm] for hm in w._metadata["version_hist"]]
            except AttributeError:
                version_hist = w._metadata["version_hist"]
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
                    data=[[t, d.real, d.imag] for t, d in zip(w.time, w.data[:, i_m])],
                    compression="gzip",
                    shuffle=True,
                )
                Data_m.attrs["ell"] = ell
                Data_m.attrs["m"] = m
        else:
            g.create_dataset("Time", data=w.time.tolist(), compression="gzip", shuffle=True)
            frame = w._metadata.get("frame", [])
            if len(frame) > 0:
                g.create_dataset("Frame", data=[[r.w, r.x, r.y, r.z] for r in frame])
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
