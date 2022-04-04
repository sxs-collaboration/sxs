"""Functions to load and save waveforms in RPXMB format"""

import sys
import warnings
import tempfile
import contextlib
import pathlib
import bz2
import json
import numpy as np
import scipy
import h5py
import quaternionic
import spherical

from . import WaveformModes
from .. import __version__
from ..utilities import default_shuffle_widths, md5checksum, xor, multishuffle


sxs_formats = ["rotating_paired_xor_multishuffle_bzip2", "rpxmb", "rpxm", "RPXMB", "RPXM"]


def save(w, file_name=None, file_write_mode="w", L2norm_fractional_tolerance=1e-10, log_frame=None, shuffle_widths=default_shuffle_widths):
    """Save a waveform in RPXMB format

    This function converts the data to "rotating paired XOR multishuffle bzip2"
    format.  In particular, it uses the corotating frame, and zeroes out bits at
    high precision to allow for optimal compression while maintaining the requested
    tolerance.

    Parameters
    ----------
    w : WaveformModes
        A waveform in either the inertial or corotating frame
    file_name : str
        Relative or absolute path to the output HDF5 file.  If this string contains
        `'.h5'` but does not *end* with that, the remainder of the string is taken
        to be the group within the HDF5 file in which the data should be stored.
        Also note that a JSON file is created in the same location, with `.h5`
        replaced by `.json` (and the corresponding data is stored under the `group`
        key if relevant).  For testing purposes, this argument may be `None`, in
        which case a temporary directory is used, just to test how large the output
        will be; it is deleted immediately upon returning.
    file_write_mode : str, optional
        One of the valid [file modes for
        h5py](https://docs.h5py.org/en/stable/high/file.html#opening-creating-files).
        Default value is `"w"`, which overwrites any existing file.  If writing
        into a group, you may prefer `"a"`, which will just ensure the file exists
        without erasing it first.
    L2norm_fractional_tolerance : float, optional
        Tolerance passed to `WaveformModes.truncate`; see that function's docstring
        for details.  Default value is 1e-10.
    log_frame : array of quaternions, optional
        If this argument is given the waveform must be in the corotating frame, and
        the given data will be used as the logarithmic frame data.  If this
        argument is `None` (the default), this will be calculated when the waveform
        is transformed to the corotating frame, or simply taken directly from the
        waveform if it is already corotating.
    shuffle_widths : iterable of ints, optional
        See `sxs.utilities.multishuffle` for details.  The default value is
        `default_shuffle_widths`.  Note that if `L2norm_fractional_tolerance` is
        0.0, this will be ignored and the standard HDF5 shuffle option will be used
        instead.

    Returns
    -------
    w_out : WaveformModes
        The output data, after conversion to the corotating frame, pairing of
        opposite `m` modes, and XOR-ing (but not shuffling).
    log_frame : array of quaternions
        The actual `log_frame` data stored in the file, and used to transform to
        the corotating frame if that was done inside this function.

    Note that the returned data are *as stored in the file*.  Specifically, they
    are presented as various types of `float` data, but have been XOR-ed, which
    makes them invalid as floats; you will see many NaNs and other nonsensical
    values unless you reverse the process.

    """
    # Make sure that we can understand the file_name and create the directory
    group = None
    if file_name is None:
        # We'll just be creating a temp directory below, to check
        warning = (
            "\nInput `file_name` is None.  Running in temporary directory.\n"
            "Note that this option is mostly for debugging purposes."
        )
        warnings.warn(warning)
    else:
        file_name = str(file_name)
        if ".h5" in file_name and not file_name.endswith(".h5"):
            file_name, group = file_name.split(".h5")
        h5_path = pathlib.Path(file_name).expanduser().resolve().with_suffix(".h5")
        h5_path.parent.mkdir(parents=True, exist_ok=True)
    if group == "/":
        group = None

    shuffle = multishuffle(tuple(shuffle_widths))

    if L2norm_fractional_tolerance == 0.0:
        log_frame = np.log(w.frame).ndarray[:, 1:]
    else:
        # We need this storage anyway, so let's just make a copy and work in-place
        w = w.copy()
        if log_frame is not None:
            log_frame = log_frame.copy()

        # Ensure waveform is in corotating frame
        if w.frame_type == "inertial":
            try:
                initial_time = w.t[0]
                relaxation_time = w.metadata.relaxation_time
                max_norm_time = w.max_norm_time()
                z_alignment_region = ((relaxation_time - initial_time) / (max_norm_time - initial_time), 0.95)
            except Exception:
                z_alignment_region = (0.1, 0.95)
            w, log_frame = w.to_corotating_frame(
                tolerance=1e-10, z_alignment_region=z_alignment_region, truncate_log_frame=True
            )
            log_frame = log_frame.ndarray[:, 1:]
        if w.frame_type != "corotating":
            raise ValueError(
                f"Frame type of input waveform must be 'corotating' or 'inertial'; it is {w.frame_type}"
            )

        # Convert mode structure to conjugate pairs
        w.convert_to_conjugate_pairs()

        # Set bits below the desired significance level to 0
        w.truncate(tol=L2norm_fractional_tolerance)

        # Compute log(frame)
        if log_frame is None:
            log_frame = np.log(w.frame).ndarray[:, 1:]

        # Change -0.0 to 0.0 (~.5% compression for non-precessing systems)
        w.t += 0.0
        np.add(w.data, 0.0, out=w.data)  # Use `.data` so WaveformModes doesn't override scalar addition
        log_frame += 0.0

        # XOR successive instants in time
        t = xor(w.t)
        data = xor(w.data.view(float), axis=w.time_axis)
        log_frame = xor(log_frame, axis=0)

    # Make sure we have a place to keep all this
    with contextlib.ExitStack() as context:
        if file_name is None:
            temp_dir = context.enter_context(tempfile.TemporaryDirectory())
            h5_path = pathlib.Path(f"{temp_dir}") / "test.h5"
        else:
            print(f'Saving H5 to "{h5_path}"')

        # Write the H5 file
        with h5py.File(h5_path, file_write_mode) as f:
            # If we are writing to a group within the file, create it
            if group is not None:
                g = f.create_group(group)
            else:
                g = f
            if L2norm_fractional_tolerance != 0.0:
                g.attrs["sxs_format"] = sxs_formats[0]
                g.attrs["n_times"] = w.n_times
                g.attrs["ell_min"] = w.ell_min
                g.attrs["ell_max"] = w.ell_max
                g.attrs["shuffle_widths"] = np.array(shuffle_widths, dtype=np.uint8)
                # warnings.warn(f'sxs_format is being set to "{sxs_formats[0]}"')
                data = np.void(
                    bz2.compress(
                        shuffle(t).tobytes()
                        + shuffle(data.flatten("F")).tobytes()
                        + shuffle(log_frame.flatten("F")).tobytes()
                    )
                )
                g.create_dataset("data", data=data)
            else:
                compression_options = {
                    "compression": "gzip",
                    "compression_opts": 9,
                    "shuffle": True,
                }
                g.attrs["sxs_format"] = f"{sxs_formats[0]}"
                g.create_dataset("time", data=w.t.view(np.uint64), chunks=(w.n_times,), **compression_options)
                g.create_dataset("modes", data=w.data.view(np.uint64), chunks=(w.n_times, 1), **compression_options)
                g["modes"].attrs["ell_min"] = w.ell_min
                g["modes"].attrs["ell_max"] = w.ell_max
                g["modes"].attrs["spin_weight"] = w.spin_weight
                if log_frame.size > 1:
                    g.create_dataset(
                        "log_frame", data=log_frame.view(np.uint64), chunks=(w.n_times, 1), **compression_options
                    )

        # Get some numbers for the JSON file
        h5_size = h5_path.stat().st_size
        if file_name is None:
            print(f"Output H5 file size: {h5_size:_} B")
        md5sum = md5checksum(h5_path)

        if file_name is not None:
            # Set up the corresponding JSON information
            json_data = {
                "sxs_format": sxs_formats[0],
                "data_info": {
                    "data_type": w.data_type,
                    "m_is_scaled_out": w._metadata.get("m_is_scaled_out", True),
                    "r_is_scaled_out": w._metadata.get("r_is_scaled_out", True),
                    "spin_weight": int(w.spin_weight),
                    "ell_min": int(w.ell_min),
                    "ell_max": int(w.ell_max),
                },
                "version_info": {
                    "python": sys.version,
                    "numpy": np.__version__,
                    "scipy": scipy.__version__,
                    "h5py": h5py.__version__,
                    "quaternionic": quaternionic.__version__,
                    "spherical": spherical.__version__,
                    "sxs": __version__,
                    # see below for "spec_version_hist"
                },
                # see below for "validation"
                # see below for "modifications"
            }
            if hasattr(w, "version_hist"):
                json_data["version_info"]["spec_version_history"] = w.version_hist
            if group is not None:
                json_data["validation"] = {
                    "n_times": w.n_times,
                }
            else:
                json_data["validation"] = {
                    "h5_file_size": h5_size,
                    "n_times": w.n_times,
                    "md5sum": md5sum
                }
            if "modifications" in w._metadata:
                json_data["modifications"] = w._metadata["modifications"]

            # Write the corresponding JSON file
            json_path = h5_path.with_suffix(".json")
            print(f'Saving JSON to "{json_path}"')
            if group is not None:
                if json_path.exists() and file_write_mode!="w":
                    with json_path.open("r") as f:
                        original_json = json.load(f)
                else:
                    original_json = {}
                original_json[group] = json_data
                json_data = original_json
            with json_path.open("w") as f:
                json.dump(json_data, f, indent=2, separators=(",", ": "), ensure_ascii=True)

    return w, log_frame


def load(file_name, ignore_validation=None, check_md5=True, transform_to_inertial=True, **kwargs):
    """Load a waveform in RPXMB format

    Parameters
    ----------
    file_name : str
        Relative or absolute path to the input HDF5 file.  If this string contains
        but does not *end* with `'.h5'`, the remainder of the string is taken to be
        the group within the HDF5 file in which the data is stored.  Also note that
        a JSON file is expected in the same location, with `.h5` replaced by
        `.json` (and the corresponding data must be stored under the `group` key if
        relevant).
    ignore_validation : bool or None, optional
        Validation checks the corresponding JSON file for (1) existence, (2) number
        of time steps, (3) H5 file size, and (4) H5 file MD5 checksum (if
        `check_md5` is `True`).  If this key is `False`, all of this will be
        ignored; if `True`, a `ValueError` will be raised if any of these checks
        fails; if `None`, warnings will be issued, but the function will continue
        as usual.
    check_md5 : bool, optional
        Default is `True`.  See `ignore_validation` for explanation.
    transform_to_inertial : bool, optional
        If `True`, the output waveform is transformed to the inertial frame;
        otherwise it is left in the frame used in the file.

    Keyword parameters
    ------------------
    data_type : str, optional
        Describes the type of data — such as "h" or "psi4".  Default is "unknown".
    m_is_scaled_out : bool, optional
        Default is True
    r_is_scaled_out : bool, optional
        Default is True
    spin_weight : int, optional
        Default is `None`.

    Note that the keyword parameters will be overridden by corresponding entries in
    the JSON file, if they exist.  If the JSON file does not exist, any keyword
    parameters not listed above will be passed through as the `json_data` field of
    the returned waveform.

    """
    def invalid(message):
        if ignore_validation:
            pass
        elif ignore_validation is None:
            warnings.warn(message)
        else:
            raise ValueError(message)

    file_name_str = str(file_name)
    group = None
    if ".h5" in file_name_str and not file_name_str.endswith(".h5"):
        file_name_str, group = file_name_str.split(".h5")
    if group == "/":
        group = None

    h5_path = pathlib.Path(file_name_str).expanduser().resolve().with_suffix(".h5")
    json_path = h5_path.with_suffix(".json")

    # This will be used for validation
    h5_size = h5_path.stat().st_size

    data_type = kwargs.pop("data_type", "unknown")
    m_is_scaled_out = kwargs.pop("m_is_scaled_out", True)
    r_is_scaled_out = kwargs.pop("r_is_scaled_out", True)
    spin_weight = kwargs.pop("spin_weight", None)

    if not json_path.exists():
        invalid(f'\nJSON file "{json_path}" cannot be found, but is expected for this data format.')
        json_data = kwargs.copy()
    else:
        with open(json_path) as f:
            json_data = json.load(f)
        if group is not None:
            json_data = json_data[group]

        data_type = json_data.get("data_info", {}).get("data_type", data_type)
        m_is_scaled_out = json_data.get("data_info", {}).get("m_is_scaled_out", m_is_scaled_out)
        r_is_scaled_out = json_data.get("data_info", {}).get("r_is_scaled_out", r_is_scaled_out)
        spin_weight = json_data.get("data_info", {}).get("spin_weight", spin_weight)

        # Make sure this is our format
        sxs_format = json_data.get("sxs_format", "")
        if sxs_format not in sxs_formats:
            invalid(
                f"\nThe `sxs_format` found in JSON file is '{sxs_format}';\n"
                f"it should be one of\n"
                f"    {sxs_formats}."
            )

        if group is None:
            # Make sure the expected H5 file size matches the observed value
            json_h5_file_size = json_data.get("validation", {}).get("h5_file_size", 0)
            if json_h5_file_size != h5_size:
                invalid(
                    f"\nMismatch between `validation/h5_file_size` key in JSON file ({json_h5_file_size}) "
                    f'and observed file size ({h5_size}) of "{h5_path}".'
                )

            # Make sure the expected H5 file hash matches the observed value
            if check_md5:
                md5sum = md5checksum(h5_path)
                json_md5sum = json_data.get("validation", {}).get("md5sum", "")
                if json_md5sum != md5sum:
                    invalid(f"\nMismatch between `validation/md5sum` key in JSON file and observed MD5 checksum.")

    with h5py.File(h5_path, "r") as f:
        if group is not None:
            g = f[group]
        else:
            g = f
        # Make sure this is our format
        sxs_format = g.attrs["sxs_format"]
        if sxs_format not in sxs_formats:
            raise ValueError(
                f'The `sxs_format` found in H5 file is "{sxs_format}"; it should be one of\n'
                f"    {sxs_formats}."
            )

        # Ensure that the 'validation' keys from the JSON file are the same as in this file
        n_times = g.attrs["n_times"]
        json_n_times = json_data.get("validation", {}).get("n_times", 0)
        if json_n_times != n_times:
            invalid(
                f"\nNumber of time steps in H5 file ({n_times}) "
                f"does not match expected value from JSON ({json_n_times})."
            )

        # Read the raw data
        sizeof_float = 8
        sizeof_complex = 2 * sizeof_float
        ell_min = g.attrs["ell_min"]
        ell_max = g.attrs["ell_max"]
        spin_weight = g.attrs.get("spin_weight", spin_weight)
        data_type = g.attrs.get("data_type", data_type)
        shuffle_widths = tuple(g.attrs["shuffle_widths"])
        unshuffle = multishuffle(shuffle_widths, forward=False)
        n_modes = ell_max * (ell_max + 2) - ell_min ** 2 + 1
        i1 = n_times * sizeof_float
        i2 = i1 + n_times * sizeof_complex * n_modes
        uncompressed_data = bz2.decompress(g["data"][...])
        t = np.frombuffer(uncompressed_data[:i1], dtype=np.uint64)
        data_tmp = np.frombuffer(uncompressed_data[i1:i2], dtype=np.uint64)
        log_frame = np.frombuffer(uncompressed_data[i2:], dtype=np.uint64)

    # Unshuffle the raw data
    t = unshuffle(t)
    data_tmp = unshuffle(data_tmp)
    log_frame = unshuffle(log_frame)

    # Reshape and re-interpret the data
    t = t.view(np.float64)
    data_tmp = data_tmp.reshape((-1, n_times)).T
    log_frame = log_frame.reshape((-1, n_times)).T.copy().view(np.float64)

    # Because of the weirdness of complex types, reshaping, and C/F order, we
    # need to create this with the layout we will eventually want, and then
    # read data into it.
    data = np.empty((n_times, n_modes), dtype=complex).view(np.uint64)

    # Un-XOR the data
    xor(t, reverse=True, preserve_dtype=True, out=t)
    xor(data_tmp, reverse=True, axis=0, out=data)
    xor(log_frame, reverse=True, preserve_dtype=True, axis=0, out=log_frame)

    frame = np.exp(quaternionic.array(np.insert(log_frame, 0, 0.0, axis=1)))

    if spin_weight is None:
        warning = f"Spin weight has not been provided for {h5_path}.  This may result in errors with some functions."
        warnings.warn(warning)

    w = WaveformModes(
        data.view(np.complex128),
        time=t,
        time_axis=0,
        modes_axis=1,
        frame=frame,
        frame_type="corotating",
        data_type=data_type,
        m_is_scaled_out=True,
        r_is_scaled_out=True,
        ell_min=ell_min,
        ell_max=ell_max,
        spin_weight=spin_weight,
    )
    w.convert_from_conjugate_pairs()

    if transform_to_inertial:
        w = w.to_inertial_frame()

    w.json_data = json_data
    w.log_frame = log_frame

    return w
