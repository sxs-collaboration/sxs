"""Functions to load and save waveforms in RPDMB format"""

import warnings
import tempfile
import contextlib
import pathlib
import numbers
import bz2
import json
import numpy as np
import h5py
import quaternionic

from ...utilities import default_shuffle_widths, md5checksum, diff, xor, multishuffle, version_info
from ... import Metadata, __version__
from .. import WaveformModes


sxs_formats = ["rotating_paired_diff_multishuffle_bzip2", "rpdmb", "RPDMB"]


def save(
        w, file_name=None, file_write_mode="w",
        L2norm_fractional_tolerance=1e-10, log_frame=None,
        shuffle_widths=default_shuffle_widths, convert_to_conjugate_pairs=True,
        compression=bz2, diff=diff, formats=None, verbose=True, allow_existing_group=False,
        version_info_update=None, max_phase_per_timestep=None
):
    """Save a waveform in RPDMB format

    This function converts the data to "rotating paired diff
    multishuffle bzip2" format.  In particular, it uses the corotating
    frame, and zeroes out bits at high precision to allow for optimal
    compression while maintaining the requested tolerance.

    Optionally, with appropriate inputs, can provide different
    formats.

    Parameters
    ----------
    w : WaveformModes
        A waveform in either the inertial or corotating frame
    file_name : str
        Relative or absolute path to the output HDF5 file.  If this
        string contains `'.h5'` but does not *end* with that, the
        remainder of the string is taken to be the group within the
        HDF5 file in which the data should be stored.  Also note that
        a JSON file is created in the same location, with `.h5`
        replaced by `.json` (and the corresponding data is stored
        under the `group` key if relevant).  To retrieve just the
        return values, or for testing purposes, this argument may be
        `None`, in which case a temporary directory is used, just to
        test how large the output will be; it is deleted immediately
        upon returning.
    file_write_mode : str, optional
        One of the valid [file modes for
        h5py](https://docs.h5py.org/en/stable/high/file.html#opening-creating-files).
        Default value is `"w"`, which overwrites any existing file.
        If writing into a group, you may prefer `"a"`, which will just
        ensure the file exists without erasing it first.
    L2norm_fractional_tolerance : float, optional
        Tolerance passed to `WaveformModes.truncate`; see that
        function's docstring for details.  Default value is 1e-10.
    log_frame : array of quaternions, optional
        If this argument is given the waveform must be in the
        corotating frame, and the given data will be used as the
        logarithmic frame data.  If this argument is `None` (the
        default), this will be calculated when the waveform is
        transformed to the corotating frame, or simply taken directly
        from the waveform if it is already corotating.
    shuffle_widths : iterable of ints, optional
        See `sxs.utilities.multishuffle` for details.  The default
        value is `default_shuffle_widths`.  Note that if
        `L2norm_fractional_tolerance` is 0.0, this will be ignored and
        the standard HDF5 shuffle option will be used instead.
    convert_to_conjugate_pairs : bool, optional
        If `True` (the default), the data is converted to
        conjugate-pair form.
    compression : module, optional
        Compression module (or any other object) with `compress` and
        `decompress` methods.  Default value is `bz2`, which is part
        of python's standard library, and performs bzip2
        de/compression.  Note that whatever this is must be available
        to the `load` code, so this should probably be oneof the
        built-in packages: `bz2`, `lmza`, `gzip`, or `zlib`.
    diff : function, optional
        Function to compare successive values.  Defaults to `diff`,
        which is the floating-point difference.  Other possibilities
        include `diffInt` (which reinterprets the data as integers
        before taking the difference) and `xor` (which is the previous
        default); neither works as well with our data at the time of
        this writing.
    formats : array[str], optional
        Possible names of the format.  Defaults to variations on
        `RPDMB`.
    allow_existing_group : bool, optional
        If `True`, allow the group being written into to exist
        already; otherwise (and by default), if the group exists, an
        error will be raised.
    version_info_update : dict, optional
        If present, the `version_info` field in the output JSON file
        will be updated with this information.  Any existing version
        info that is not found in this input will remain unchanged,
        but any that are in this input will be either altered or
        added.
    max_phase_per_timestep : float, optional
        Maximum phase change per time step.  If given, the approximate
        phase change per time step is calculated when transforming to
        the corotating frame.  If it exceeds this value, a ValueError
        is raised.  A sensible value here is probably 10 or so — high
        enough that random spikes from junk radiation or late ringdown
        noise won't trigger it, but low enough that the integration
        will still proceed reasonably quickly for anything below this
        threshold.

    Returns
    -------
    w_out : WaveformModes
        The output data, after conversion to the corotating frame,
        pairing of opposite `m` modes, and differencing (but not
        shuffling).
    log_frame : array of quaternions
        The actual `log_frame` data stored in the file, and used to
        transform to the corotating frame if that was done inside this
        function.

    Note that the returned data are *as stored in the file*.
    Specifically, they are presented as various types of `float` data,
    but have been differenced, which makes them invalid as floats; you
    will see many NaNs and other nonsensical values unless you reverse
    the process.

    """
    if formats is None:
        formats = sxs_formats

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

    if L2norm_fractional_tolerance is None:
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
            try:
                w, log_frame = w.to_corotating_frame(
                    tolerance=L2norm_fractional_tolerance,
                    z_alignment_region=z_alignment_region,
                    truncate_log_frame=True,
                    max_phase_per_timestep=max_phase_per_timestep
                )
            except ValueError as e:
                print("\nError transforming to corotating frame:")
                raise e
            log_frame = log_frame.ndarray[:, 1:]
        if w.frame_type != "corotating":
            raise ValueError(
                f"Frame type of input waveform must be 'corotating' or 'inertial'; it is {w.frame_type}"
            )

        # Convert mode structure to conjugate pairs
        if convert_to_conjugate_pairs:
            w.convert_to_conjugate_pairs()

        # Set bits below the desired significance level to 0
        if L2norm_fractional_tolerance > 0.0:
            w.truncate(tol=L2norm_fractional_tolerance)

        # Compute log(frame)
        if log_frame is None:
            log_frame = np.log(w.frame).ndarray[:, 1:]

        # Change -0.0 to 0.0 (~.5% compression for non-precessing systems)
        w.t += 0.0
        np.add(w.data, 0.0, out=w.data)  # Use `.data` so WaveformModes doesn't override scalar addition
        log_frame += 0.0

        # diff successive instants in time
        t = xor(w.t)
        data = diff(w.data.view(float), axis=w.time_axis)
        log_frame = diff(log_frame, axis=0)

    # Make sure we have a place to keep all this
    with contextlib.ExitStack() as context:
        if file_name is None:
            temp_dir = context.enter_context(tempfile.TemporaryDirectory())
            h5_path = pathlib.Path(f"{temp_dir}") / "test.h5"
        elif verbose:
            print(f'Saving H5 to "{h5_path}"')

        # Write the H5 file
        with h5py.File(h5_path, file_write_mode) as f:
            # If we are writing to a group within the file, create it
            if group is not None:
                if allow_existing_group and group in f:
                    g = f[group]
                else:
                    g = f.create_group(group)
            else:
                g = f
            if L2norm_fractional_tolerance is not None:
                g.attrs["sxs_format"] = formats[0]
                g.attrs["n_times"] = w.n_times
                g.attrs["ell_min"] = w.ell_min
                g.attrs["ell_max"] = w.ell_max
                g.attrs["shuffle_widths"] = np.array(shuffle_widths, dtype=np.uint8)
                data = np.void(
                    compression.compress(
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
                g.attrs["sxs_format"] = f"{formats[0]}"
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
                "sxs_format": formats[0],
                "data_info": {
                    "data_type": w.data_type,
                    "m_is_scaled_out": w._metadata.get("m_is_scaled_out", True),
                    "r_is_scaled_out": w._metadata.get("r_is_scaled_out", True),
                    "spin_weight": int(w.spin_weight),
                    "ell_min": int(w.ell_min),
                    "ell_max": int(w.ell_max),
                },
                "version_info": version_info(),  # see below for "spec_version_hist"
                # see below for "validation"
                # see below for "modifications"
            }

            version_hist = getattr(w, "version_hist", w._metadata.get("version_hist", None))
            if version_hist is not None:
                json_data["version_info"]["spec_version_hist"] = version_hist
            if version_info_update is not None:
                json_data["version_info"].update(version_info_update)

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
            if verbose:
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


def load(
        file_name, ignore_validation=None, check_md5=True,
        transform_to_inertial=True, convert_from_conjugate_pairs=True,
        compression=bz2, diff=diff, formats=None, metadata=None,
        **kwargs
):
    """Load a waveform in RPDMB format

    Parameters
    ----------
    file_name : str
        Relative or absolute path to the input HDF5 file.  If this
        string contains but does not *end* with `'.h5'`, the remainder
        of the string is taken to be the group within the HDF5 file in
        which the data is stored.  Also note that a JSON file is
        expected in the same location, with `.h5` replaced by `.json`
        (and the corresponding data must be stored under the `group`
        key if relevant).
    ignore_validation : bool or None, optional
        Validation checks the corresponding JSON file for (1)
        existence, (2) number of time steps, (3) H5 file size, and (4)
        H5 file MD5 checksum (if `check_md5` is `True`).  If this key
        is `False`, all of this will be ignored; if `True`, a
        `ValueError` will be raised if any of these checks fails; if
        `None`, warnings will be issued, but the function will
        continue as usual.
    check_md5 : bool, optional
        Default is `True`.  See `ignore_validation` for explanation.
    transform_to_inertial : bool, optional
        If `True`, the output waveform is transformed to the inertial
        frame; otherwise it is left in the frame used in the file.
    convert_from_conjugate_pairs : bool, optional
        If `False`, the data is left in its conjugate-pair form;
        default is `True`.
    compression : module, optional
        Compression module (or any other object) with `compress` and
        `decompress` methods.  Default value is `bz2`, which is part
        of python's standard library, and performs bzip2
        de/compression.  Note that whatever this is must be available
        to the `load` code, so this should probably be oneof the
        built-in packages: `bz2`, `lmza`, `gzip`, or `zlib`.
    diff : function, optional
        Function to compare successive values.  Defaults to `diff`,
        which is the floating-point difference.  Other possibilities
        include `diffInt` (which reinterprets the data as integers
        before taking the difference) and `xor` (which is the previous
        default); neither works as well with our data at the time of
        this writing.
    formats : array[str], optional
        Possible names of the format.  Defaults to variations on
        `RPDMB`.

    Keyword parameters
    ------------------
    data_type : str, optional
        Describes the type of data — such as "h" or "psi4".  Default
        is "unknown".
    m_is_scaled_out : bool, optional
        Default is True
    r_is_scaled_out : bool, optional
        Default is True
    spin_weight : int, optional
        Default is `None`.
    drop_times_before : {float,"begin","reference","merger"}, optional
        If given, all times before the given time will be dropped.  If
        a string is given, the corresponding time will be used.  The
        "begin" option represents the earliest time in the data, so no
        data will be dropped; "reference" is the `reference_time`
        reported in the metadata.  If `None`, the default, t=0 will be
        used if present in the data, and the first time step
        otherwise.
    metadata : Metadata, optional
        If given, this metadata will be used instead of attempting to
        load the metadata from an accompanying file.

    Note that the keyword parameters will be overridden by
    corresponding entries in the JSON file, if they exist.  If the
    JSON file does not exist, any keyword parameters not listed above
    will be passed through as the `json_data` field of the returned
    waveform.

    """
    if formats is None:
        formats = sxs_formats

    def invalid(message):
        if ignore_validation:
            pass
        elif ignore_validation is None:
            warnings.warn(message)
        else:
            raise ValueError(message)

    file_name_str = str(file_name)
    group = kwargs.get("group", None)
    if ".h5" in file_name_str and not file_name_str.endswith(".h5"):
        file_name_str, group = file_name_str.split(".h5")
    if group == "/":
        group = None

    # At this point, file_name_str may or may not have an h5 suffix
    # (because the user may have passed in a file without the suffix,
    # or because there was an .h5 in the middle of the filename.)
    # But that's ok, we will use with_suffix below.

    # In a common use case, the h5 and json files will be named
    # bla.h5 and bla.json as expected, but they will be git-annex symlinks.
    # This means that the files pointed to will have strange names, like:
    # bla.h5   -> /some/crazy/path/some_crazy_hash.h5
    # bla.json -> /some/different/crazy/path/a_different_crazy_hash.json
    # where the paths are determined by git-annex and the hashes are basically
    # SHA256 hashes of the contents.
    #
    # This means we need to change the suffix *before* the resolve() call.
    h5_path = pathlib.Path(file_name_str).with_suffix(".h5").expanduser().resolve()
    json_path = pathlib.Path(file_name_str).with_suffix(".json").expanduser().resolve()
    metadata_path = (pathlib.Path(file_name_str).parent / "metadata").expanduser().resolve()

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
        if sxs_format not in formats:
            invalid(
                f"\nThe `sxs_format` found in JSON file is '{sxs_format}';\n"
                f"it should be one of\n"
                f"    {formats}."
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
        if sxs_format not in formats:
            raise ValueError(
                f'The `sxs_format` found in H5 file is "{sxs_format}"; it should be one of\n'
                f"    {formats}."
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
        uncompressed_data = compression.decompress(g["data"][...])
        t = np.frombuffer(uncompressed_data[:i1], dtype=np.uint64)
        data_tmp = np.frombuffer(uncompressed_data[i1:i2], dtype=np.uint64)
        log_frame = np.frombuffer(uncompressed_data[i2:], dtype=np.uint64)

    # Unshuffle the raw data
    t = unshuffle(t)
    data_tmp = unshuffle(data_tmp)
    log_frame = unshuffle(log_frame)

    # Reshape and re-interpret the data
    t = t.view(np.float64)
    data_tmp = data_tmp.reshape((-1, n_times)).T.copy()
    log_frame = log_frame.reshape((-1, n_times)).T.copy().view(np.float64)

    # Because of the weirdness of complex types, reshaping, and C/F order, we
    # need to create this with the layout we will eventually want, and then
    # read data into it.
    data = np.empty((n_times, n_modes), dtype=complex)

    # Un-diff the data
    xor(t, reverse=True, preserve_dtype=True, out=t)
    diff(data_tmp, reverse=True, axis=0, out=data)
    diff(log_frame, reverse=True, preserve_dtype=True, axis=0, out=log_frame)

    frame = np.exp(quaternionic.array.from_vector_part(log_frame))

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
        m_is_scaled_out=m_is_scaled_out,
        r_is_scaled_out=r_is_scaled_out,
        ell_min=ell_min,
        ell_max=ell_max,
        spin_weight=spin_weight,
    )

    if convert_from_conjugate_pairs:
        w.convert_from_conjugate_pairs()

    if transform_to_inertial:
        w = w.to_inertial_frame()

    if metadata is None:
        try:
            metadata = Metadata.from_file(metadata_path)
        except ValueError as e:
            invalid(f"\n{e},\nbut one is expected for this data format.")

    dtb = kwargs.pop("drop_times_before", 0)
    if dtb=="begin":
        i0 = 0
    elif dtb==0:
        i0 = np.argmin(w.t < 0)
    elif dtb=="reference":
        if metadata is None:
            raise ValueError(f"Metadata is required if `drop_times_before` is set to 'reference'")
        i0 = np.argmin(w.t < metadata.reference_time)
    elif isinstance(dtb, numbers.Real):
        i0 = np.argmin(w.t < dtb)
    elif dtb=="merger":
        i0 = w.max_norm_index()
    else:
        raise ValueError(f"Invalid value for `drop_times_before`: {dtb}")
    if i0 != 0:
        w = w[i0:]

    w.json_data = json_data
    w.log_frame = log_frame
    if metadata is not None:
        w.metadata = metadata

    return w


def load_time(
        file_name, ignore_validation=None, check_md5=True,
        transform_to_inertial=True, convert_from_conjugate_pairs=True,
        compression=bz2, diff=diff, formats=None,
        **kwargs
):
    """Load time data from a waveform in RPDMB format

    This function is just an abbreviated form of the `load` function, without any
    type of validation, or conversion of the mode or frame data; only the `time`
    dataset is returned.

    For compatibility, this function has the same calling signature as `load`,
    though many of the parameters are meaningless here.  Most importantly, the
    `file_name` is processed in the same way, meaning that groups can also be
    specified.  Otherwise, the only possibly relevant parameter is `compression`.

    """
    file_name_str = str(file_name)
    group = None
    if ".h5" in file_name_str and not file_name_str.endswith(".h5"):
        file_name_str, group = file_name_str.split(".h5")
    if group == "/":
        group = None

    h5_path = pathlib.Path(file_name_str).with_suffix(".h5").expanduser().resolve()

    with h5py.File(h5_path, "r") as f:
        if group is not None:
            g = f[group]
        else:
            g = f

        # Read the raw data
        n_times = g.attrs["n_times"]
        sizeof_float = 8
        i1 = n_times * sizeof_float
        shuffle_widths = tuple(g.attrs["shuffle_widths"])
        unshuffle = multishuffle(shuffle_widths, forward=False)
        uncompressed_data = compression.decompress(g["data"][...])
        t = np.frombuffer(uncompressed_data[:i1], dtype=np.uint64)

    # Unshuffle the raw data
    t = unshuffle(t)

    # Reshape and re-interpret the data
    t = t.view(np.float64)

    # Un-diff the data
    xor(t, reverse=True, preserve_dtype=True, out=t)

    return t
