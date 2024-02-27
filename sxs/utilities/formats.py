"""Function to inspect file formats"""


def file_format(file, group=None):
    """Return format stored in file

    If a format is specified, we assume that it is a top-level attribute or member
    named "sxs_format" or just "format".  If neither exists, return None; the
    calling function should check for this possibility.

    Parameters
    ----------
    file : file-like object, string, or pathlib.Path
    group : str, optional
        If the file is an HDF5 or JSON file, this is the group to inspect.

    Returns
    -------
    format : {None, str}
        If the format is not specified, return None.  Otherwise, return the format
        as a string.

    """
    from pathlib import Path
    import json
    import h5py
    if isinstance(file, (str, Path)):
        path = Path(file).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"Could not find {path}")
        if h5py.is_hdf5(str(path)):
            with h5py.File(str(path), "r") as f:
                if group is not None:
                    f = f[group]
                return f.attrs.get("sxs_format", f.attrs.get("format", None))
        else:
            try:
                with path.open("r") as f:
                    f_json = json.load(f)
                    if group is not None:
                        f_json = f_json[group]
                    return f_json.get("sxs_format", f_json.get("format", None))
            except Exception as e:
                raise ValueError(f"Failed to interpret '{path}' as HDF5 or JSON file") from e
    elif hasattr(file, "read"):
        try:
            file.seek(0)
            mode = "binary" if isinstance(file.read(0), bytes) else "text"
            if mode == "binary":
                # This file is open in binary mode, so it must be an HDF5 file
                with h5py.File(file, "r") as f:
                    if group is not None:
                        f = f[group]
                    return f.attrs.get("sxs_format", f.attrs.get("format", None))
            else:
                # This file is open in text mode, so it must be a JSON file
                f_json = json.load(file)
                if group is not None:
                    f_json = f_json[group]
                return f_json.get("sxs_format", f_json.get("format", None))
        except Exception as e:
            if mode == "binary":
                if group is None:
                    raise ValueError(
                        "\nInput `file` is open in binary mode, but is not an open HDF5 file."
                        "\nIf this is a JSON file, open it in text mode: `open(file, 'r')`"
                    ) from e
                else:
                    raise ValueError(
                        "\nInput `file` is open in binary mode, but is not an open HDF5 file,"
                        f"\nor the group '{group}' does not exist in the file."
                        "\nIf this is a JSON file, open it in text mode: `open(file, 'r')`"
                    ) from e
            else:
                if group is None:
                    raise ValueError(
                        "\nInput `file` is open in text mode, but is not an open JSON file."
                        "\nIf this is an HDF5 file, open it in binary mode: `open(file, 'rb')`"
                    ) from e
                else:
                    raise ValueError(
                        "\nInput `file` is open in text mode, but is not an open JSON file,"
                        f"\nor the group '{group}' does not exist in the file."
                        "\nIf this is an HDF5 file, open it in binary mode: `open(file, 'rb')`"
                    ) from e
        finally:
            file.seek(0)
    else:
        raise TypeError(f"Input `file` has type {type(file)}; it should be str, pathlib.Path, or a file-like object")
