"""A core utility function for downloading efficiently and robustly"""

def download_file(url, path, progress=False, if_newer=True):
    """Download large file efficiently from url into path

    Parameters
    ----------
    url : str
        The URL to download from.  Redirects are followed.
    path : {str, pathlib.Path}
        Path to the file in which the download will be stored.  If this is an
        existing directory or ends in a path separator, the "path" component of the
        URL will be used as the file name, and the full directory path will be
        created.
    progress : bool, optional
        If True, and a nonzero Content-Length header is returned, a progress bar
        will be shown during the download.
    if_newer : {bool, datetime, pathlib.Path}, optional
        If True (the default), the file will only be downloaded if the version on
        the server is newer than the "mtime" of the local version.  If this flag is
        False, or there is no local version, or the server does not reply with a
        'Last-Modified' header, the file is downloaded as usual.  If a datetime
        object is passed, it is used instead of the local file's mtime.  If a Path
        object is passed, its mtime is used instead of the output path's, and this
        path is returned if it is newer than the server's file.

    Returns
    -------
    local_filename : pathlib.Path

    """
    import functools
    import pathlib
    import os
    import shutil
    import urllib.parse
    import requests
    from tqdm.auto import tqdm
    from datetime import datetime, timezone

    url_path = urllib.parse.urlparse(url).path
    path = pathlib.Path(path).expanduser().resolve()
    if path.is_dir():
        path = path / url_path[1:]  # May have some new directories
    directory = path.parent
    filename = path.name
    directory.mkdir(parents=True, exist_ok=True)
    if not os.access(str(directory), os.W_OK) or not directory.is_dir():
        raise ValueError(f"Path parent '{directory}' is not writable or is not a directory")
    local_filename = directory / filename

    r = requests.get(url, stream=True, allow_redirects=True)
    if r.status_code != 200:
        print(f"An error occurred when trying to access <{url}>.")
        try:
            print(r.json())
        except Exception:
            pass
        r.raise_for_status()
        raise RuntimeError()  # Will only happen if the response was not strictly an error

    if if_newer and "Last-Modified" in r.headers:
        remote_timestamp = datetime.strptime(
            r.headers["Last-Modified"], "%a, %d %b %Y %H:%M:%S GMT"
        ).replace(tzinfo=timezone.utc)
        if isinstance(if_newer, datetime):
            local_timestamp = if_newer
        elif isinstance(if_newer, pathlib.Path) and if_newer.exists():
            local_timestamp = datetime.fromtimestamp(if_newer.stat().st_mtime, timezone.utc)
        elif local_filename.exists():
            local_timestamp = datetime.fromtimestamp(local_filename.stat().st_mtime, timezone.utc)
        else:
            local_timestamp = remote_timestamp  # Just to make the next condition evaluate to False
        if local_timestamp > remote_timestamp:
            if progress:
                print(f"Skipping download from '{url}' because local file is newer")
            if isinstance(if_newer, pathlib.Path) and if_newer.exists():
                return if_newer
            return local_filename

    file_size = int(r.headers.get('Content-Length', 0))
    r.raw.read = functools.partial(r.raw.read, decode_content=True)

    output_path = local_filename.parent / (local_filename.name + '.part')
    try:
        with output_path.open("wb") as f:
            if progress and file_size:
                desc = "(Unknown total file size)" if file_size == 0 else ""
                print(f"Downloading to {path}:", flush=True)
                with tqdm.wrapattr(r.raw, "read", total=file_size, desc=desc, dynamic_ncols=True) as r_raw:
                    shutil.copyfileobj(r_raw, f)
            else:
                shutil.copyfileobj(r.raw, f)
    except Exception as e:
        raise RuntimeError(f"Failed to download {url} to {local_filename}; original file remains") from e
    else:
        output_path.replace(local_filename)
    finally:
        try:
            output_path.unlink()  # missing_ok is only available in python 3.8
        except FileNotFoundError:
            pass

    return local_filename
