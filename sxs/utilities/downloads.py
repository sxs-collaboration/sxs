def download_file(url, path, progress=False):
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

    """
    import functools
    import pathlib
    import os
    import shutil
    import urllib.parse
    import requests
    import tqdm

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
        except:
            pass
        r.raise_for_status()
        raise RuntimeError()  # Will only happen if the response was not strictly an error

    file_size = int(r.headers.get('Content-Length', 0))
    r.raw.read = functools.partial(r.raw.read, decode_content=True)
    with local_filename.open("wb") as f:
        if progress:
            desc = "(Unknown total file size)" if file_size == 0 else ""
            print(f"Downloading to {path}:", flush=True)
            with tqdm.tqdm.wrapattr(r.raw, "read", total=file_size, desc=desc, ncols=79) as r_raw:
                shutil.copyfileobj(r_raw, f)
        else:
            shutil.copyfileobj(r.raw, f)

    return local_filename
