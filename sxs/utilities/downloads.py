def download_file(url, path, verbosity=0):
    """Download large file efficiently from url into path

    Parameters
    ----------
    url : str
        The URL to download from.  Redirects are followed.
    path : str
        Relative or absolute path to the file in which the download will be stored.
        If this is an existing directory or ends in a path separator, the "path"
        component of the URL will be used as the file name, and the full directory
        path will be created.
    verbosity : int
        If greater than zero, show the response to this request.

    """
    from pathlib import Path
    from urllib.parse import urlparse
    from functools import partial
    from shutil import copyfileobj
    import requests
    from requests_toolbelt.utils import dump

    url_path = urlparse(url).path
    path = Path(path).expanduser().resolve()
    if path.is_dir():
        path = path / url_path[1:]  # May have some new directories
        directory = path.parent
        filename = path.name
        if not directory.exists():
            directory.mkdir(parents=True, exist_ok=True)
        local_filename = directory / filename
    else:
        directory = path.parent
        filename = path.name
        if not directory.exists():
            directory.mkdir(parents=True, exist_ok=True)
        if not filename:
            filename = url_path
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

    r.raw.read = partial(r.raw.read, decode_content=True)
    with open(local_filename, "wb") as f:
        copyfileobj(r.raw, f)
    if verbosity > 0:
        data = dump.dump_all(r)
        print(data.decode("utf-8"))
    return local_filename
