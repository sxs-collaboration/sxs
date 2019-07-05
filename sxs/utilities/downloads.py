def download(url, path, verbosity=0):
    """Download large file efficiently from url into path

    Parameters
    ==========
    url: string
        The URL to download from.  Redirects are followed.
    path: string
        Relative or absolute path to the file in which the download will be stored.  If this is
        an existing directory or ends in a path separator, the "path" component of the URL will
        be used as the file name, and the full directory path will be created.
    verbosity: integer
        If greater than zero, dump the response to this request.

    """
    from shutil import copyfileobj
    from os import makedirs
    from os.path import split, exists, join, isdir
    from functools import partial
    import requests
    try:
        from urllib.parse import urlparse
    except ImportError:
        from urlparse import urlparse
    url_path = urlparse(url).path
    if isdir(path):
        path = join(path, url_path[1:])
        directory, filename = split(path)
        if not exists(directory):
            makedirs(directory)
        local_filename = join(directory, filename)
    else:
        directory, filename = split(path)
        if not exists(directory):
            makedirs(directory)
        if not filename:
            filename = url_path
        local_filename = join(directory, filename)
    r = requests.get(url, stream=True, allow_redirects=True)
    if r.status_code != 200:
        print('An error occurred when trying to access <{0}>.'.format(url))
        try:
            print(r.json())
        except:
            pass
        r.raise_for_status()
        raise RuntimeError()  # Will only happen if the response was not strictly an error
    r.raw.read = partial(r.raw.read, decode_content=True)
    with open(local_filename, 'wb') as f:
        copyfileobj(r.raw, f)
    if verbosity>0:
        from requests_toolbelt.utils import dump
        data = dump.dump_all(r)
        print(data.decode('utf-8'))
    return local_filename
