"""Several convenient functions for interacting with files"""

from contextlib import contextmanager


@contextmanager
def lock_file_manager(file_name, mode='a+'):
    """Create and lock `file_name` as context manager
    
    Use this manager as

    >>> with lock_file_manager('my_file.lock'):
    ...     # do stuff here

    The given file will be created, if needed.  While the code inside the `with`
    block is executing, the file is locked using `fnctl.flock`, preventing any
    other process from accessing it.  As soon as that code finishes executing, the
    file is unlocked and closed.  This is true even if the code raises an
    exception.

    Note that the file descriptor is returned, so it is possible to write to the
    locked file, as in

    >>> with lock_file_manager('my_file.lock') as lf:
    ...     print('Hello, world!', file=lf)
    ...     # do other stuff here

    """
    import fcntl
    with open(file_name, mode) as file_descriptor:
        try:
            fcntl.flock(file_descriptor, fcntl.LOCK_EX)
            yield file_descriptor
        finally:
            fcntl.flock(file_descriptor, fcntl.LOCK_UN)


def find_simulation_directories(root):
    """Search for SpEC-simulation output directories under `root`"""
    import os

    # We can't just do this as a list comprehension because of how os.walk works
    simulation_directories = []

    # Now, recursively walk the directories below `root`
    for dirpath, dirnames, filenames in os.walk(root):
        # We look for files named common-metadata.txt to indicate a simulation //may// be present
        if 'common-metadata.txt' in filenames:
            # Confirm that this is a simulation only if Lev* is present
            if any(s.startswith('Lev') for s in dirnames):
                simulation_directories.append(os.path.join(root, dirpath))
                dirnames[:] = []  # Skip any subdirectories

    return simulation_directories


def find_files(top_directory, exclude=[], include_top_directory_in_name=True):
    """Recursively find all files in `top_directory` and give them relative names

    This function returns a list of pairs.  Each pair gives (first) the full path
    to the file, and (second) its name relative to the parent directory of
    `top_directory`.  These are the parameters needed to upload a file to Zenodo:
    first where to find the file, and second where to put it in the Zenodo deposit.

    Parameters
    ----------
    top_directory : string
        Absolute or relative path to the top directory from which the recursive
        search for files begins.
    exclude : list of strings
        Each string is compiled as a regular expression.  The path to each
        directory and file relative to `top_directory` is searched for a match, and
        if found that item is excluded.  In particular, if a directory matches, no
        files from that directory will be uploaded.
    include_top_directory_in_name : bool [defaults to True]
        If True, the name of the top_directory (relative to its parent) will be
        included in the output names.

    """
    import os
    import re
    paths_and_names = []
    exclude = [re.compile(exclusion) for exclusion in exclude]
    top_directory = os.path.abspath(os.path.expanduser(top_directory))
    parent_directory = os.path.dirname(top_directory)
    for root, dirs, files in os.walk(top_directory, topdown=True):
        dirs.sort(key=str.lower)  # Go in case-insensitive alphabetical order
        files.sort(key=str.lower)  # Go in case-insensitive alphabetical order
        for exclusion in exclude:
            for d in dirs:
                if exclusion.search(os.path.relpath(d, top_directory)):
                    dirs.remove(d)
            for f in files:
                if exclusion.search(os.path.relpath(f, top_directory)):
                    files.remove(f)
        for f in files:
            path = os.path.join(root, f)
            if include_top_directory_in_name:
                name = os.path.relpath(path, parent_directory)
            else:
                name = os.path.relpath(path, top_directory)
            paths_and_names.append([path, name])
    return paths_and_names


def md5checksum(file_name):
    """Compute MD5 checksum on a file, even if it is quite large"""
    from hashlib import md5
    hash_md5 = md5()
    with open(file_name, "rb") as f:
        for chunk in iter(lambda: f.read(32768), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
