from contextlib import contextmanager


@contextmanager
def lock_file_manager(file_name, mode='a+'):
    """Create and lock `file_name` as context manager
    
    Use this manager as

    >>> with lock_file_manager('my_file.lock'):
    ...     # do stuff here

    The given file will be created, if needed.  While the code inside the `with` block is executing,
    the file is locked using `fnctl.flock`, preventing any other process from accessing it.  As soon
    as that code finishes executing, the file is unlocked and closed.  This is true even if the code
    raises an exception.

    Note that the file descriptor is returned, so it is possible to write to the locked file, as in

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

    This function returns a list of pairs.  Each pair gives (first) the full path to the file, and
    (second) its name relative to the parent directory of `top_directory`.  These are the parameters
    needed to upload a file to Zenodo: first where to find the file, and second where to put it in
    the Zenodo deposit.

    Parameters
    ==========
    top_directory: string
        Absolute or relative path to the top directory from which the recursive search for files
        begins.
    exclude: list of strings
        Each string is compiled as a regular expression.  The path to each directory and file
        relative to `top_directory` is searched for a match, and if found that item is excluded.
        In particular, if a directory matches, no files from that directory will be uploaded.
    include_top_directory_in_name: bool [defaults to True]
        If True, the name of the top_directory (relative to its parent) will be included in the
        output names.

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


def update_checksum_map(checksum_map_file_name='checksums.json', verify_checksums=False,
                        top_directory='.', exclude=[], raise_errors=False, *args, **kwargs):
    """Update (or construct) a map from checksums to file paths

    The map's format is a dictionary with checksums as the keys, and corresponding values are
    dictionaries with at least the keys 'size' (giving the file size in bytes) and 'path' (giving
    the path to a non-link copy of the file).  There may be additional keys in this second
    dictionary, such as a list of files that are simply links to the original, or the modification
    time of the original.  As much as possible, those additional entries will be preserved.

    Parameters
    ==========
    checksum_map_file_name: str [defaults to 'checksums.json']
        Absolute or relative path to JSON file containing the checksums of all the files within
        `top_directory`.  Note that a temporary file with this name plus '_tmp.json' is also written
        to preserve all the changes as the function goes along.
    verify_checksums: bool [defaults to False]
        If True, each checksum will be rechecked.  By default, only file sizes are checked.
    top_directory: string [defaults to '.']
        Absolute or relative path to the top directory from which the recursive search for files
        begins.
    exclude: list of strings [defaults to empty list]
        Each string is compiled as a regular expression.  The path to each directory and file
        relative to `top_directory` is searched for a match, and if found that item is excluded.
        In particular, if a directory matches, no files from that directory will be uploaded.
    raise_errors: bool [defaults to False]
        If True, rather than just issuing a warning, an error will raise an exception, which will
        stop the execution of this function.  Since this function may take a very long time to run,
        it may be preferable to just issue warnings.

    """
    import os.path
    import json
    import warnings
    
    verbosity = kwargs.pop('verbosity', 0)
    if not os.path.isfile(checksum_map_file_name) or not os.path.getsize(checksum_map_file_name) > 1:
        with open(checksum_map_file_name, 'w') as f:
            f.write('{}')
    with open(checksum_map_file_name, 'r') as f:
        checksum_map = json.load(f)
    path_map = {checksum_map[c]['path']: (c, checksum_map[c]) for c in checksum_map}
    changes = {}
    paths_and_names = find_files(top_directory, exclude)
    for path, name in paths_and_names:
        if verbosity > 1:
            print('Processing {0}'.format(path))
        changed = False
        if os.path.isfile(path) and not os.path.islink(path):
            if verbosity > 0:
                print('Checking {0}'.format(path))
            size = os.path.getsize(path)
            if path in path_map:
                checksum, map_entry = path_map[path]
                if map_entry['size'] != size:
                    if map_entry['size'] > 0:
                        error_message = 'Expected size of {0} does not match actual size of {1} for {2}'.format(map_entry['size'], size, path)
                        if raise_errors:
                            raise ValueError(error_message)
                        else:
                            warnings.warn(error_message)
                    else:
                        if verbosity > 0:
                            print('\tThe expected file size has changed.  Re-computing checksum.')
                        old_checksum = checksum
                        checksum = md5checksum(path)
                        map_entry.update({
                            'path': path,
                            'size': size,
                        })
                        checksum_map[checksum] = map_entry
                        changes[checksum] = map_entry
                        changed = True
                        checksum_map.pop(old_checksum, None)  # Remove the old entry from the map
                if verify_checksums and not changed:
                    new_checksum = md5checksum(path)
                    if checksum != new_checksum:
                        error_message = 'Expected checksum {0} does not match actual checksum {1} for {2}'.format(checksum, new_checksum, path)
                        if raise_errors:
                            raise ValueError(error_message)
                        else:
                            warnings.warn(error_message)
            else:
                checksum = md5checksum(path)
                map_entry = {
                    'path': path,
                    'size': size,
                }
                checksum_map[checksum] = map_entry
                changes[checksum] = map_entry
                changed = True
            if changed:
                with open(checksum_map_file_name+'_tmp.json', 'w') as f:
                    json.dump(changes, f, indent=4, separators=(',', ': '))
        else:
            if verbosity > 1:
                print('\tThis is either not a file or just a link; skipping')
    if changes:
        with open(checksum_map_file_name, 'w') as f:
            json.dump(checksum_map, f, indent=4, separators=(',', ': '))
    if verbosity > 1:
        return checksum_map
    else:
        return
