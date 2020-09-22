"""Utilities that inherit from dict"""

class KeyPassingDict(dict):
    """Subclass of dict that enables extended keys

    When indexing this object, the key is split into two pieces: the part preceding
    the first "/" and everything after it.  The second part is passed to the
    __getitem__ method of whatever is returned by indexing this object with the
    first part of the key.  This is helpful when trying to emulate the behavior of
    HDF5 files, for example.

    """
    def __getitem__(self, key):
        parts = key.split("/", 1)
        value = super().__getitem__(parts[0])
        if len(parts) > 1:
            return value[parts[1]]
        else:
            return value
