"""Sub-submodule to validate BHNS systems"""


def horizons(path, short_circuit=False, validity=None):
    """Validate Horizons.h5 file for BHNS simulation

    """
    import os.path
    import numpy as np
    import h5py
    from . import _Validity

    v = validity or _Validity(short_circuit)  # If validity is None, this line constructs it

    raise NotImplementedError()

    return v.valid


def waveforms(path, *args, short_circuit=False, validity=None, extrapolated=None):
    """Validate waveform files for BHNS simulation

    """
    import re
    import os.path
    import numpy as np
    import h5py
    from . import _Validity

    waveforms = list(set(args))  # Make sure the input file names are all unique

    v = validity or _Validity(short_circuit)  # If validity is None, this line constructs it

    raise NotImplementedError()
    
    return v.valid


def simulation(path='.', short_circuit=False):
    """Check that a directory contains a valid BHNS simulation

    Parameters
    ==========
    path: string, optional [defaults to working directory]
        Absolute or relative path to directory to be validated (or a file within that directory)
    short_circuit: bool, optional [defaults to False]
        If True, raise a ValueError as soon as any invalid condition is encountered.  Otherwise,
        continue checking each requirement and printing messages about any invalid conditions, but
        return False.

    Returns
    =======
    validity: bool
        True or False, depending on whether the directory does or does not (respectively) contain a
        valid BHNS simulation.

    Raises
    ======
    ValueError
        If input argument `short_circuit` is True and any invalidity is found.

    """
    import re
    import os.path
    import numpy as np
    from .. import sxs_identifier_regex
    from ..metadata import Metadata
    from . import _Validity

    v = _Validity(short_circuit)

    # Ensure that the input `path` is a directory that exists
    path = os.path.normpath(path)
    if not os.path.isdir(path) and os.path.isfile(path):
        path = os.path.dirname(path)
    if not os.path.isdir(path):
        v.invalid('The input path name "{0}" does not appear to be a directory.'.format(path))
        return False

    # Check that common-metadata.txt exists, can be parsed, and contains an SXS ID
    common_metadata_path = os.path.join(path, 'common-metadata.txt')
    if not os.path.isfile(common_metadata_path):
        v.invalid('Could not find "common-metadata.txt" file in "{0}".'.format(path))
    try:
        common_metadata = Metadata.from_txt_file(common_metadata_path)
    except:
        v.invalid('Could not parse "{0}".'.format(common_metadata_path))
        common_metadata = {}
    if not 'alternative_names' in common_metadata:
        v.invalid('No "alternative_names" key in "{0}".'.format(common_metadata_path))
    else:
        alternative_names = common_metadata['alternative_names']
        if isinstance(alternative_names, list):
            alternative_names = ' '.join(alternative_names)
        match = re.search(sxs_identifier_regex, alternative_names)
        if not match or not match['sxs_identifier']:
            v.invalid('No valid SXS ID found in "{0}".'.format(common_metadata_path))

    raise NotImplementedError()

    return v.valid


_main = simulation  # This will be the default executable when run via command line
