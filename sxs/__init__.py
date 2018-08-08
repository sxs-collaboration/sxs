from __future__ import division, print_function, absolute_import

__doc_title__ = "SXS python code"
__doc__ = "A collection of python code used by the SXS collaboration"

from ._version import __version__

import sxs.doxygen
import sxs.metadata
import sxs.references
import sxs.zenodo

sxs_identifier_regex = r'(?P<sxs_identifier>SXS:(?P<simulation_type>BBH|BHNS|NSNS):(?P<sxs_number>[0-9]*))'


def sxs_id(s):
    """Return the SXS ID contained in the input string

    An SXS ID is anything that matches the following regular expression:

        SXS:(BBH|BHNS|NSNS):[0-9]*

    If no match is found, the empty string is returned.  If multiple matches are found, only the
    first is returned.

    """
    import re
    m = re.search(sxs_identifier_regex, s)
    if m:
        return m['sxs_identifier']
    else:
        return ''


def uses_new_initial_data(system, raise_on_missing_file=False):
    """Look for evidence that new initial-data solver was used for this system
    
    This function looks in EvID/ID_Files.tgz for a file named ID_Init_FuncTrans.txt.
    If that file is present, it looks for keys named dTx, etc.  If those keys are
    present and any are nonzero, then this system was created with initial data that
    was adjusted to eliminate residual linear momentum, as described in this paper:
    <https://arxiv.org/abs/1506.01689>.
    
    Parameters
    ==========
    system: string
        Absolute or relative path to the top directory of a simulation.  This directory
        is expected to contain the EvID subdirectory, which in turn is expected to
        contain the 'ID_Files.tgz' file.
    raise_on_missing_file: bool [defaults to False]
        If True, raise an error if 'EvID/ID_Files.tgz' is not found; otherwise just
        print a warning and return False.
        
    Returns
    =======
    bool:
        Returns True if there is evidence that this system was created using the new
        initial data, and False otherwise.
    
    """
    import os.path
    import tarfile
    import ast
    import warnings
    
    tar_file_path = os.path.join(system, 'EvID', 'ID_Files.tgz')

    # Open the tarfile for reading
    try:
        tgz = tarfile.open(tar_file_path)
    except FileNotFoundError:
        if raise_on_missing_file:
            raise
        else:
            warnings.warn('\n\tCould not find "{0}".'.format(tar_file_path))
            return False
    try:
        # Try to extract ID_Init_FuncTrans.txt from the tarfile
        extracted_file = tgz.extractfile('ID_Init_FuncTrans.txt')
    except KeyError:
        # The file is not even present for very old systems
        return False
    # Convert extracted file to a single string
    file = ''.join([line.decode() for line in extracted_file.readlines()]).strip()
    # Separate the file into a list of strings like 'Time= 0.0'
    assignments = file.split(';')
    # Now, make each such string into a key-value pair in this dictionary
    d = {
        kv[0].strip(): ast.literal_eval(kv[1].strip())
        for assignment in assignments
        for kv in [assignment.split('=')]
        if len(kv)>1
    }
    # Check for nonzero initial momentum correction
    return abs(d.get('dTx', 0.0))+abs(d.get('dTy', 0.0))+abs(d.get('dTz', 0.0)) > 0.0
