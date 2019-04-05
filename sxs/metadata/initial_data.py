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
    untarred_file_path = os.path.join(system, 'EvID', 'ID_Init_FuncTrans.txt')

    if os.path.isfile(tar_file_path):
        # Open the tarfile for reading
        with tarfile.open(tar_file_path) as tgz:
            if 'ID_Init_FuncTrans.txt' not in tgz.getnames():
                # The file is not even present for very old systems
                return False
            with tgz.extractfile('ID_Init_FuncTrans.txt') as extracted_file:
                file = ''.join([line.decode() for line in extracted_file.readlines()]).strip()
    elif os.path.isfile(untarred_file_path):
        with open(untarred_file_path, 'r') as f:
            file = ''.join([line.decode() for line in f.readlines()]).strip()
    else:
        if raise_on_missing_file:
            raise FileNotFoundError('Could not find\n\t"{0}"\nor\n\t"{1}".'.format(tar_file_path, untarred_file_path))
        return False

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


def get_initial_data_details(system, raise_on_missing_info=False):
    """Look for initial-data type and construction date in EvID/ID_Params.perl"""
    import re
    import os.path
    from warnings import warn

    pattern1 = re.compile(r'"(?P<id_type>.*?) constructed on(?: or before)? (?P<datetime>.*?),? using')
    pattern2 = re.compile(r'"(?P<id_type>.*?) converted to modern format by .*? on (?P<datetime>.*?)\.')
    id_params_path = os.path.join(system, 'EvID', 'ID_Params.perl')

    if not os.path.exists(id_params_path):
        message = "Could not find ID file {0}".format(id_params_path)
        if raise_on_missing_info:
            raise ValueError(message)
        else:
            warn(message)
        return '', ''
    
    with open(id_params_path, 'r') as f:
        for line in f:
            if "$ID_Origin" in line:
                results = pattern1.search(line)
                if not results:
                    results = pattern2.search(line)
                if results:
                    datetime = results['datetime'] if '???' not in results['datetime'] else 'NaT'
                    return results['id_type'], datetime
    
    message = "Could not find ID type and datetime in {0}".format(id_params_path)
    if raise_on_missing_info:
        raise ValueError(message)
    else:
        warn(message)
        return '', ''
