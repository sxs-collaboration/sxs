"""Sub-submodule to validate BBH systems"""


def horizons(path, short_circuit=False, validity=None):
    """Validate Horizons.h5 file for BBH simulation

    * The file given by `path` must exist and be a readable HDF5 file
    * Its root must contain the groups 'AhA.dir', 'AhB.dir', and 'AhC.dir'
    * Each of those groups must contains these float-valued subgroups with the given shapes
      * ArealMass.dat (N, 2)
      * ChristodoulouMass.dat (N, 2)
      * CoordCenterInertial.dat (N, 4)
      * DimensionfulInertialSpin.dat (N, 4)
      * DimensionfulInertialSpinMag.dat (N, 2)
      * chiInertial.dat (N, 4)
      * chiMagInertial.dat (N, 2)
    * The shape parameter N must be identical within each group, and between 'AhA.dir' and 'AhB.dir'.
    * The first (0th) column of each data is supposed to correspond to time, and therefore must be
      identical within each group, and between 'AhA.dir' and 'AhB.dir'.

    Parameters
    ==========
    path: string
        Absolute or relative path to Horizons.h5 file to be validated
    short_circuit: bool, optional [defaults to False]
        If True, raise a ValueError as soon as any invalid condition is encountered.  Otherwise,
        continue checking each requirement and printing messages about any invalid conditions, but
        return False.  Note that if the `validity` argument is given, this argument is ignored.
    validity: _Validity, optional
        This is the utility object to track validity and print or raise messages on invalidation.
        This can be passed from other functions so that they share validity states.  If not given,
        one will be constructed for use inside this function, using the `short_circuit` argument.

    Returns
    =======
    validity: bool
        True or False, depending on whether the directory does or does not (respectively) contain a
        valid BBH simulation.

    Raises
    ======
    ValueError
        If any invalidity is found and either the input argument `short_circuit` is True or the
        input `validity` object's `short_circuit` attribute is True.

    """
    import os.path
    import numpy as np
    import h5py
    from . import _Validity

    v = validity or _Validity(short_circuit)  # If validity is None, this line constructs it
    groups = ['AhA.dir', 'AhA.dir', 'AhC.dir']
    subgroups = ['ArealMass.dat', 'ChristodoulouMass.dat', 'CoordCenterInertial.dat',
                 'DimensionfulInertialSpin.dat', 'DimensionfulInertialSpinMag.dat',
                 'chiInertial.dat', 'chiMagInertial.dat']
    subgroup_shape1 = [2, 2, 4, 4, 2, 4, 2]

    # Check that the file exists
    if not os.path.isfile(path):
        v.invalid('No file exists in "{0}".'.format(path))
        return False

    # Check that it can be read by h5py
    try:
        h5py.File(path, 'r').close()
    except:
        v.invalid('Could not open "{0}" for reading as an HDF5 file.'.format(path))
        return False

    with h5py.File(path, 'r') as f:
        # Check that 'AhA.dir', 'AhB.dir', and 'AhC.dir' are in the h5 root
        for group in groups:
            if group not in f:
                v.invalid('Group "{0}" not found in "{1}".'.format(group, path))
                return False

        # Check that each subgroup is present in each group, and has dtype==float
        for group in groups:
            for subgroup in subgroups:
                if subgroup not in f[group]:
                    v.invalid('Subgroup "{0}" not found in group "{1}" of "{2}".'.format(
                        subgroup, group, path))
                elif f[group][subgroup].dtype != np.float:
                    v.invalid('Subgroup "{0}" has wrong dtype "{1}" in group "{2}" of "{3}".'.format(
                        subgroup, f[group][subgroup].dtype, group, path))

        # Check that subgroup shapes and time data match expected shapes and time data
        t_AB = f['AhA.dir'][list(f['AhA.dir'])[0]][:, 0]  # We just pick the first time data as the standard
        t_C = f['AhC.dir'][list(f['AhC.dir'])[0]][:, 0]  # We just pick the first time data as the standard
        for group, t in zip(groups, [t_AB, t_AB, t_C]):
            shape0 = t.shape[0]
            for subgroup, shape1 in zip(subgroups, subgroup_shape1):
                if subgroup not in f[group]:
                    continue  # This failure was already checked and noted above
                shape = f[group][subgroup].shape
                if shape != (shape0, shape1):
                    v.invalid('Wrong shape {0} found in "{1}/{2}" of "{3}".  [Expected {4}.]'.format(
                        shape, group, subgroup, path, (shape0, shape1)))
                elif not np.array_equal(f[group][subgroup][:, 0], t):
                    v.invalid('Inconsistent time data found in group "{0}/{1}" of "{2}".'.format(
                        group, subgroup, path))

    return v.valid


def waveforms(path, *args, short_circuit=False, validity=None, extrapolated=None):
    """Validate waveform files for a BBH simulation

    * Each of the waveform data files must be readable as an HDF5 file
    * Have a series of groups at the top level of the file, each containing a waveform
    * Each group must contain mode subgroups named as 'Y_l{ell}_m{m}.dat', where
      * ell ranges from 2 to 8 (inclusive)
      * m ranges from -ell to ell (inclusive)
    * Each subgroup must be a dataset of shape (N, 3) and dtype `float`
    * The first (0th) column in every dataset seen by this function must be identical.  That is
      * Each dataset must have identical time series /within/ each waveform
      * Each dataset must have identical time series /between/ all input waveforms

    In all cases, a top-level group named 'VersionHist.ver' is simply ignored.  Other top-level
    groups are unconstrained if the argument `extrapolated` is None.  If that argument is True, the
    other groups must be exactly

        ['Extrapolated_N2.dir', 'Extrapolated_N3.dir', 'Extrapolated_N4.dir', 'OutermostExtraction.dir']

    If that argument is False, all other groups are assumed to have a name matching 'R[0-9]*.dir'.


    Parameters
    ==========
    path: string
        Absolute or relative path to the directory containing the waveform files to be validated.
    *args: one or more strings
        Names of waveform files to be validated, relative to `path`.
    short_circuit: bool, optional [defaults to False]
        If True, raise a ValueError as soon as any invalid condition is encountered.  Otherwise,
        continue checking each requirement and printing messages about any invalid conditions, but
        return False.  Note that if the `validity` argument is given, this argument is ignored.
    validity: _Validity, optional
        This is the utility object to track validity and print or raise messages on invalidation.
        This can be passed from other functions so that they share validity states.  If not given,
        one will be constructed for use inside this function, using the `short_circuit` argument.
    extrapolated: None or bool, optional [default is None]
        If True, the waveform must have the four extrapolation groups noted above; if False, it must
        have multiple top-level groups with names matching 'R[0-9]*.dir'; if None, no group names
        are required, but all groups are checked.

    Returns
    =======
    validity: bool
        True or False, depending on whether the directory does or does not (respectively) contain a
        valid BBH simulation.

    Raises
    ======
    ValueError
        If any invalidity is found and either the input argument `short_circuit` is True or the
        input `validity` object's `short_circuit` attribute is True.

    """
    import re
    import os.path
    import numpy as np
    import h5py
    from . import _Validity

    waveforms = list(set(args))  # Make sure the input file names are all unique

    v = validity or _Validity(short_circuit)  # If validity is None, this line constructs it

    # Check that each file exists and can be read by h5py
    missing_waveforms = []
    for waveform in waveforms:
        waveform_path = os.path.join(path, waveform)
        if not os.path.isfile(waveform_path):
            v.invalid('No file exists in "{0}".'.format(waveform_path))
            missing_waveforms.append(waveform)
        else:
            try:
                h5py.File(waveform_path, 'r').close()
            except:
                v.invalid('Could not open "{0}" for reading as an HDF5 file.'.format(waveform_path))
                missing_waveforms.append(waveform)

    # Remove any missing waveforms, and check that we still have at least one
    waveforms = [w for w in waveforms if w not in missing_waveforms]
    if len(waveforms) == 0:
        return False

    # Check for top-level groups
    for waveform in waveforms:
        waveform_path = os.path.join(path, waveform)
        with h5py.File(waveform_path, 'r') as f:
            groups = [group for group in f if group != 'VersionHist.ver']
            if not groups:
                v.invalid('No data groups found in "{0}".'.format(waveform_path))
            else:
                if extrapolated is None:
                    pass  # Don't require any particular names
                elif extrapolated:
                    if set(groups) != set(['Extrapolated_N2.dir', 'Extrapolated_N3.dir',
                                           'Extrapolated_N4.dir', 'OutermostExtraction.dir']):
                        v.invalid('Wrong top-level groups for extrapolated file in "{0}".'.format(waveform_path))
                else:
                    if not all(bool(re.match(r'R[0-9]*\.dir', group)) for group in groups):
                        v.invalid('Found incorrectly named top-level groups in finite-radius file "{0}".'.format(waveform_path))

    # Check for correct subgroups and consistent time data
    t = None
    for waveform in waveforms:
        waveform_path = os.path.join(path, waveform)
        with h5py.File(waveform_path, 'r') as f:
            groups = [group for group in f if group != 'VersionHist.ver']
            if not groups:
                continue
            for group in groups:
                for ell in range(2, 8+1):
                    for m in range(-ell, ell+1):
                        subgroup = 'Y_l{0}_m{1}.dat'.format(ell, m)
                        if subgroup not in f[group]:
                            v.invalid('Subgroup "{0}/{1}" not found in "{2}".'.format(group, subgroup, waveform_path))
                            continue
                        if f[group][subgroup].dtype != np.float:
                            v.invalid('Subgroup "{0}/{1}" of "{2}" has wrong dtype: {3}.'.format(group, subgroup, waveform_path,
                                                                                                 f[group][subgroup].dtype))
                            continue
                        if f[group][subgroup].ndim != 2 or f[group][subgroup].shape[1] != 3:
                            v.invalid('Subgroup "{0}/{1}" of "{2}" has wrong shape: {3}.'.format(group, subgroup, waveform_path,
                                                                                                 f[group][subgroup].shape))
                            continue
                        if t is None:
                            t = f[group][subgroup][:, 0]
                        elif not np.array_equal(t, f[group][subgroup][:, 0]):
                            v.invalid('Inconsistent time data found in subgroup "{0}/{1}" of "{2}".'.format(group, subgroup, waveform_path))
    
    return v.valid


def simulation(path='.', short_circuit=False):
    """Check that a directory contains a valid BBH simulation

    The top-level directory must contain:
    * A file called 'common-metadata.txt' that can be parsed by this module's 'Metadata.from_file'
      method and contains a key 'alternative-names' which includes an entry that matches
      'SXS:BBH:[0-9]*'
    * A directory named 'EvID' containing
      * ID_Files.tgz
      * ID_Params.perl
    * At least one directory matching 'Lev[-0-9]*', each of which must contain
      * A file called 'metadata.txt' that can be parsed by this module's 'Metadata.from_file' method with
        * object1 = 'bh'
        * object2 = 'bh'
        * initial_mass1 > 0
        * initial_mass2 > 0
        * initial_dimensionless_spin1 or initial_spin1 (as a 3-vector)
        * initial_dimensionless_spin2 or initial_spin2 (as a 3-vector)
        * initial_orbital_frequency > 0
      * Horizons.h5 (see `horizons` function for requirements)
      * The standard waveform files (see `waveform` function for requirements):
        * rh_FiniteRadii_CodeUnits.h5 and rPsi4_FiniteRadii_CodeUnits.h5
        * rhOverM_Asymptotic_GeometricUnits.h5 and rMPsi4_Asymptotic_GeometricUnits.h5
        * rhOverM_Asymptotic_GeometricUnits_CoM.h5 and rMPsi4_Asymptotic_GeometricUnits_CoM.h5

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
        valid BBH simulation.

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

    # Check that EvID exists and contains its two files
    evid_path = os.path.join(path, 'EvID')
    if not os.path.isdir(evid_path):
        v.invalid('No "EvID" subdirectory found in "{0}".'.format(path))
    else:
        if not os.path.isfile(os.path.join(evid_path, 'ID_Files.tgz')):
            v.invalid('No "ID_Files.tgz" file found in "{0}".'.format(evid_path))
        if not os.path.isfile(os.path.join(evid_path, 'ID_Params.perl')):
            v.invalid('No "ID_Params.perl" file found in "{0}".'.format(evid_path))

    lev_re = re.compile('Lev[-0-9]+')
    levs = [d for d in os.listdir(path) if lev_re.match(d)]
    if not levs:  # Check that Lev subdirectories exist
        v.invalid('No "Lev" directories found in "{0}".'.format(path))
    else:
        # Loop through, checking each Lev subdirectory
        for lev in levs:
            lev_path = os.path.join(path, lev)

            # Check that metadata.txt exists, can be parsed, and contains the essential data
            metadata_path = os.path.join(lev_path, 'metadata.txt')
            if not os.path.isfile(metadata_path):
                v.invalid('No "metadata.txt" file found in "{0}".'.format(lev_path))
            else:
                try:
                    metadata = Metadata.from_txt_file(metadata_path)
                except:
                    v.invalid('Could not parse "{0}".'.format(metadata_path))
                    metadata = {}
                if metadata.get('object1', '')!='bh' or metadata.get('object2', '')!='bh':
                    v.invalid('Binary object types not "bh" in "{0}".'.format(metadata_path))
                if metadata.get('initial-mass1', 0.0) <= 0.0 or metadata.get('initial-mass2', 0.0) <= 0.0:
                    v.invalid('Masses are not greater than 0 in "{0}".'.format(metadata_path))
                s1 = np.array(metadata.get('initial_dimensionless_spin1', metadata.get('initial_spin1', [])), dtype=float)
                s2 = np.array(metadata.get('initial_dimensionless_spin2', metadata.get('initial_spin2', [])), dtype=float)
                if s1.shape != (3,) or s2.shape != (3,):
                    v.invalid('Initial spins are not found as 3-vectors in "{0}".'.format(metadata_path))
                if metadata.get('initial_orbital_frequency', 0.0) <= 0.0:
                    v.invalid('Initial orbital frequency not greater than 0 in "{0}".'.format(metadata_path))

            # Check that Horizons.h5 exists and contains the expected data
            horizons(os.path.join(lev_path, 'Horizons.h5'), validity=v)

            # Check each pair of waveforms
            waveforms(lev_path, 'rh_FiniteRadii_CodeUnits.h5', 'rPsi4_FiniteRadii_CodeUnits.h5', validity=v, extrapolated=False)
            waveforms(lev_path, 'rhOverM_Asymptotic_GeometricUnits.h5', 'rMPsi4_Asymptotic_GeometricUnits.h5', validity=v, extrapolated=True)
            waveforms(lev_path, 'rhOverM_Asymptotic_GeometricUnits_CoM.h5', 'rMPsi4_Asymptotic_GeometricUnits_CoM.h5', validity=v, extrapolated=True)

    return v.valid


_main = simulation  # This will be the default executable when run via command line
