
def bbh(path):
    """Check that a directory contains a valid BBH simulation

    The top-level directory must contain:
    * A file called 'common-metadata.txt' that can be parsed by this module's 'Metadata.from_file'
      method and contains a key 'alternative-names' which includes an entry that matches
      'SXS:BBH:[0-9]*'
    * At least one directory matching 'Lev[0-9]', each of which must contain
      * A file called 'metadata.txt' that can be parsed by this module's 'Metadata.from_file' method
      * The standard data files:
        * Horizons.h5
        * rh_FiniteRadii_CodeUnits.h5
        * rPsi4_FiniteRadii_CodeUnits.h5
        * rhOverM_Asymptotic_GeometricUnits.h5
        * rMPsi4_Asymptotic_GeometricUnits.h5
        * rhOverM_Asymptotic_GeometricUnits_CoM.h5
        * rMPsi4_Asymptotic_GeometricUnits_CoM.h5

    """
    raise NotImplementedError()
