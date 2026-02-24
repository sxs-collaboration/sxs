from .. import WaveformModes


def load(file_name, radius="outermost"):
    """Load waveforms from the new MAYA file format.

    Parameters
    ==========
    file_name : str or Path
        Relative or absolute path to the input HDF5 file.  If this
        string contains but does not *end* with `'.h5'`, the remainder
        of the string is taken to be the group within the HDF5 file in
        which the data is stored.

    radius : str
        The extraction radius to use.
        Default is "outermost" which corresponds to largest available extraction
        radius. Another option is "innermost", corresponding to the smallest
        extraction radius. The user can also specify a string correspond to any
        extraction radius present in the HDF5 file. See notes.

    Notes
    =====
    Waveforms from MAYA are distributed as HDF5 files, containing finite-radius
    data for Psi4. Each file includes finite-radii waveforms corresponding to a series of extraction radii, which may include:

    * "50.00"
    * "60.00"
    * "70.00"
    * "80.00"
    * "90.00"
    * "100.00"
    * "120.00"
    * "140.00"

    Note the two 0s after the decimal point. These radii correspond to the keys
    of the group located at radiative/psi4 path in the HDF5 file.

    Returns
    =======
    WaveformModes object.
    """
    import numpy as np
    import h5py
    import spherical
    from pathlib import Path
    import re

    # Make sure the file exists
    path = Path(file_name).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Could not find {path}")

    if not isinstance(radius, str):
        raise ValueError(f"Radius can only be specified from a str instance; this object is of type `{type(radius)}`.")

    radius_re = re.compile(r"^radius=(?P<radius>.+)$")

    with h5py.File(path) as f:
        radii_group = f["radiative/psi4"]

        radii = [m.group("radius") for key in radii_group.keys() if (m := radius_re.match(key))]

        if radius == "outermost":
            input_radius = max(radii, key=float)
        elif radius == "innermost":
            input_radius = min(radii, key=float)
        elif radius in radii:
            input_radius = radius
        else:
            raise ValueError(f"Invalid radius string: {radius}. Must be one of {radii}.")

        radius_key = f"radius={input_radius}"

        time = radii_group[f"{radius_key}/time"][()]
        modes_group = radii_group[f"{radius_key}/modes"]

        ells_re = re.compile(r"^l=(?P<ell_number>[0-9]+)$")

        ells = [int(m.group("ell_number")) for key in modes_group.keys() if (m := ells_re.match(key))]
        ell_min = min(ells)
        ell_max = max(ells)

        emms_re = re.compile(r"^m=(?P<emm_number>-?[0-9]+)$")

        data = np.zeros((len(time), spherical.Ysize(int(ell_min), int(ell_max))), dtype=complex)

        for ell in ells:
            ell_group = modes_group[f"l={ell}"]
            emms = [int(m.group("emm_number")) for key in ell_group.keys() if (m := emms_re.match(key))]
            for emm in emms:
                emm_group = ell_group[f"m={emm}"]
                data[:, spherical.Yindex(ell, emm, ell_min=ell_min)] = (
                    emm_group["real"][()] + 1j * emm_group["imaginary"][()]
                )

        return WaveformModes(
            data,
            time=time,
            time_axis=0,
            modes_axis=1,
            frame_type="inertial",
            data_type="psi4",
            m_is_scaled_out=True,
            r_is_scaled_out=False,
            ell_min=ell_min,
            ell_max=ell_max,
            spin_weight=-2,
        )
