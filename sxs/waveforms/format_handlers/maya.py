from .. import WaveformModes
from enum import Enum


class RadiusOption(Enum):
    INNERMOST = 1
    OUTERMOST = 2

def load(file_name, radius=RadiusOption.OUTERMOST):
    """Load waveforms from the new MAYA file format.

    Parameters
    ==========
    file : str or Path
        Relative or absolute path to the input HDF5 file.  If this
        string contains but does not *end* with `'.h5'`, the remainder
        of the string is taken to be the group within the HDF5 file in
        which the data is stored.

    radius : str or Enum
        The extraction radius to use. Default is the outermost radii RadiusOption.OUTERMOST.

    Notes
    =====
    Waveforms from MAYA are distributed as HDF5 files, containing finite-radius
    data, for Psi4. Within each of those h5 files, we have finite-radii
    waveforms corresponding to a series of extraction radii.
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

    radius_re = re.compile(r"^radius=(?P<radius>.+)$")

    with h5py.File(path) as f:
        radii_group = f["radiative/psi4"]

        radii = [m.group("radius")
                 for key in radii_group.keys() if (m := radius_re.match(key))]

        if isinstance(radius, RadiusOption):
            if radius == RadiusOption.OUTERMOST:
                input_radius = max(radii, key=float)
            elif radius == RadiusOption.INNERMOST:
                input_radius = min(radii, key=float)
        elif isinstance(radius, str):
            if radius in radii:
                input_radius = radius
            else:
                raise ValueError(f"Invalid radius string: {radius}. Must be one of {radii}.")
        else:
            raise ValueError(
                f"Radius can only be specified from a RadiusOption or str instance; this object is of type `{type(radius)}`."
            )

        radius_key = f"radius={input_radius}"

        time = radii_group[f"{radius_key}/time"][()]
        modes_group = radii_group[f"{radius_key}/modes"]

        ells_re = re.compile(r"^l=(?P<ell_number>\d+)$")

        ells = [int(m.group("ell_number"))
                for key in modes_group.keys() if (m := ells_re.match(key))]
        ell_min = min(ells)
        ell_max = max(ells)

        emms_re = re.compile(r"^m=(?P<emm_number>-?\d+)$")

        data = np.zeros((len(time), spherical.Ysize(int(ell_min), int(ell_max))), dtype=complex)

        for ell in ells:
            ell_group = modes_group[f"l={ell}"]
            emms = [int(m.group("emm_number"))
                    for key in ell_group.keys() if (m := emms_re.match(key))]
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
