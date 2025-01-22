"""Functions to load and save waveforms in SpECTRE CCE v1 format"""

from .. import WaveformModes


def save(*args, **kwargs):
    raise NotImplementedError("Saving waveforms in SpECTRE CCE format (v1) is not supported")


def load(file_name, **kwargs):
    """Load a waveform in SpECTRE CCE format (version 1)

    Parameters
    ----------
    file_name : str or Path
        Relative or absolute path to the input HDF5 file.  If this
        string contains but does not *end* with `'.h5'`, the remainder
        of the string is taken to be the group within the HDF5 file in
        which the data is stored.  Also note that a JSON file is
        expected in the same location, with `.h5` replaced by `.json`
        (and the corresponding data must be stored under the `group`
        key if relevant).

    Required keyword argument
    -------------------------
    group : str
        The group within the HDF5 file in which the data is stored.

    """
    import re
    from pathlib import Path
    import numpy as np
    import h5py
    import spherical

    # Make sure the file exists
    path = Path(file_name).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Could not find {path}")

    # Get the group name
    group = kwargs.pop("group", None)
    if group is None:
        raise ValueError("The 'group' keyword argument is required")

    # Determine the data type and spin weight from the group name
    if "news" in group.lower():
        data_type = "news"
        spin_weight = -2
    elif "psi0" in group.lower():
        data_type = "psi0"
        spin_weight = 2
    elif "psi1" in group.lower():
        data_type = "psi1"
        spin_weight = 1
    elif "psi2" in group.lower():
        data_type = "psi2"
        spin_weight = 0
    elif "psi3" in group.lower():
        data_type = "psi3"
        spin_weight = -1
    elif "psi4" in group.lower():
        data_type = "psi4"
        spin_weight = -2
    elif "strain" in group.lower():
        data_type = "h"
        spin_weight = -2
    else:
        raise ValueError(f"Unrecognized data type in group name '{group}'")

    m_is_scaled_out = True
    r_is_scaled_out = True

    with h5py.File(str(path), "r") as h5:
        if group is not None:
            h5 = h5[group]
        legend = list(h5.attrs["Legend"])
        if len(legend) != h5.shape[1]:
            raise ValueError(
                f"Number of columns in the data ({h5.shape[1]}) "
                f"does not match the number of entries in the legend ({len(legend)})"
            )
        regex = re.compile(r"(Real|Imag) Y_([0-9]+),([-0-9]+)")

        time_index = legend.index("time")
        time = h5[:, time_index]
        indices = np.argsort(time)
        time = time[indices]
        lm = np.array(
            [[m[2], m[3]] for l in legend if (m:=regex.match(l))],
            dtype=int
        )
        ell_min = min(lm[:, 0])
        ell_max = max(lm[:, 0])

        data = np.zeros((len(time), spherical.Ysize(ell_min, ell_max)), dtype=complex)
        for i, legend_entry in enumerate(legend):
            if (match:=regex.match(legend_entry)):
                ell, m = int(match[2]), int(match[3])
                if match[1] == "Real":
                    data[:, spherical.Yindex(ell, m, ell_min=ell_min)] += h5[:, i][indices]
                elif match[1] == "Imag":
                    data[:, spherical.Yindex(ell, m, ell_min=ell_min)] += 1j * h5[:, i][indices]
                else:
                    raise ValueError(f"Unrecognized legend entry '{legend_entry}'")

        return WaveformModes(
            data,
            time=time,
            time_axis=0,
            modes_axis=1,
            frame_type="inertial",
            data_type=data_type,
            m_is_scaled_out=m_is_scaled_out,
            r_is_scaled_out=r_is_scaled_out,
            ell_min=ell_min,
            ell_max=ell_max,
            spin_weight=spin_weight,
        )
