import re
import numpy as np
import h5py
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
import spherical as sf
from ... import waveforms

def load(file_name, data_type="h"):
    """Load a waveform from an HDF5 file in the LVCNR format.

    This function reads data from HDF5 files in the LVCNR format.
    It extracts amplitude and phase information for spherical harmonic modes
    and reconstructs the complex-valued data. This function can be extended
    to further data types.

    Parameters
    ----------
    file_name : str
        Path to the HDF5 file to be loaded.
        
    data_type : str, optional
        The type of data to be loaded. Default is 'h' for strain.

    Returns
    -------
    waveforms.WaveformModes
        The loaded waveform data in WaveformModes format with spin weight
        -2 and data type 'h' (strain). The waveform includes mode
        decomposition information with ell_min and ell_max determined from
        the file contents. The data type can be overriden if necessary, to
        accomodate other types of data stored in LVCNR format.
    """

    phase_re = re.compile("phase_l(?P<ell>.*)_m(?P<m>.*)")

    with h5py.File(file_name, "r") as f:

        # Find the NRTimes array. RIT and MAYA have different 
        # capitalizations.
        nrtimes_matches = [
            key for key in f if key.lower() == "nrtimes"
        ]
        if not (len(nrtimes_matches) == 1 and
                isinstance(f[nrtimes_matches[0]], h5py.Dataset)):
            raise KeyError(
                f"File {file_name} must contain exactly 1 Dataset "
                "named NRTimes (case insensitive)"
            )

        t = f[nrtimes_matches[0]]

        # Extract l and m values from the dataset keys
        ell_m = np.array(
            [
                [int(match["ell"]), int(match["m"])]
                for key in f
                for match in [phase_re.match(key)]
                if match
            ]
        )
        ell_min = np.min(ell_m[:, 0])
        ell_max = np.max(ell_m[:, 0])
        data = np.empty(
            (t.size, sf.LM_total_size(ell_min, ell_max)),
            dtype=complex
        )

        # Compose the data array using amplitude and phase datasets
        for ell in range(ell_min, ell_max + 1):
            for m in range(-ell, ell + 1):
                amp = Spline(
                    f[f"amp_l{ell}_m{m}/X"][:],
                    f[f"amp_l{ell}_m{m}/Y"][:],
                    k=int(f[f"amp_l{ell}_m{m}/deg"][()])
                )(t)
                phase = Spline(
                    f[f"phase_l{ell}_m{m}/X"][:],
                    f[f"phase_l{ell}_m{m}/Y"][:],
                    k=int(f[f"phase_l{ell}_m{m}/deg"][()])
                )(t)
                idx = sf.LM_index(ell, m, ell_min)
                data[:, idx] = amp * np.exp(1j * phase)

        if "auxiliary-info" in f and \
           "history.txt" in f["auxiliary-info"]:
            history_txt = f["auxiliary-info/history.txt"][()].decode()
            history = (
                "### " + history_txt.replace("\n", "\n### ")
            ).split("\n")
        else:
            history = [""]

        kwargs = dict(
            time=t,
            time_axis=0,
            modes_axis=1,
            spin_weight=-2,
            data_type=data_type,
            frame_type="inertial",
            history=history,
            version_hist=[],
            r_is_scaled_out=True,
            m_is_scaled_out=True,
            ell_min=ell_min,
            ell_max=ell_max,
        )

        w = waveforms.WaveformModes(data, **kwargs)

    return w
