import numpy as np
import quaternionic
from .. import WaveformModes

import juliacall


PostNewtonian = juliacall.newmodule("PN")
PostNewtonian.seval("using PostNewtonian")

from . import GWFrames

def pkg_update():
    """Update all installed Julia packages (but not julia itself)."""
    PostNewtonian.seval("import Pkg")
    return PostNewtonian.seval("Pkg.update()")


def PNWaveform(
        M1, M2, chi1, chi2, Omega_i, *,
        inertial=True,
        ell_min=2, ell_max=8, waveform_pn_order=None,
        **orbital_evolution_kwargs
):
    """Generate a PN waveform.

    The return value is an `sxs.WaveformModes` object with the
    following additional fields:

        - `M1` (array): The primary mass as a function of time.
        - `M2` (array): The secondary mass as a function of time.
        - `chi1` (array): The primary spin as a function of time.
        - `chi2` (array): The secondary spin as a function of time.
        - `frame` (array): The quaternionic frame as a function of
          time.
        - `v` (array): The orbital velocity as a function of time.
        - `orbital_phase` (array): The orbital phase as a function of
          time.

    This is a wrapper around the Julia functions
    `PostNewtonian.orbital_evolution` and
    `PostNewtonian.inertial_waveform`.  See [the Julia
    documentation](https://moble.github.io/PostNewtonian.jl/dev/internals/dynamics/#PostNewtonian.orbital_evolution)
    for details on the optional keyword arguments.

    Note that the full Julia interface is also accessible from this
    Python module.  See [the Julia
    docs](https://moble.github.io/PostNewtonian.jl/dev/interface/python/)
    for details.  Also, the GWFrames submodule of this module provides
    another interface.

    """
    # Integrate the orbital dynamics
    inspiral = PostNewtonian.orbital_evolution(
        M1, M2, chi1, chi2, Omega_i, **orbital_evolution_kwargs
    )
    values = PostNewtonian.stack(inspiral.u)

    # Compute the waveform in the inertial frame
    waveform_pn_order = waveform_pn_order or PostNewtonian.typemax(PostNewtonian.Int)
    coorbital_frame = quaternionic.array(values[8:12, :].to_numpy().T)
    if inertial:
        frame = np.array([quaternionic.one])
        frame_type = "inertial"
        h = PostNewtonian.inertial_waveform(
            inspiral, ell_min=ell_min, ell_max=ell_max, PNOrder=waveform_pn_order
        ).to_numpy().T
    else:
        frame = coorbital_frame
        frame_type = "coorbital"
        h = PostNewtonian.coorbital_waveform(
            inspiral, ell_min=ell_min, ell_max=ell_max, PNOrder=waveform_pn_order
        ).to_numpy().T

    w = WaveformModes(
        h,
        time=inspiral.t,
        modes_axis=1,
        ell_min=ell_min,
        ell_max=ell_max,
        spin_weight=-2,
        frame=frame,
        frame_type=frame_type,
        data_type="h",
        M1=values[0, :].to_numpy(),
        M2=values[1, :].to_numpy(),
        chi1=values[2:5, :].to_numpy().T,
        chi2=values[5:8, :].to_numpy().T,
        coorbital_frame=coorbital_frame,
        v=values[12, :].to_numpy(),
        orbital_phase=values[13, :].to_numpy(),
    )

    # Allow the extra fields to be accessed naturally
    w.M1 = w._metadata["M1"]
    w.M2 = w._metadata["M2"]
    w.chi1 = w._metadata["chi1"]
    w.chi2 = w._metadata["chi2"]
    w.coorbital_frame = w._metadata["coorbital_frame"]
    w.v = w._metadata["v"]
    w.orbital_phase = w._metadata["orbital_phase"]

    return w
