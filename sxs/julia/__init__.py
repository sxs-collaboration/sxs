import quaternionic
from .. import WaveformModes

import juliacall


PostNewtonian = juliacall.newmodule("PN")
PostNewtonian.seval("using PostNewtonian")


def pkg_update():
    """Update all installed Julia packages (but not julia itself)."""
    PostNewtonian.seval("import Pkg")
    return PostNewtonian.seval("Pkg.update()")


def PNWaveform(M1, M2, chi1, chi2, Omega_i, **kwargs):
    """Generate a PN waveform.
    
    This is a wrapper around the Julia functions
    `PostNewtonian.orbital_evolution` and
    `PostNewtonian.inertial_waveform`.  See [the Julia
    documentation](https://moble.github.io/PostNewtonian.jl/dev/internals/dynamics/#PostNewtonian.orbital_evolution)
    for details on the optional keyword arguments.
    
    """
    # Integrate the orbital dynamics
    inspiral = PostNewtonian.orbital_evolution(M1, M2, chi1, chi2, Omega_i, **kwargs)

    # Compute the waveform in the inertial frame
    h = PostNewtonian.inertial_waveform(inspiral).to_numpy().T

    w = WaveformModes(
        h,
        time=inspiral.t,
        modes_axis=1,
        ell_min=2,
        ell_max=8,
        M1=inspiral[0, :].to_numpy(),
        M2=inspiral[1, :].to_numpy(),
        chi1=inspiral[2:5, :].to_numpy().T,
        chi2=inspiral[5:8, :].to_numpy().T,
        frame=quaternionic.array(inspiral[8:12, :].to_numpy().T),
        v=inspiral[12, :].to_numpy(),
        orbital_phase=inspiral[13, :].to_numpy(),
    )

    # Allow the extra fields to be accessed naturally
    w.M1 = w._metadata["M1"]
    w.M2 = w._metadata["M2"]
    w.chi1 = w._metadata["chi1"]
    w.chi2 = w._metadata["chi2"]
    #w.frame = w._metadata["frame"]  ## Already done
    w.v = w._metadata["v"]
    w.orbital_phase = w._metadata["orbital_phase"]

    return w
