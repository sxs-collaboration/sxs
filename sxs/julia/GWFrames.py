import quaternionic
from .. import WaveformModes
from . import PostNewtonian

def PNWaveform(
    Approximant, delta, chi1_i, chi2_i, Omega_orb_i, *,
    Omega_orb_0=None, R_frame_i=None, MinStepsPerOrbit=32,
    PNWaveformModeOrder=4.0, PNOrbitalEvolutionOrder=4.0,
    inertial=False, dt=0.0, quiet=True,
    ell_min=2, ell_max=8, Lambda1=0, Lambda2=0,
    **kwargs
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
    
    This is a wrapper around the Julia function
    `PostNewtonian.GWFrames.PNWaveform`.  See [the Julia
    documentation](https://moble.github.io/PostNewtonian.jl/dev/interface/gwframes/#PostNewtonian.GWFrames.PNWaveform)
    for details on the optional keyword arguments.
    
    """
    Omega_orb_0 = Omega_orb_0 or Omega_orb_i
    R_frame_i = R_frame_i or [1.0]

    # Integrate the orbital dynamics
    w1 = PostNewtonian.GWFrames.PNWaveform(
        Approximant, delta, chi1_i, chi2_i, Omega_orb_i,
        Omega_orb_0=Omega_orb_0, R_frame_i=R_frame_i,
        MinStepsPerOrbit=MinStepsPerOrbit,
        PNWaveformModeOrder=PNWaveformModeOrder,
        PNOrbitalEvolutionOrder=PNOrbitalEvolutionOrder,
        inertial=inertial, dt=dt, quiet=quiet,
        ell_min=ell_min, ell_max=ell_max, Lambda1=Lambda1, Lambda2=Lambda2,
        **kwargs
    )

    w = WaveformModes(
        w1.data,
        time=w1.t,
        modes_axis=1,
        ell_min=2,
        ell_max=8,
        M1=w1.M1.to_numpy(),
        M2=w1.M2.to_numpy(),
        chi1=w1.chi1.to_numpy(),
        chi2=w1.chi2.to_numpy(),
        frame=quaternionic.array(w1.frame.to_numpy()),
        v=w1.v.to_numpy(),
        orbital_phase=w1.Phi.to_numpy(),
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
