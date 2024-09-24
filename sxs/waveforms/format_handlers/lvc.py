import numpy as np
from scipy.optimize import root_scalar
from scipy.interpolate import CubicSpline, UnivariateSpline
import quaternionic

from ... import jit, TimeSeries
from ..waveform_modes import WaveformModesDict

from numpy import exp
from numpy import pi as π
from quaternionic import z


@jit
def find_i1(f, i2):
    """Convenience function to find previous non-increasing index.
    
    Given index `i2` find the nearest preceding `i1` such that
    `f[i1-1] > f[i1]`; if there is no such `i1`, return 0.
    """
    for i in range(i2, 0, -1):
        if f[i] <= f[i-1]:
            return i
    return 0


def to_lvc_conventions(
        h, horizons, *,
        t_ref=None, f_ref=None,
        dt=None,
        f_low=None,
        ell_max=None,
        phi_ref=None, inclination=None,
        ell_max_epoch=None,
        **kwargs
):
    r"""Return SXS waveform in LVC conventions.

    This is the function underpinning `sxs.load_lvc`.  Essentially,
    that function is a fairly thin wrapper around this one; it just
    loads the `h` and `horizons` data and passes it to this function.
    For explanation of what this function does, the inputs, and
    outputs, see that function's docstring.

    """
    if t_ref is None and f_ref is None:
        raise ValueError("One of `t_ref` or `f_ref` must be specified")
    if t_ref is not None and f_ref is not None:
        raise ValueError("Only one of `t_ref` or `f_ref` may be specified")
    
    # Transform to corotating frame, if necessary
    if h.frame_type != "corotating":
        h = h.to_corotating_frame()
    
    # Find the spin data
    chi1 = horizons.A.chi_inertial
    chi2 = horizons.B.chi_inertial
    t_chi = chi1.t.copy()
    t_remnant = horizons.C.time.copy()
    mass_remnant = horizons.C.christodoulou_mass.ndarray

    # If requested, restrict `ell_max`
    if ell_max is not None:
        h = h[:, h.index(2,-2):h.index(ell_max, ell_max)+1]
        ell_max_epoch = min(ell_max_epoch or ell_max, ell_max)
    else:
        ell_max_epoch = ell_max_epoch or h.ell_max

    # Set `t` to equal 0 at the moment of maximum ell<=ell_max_epoch norm
    epoch_time = h[:, h.index(2,-2):h.index(ell_max_epoch, ell_max_epoch)+1].max_norm_time(interpolate=True)
    h.t -= epoch_time
    assert h.t[0] < 0.0 and h.t[-1] > 0.0, "Epoch time is not within waveform data"
    t_chi -= epoch_time
    t_remnant -= epoch_time

    # If `t_ref` is not specified, find `t_ref` from `f_ref`.
    # Otherwise, use `t_ref` to find `f_ref`.
    omega = h.frame.to_angular_velocity(h.t)  # Note that this is a vector-valued function of time
    f = np.linalg.norm(omega, axis=1) / (2 * np.pi)
    f = UnivariateSpline(h.t, f, s=sum((1e-4*f)**2))(h.t)  # Smooth the data a little
    i2 = h.max_norm_index()
    i1 = find_i1(f, i2)
    if f_low is not None and (f_low < f[i1] or f_low > f[i2]):
        raise ValueError(
            f"Requested frequency {f_low=} is outside this waveform's "
            f"frequency band ({f[i1]}, {f[i2]})"
        )
    if f_ref is not None and (f_ref < f[i1] or f_ref > f[i2]):
        raise ValueError(
            f"Requested frequency {f_ref=} is outside this waveform's "
            f"frequency band ({f[i1]}, {f[i2]})"
        )
    f_spline = CubicSpline(h.t[i1:i2], f[i1:i2])
    if f_ref is not None:
        t_ref = root_scalar(lambda t: f_spline(t)-f_ref, bracket=(h.t[i1], h.t[i2])).root
    else:
        f_ref = f_spline(t_ref)
    if f_low is not None:
        t_low = root_scalar(lambda t: f_spline(t)-f_low, bracket=(h.t[i1], h.t[i2])).root
    else:
        t_low = h.t[i1]
        f_low = f[i1]
        h = h[i1:, :]

    # Interpolate to find frame at `t_ref` and account for ϕ at that time
    h_ref = h.interpolate([t_ref])
    frame_ref = h_ref.frame[0]

    # Use omega at that time as a guess for the rough direction of the radiation axis,
    # and find the dominant principal axis of the <LL> matrix at that time
    omega_ref = omega[h.index_closest_to(t_ref)]
    rough_direction = omega_ref / np.linalg.norm(omega_ref)
    LL_axis = quaternionic.array.from_vector_part(
        h_ref.dominant_eigenvector_LL(rough_direction=rough_direction)[0]
    )

    # Now find the rotation *in the corotating frame* that aligns the
    # z axis with the radiation axis.
    R_ref = np.sqrt(-LL_axis * quaternionic.z)  # R_ref * z * R_ref.conj() ≈ LL_axis

    # If we rotate this waveform by `R_ref`, the LL_axis should be aligned
    # with the z-axis to machine precision.
    h_ref = h_ref.rotate(R_ref)

    # Now we augment the frame transformation to account for the
    # reference phase, as described in the docstring:
    #   * Im{h_{2,2}^sym} = 0  # There will be 4 values of phi that satisfy this
    #   * Re{h_{2,2}^sym} < 0  # This will reduce the number of solutions to 2
    #   * Im{h_{2,1}^sym} < 0  # An odd-m mode is required to break the degeneracy
    α = -np.angle(
        h_ref[0, h_ref.index(2,2)]
        + h_ref[0, h_ref.index(2,-2)].conjugate()
    ) / 2
    h_ref2 = h_ref.rotate(exp(z * α/2))
    β = (
        h_ref2[0, h_ref2.index(2,2)]
        + h_ref2[0, h_ref2.index(2,-2)].conjugate()
    )
    if β.real > 0:
        α += π / 2
    h_ref3 = h_ref.rotate(exp(z * α/2))
    γ = (
        h_ref3[0, h_ref3.index(2,1)]
        + h_ref3[0, h_ref3.index(2,-1)].conjugate()
    )
    if γ.imag > 0:
        α += π
    R_ref = R_ref * exp(z * α/2)
    
    # We want to re-interpret the corotating waveform so that the
    # waveform *in this frame* remains the same, but appears to
    # inertial observers to be moving with an angular velocity (as a
    # function of time) that is just rotated by a constant.  This
    # means that we need to multiply `h.frame` *on the left* by a
    # constant.  And we want the resulting frame as a function of time
    # to be `R_ref.conj()` at time `t_ref`.
    R_adjustment = R_ref.conj() * frame_ref.conj()
    frame_quat = R_adjustment * h.frame
    h._metadata["frame"] = frame_quat

    # Transform to inertial-frame waveform, spins, and omega to `t_ref` frame
    h = h.to_inertial_frame()
    omega = R_adjustment.rotate(omega)
    chi1 = TimeSeries(R_adjustment.rotate(chi1), t_chi)
    chi2 = TimeSeries(R_adjustment.rotate(chi2), t_chi)
    chi_remnant = TimeSeries(R_adjustment.rotate(horizons.C.chi_inertial), t_remnant)

    # If requested, interpolate to time step `dt`, ensuring that `t=0` is preserved
    if dt is not None:
        t = np.concatenate((
            np.arange(0.0, h.t[0], -dt)[::-1],
            np.arange(dt, h.t[-1], dt)
        ))
        h = h.interpolate(t)

    # If either `phi_ref` or `inclination` is not None, set default for the other
    if phi_ref is not None:
        inclination = inclination or 0.0
    if inclination is not None:
        phi_ref = phi_ref or 0.0

    chi1_ref = chi1.interpolate([t_ref+epoch_time]).ndarray[0]
    chi2_ref = chi2.interpolate([t_ref+epoch_time]).ndarray[0]

    dynamics_dict = {
        "t_ref": t_ref,
        "f_ref": f_ref,
        "t_low": t_low,
        "f_low": f_low,
        "chi1_ref": chi1_ref,
        "chi2_ref": chi2_ref,
        "frame_quat": frame_quat,
        "frame_omega": omega,
        "times_spins": t_chi,
        "chi1": chi1,
        "chi2": chi2,
        "times_remnant": t_remnant,
        "chi_remnant": chi_remnant,
        "mass_remnant": mass_remnant,
    }

    # If `phi_ref` and `inclination` are not None, return polarizations
    if phi_ref is not None:
        hp, hc = h.evaluate(inclination, π/2 - phi_ref).ndarray.view((float, 2)).T
        return h.t, hp, hc, dynamics_dict
    else:
        # Could do `dict(WaveformModesDict(h))` to convert to a plain dict
        hlm_dict = WaveformModesDict(h)
        return h.t, hlm_dict, dynamics_dict
