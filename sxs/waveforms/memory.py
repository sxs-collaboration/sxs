"""Compute contributions to memory in asymptotic waveforms

This code is based on the paper "Adding Gravitational Memory to Waveform
Catalogs using BMS Balance Laws" by Mitman et al.  The main result of that
paper is encapsulated in the `add_memory` function.  Basic usage looks like
this:

```python
h = sxs.load("SXS:BBH:0123/Lev/rhOverM", extrapolation_order=3)
h_with_memory = sxs.waveforms.memory.add_memory(h, start_time=1000.0)
```

"""

import numpy as np
import spherical
from . import WaveformModes


class _ModesTimeSeries(WaveformModes):

    @property
    def s(self):
        return self.spin_weight

    @property
    def u(self):
        return self.time

    from spherical.modes.algebra import (
        conjugate, bar, _real_func, real, _imag_func, imag, norm,
        add, subtract, multiply, divide
    )

    conj = conjugate

    from spherical.modes.utilities import (
        truncate_ell, _check_broadcasting
    )

    from spherical.modes.ufuncs import __array_ufunc__


def MTS(*args, **kwargs):
    kwargs.setdefault('multiplication_truncator', max)
    return _ModesTimeSeries(*args, **kwargs)


def 𝔇(h_mts):
    """Differential operator 𝔇 acting on spin-weight s=0 function

    As defined in Eq. (7b) of 'Adding Gravitational Memory to Waveform Catalogs
    using BMS Balance Laws'

    """
    h = h_mts.copy()
    s = h.ndarray

    for ell in range(h.ell_min, h.ell_max + 1):
        if ell < 2:
            𝔇_value = 0
        else:
            𝔇_value = 0.125 * (ell + 2) * (ell + 1) * (ell) * (ell - 1)
        s[..., h.index(ell, -ell) : h.index(ell, ell) + 1] *= 𝔇_value

    return h


def 𝔇inverse(h_mts):
    """Inverse of differential operator 𝔇 acting on spin-weight s=0 function

    As defined in Eq. (7b) of 'Adding Gravitational Memory to Waveform Catalogs
    using BMS Balance Laws'

    """
    h = h_mts.copy()
    s = h.ndarray

    for ell in range(h.ell_min, h.ell_max + 1):
        if ell < 2:
            𝔇inverse_value = 0
        else:
            𝔇inverse_value = 1.0 / (0.125 * (ell + 2) * (ell + 1) * (ell) * (ell - 1))
        s[..., h.index(ell, -ell) : h.index(ell, ell) + 1] *= 𝔇inverse_value

    return h


def 𝔇inverseLaplacianinverse(h_mts):
    """Inverse of differential operator D² 𝔇 acting on spin-weight s=0 function

    As defined in Eqs. (7) of 'Adding Gravitational Memory to Waveform Catalogs
    using BMS Balance Laws'

    """
    h = h_mts.copy()
    s = h.ndarray

    for ell in range(h.ell_min, h.ell_max + 1):
        if ell < 2:
            𝔇inverse_value = 0
        else:
            𝔇inverse_value = -1.0 / (0.125 * (ell + 2) * (ell + 1)**2 * (ell)**2 * (ell - 1))
        s[..., h.index(ell, -ell) : h.index(ell, ell) + 1] *= 𝔇inverse_value

    return h


def mass_aspect(Psi2, h):
    h = MTS(h)
    Psi2 = MTS(Psi2)
    return - (Psi2 + 0.25 * h.dot * h.bar).re


def J_m(h, Psi2):
    """Bondi mass aspect contribution to electric part of strain

    Calculated according to Eq. (17a) of 'Adding Gravitational Memory to Waveform
    Catalogs using BMS Balance Laws'

    Parameters
    ----------
    h : WaveformModes
        WaveformModes object corresponding to the strain
    Psi2 : WaveformModes
        WaveformModes object corresponding to Psi2

    Returns
    -------
    J_m : WaveformModes
        Bondi mass aspect contribution to the strain

    """
    h = MTS(h)
    Psi2 = MTS(Psi2)
    m = mass_aspect(Psi2, h)
    J_m = 0.5 * 𝔇inverse(m).ethbar.ethbar

    return J_m


def J_E(h, start_time=None):
    """Energy flux contribution to electric part of strain

    Calculated according to Eq. (17b) of 'Adding Gravitational Memory to Waveform
    Catalogs using BMS Balance Laws'

    Parameters
    ----------
    h : WaveformModes
        WaveformModes object corresponding to the strain
    start_time : float, optional
        The time at which the energy flux integral should begin.  The default is
        `h.t[0]`.

    Returns
    -------
    J_E : WaveformModes
        Energy flux contribution to the strain

    """

    hdot = MTS(h).dot

    J_ℰ = 0.5 * 𝔇inverse(0.25 * (hdot * hdot.bar).int).ethbar.ethbar

    if start_time is not None:
        start_time_index = np.argmin(abs(h.t - start_time))
        J_ℰ -= J_ℰ[start_time_index, :]

    return J_ℰ


def J_Nhat(h, Psi2):
    """Angular momentum aspect contribution to magnetic part of strain

    Calculated according to Eq. (17c) of 'Adding Gravitational Memory to Waveform
    Catalogs using BMS Balance Laws'

    Parameters
    ----------
    h : WaveformModes
        WaveformModes object corresponding to the strain
    Psi2 : WaveformModes
        WaveformModes object corresponding to Psi2

    Returns
    -------
    J_Nhat : WaveformModes
        Bondi angular momentum aspect contribution to the strain

    """

    h = MTS(h)
    Psi2 = MTS(Psi2)

    # # Note that the contributions from the last two terms drop out as soon as we
    # # take the imaginary part below.
    # m = mass_aspect(Psi2, h)
    # N̂ = 2 * Psi1 - 0.25 * (h.bar * h.eth) - h.t * m.eth - 0.125 * (h * h.bar).eth

    # Ψ̇₁ = - ½ (ðΨ₂ - ½ h̄ ðḣ)
    # N̂̇ = 2 Ψ̇₁ - ¼ h̄ ðḣ - ¼ h̄̇ ðh
    #   = -ðΨ₂ + ½ h̄ ðḣ - ¼ h̄ ðḣ - ¼ h̄̇ ðh
    #   = -ðΨ₂ + ¼ h̄ ðḣ - ¼ h̄̇ ðh

    N̂dot = -Psi2.eth + 0.25 * (h.bar * h.dot.eth - h.bar.dot * h.eth)
    J_N̂ = 0.5j * 𝔇inverseLaplacianinverse(N̂dot.ethbar.im).ethbar.ethbar

    return J_N̂


def J_J(h):
    """Angular momentum flux contribution to magnetic part of strain

    Calculated according to Eq. (17d) of 'Adding Gravitational Memory to Waveform
    Catalogs using BMS Balance Laws'

    Parameters
    ----------
    h : WaveformModes
        WaveformModes object corresponding to the strain

    Returns
    -------
    J_J : WaveformModes
        Angular momentum flux contribution to the strain

    """

    h = MTS(h)
    hdot = h.dot
    J_𝒥 = 0.5j * 𝔇inverseLaplacianinverse(
        0.125 * (3 * h * hdot.bar.ethbar - 3 * hdot * h.bar.ethbar + hdot.bar * h.ethbar - h.bar * hdot.ethbar).eth.im
    ).ethbar.ethbar

    return J_𝒥


def add_memory(h, start_time=None):
    """Add electric component of null memory to strain

    This adds the contribution from the energy flux to the strain.

    Parameters
    ----------
    h : WaveformModes
        WaveformModes object corresponding to the strain
    start_time : float, optional
        Time at which the energy flux integral should begin.  The default is
        `h.t[0]`.

    Returns
    -------
    h_with_memory : WaveformModes
        WaveformModes object corresponding to the strain with electric memory

    """
    h_with_memory = MTS(h) + J_E(h, start_time=start_time)
    return WaveformModes(h_with_memory)
