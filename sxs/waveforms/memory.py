import numpy as np
import spherical
from . import WaveformModes

class ModesTimeSeries(WaveformModes):

    @property
    def s(self):
        return self.spin_weight

    @property
    def u(self):
        return self.time

    from spherical.SWSH_modes.algebra import (
        conjugate, bar, _real_func, real, _imag_func, imag, norm,
        add, subtract, multiply, divide
    )

    conj = conjugate
    re = real
    im = imag

    from spherical.SWSH_modes.derivatives import (
        Lsquared, Lz, Lplus, Lminus,
        Rsquared, Rz, Rplus, Rminus,
        eth, ethbar
    )

    @property
    def eth_GHP(self):
        """Raise spin-weight with GHP convention"""
        return self.eth / np.sqrt(2)

    @property
    def ethbar_GHP(self):
        """Lower spin-weight with GHP convention"""
        return self.ethbar / np.sqrt(2)

    from spherical.SWSH_modes.utilities import (
        truncate_ell, _check_broadcasting
    )

    from spherical.SWSH_modes.ufuncs import __array_ufunc__


def MTS(*args, **kwargs):
    kwargs.setdefault('multiplication_truncator', max)
    return ModesTimeSeries(*args, **kwargs)


def 𝔇(h_mts):
    """Differential operator 𝔇 acting on spin-weight s=0 function

    As defined in Eq. (7b) of 'Adding Gravitational Memory to Waveform Catalogs
    using BMS Balance Laws'

    """
    h = h_mts.copy()
    s = h_mts.ndarray

    for ell in range(h.ell_min, h.ell_max + 1):
        if ell < 2:
            𝔇_value = 0
        else:
            𝔇_value = 0.125 * (ell + 2) * (ell + 1) * (ell) * (ell - 1)
            # 𝔇_value = -0.125 * (ell + 2) * (ell + 1)**2 * (ell)**2 * (ell - 1)
            # # 𝔇_value = 0.125 * (-ell * (ell + 1)) * (-ell * (ell + 1) + 2)
            # # 𝔇_value = 0.125 * (-ell * (ell + 1)) ** 2 * (-ell * (ell + 1) + 2)
        s[..., h.index(ell, -ell) : h.index(ell, ell) + 1] *= 𝔇_value

    return h


def 𝔇inverse(h_mts):
    """Inverse of differential operator 𝔇 acting on spin-weight s=0 function

    As defined in Eq. (7b) of 'Adding Gravitational Memory to Waveform Catalogs
    using BMS Balance Laws'

    """
    h = h_mts.copy()
    s = h_mts.ndarray

    for ell in range(h.ell_min, h.ell_max + 1):
        if ell < 2:
            𝔇inverse_value = 0
        else:
            𝔇inverse_value = 1.0 / (0.125 * (ell + 2) * (ell + 1) * (ell) * (ell - 1))
            # 𝔇inverse_value = 1.0 / (0.125 * (ell + 2) * (ell + 1)**2 * (ell)**2 * (ell - 1))
            # # 𝔇inverse_value = 1.0 / (0.125 * (-ell * (ell + 1)) * (-ell * (ell + 1) + 2))
            # # 𝔇inverse_value = 1.0 / (0.125 * (-ell * (ell + 1)) ** 2 * (-ell * (ell + 1) + 2))
        s[..., h.index(ell, -ell) : h.index(ell, ell) + 1] *= 𝔇inverse_value

    return h


def 𝔇inverseLaplacianinverse(h_mts):
    """Inverse of differential operator D² 𝔇 acting on spin-weight s=0 function

    As defined in Eqs. (7) of 'Adding Gravitational Memory to Waveform Catalogs
    using BMS Balance Laws'

    """
    h = h_mts.copy()
    s = h_mts.ndarray

    for ell in range(h.ell_min, h.ell_max + 1):
        if ell < 2:
            𝔇inverse_value = 0
        else:
            # 𝔇inverse_value = 1.0 / (0.125 * (ell + 2) * (ell + 1) * (ell) * (ell - 1))
            𝔇inverse_value = 1.0 / (0.125 * (ell + 2) * (ell + 1)**2 * (ell)**2 * (ell - 1))
            # # 𝔇inverse_value = 1.0 / (0.125 * (-ell * (ell + 1)) * (-ell * (ell + 1) + 2))
            # # 𝔇inverse_value = 1.0 / (0.125 * (-ell * (ell + 1)) ** 2 * (-ell * (ell + 1) + 2))
        s[..., h.index(ell, -ell) : h.index(ell, ell) + 1] *= 𝔇inverse_value

    return h


# def 𝔇inverse(h_mts, use_laplacian=False):
#     """Inverse of differential operator 𝔇

#     As defined in Eq. (7b) of 'Adding Gravitational Memory to Waveform Catalogs
#     using BMS Balance Laws'

#     """
#     h = h_mts.copy()
#     s = h_mts.ndarray

#     for ell in range(h.ell_min, h.ell_max + 1):
#         if ell < 2:
#             𝔇inverse_value = 0
#         else:
#             if not use_laplacian:
#                 𝔇inverse_value = 1.0 / (0.125 * (-ell * (ell + 1)) * (-ell * (ell + 1) + 2))
#             else:
#                 𝔇inverse_value = 1.0 / (0.125 * (-ell * (ell + 1)) ** 2 * (-ell * (ell + 1) + 2))
#         s[..., h.index(ell, -ell) : h.index(ell, ell) + 1] *= 𝔇inverse_value

#     return h


def mass_aspect(Psi2, h):
    h = MTS(h)
    Psi2 = MTS(Psi2)
    return - (Psi2 + 0.25 * h.dot * h.bar).re


def J_m(h, Psi2):
    """Bondi mass aspect contribution to electric part of strain

    Calculated according to Eq. (16a) of 'Adding Gravitational Memory to Waveform
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

    Calculated according to Eq. (16b) of 'Adding Gravitational Memory to Waveform
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
        start_time_idx = np.argmin(abs(h.t - start_time))
        J_ℰ -= J_ℰ[start_time_idx, :]

    return J_ℰ


def J_Nhat(h, Psi1):
    """Angular momentum aspect contribution to magnetic part of strain

    Calculated according to Eq. (16c) of 'Adding Gravitational Memory to Waveform
    Catalogs using BMS Balance Laws'

    Parameters
    ----------
    h : WaveformModes
        WaveformModes object corresponding to the strain
    Psi1 : WaveformModes
        WaveformModes object corresponding to Psi1

    Returns
    -------
    J_Nhat : WaveformModes
        Bondi angular momentum aspect contribution to the strain

    """

    h = MTS(h)
    Psi1 = MTS(Psi1)

    # # Note that the contributions from the last two terms drop out as soon as we
    # # take the imaginary part below.
    # m = mass_aspect(Psi2, h)
    # N̂ = 2 * Psi1 - 0.25 * (h.bar * h.eth) - h.t * m.eth - 0.125 * (h * h.bar).eth

    N̂ = 2 * Psi1 - 0.25 * (h.bar * h.eth)
    J_N̂ = 0.5j * 𝔇inverseLaplacianinverse(N̂.dot.ethbar.im).ethbar.ethbar

    return J_N̂


def J_J(h):
    """Angular momentum flux contribution to magnetic part of strain

    Calculated according to Eq. (16d) of 'Adding Gravitational Memory to Waveform
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
    h_with_mem : WaveformModes
        WaveformModes object corresponding to the strain with electric memory

    """
    h_with_mem = h.copy()
    h_with_mem += J_E(h, start_time=start_time).data

    return h_with_mem


# def BMS_strain(h, Psi2, news=None, start_time=None, match_time=None):
#     """Strain inferred from the BMS balance laws

#     Calculated according to Eqs. () of 'Adding Gravitational Memory to Waveform
#     Catalogs using BMS Balance Laws'

#     Parameters
#     ----------
#     h : WaveformModes
#         WaveformModes object corresponding to the strain
#     Psi2 : WaveformModes
#         WaveformModes object corresponding to Psi2
#     news : WaveformModes, optional
#         WaveformModes object corresponding to the news.  The default is to compute
#         this as the time-derivative of the strain, `h.dot`.
#     start_time : float, optional
#         Time at which the energy flux integral should begin.  Default is `h.t[0]`.
#     match_time : float, optional
#         Time at which the strain and BMS strain should match.

#     Returns
#     -------
#     (h_BMS, constraint) : tuple of WaveformModes

#     """

#     if news is None:
#         news = h.copy()
#         news.data = h.dot
#     h_mts = MTS(h)
#     news_mts = MTS(news)
#     Psi2_mts = MTS(Psi2)

#     Psi2_term = 0.5 * 𝔇inverse(-(Psi2_mts + 0.25 * news_mts * h_mts.bar)).ethbar.ethbar

#     h_BMS = h.copy()
#     E = mem_energy_flux_contribution(h, news=news, start_time=start_time)
#     h_BMS.data = E.data + np.array(Psi2_term[:, LM_index(2, -2, 0) :])

#     if not match_time is None:
#         match_time_idx = np.argmin(abs(h.t - match_time))
#         h_BMS.data = h_BMS.data + (h.data[match_time_idx, LM_index(2, -2, h.ell_min) :] - h_BMS.data[match_time_idx, :])

#     Constraint = h_BMS.copy()
#     Constraint.data = h.data[:, LM_index(2, -2, h.ell_min) :] - h_BMS.data

#     return (h_BMS, Constraint)


# def BMS_strain_other(h, Psi2, Psi1, news=None, start_time=None, match_time=None):
#     """Strain inferred from the BMS balance laws

#     Calculated according to Eqs. () of 'Adding Gravitational Memory to Waveform
#     Catalogs using BMS Balance Laws'

#     Parameters
#     ----------
#     h : WaveformModes
#         WaveformModes object corresponding to the strain
#     Psi2 : WaveformModes
#         WaveformModes object corresponding to Psi2
#     Psi1 : WaveformModes
#         WaveformModes object corresponding to Psi1
#     news : WaveformModes, optional
#         WaveformModes object corresponding to the news.  The default is to compute
#         this as the time-derivative of the strain, `h.dot`.
#     start_time : float, optional
#         Time at which the energy flux integral should begin.  Default is `h.t[0]`.
#     match_time : float, optional
#         Time at which the strain and BMS strain should match.

#     Returns
#     -------
#     (h_BMS, Constraint, M, E, Ndot, Jdot) : tuple of WaveformModes

#     """

#     h_BMS = h.copy()
#     M = mem_mass_aspect_contribution(h, Psi2, news=news)
#     E = mem_energy_flux_contribution(h, news=news, start_time=start_time)
#     Ndot = mem_angular_momentum_aspect_contribution(h, Psi1, news=news)
#     Jdot = mem_angular_momentum_flux_contribution(h, news=news)
#     h_BMS.data = M.data + E.data + Ndot.data + Jdot.data

#     if not match_time is None:
#         match_time_idx = np.argmin(abs(h.t - match_time))
#         h_BMS.data = h_BMS.data + (h.data[match_time_idx, LM_index(2, -2, h.ell_min) :] - h_BMS.data[match_time_idx, :])

#     Constraint = h_BMS.copy()
#     Constraint.data = h.data[:, LM_index(2, -2, h.ell_min) :] - h_BMS.data

#     return (h_BMS, Constraint, M, E, Ndot, Jdot)
