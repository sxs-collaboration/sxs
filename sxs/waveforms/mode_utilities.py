"""Utility functions for WaveformModes objects"""

import numpy as np
from .. import jit


@jit
def _expectation_value_Ldt(data, datadot, lm, Ldt):
    """Helper function for the WaveformModes.expectation_value_Ldt method"""
    # L+ = Lx + i Ly      Lx =    (L+ + L-) / 2
    # L- = Lx - i Ly      Ly = -i (L+ - L-) / 2

    for i_mode in range(lm.shape[0]):
        L = lm[i_mode, 0]
        M = lm[i_mode, 1]
        for i_time in range(data.shape[0]):
            # Compute first in (+,-,z) basis
            Lp = (
                np.conjugate(data[i_time, i_mode + 1]) * datadot[i_time, i_mode] * np.sqrt(L * (L + 1) - M * (M + 1))
                if M + 1 <= L
                else 0.0 + 0.0j
            )
            Lm = (
                np.conjugate(data[i_time, i_mode - 1]) * datadot[i_time, i_mode] * np.sqrt(L * (L + 1) + M * (-M + 1))
                if M - 1 >= -L
                else 0.0 + 0.0j
            )
            Lz = np.conjugate(data[i_time, i_mode]) * datadot[i_time, i_mode] * M

            # Convert into (x,y,z) basis
            Ldt[i_time, 0] += 0.5 * (Lp.imag + Lm.imag)
            Ldt[i_time, 1] += -0.5 * (Lp.real - Lm.real)
            Ldt[i_time, 2] += Lz.imag
    return


@jit
def _expectation_value_LL(data, lm, LL):
    """Helper function for the WaveformModes.expectation_value_LL method"""
    # Big, bad, ugly, obvious way to do the calculation
    # =================================================
    # L+ = Lx + i Ly      Lx =    (L+ + L-) / 2     Im(Lx) =  ( Im(L+) + Im(L-) ) / 2
    # L- = Lx - i Ly      Ly = -i (L+ - L-) / 2     Im(Ly) = -( Re(L+) - Re(L-) ) / 2
    # Lz = Lz             Lz = Lz                   Im(Lz) = Im(Lz)
    # LxLx =   (L+ + L-)(L+ + L-) / 4
    # LxLy = -i(L+ + L-)(L+ - L-) / 4
    # LxLz =   (L+ + L-)(  Lz   ) / 2
    # LyLx = -i(L+ - L-)(L+ + L-) / 4
    # LyLy =  -(L+ - L-)(L+ - L-) / 4
    # LyLz = -i(L+ - L-)(  Lz   ) / 2
    # LzLx =   (  Lz   )(L+ + L-) / 2
    # LzLy = -i(  Lz   )(L+ - L-) / 2
    # LzLz =   (  Lz   )(  Lz   )

    for i_mode in range(lm.shape[0]):
        L = lm[i_mode, 0]
        M = lm[i_mode, 1]
        for i_time in range(data.shape[0]):
            # Compute first in (+,-,z) basis
            LpLp = (
                np.conjugate(data[i_time, i_mode + 2]) * data[i_time, i_mode]
                * np.sqrt((L * (L + 1) - (M + 1) * (M + 2)) * (L * (L + 1) - M * (M + 1)))
                if M + 2 <= L
                else 0.0 + 0.0j
            )
            LpLm = (
                np.conjugate(data[i_time, i_mode]) * data[i_time, i_mode]
                * np.sqrt((L * (L + 1) - (M - 1) * M) * (L * (L + 1) + M * (-M + 1)))
                if M - 1 >= -L
                else 0.0 + 0.0j
            )
            LmLp = (
                np.conjugate(data[i_time, i_mode]) * data[i_time, i_mode]
                * np.sqrt((L * (L + 1) - (M + 1) * M) * (L * (L + 1) - M * (M + 1)))
                if M + 1 <= L
                else 0.0 + 0.0j
            )
            LmLm = (
                np.conjugate(data[i_time, i_mode - 2]) * data[i_time, i_mode]
                * np.sqrt((L * (L + 1) - (-M + 1) * (-M + 2)) * (L * (L + 1) + M * (-M + 1)))
                if M - 2 >= -L
                else 0.0 + 0.0j
            )
            LpLz = (
                np.conjugate(data[i_time, i_mode + 1]) * data[i_time, i_mode]
                * (np.sqrt(L * (L + 1) - M * (M + 1)) * M)
                if M + 1 <= L
                else 0.0 + 0.0j
            )
            LzLp = (
                np.conjugate(data[i_time, i_mode + 1]) * data[i_time, i_mode]
                * ((M + 1) * np.sqrt(L * (L + 1) - M * (M + 1)))
                if M + 1 <= L
                else 0.0 + 0.0j
            )
            LmLz = (
                np.conjugate(data[i_time, i_mode - 1]) * data[i_time, i_mode]
                * (np.sqrt(L * (L + 1) + M * (-M + 1)) * M)
                if M - 1 >= -L
                else 0.0 + 0.0j
            )
            LzLm = (
                np.conjugate(data[i_time, i_mode - 1]) * data[i_time, i_mode]
                * ((M - 1) * np.sqrt(L * (L + 1) + M * (-M + 1)))
                if M - 1 >= -L
                else 0.0 + 0.0j
            )
            LzLz = np.conjugate(data[i_time, i_mode]) * data[i_time, i_mode] * M ** 2

            # Convert into (x,y,z) basis
            LxLx = 0.25 * (LpLp + LmLm + LmLp + LpLm)
            LxLy = -0.25j * (LpLp - LmLm + LmLp - LpLm)
            LxLz = 0.5 * (LpLz + LmLz)
            LyLx = -0.25j * (LpLp - LmLp + LpLm - LmLm)
            LyLy = -0.25 * (LpLp - LmLp - LpLm + LmLm)
            LyLz = -0.5j * (LpLz - LmLz)
            LzLx = 0.5 * (LzLp + LzLm)
            LzLy = -0.5j * (LzLp - LzLm)
            # LzLz = (LzLz)

            # Symmetrize
            LL[i_time, 0, 0] += LxLx.real
            LL[i_time, 0, 1] += (LxLy + LyLx).real / 2.0
            LL[i_time, 0, 2] += (LxLz + LzLx).real / 2.0
            LL[i_time, 1, 0] += (LyLx + LxLy).real / 2.0
            LL[i_time, 1, 1] += LyLy.real
            LL[i_time, 1, 2] += (LyLz + LzLy).real / 2.0
            LL[i_time, 2, 0] += (LzLx + LxLz).real / 2.0
            LL[i_time, 2, 1] += (LzLy + LyLz).real / 2.0
            LL[i_time, 2, 2] += LzLz.real
    return
