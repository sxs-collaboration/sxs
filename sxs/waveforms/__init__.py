"""Containers for various types of waveforms

Currently, the most interesting object in this submodule is
[`WaveformModes`](../sxs.waveforms.waveform_modes/).

"""

from .waveform_mixin import WaveformMixin
from .waveform_modes import WaveformModes, WaveformModesDict
# from .waveform_grid import WaveformGrid
# from .waveform_signal import WaveformSignal

from .format_handlers import (
    nrar,
    rotating_paired_diff_multishuffle_bzip2,
    rotating_paired_xor_multishuffle_bzip2,
    spectre_cce_v1,
    grathena,
    lvcnr,
    maya,
)
from .format_handlers.lvc import to_lvc_conventions
from . import memory, transformations, alignment, norms

from .flux import (
    energy_flux,
    momentum_flux,
    angular_momentum_flux,
    poincare_fluxes,
)

WaveformModes.energy_flux = energy_flux
WaveformModes.momentum_flux = momentum_flux
WaveformModes.angular_momentum_flux = angular_momentum_flux
WaveformModes.poincare_fluxes = poincare_fluxes

# Map a format string to the module that should be used to load a waveform with that format
formats = {
    None: nrar,
    "": nrar,
    "maya": maya,
    "nrar": nrar,
    "lvcnr": lvcnr,
    "rotating_paired_xor_multishuffle_bzip2": rotating_paired_xor_multishuffle_bzip2,
    "rpxm": rotating_paired_xor_multishuffle_bzip2,
    "rpxmb": rotating_paired_xor_multishuffle_bzip2,
    "rotating_paired_diff_multishuffle_bzip2": rotating_paired_diff_multishuffle_bzip2,
    "rpdmb": rotating_paired_diff_multishuffle_bzip2,
    "SpECTRE_CCE_v1": spectre_cce_v1,
}
