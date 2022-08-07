"""Containers for various types of waveforms

Currently, the most interesting object in this submodule is [`WaveformModes`](../sxs.waveforms.waveform_modes/).

"""

from .waveform_mixin import WaveformMixin
from .waveform_modes import WaveformModes
# from .waveform_grid import WaveformGrid
# from .waveform_signal import WaveformSignal

from . import nrar, rotating_paired_xor_multishuffle_bzip2, rotating_paired_diff_multishuffle_bzip2, memory, transformations, alignment

# Map a format string to the module that should be used to load a waveform with that format
formats = {
    None: nrar,
    "": nrar,
    "nrar": nrar,
    "rotating_paired_xor_multishuffle_bzip2": rotating_paired_xor_multishuffle_bzip2,
    "rpxm": rotating_paired_xor_multishuffle_bzip2,
    "rpxmb": rotating_paired_xor_multishuffle_bzip2,
    "rotating_paired_diff_multishuffle_bzip2": rotating_paired_diff_multishuffle_bzip2,
    "rpdmb": rotating_paired_diff_multishuffle_bzip2,
}
