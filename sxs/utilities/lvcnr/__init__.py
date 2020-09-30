"""Convert to/from LVC-NR format

The LVC-NR format was first described in https://arxiv.org/abs/1703.01076.  The de
facto format definition is in the code that checks it, which seems to be the
`lvcnrpy` module — more specifically, the classes in lvcnrpy.format.specs at
https://git.ligo.org/waveforms/lvcnrpy

Note that most of the code in this submodule is adapted from code written by Alyssa
Garcia, Geoffrey Lovelace, and Patricia Schmidt.

"""

from .dataset import Dataset
from .waveform_amp_phase import WaveformAmpPhase
from . import waveforms, horizons, metadata, comparisons, conversion
from .comparisons import compare
from .conversion import SimulationConverter, convert_simulation


def bbh_keys_from_simulation_keys(simulation_keys):
    """Extract BBH simulations from a list of all simulations

    Note that this function is maintained here for precise backwards-compatibility.
    More useful functions may be found in `sxs.utilities`.

    """
    return [simulation_key for simulation_key in simulation_keys
            if simulation_key.split(':')[-2] == "BBH"]
