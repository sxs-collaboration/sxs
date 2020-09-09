# Copyright (c) 2020, Michael Boyle
# See LICENSE file for details:
# <https://github.com/sxs-collaboration/sxs/blob/master/LICENSE>

"""Interface to SXS data and utilities"""

try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:  # pragma: no cover
    import importlib_metadata

__version__ = importlib_metadata.version(__name__)

from . import utilities
from .utilities import file_format, sxs_directory, read_config, write_config
from .time_series import TimeSeries
from .catalog import Catalog
from .metadata import Metadata
from .horizons import Horizons, HorizonQuantities
from .waveforms import WaveformModes, WaveformGrid, WaveformSignal
from .loader import load
from . import catalog, metadata, horizons, waveforms
