# Copyright (c) 2020, Michael Boyle
# See LICENSE file for details:
# <https://github.com/sxs-collaboration/sxs/blob/main/LICENSE>

"""Interface to SXS data and utilities"""

try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:  # pragma: no cover
    import importlib_metadata

__version__ = importlib_metadata.version(__name__)

from . import utilities
from .utilities import (
    file_format, sxs_directory, read_config, write_config, sxs_id, lev_number,
    jit, vectorize, guvectorize, version_info
)
from .time_series import TimeSeries
from .metadata import Metadata
from .catalog import Catalog
from .horizons import Horizons, HorizonQuantities
from .waveforms import WaveformModes #, WaveformGrid, WaveformSignal
from .waveforms import rotating_paired_xor_multishuffle_bzip2 as rpxmb
from .waveforms import rotating_paired_diff_multishuffle_bzip2 as rpdmb
from . import catalog, metadata, horizons, waveforms, zenodo, caltechdata
from .handlers import load, loadcontext

# The speed of light is, of course, defined to be exactly
speed_of_light = 299_792_458.0  # m/s

# By "IAU 2012 Resolution B2", the astronomical unit is defined to be exactly
astronomical_unit = 149_597_870_700.0  # m

# The parsec is, in turn, defined as "The distance at which 1 au subtends 1 arc sec: 1 au divided by
# pi/648000."  Thus, the extended-precision value of the parsec in meters is
parsec_in_meters = 3.0856775814913672789139379577965e16  # m

# The value of the solar mass parameter G*M_sun is known to higher accuracy than either of its
# factors.  The value here is taken from the publication "2021 Selected Astronomical Constants",
# which can be found at <https://asa.hmnao.com/static/files/2021/Astronomical_Constants_2021.pdf>.
# In the TDB (Barycentric Dynamical Time) time scale — which seems to be the more relevant one and
# looks like the more standard one for LIGO — we have
solar_mass_parameter = 1.32712440041e20  # m^3/s^2

# Dividing by the speed of light squared, we get the mass of the sun in meters; dividing again, we
# get the mass of the sun in seconds. For consistency with the nominal value above, we retain more
# precision than is warranted by the measurement.
m_sun_in_meters = 1476.6250385063112526099633973363  # m
m_sun_in_seconds = 4.925490949162941425997992909269e-06  # s
