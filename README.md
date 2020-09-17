[![Test Status](https://github.com/sxs-collaboration/sxs/workflows/tests/badge.svg)](https://github.com/sxs-collaboration/sxs/actions)
[![Documentation Status](https://readthedocs.org/projects/sxs/badge/?version=main)](https://sxs.readthedocs.io/en/main/?badge=main)
[![PyPI Version](https://img.shields.io/pypi/v/sxs?color=)](https://pypi.org/project/sxs/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/sxs.svg?color=)](https://anaconda.org/conda-forge/sxs)
[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/sxs-collaboration/sxs/blob/main/LICENSE)


# Simulating eXtreme Spacetimes python package

The `sxs` python package provides a high-level interface for using data
produced by the SXS collaboration.  In particular, the function `sxs.load` can
automatically find, download, and load data, returning objects that provide
common interfaces to the various types of data, without forcing the user to
worry about details like data formats or where to find the data.  It can also
automatically select the newest or highest-resolution dataset for a given
simulation, or return a range of versions or resolutions.  Currently, the
high-level objects encapsulate

  * Catalog — a listing of all data produced by the SXS collaboration
  * Metadata — data describing the simulation parameters
  * Horizons — time-series data describing the apparent horizons
  * Waveforms — time-series data describing the extrapolated gravitational-wave
    modes


## Installation

Because this package is pure python code, installation is very simple.  In
particular, with a reasonably modern installation, you can just run a command
like

```bash
conda install -c conda-forge sxs
```

or

```bash
python -m pip install sxs
```

Either of these will download and install the package.  You may also want to
set some sensible defaults to automatically download and cache data:

```bash
python -c "import sxs; sxs.write_config(download=True, cache=True)"
```

This will write the configuration file in the directory returned by
`sxs.sxs_directory("config")`, and data will be cached in the directory
returned by `sxs.sxs_directory("cache")`.


## Usage

There are four important objects to understand in this package:

```python
import sxs

catalog = sxs.load("catalog")
metadata = sxs.load("SXS:BBH:0123/Lev/metadata.json")
horizons = sxs.load("SXS:BBH:0123/Lev/Horizons.h5")
waveform = sxs.load("SXS:BBH:0123/Lev/rhOverM", extrapolation_order=2)
```

[The `catalog` object](api/sxs.catalog.catalog/#sxs.catalog.catalog.Catalog)
contains information about every simulation in the catalog, including all
available data files, and information about how to get them.  You probably
don't need to actually know about details like where to get the data, but
`catalog` can help you find the simulations you care about.  Most importantly,
`catalog.simulations` is a `dict` object, where the keys are names of
simulations (like "SXS:BBH:0123") and the values are the same types as [the
`metadata` object](api/sxs.metadata.metadata/#sxs.metadata.metadata.Metadata),
which contains metadata about that simulation — things like mass ratio, spins,
etc.  This `metadata` reflects the actual output of the simulations, which
leads to some inconsistencies in their formats.  A more consistent interface
(though it is biased toward returning NaNs where a human might glean more
information) is provided by `catalog.simulations_dataframe`, which returns a
[`pandas`](https://pandas.pydata.org/docs/) `DataFrame` with specific data
types for each column.

The actual data itself is primarily contained in the next two objects.  [The
`horizons` object](api/sxs.horizons/#sxs.horizons.Horizons) has three
attributes — `horizons.A`, `horizons.B`, and `horizons.C` — typically
representing the original two horizons of the black-hole binary and the common
horizon that forms at merger.  In matter simulations, one or more of these may
be `None`.  Otherwise, each of these three is a
[`HorizonQuantities`](api/sxs.horizons/#sxs.horizons.HorizonQuantities) object,
containing several timeseries relating to mass, spin, and position.

Finally, the
[`waveform`](api/sxs.waveforms.waveform_modes/#sxs.waveforms.waveform_modes.WaveformModes)
encapsulates the modes of the waveform and the corresponding time information,
along with relevant metadata like data type, spin weight, etc., and useful
features like numpy-array-style slicing.
