[![Test Status](https://github.com/sxs-collaboration/sxs/workflows/tests/badge.svg)](https://github.com/sxs-collaboration/sxs/actions)
[![Test Coverage](https://codecov.io/gh/sxs-collaboration/sxs/branch/master/graph/badge.svg)](https://codecov.io/gh/sxs-collaboration/sxs)
[![Documentation Status](https://readthedocs.org/projects/sxs/badge/?version=latest)](https://sxs.readthedocs.io/en/latest/?badge=latest)
[![PyPI Version](https://img.shields.io/pypi/v/sxs?color=)](https://pypi.org/project/sxs/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/sxs.svg?color=)](https://anaconda.org/conda-forge/sxs)
[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/sxs-collaboration/sxs/blob/master/LICENSE)


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

These will download and install the package.  You may also want to set some
sensible defaults to automatically download and cache data:

```bash
python -c "import sxs; sxs.write_config(download=True, cache=True)"
```

This will write the configuration file in the directory returned by
`sxs.sxs_directory("config")`, and data will be cached in the directory
returned by `sxs.sxs_directory("cache")`.
