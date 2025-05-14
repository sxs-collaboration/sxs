[![Test Status](https://github.com/sxs-collaboration/sxs/workflows/tests/badge.svg)](https://github.com/sxs-collaboration/sxs/actions)
[![Documentation Status](https://readthedocs.org/projects/sxs/badge/?version=main)](https://sxs.readthedocs.io/en/main/?badge=main)
[![PyPI Version](https://img.shields.io/pypi/v/sxs?color=)](https://pypi.org/project/sxs/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/sxs.svg?color=)](https://anaconda.org/conda-forge/sxs)
[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/sxs-collaboration/sxs/blob/main/LICENSE)
[![Open with marimo](https://marimo.io/shield.svg)](https://data.black-holes.org/)


# Simulating eXtreme Spacetimes package

The `sxs` python package provides a high-level interface for using data
produced by the SXS collaboration.  In particular, the function `sxs.load` can
automatically find, download, and load data, returning objects that provide
common interfaces to the various types of data, without forcing the user to
worry about details like data formats or where to find the data.  It can also
automatically select the newest or highest-resolution dataset for a given
simulation, or return a range of versions or resolutions.  Currently, the
high-level objects encapsulate

  * Dataframe — a catalog of all simulations produced by the SXS collaboration
  * Simulation — an object encapsulating all data for a single simulation
  * Metadata — data describing the simulation parameters
  * Horizons — time-series data describing the apparent horizons
  * Waveform — time-series data describing the extrapolated gravitational-wave
    modes


## Installation

Because this package is pure python code, installation is very simple.  In
particular, with a reasonably modern installation, you can just run a command
like

```bash
python -m pip install sxs
```

or

```bash
mamba install -c conda-forge sxs
```

Here, the first command assumes that you have an appropriate python
environment set up in some other way;
[`mamba`](https://mamba.readthedocs.io/en/latest/index.html) is the
newer replacement for `conda`, and is a convenient way to install
python and manage environments.  Either of these commands will
download and install the `sxs` package and its most vital
requirements.

You may also want to set some convenient defaults to automatically
download and cache data:

```bash
python -c "import sxs; sxs.write_config(download=True, cache=True)"
```

This will write a configuration file in the directory returned by
`sxs.sxs_directory("config")`, and downloaded data will be cached in the
directory returned by `sxs.sxs_directory("cache")`.  See [that function's
documentation](https://sxs.readthedocs.io/en/main/api/sxs.utilities.sxs_directories/#sxsutilitiessxs_directoriessxs_directory)
for details.

## Citing this package and/or data

If you use this package and/or the data it provides in your research,
please cite them, including the *specific version of the data* that
you use (see below).  To help with this, we provide the function
`sxs.cite`.  Use `print(sxs.cite())` to see BibTeX citations for the
version of this package you are using, the most recent paper
describing the catalog, and the catalog data itself.  Use, e.g.,
`print(sxs.cite("SXS:BBH:0001", "SXS:BBH:4001"))` to include citations
for those specific simulations *and* the papers that introduced them.

## Usage

An extensive demonstration of this package's capabilities is available
[here](https://mybinder.org/v2/gh/moble/sxs_notebooks/main), in the form of
interactive jupyter notebooks that are actually running this code and some
pre-downloaded data.  The following is just a very brief overview of the `sxs`
package's main components.

### Loading a specific version of the catalog

For the purposes of reproducibility — both reproducing your own
results and allowing others to reproduce them — it is important to be
aware of which version of the catalog you are using, and to cite it
when you publish results using SXS data.  Whenever `sxs` tries to load
data, it most first load some version of the catalog.  If you do not
specify a version, it will automatically find and load the most recent
version available via github, and print out a message telling you
which version it is using, like

```python
Loading SXS simulations using latest tag 'v3.0.0', published at 2025-05-12T10:00:00Z.
```

For the rest of that Python session, all data loaded will be from that
version of the catalog.  If you want to use a different version, you
can specify it explicitly while loading the catalog — preferably as

```python
sxs.load("dataframe", tag="3.0.0")
```

Even if you do not use the returned object from this command, it will
ensure that all data will be loaded from the specified version of the
catalog.  Thus, it is best practice to make this call as soon as you
import the `sxs` package.


### Interacting with the data


There are four important objects to understand in this package:

```python
import sxs

# Load a specific version of the catalog for reproducibility
df = sxs.load("dataframe", tag="3.0.0")

# Load a specific simulation
sim = sxs.load("SXS:BBH:4001")

# Obtain data about the horizons
horizons = sim.horizons

# Obtain data about the gravitational-wave strain
h = sim.h
```

Note that `tag` is optional, but is good to include because it sets
the version of the catalog from which data is loaded, which ensures
reproducibility.  Leave it out to see the most recent version
available, and then use that version consistently in any analysis.  Be
sure to cite the specific version of the catalog you used in any
publications.

[The "dataframe"](https://sxs.readthedocs.io/en/main/api/dataframe/)
`df` contains information about every simulation in the catalog,
including all available data files, and information about how to get
them.  You probably don't need to actually know about details like
where to get the data, but `df` can help you find the simulations you
care about.  It is a
[`pandas.DataFrame`](https://pandas.pydata.org/docs/) object, where
the rows are names of simulations (like "SXS:BBH:0123") and the
columns include [the
`metadata`](https://sxs.readthedocs.io/en/main/api/sxs.metadata.metadata/#sxs.metadata.metadata.Metadata)
for the simulations — things like mass ratio, spins, eccentricity,
etc. — in addition to extra refinements like spin magnitudes, etc.

Once you have found a simulation you want to work with, you can load
it with, e.g., `sxs.load("SXS:BBH:4001")`, which will return a
[`Simulation`](https://sxs.readthedocs.io/en/main/api/sxs.simulation/#sxs.simulation.Simulation)
object, which contains metadata about the simulation, and allows you
to load data from the simulation.  By default, it uses the
highest-resolution run of the simulation, though this lower
resolutions can be specified.

The actual data itself is primarily contained in the next two objects.  [The
`horizons`
object](https://sxs.readthedocs.io/en/main/api/sxs.horizons/#sxs.horizons.Horizons)
has three attributes — `horizons.A`, `horizons.B`, and `horizons.C` — typically
representing the original two horizons of the black-hole binary and the common
horizon that forms at merger.  In matter simulations, one or more of these may
be `None`.  Otherwise, each of these three is a
[`HorizonQuantities`](https://sxs.readthedocs.io/en/main/api/sxs.horizons/#sxs.horizons.HorizonQuantities)
object, containing several timeseries relating to mass, spin, and position.

Finally, the [`h`
waveform](https://sxs.readthedocs.io/en/main/api/sxs.waveforms.waveform_modes/#sxs.waveforms.waveform_modes.WaveformModes)
encapsulates the modes of the strain waveform and the corresponding
time information, along with relevant metadata like data type, spin
weight, etc., with useful features like numpy-array-style slicing.

There is also `psi4` data available, which is computed with entirely
different methods; `h` and `psi4` are not just computed one from the
other by a double integral or differentiation.  As a result, we
generally recommend using `h` instead of `psi4` unless you have very
specific requirements.

## Contributing

Contributions are welcome!  There are at least two ways to contribute
to this codebase:

1. If you find a bug or want to suggest an enhancement, use the [issue
   tracker](https://github.com/sxs-collaboration/sxs/issues) on
   GitHub.  It's a good idea to look through past issues, too, to see
   if anybody has run into the same problem or made the same
   suggestion before.
2. If you will write or edit the python code, we use the [fork and
   pull
   request](https://help.github.com/articles/creating-a-pull-request-from-a-fork/)
   model.

You are also allowed to make use of this code for other purposes, as
detailed in the [MIT license](LICENSE).  For any type of contribution,
please follow the [code of
conduct](https://github.com/sxs-collaboration/.github/blob/master/CODE_OF_CONDUCT.md).

## Reporting catalog data issues

If you find an issue with our data or metadata, please let us know!
[Fill out an issue with the catalog data
template](https://github.com/sxs-collaboration/sxs/issues/new?template=catalog-data-issue-template.md)
and we will take a look as soon as possible.
