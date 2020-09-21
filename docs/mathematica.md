# Calling `sxs` from Mathematica

Newer versions of Mathematica (version 11.2 and greater) include a convenient
way to call Python code from a Mathematica session: you just tell Mathematica
where to find the `python` executable, and then evaluate Mathematica commands
containing the Python code you want to run.

## Getting Mathematica working with Python

There are two main points to consider when making sure that Mathematica will
interact correctly with Python:

  1. As of this writing, Mathematica does not support Python 3.8 or greater.
     Since `sxs` only supports python 3.6 or greater, this means the only
     possible versions are 3.6 or 3.7.
  2. In addition to the `sxs` package and its requirements, Mathematica
     requires the `pyzmq` package.

If you happen to be running Python 3.6 or 3.7, you can just install `pyzmq` and
move on.  (For example, replace `sxs` with `pyzmq` in [the original
installation directions](../#installation).)

Otherwise, the easiest thing to do is to create a separate environment with
just these requirements.  With
[`conda`](https://docs.anaconda.com/anaconda/install/), you can run

```bash
conda create -n sxs_mathematica python=3.7 sxs pyzmq
conda activate sxs_mathematica
```

The second command here just "activates" your new environment, which basically
just changes the `PATH` in your current terminal so that it uses this version
of python by default; that change will disappear the next time you open a
terminal or if you call `conda deactivate`.

Finally, you'll need to tell Mathematica the full path to your python
executable.  You can print that out with

```bash
python -c 'import sys; print(sys.executable)'
```

Just make sure that you run this command using the version of Python that you
actually intend to use — which was the purpose of the `conda activate` command
above.  You may never need to activate the environment again except to update;
just don't delete it if you intend to continue using it from Mathematica.


## Usage

In Mathematica, you can call something like the following:

```mathematica
session = StartExternalSession[<|
    "System" -> "Python",
    "Executable" -> "/full/path/to/bin/python"
    |>];
py = ExternalEvaluate[session];
```

Remember to set the full path to your python executable in the first command.

Running this code will start up a python session, which will remain running
until you quit this Mathematica kernel or call `DeleteObject[session]`.  And
you pass all your python code to the `py` function as strings.  To get started,
import the `sxs` module like this:

```mathematica
py["import sxs"];
```

To load the catalog, you can run

```mathematica
py["catalog = sxs.load('catalog')"]
```

The python session will hang onto this `catalog` object, and you can use it
again later.  For example, you can extract the table of metadata for all
simulations like this:

```mathematica
sims = py["catalog.table"];
```

The result is a Mathematica `Dataset`, which you can search and slice in ways
comparable to the `pandas` interface in python.  For example, to select systems
with mass ratios greater than 8 and both spin magnitudes less than 0.1, then
list only selected columns (initial mass ratio, and so on), use code like this:

```mathematica
sims[
 Select[#["reference_mass_ratio"] > 8
    && #["reference_dimensionless_spin1_mag"] < 0.1
    && #["reference_dimensionless_spin2_mag"] < 0.1 &],
 {"initial_mass_ratio", "initial_separation", "reference_chi1_mag",
  "reference_chi2_mag", "reference_eccentricity_bound"}]
```

Of course, more complex python code can be used.  For example, we could write
the selections above directly in python using the usual `pandas` interface to
achieve the same result.  But this approach is presumably less natural for
Mathematica users.

We can also extract horizon data:

```mathematica
py["horizons = sxs.load('SXS:BBH:1107/Lev/Horizons.h5')"]
time = Normal[py["horizons.A.time"]];
arealmass = Normal[py["horizons.A.areal_mass"]];
coordcenterinertial = Normal[py["horizons.A.coord_center_inertial"]];
```

Here, we have wrapped the `py` calls in `Normal`, which converts the numeric
arrays into standard Mathematica Lists, so that we can use them more naturally.
(If you can use the `NumericArray` object directly, just omit this call.)  For
example, to plot the coordinate trajectory of horizon A as a function of time,
we can run

```mathematica
ListLinePlot[
 Table[Transpose[{time, coordcenterinertial[[All, i]]}], {i, 1, 3}]]
```

to plot the components of the coordinate trajectory as functions of time.

We can similarly extract the data from waveforms:

```mathematica
py["waveform = sxs.load('SXS:BBH:1107/Lev/rhOverM', extrapolation_order=2)"]
time = Normal[py["waveform.time"]];
mode22 = Normal[py["waveform.data[:, waveform.index(2, 2)]"]];
ListLinePlot[{Transpose[{time, Re[mode22]}], Transpose[{time, Im[mode22]}]}]
```

This plots the real and imaginary parts of the (2,2) mode of the waveform.
And, of course, you can still access all the features of the python object
containing the waveform, like evaluation in some particular direction:

```mathematica
signal = Normal[py["waveform.evaluate(0.1, 0.2)"]]  (* (θ, ϕ) = (0.1, 0.2) *)
ListLinePlot[{Transpose[{time, Re[signal]}], Transpose[{time, Im[signal]}]}]
```
