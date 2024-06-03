# Calling `sxs` from Julia

Julia makes it extraordinarily simple to call Python packages from within Julia.

!!! note

    This page describes calling the `sxs` Python package from Julia.  Instructions
    for calling Julia components from Python are found elsewhere.

## Optional preliminary step: A new environment

Julia also makes it much easier to work with "environments" than
Python does.  If you want to create one from working with `sxs` in
Julia, just create a new directory, `cd` into it, and run
```bash
julia --project=.
```
This will create a file named `Project.toml` in that directory
containing this new environment.  Every time you need to work with
this environment, you can tell Julia to use it by adding `--project=`
followed by the path to that directory.  Alternatively, if Julia is already running,
just hit `]` to enter the package manager, and run
```julia
activate /path/to/new/directory
```

## Getting Julia working with Python

Once Julia is running (with the project you want), just run
```julia
import Pkg
Pkg.add(["PythonCall", "CondaPkg"])
using PythonCall
import CondaPkg
CondaPkg.add("numba", channel="numba")
CondaPkg.add("llvmlite", channel="numba")
CondaPkg.add("sxs", channel="conda-forge")
```
Once everything is installed, the main difference from Python is that
you have to import the module as
```julia
sxs = pyimport("sxs")
```
Now, you can call functions from the `sxs` package nearly identically
to how they are called in Python.  For example, the simple call to
load a waveform from the README of this package is
```julia
waveform = sxs.load("SXS:BBH:0123/Lev/rhOverM", extrapolation_order=2)
```
That exact line works just in Julia like it does in Python.  (Note
that strings must use double quotes in Julia, and Python's `True` and
`False` are lowercased in Julia.)  And accessing members and methods
of Python objects is just as natural.  For example, `waveform.t` and
`waveform.data` return the same data as in Python.

See the [PythonCall
documentation](https://juliapy.github.io/PythonCall.jl/stable/pythoncall/)
for more details on how to call Python from Julia, and any conversions
that may need to be done.  Note that `pyconvert(Any, x)` will usually
convert `x` to its most natural Julia equivalent.
