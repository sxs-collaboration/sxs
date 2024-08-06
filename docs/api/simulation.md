# `Simulation` classes

The `Simulation` classes are the primary interface for accessing
simulation data.  They are the top-level objects that contain all the
data for a single simulation.  These classes are designed to be
lightweight and to load data only when needed.

Note that the user will not usually need to interact with these
classes (or the factory) automatically; the `sxs.load` function will
automatically select the appropriate version of the simulation and
return the appropriate object.

## `sxs.Simulation` factory function
::: sxs.Simulation

## `SimulationBase` class
::: sxs.simulations.simulation.SimulationBase

## `Simulation_v1` class
::: sxs.simulations.simulation.Simulation_v1

## `Simulation_v2` class
::: sxs.simulations.simulation.Simulation_v2
