"""Collection of functions to reduce data size

This submodule provides several related, though not quite directly compatible,
routines to reduce the number of points in a data set while retaining a given
level of accuracy:

* The original NINJA2 routine in `linear_bisection`
* A routine emulating `romspline`'s in `greedy_spline`
* A new routine in `peak_greed`

The original NINJA2 routine was limited by the fact that it had to be written
as a single easy-to-understand and easy-to-compile C file.  As a result, it
used linear interpolation, rather than the more reasonable cubic spline.  It is
included here for historical reasons, but probably should not generally be
used.

The `romspline` version is very slow, because the algorithm essentially scales
as N**2, where N is the number of points in the final output.  This is because
one point is added per iteration, each iteration requires constructing a new
spline for all the current points, and the splines dominate the timing.

The `peak_greed` routine gets around this problem by selecting multiple new
points on each iteration, guided by the peaks of the error quantity.  The
result has nearly the same number of points, but the timing scales with N,
typically providing results dozens to hundreds of times more quickly.

"""

from .suppression import suppressor, suppress
