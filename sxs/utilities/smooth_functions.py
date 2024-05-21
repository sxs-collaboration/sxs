import numpy as np
from . import jit


def transition_function(x, x0, x1, y0=0, y1=1, msquared=3):
    """Evaluate a C^∞ function transitioning from `y0` to `y1` in the range (`x0`, `x1`).

    This is a special type of [sigmoid
    function](https://en.wikipedia.org/wiki/Sigmoid_function) in which
    the transition region has finite size (between `x0` and `x1`), the
    function is constant outside of that region (equal to `y0` for
    values of `x <= x0` and `y1` for values of `x >= x1`), but the
    function is monotonic and infinitely differentiable everywhere —
    including *at* `x0` and `x1`.

    This particular function might be more descriptively called the
    "compactified tanh" function, because the nontrivial part is
    essentially the tanh function with its domain `(-∞,∞)` remapped —
    or compactified — to `(x0,x1)`.  This compactification could be
    applied to any sigmoid function; here we use the standard tanh
    function.

    The parameter `msquared` controls the abruptness of the
    transition.  The value `m` is the "fractional" slope of the
    function at the midpoint — that is, the slope is `m * (y1 - y0) /
    (x1 - x0)`.  The default value of `msquared = 3` corresponds to
    the most gradual transition with only one inflection point.
    Larger values will be more abrupt; smaller values will have
    multiple inflection points.  The least extreme second derivatives
    are obtained when `msquared` is approximately 4.0242806285.

    Parameters
    ----------
    x : array-like
        Input values.
    x0 : float
        Lower bound of the transition range.
    x1 : float
        Upper bound of the transition range.
    y0 : float, optional
        Function value at and below `x0`. Defaults to 0.
    y1 : float, optional
        Function value at and above `x1`. Defaults to 1.
    msquared : float, optional
        Controls the abruptness of the transition. Defaults to 3.

    Returns
    -------
    array-like
        The function values at `x`.
    """
    y = np.empty_like(x, dtype=np.promote_types(x.dtype, type(y1)))
    return transition_function_inplace(y, x, x0, x1, y0, y1, msquared)


@jit
def transition_function_inplace(y, x, x0, x1, y0=0, y1=1, msquared=3):
    """Evaluate `transition_function` in place on `y`
    
    See that function's docstring for details.
    """
    assert x0 <= x1
    m = np.sqrt(msquared)
    Δy = (y1 - y0) / 2
    Δx = (x1 - x0) / 2
    ȳ = (y0 + y1) / 2
    x̄ = (x0 + x1) / 2
    for i in range(x.size):
        x̂ = (x[i] - x̄) / Δx
        # Note that the conditions after the `or`s in the following
        # conditions will be *almost* always equivalent to the first
        # conditions, but they are necessary to protect against
        # roundoff errors.
        if x̂ <= -1 or x[i] <= x0:
            y[i] = y0
        elif x̂ >= 1 or x[i] >= x1:
            y[i] = y1
        else:
            y[i] = ȳ + Δy * np.tanh(m * x̂ / (1 - x̂**2))
    return y


def transition_function_derivative(x, x0, x1, y0=0, y1=1, msquared=3):
    """Evaluate the derivative of `transition_function`.

    See the `transition_function` docstring for details about the
    parameters.
    """
    y = np.empty_like(x, dtype=np.promote_types(x.dtype, type(y1)))
    return transition_function_derivative_inplace(y, x, x0, x1, y0, y1, msquared)


@jit
def transition_function_derivative_inplace(yprime, x, x0, x1, y0=0, y1=1, msquared=3):
    """Evaluate `transition_function_derivative` in place on `yprime`
    
    See that function's docstring for details.
    """
    assert x0 <= x1
    m = np.sqrt(msquared)
    Δy = (y1 - y0) / 2
    Δx = (x1 - x0) / 2
    x̄ = (x0 + x1) / 2
    for i in range(x.size):
        x̂ = (x[i] - x̄) / Δx
        if x̂ <= -1 or x[i] <= x0 or x̂ >= 1 or x[i] >= x1:
            yprime[i] = 0
        else:
            yprime[i] = Δy * m * (1 + x̂**2) / (Δx * (1 - x̂**2)**2 * np.cosh(m * x̂ / (1 - x̂**2))**2)
    return yprime


def transition_to_constant(f, t, t1, t2):
    """Smoothly transition from the function to a constant.

    This works (implicitly) by multiplying the derivative of `f` with
    the transition function going from 1 at `t<=t1` to 0 at `t>=t2`,
    and then integrating.  Using integration by parts, this simplifies
    to multiplying `f` itself by the transition function, and then
    subtracting the integral of `f` times the derivative of the
    transition function.  This integral is effectively restricted to
    the region (t1, t2).  Note that the final value (after t2) will
    depend on the precise values of `t1` and `t2`, and the behavior of
    `f` in between.

    Parameters
    ----------
    f : array_like
        Array corresponding to the following `t` parameter.  If this
        array is multi-dimensional, the first dimension must be the
        same size as `t`.
    t : array_like
        One-dimensional monotonic array of floats.
    t1 : float
        Value before which the output will equal `f`.
    t2 : float
        Value after which the output will be constant.

    """
    fprime = np.asanyarray(f).copy()
    transition_to_constant_inplace(fprime, t, t1, t2)
    return fprime


def transition_to_constant_inplace(f, t, t1, t2):
    """Evaluate `transition_to_constant` in place on `f`
    
    See that function's docstring for details.
    """
    from scipy.interpolate import CubicSpline
    assert f.shape[0] == t.size
    i1, i2 = np.searchsorted(t, [t1, t2])
    τ = transition_function(t, t1, t2, y0=1, y1=0)
    τ̇ = transition_function_derivative(t, t1, t2, y0=1, y1=0)
    # Note that numpy's broadcasting rules add extra dimensions on the
    # left, which means we need to move the time dimension to the
    # right to multiply `f` by `t`, then move it back.  This is the
    # reason for the double transposes.
    intfτ̇ = CubicSpline(
        t[i1:i2], (f[i1:i2].T * τ̇[i1:i2]).T
    ).antiderivative()(t[i1:i2])
    f[i1:i2] = (f[i1:i2].T * τ[i1:i2]).T - intfτ̇
    f[i2:] = f[i2 - 1]
    return f


def bump_function(x, x0, x1, x2, x3, y0=0, y12=1, y3=None, m01squared=3, m23squared=None):
    """Evaluate a C^∞ "bump" function.

    This is essentially the same as `transition_function`, but with
    *two* transitions.  See that function's documentation for details
    about the meanings of the parameters.

    In short, the function's value is `y0` for `x <= x0`, smoothly
    transitions from there to `y12` on `x0 < x < x1`, stays at `y12`
    for `x1 <= x <= x2`, smoothly transitions from there to `y3` for
    `x2 < x < x3`, and stays at `y3` for `x >= x3`.

    Parameters
    ----------
    x : array-like
        Input values.
    x0 : float
        Lower bound of the first transition range.
    x1 : float
        Upper bound of the first transition range.
    x2 : float
        Lower bound of the second transition range.
    x3 : float
        Upper bound of the second transition range.
    y0 : float, optional
        Function value at and below `x0`. Defaults to 0.
    y12 : float, optional
        Function value in the first transition range. Defaults to 1.
    y3 : float, optional
        Function value at and above `x3`. Defaults to 0.
    m01squared : float, optional
        Controls the abruptness of the first transition. Defaults to 3.
    m23squared : float, optional
        Controls the abruptness of the second transition. Defaults to `m01squared`.

    Returns
    -------
    array-like
        The function values at `x`.
    """
    y3 = y3 if y3 is not None else y0
    y = np.empty_like(x, dtype=np.promote_types(x.dtype, type(y3)))
    return bump_function_inplace(y, x, x0, x1, x2, x3, y0, y12, y3, m01squared, m23squared)


@jit
def bump_function_inplace(y, x, x0, x1, x2, x3, y0=0, y12=1, y3=0, m01squared=3, m23squared=None):
    """Evaluate `bump_function` in place on `y`
    
    See that function's docstring for details.
    """
    assert x0 <= x1 <= x2 <= x3
    m01 = np.sqrt(m01squared)
    m23 = np.sqrt(m23squared if m23squared is not None else m01squared)
    Δy01 = (y12 - y0) / 2
    Δy23 = (y3 - y12) / 2
    Δx01 = (x1 - x0) / 2
    Δx23 = (x3 - x2) / 2
    ȳ01 = (y0 + y12) / 2
    ȳ23 = (y12 + y3) / 2
    x̄01 = (x0 + x1) / 2
    x̄23 = (x2 + x3) / 2
    for i in range(x.size):
        x̂01 = (x[i] - x̄01) / Δx01
        x̂23 = (x[i] - x̄23) / Δx23
        # Note that the conditions after the `or`s in the following
        # conditions will be *almost* always equivalent to the first
        # conditions, but they are necessary to protect against
        # roundoff errors.
        if x̂01 <= -1 or x[i] <= x0:
            y[i] = y0
        elif (x̂01 >= 1 or x[i] >= x1) and (x̂23 <= -1 or x[i] <= x2):
            y[i] = y12
        elif x̂23 >= 1 or x[i] >= x3:
            y[i] = y3
        elif -1 < x̂01 < 1:
            y[i] = ȳ01 + Δy01 * np.tanh(m01 * x̂01 / (1 - x̂01**2))
        else:  # -1 < x̂23 < 1
            y[i] = ȳ23 + Δy23 * np.tanh(m23 * x̂23 / (1 - x̂23**2))
    return y
