import contextlib


@contextlib.contextmanager
def environment_context(*undefine, **define):
    """Temporarily alter the `os.environ` dictionary

    Parameters
    ----------
    undefine : tuple of str
        Each string passed as a positional parameter is removed (if found) from the
        environment within the context.
    define : dict
        Each keyword argument sets an environment variable within the context.

    Examples
    --------
    >>> with environment_context("PWD", HOME="/some/crazy/path"):
    >>>     print(os.environ.get("PWD", None))  # prints None
    >>>     print(os.environ["HOME"])  # prints "/some/crazy/path"

    """
    import os
    env = os.environ

    # Set of current environment variables being altered
    conflicting = (set(define) | set(undefine)) & set(os.environ)

    # Environment variables to restore on exit
    define_after = {k: os.environ[k] for k in conflicting}

    # Environment variables to remove on exit
    undefine_after = tuple(k for k in define if k not in os.environ)

    try:
        for k in undefine:
            os.environ.pop(k, None)
        os.environ.update(define)
        yield
    finally:
        for k in undefine_after:
            os.environ.pop(k)
        os.environ.update(define_after)
