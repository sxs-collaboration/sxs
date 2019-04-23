"""Validate that SXS simulations are complete and suitable for upload

This submodule contains several functions for validating SXS simulations, primarily the functions
`bbh`, `bhns`, and `nsns`, which validate the corresponding types of systems.  See each function's
docstring for details about what is required.

"""


class _Validity(object):
    """Simple utility class to track validity and print or raise messages on invalidation"""
    def __init__(self, short_circuit):
        self.valid = True
        self.short_circuit = short_circuit

    def invalid(self, message):
        self.valid = False
        if self.short_circuit:
            raise ValueError(message)
        else:
            print(message)





def bhns(path):
    """Check that a path contains a valid BHNS simulation
    """
    raise NotImplementedError()


def nsns(path):
    """Check that a path contains a valid NSNS simulation
    """
    raise NotImplementedError()


