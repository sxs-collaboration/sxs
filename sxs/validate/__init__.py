"""Check that SXS simulations are complete and suitable for upload

This submodule contains sub-submodules `bbh`, `bhns`, and `nsns`, which validate the corresponding
types of systems.

"""

from . import bbh
from . import bhns
from . import nsns


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
