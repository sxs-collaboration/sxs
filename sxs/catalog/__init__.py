"""Interface to SXS catalog of data



"""



def load():
    """Load SXS catalog from file, optionally downloading

    """
    raise NotImplementedError()


class Catalog(object):
    def __init__(self, *args, **kwargs):
        raise NotImplementedError()

    @property
    def description(self):
        return self['description']

    @property
    def modified(self):
        return self['modified']

    @property
    def records(self):
        return self['records']

    @property
    def simulations(self):
        return self['simulations']
