"""Interface to SXS catalog of data



"""

import functools


@functools.lru_cache()
def load():
    """Load SXS catalog from file, optionally downloading

    """
    raise NotImplementedError()


class Catalog(object):
    def __init__(self, *args, **kwargs):
        raise NotImplementedError()

    def __call__(self, *args, **kwargs):
        # construct a map of keys in SXS IDs (with versions) and values of
        #     maps of keys in levs and values of
        #         maps of keys in file names and values of
        #             URLs

        # Split input string by '/'
        # i = 0
        # j = 0
        # for j in range(3)
        #   search map_j for split
        #   if exact match i += 1 and select match
        #   if partial match i += 1 and select highest match
        #   if regex match i += 1 and select all matches
        #   if not select highest key

        # "SXS:BBH:0123v4/Lev4/h_Extrapolated_N2.h5"

        # return tuple of 2-tuples (full_match_string, URL)

    @property
    @functools.lru_cache()
    def nested_maps(self):
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
