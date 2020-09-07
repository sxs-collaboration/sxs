"""Interface to the catalog of SXS data

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

    def select(self, location, /, *args, **kwargs):
        import re
        # construct a map of keys in SXS IDs (with versions) and values of
        #     maps of keys in levs and values of
        #         maps of keys in file names and values of
        #             URLs

        # Split input string by '/'
        # i = 0
        # j = 0
        # for j in range(4)
        #   search map_j for split
        #   if exact match i += 1 and select match
        #   if partial match i += 1 and select highest match
        #   if regex match i += 1 and select all matches
        #   if not select highest key

        # "SXS:BBH:0123v4/Lev4/h_Extrapolated_N2.h5"
        # "SXS:BBH:0123v1/Lev4/rhOverM_Asymptotic_GeometricUnits_CoM.h5/Extrapolated_N2.dir"

        # return tuple of 3-tuples (full_match_string, h5group, URL)
        raise NotImplementedError()

        parts = location.split("/")
        if parts[0].lower().startswith(("catalog", "sxs")) and "v" in parts[0]:
            # The user has asked for a specific version
            files = self.files_with_all_versions
        else:
            files = self.files_with_latest_version

        matches = [""]
        part_index = 0
        for matching_index in range(len(file_parts)):
            if part_index >= len(parts):
                break
            part = parts[part_index]
            files_index = files[matching_index]
            if not part:
                matches = [sorted(files_index)[-1]]
                part_index += 1
            else:
                part_lower = part.lower()
                matches = [f for f in files_index if f.lower().startswith(part_lower)]
                if matches:
                    # We have matches of the raw string, so we just take the "largest"
                    matches = [sorted(matches)[-1]]
                else:
                    # We'll try to compile this part as a regex, and take *everything* that matches.
                    # Note that re.Pattern.match requires a match at the *beginning* of the string.
                    part = re.compile(part, re.IGNORECASE)
                    matches = [f for f in files_index if part.match(f)]
                if matches:
                    part_index += 1
                else:
                    matches = [sorted(files_index)[-1]]

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
