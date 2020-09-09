class Catalog(object):
    def __init__(self, *args, **kwargs):
        raise NotImplementedError()

    def select(self, path_pattern):
        """Select from all catalog files by progressively matching path components

        Parameters
        ----------
        path_pattern : str
            A pattern to search for among the catalog files.  This is first searched as
            a literal substring, which will return a set of one or more matches.  If no
            matches were found, the pattern is split into components by "/", which are
            used progressively to match corresponding path components — first as
            literal substrings, then as (python-style) regular expressions.  Each
            partial match in this step will pass all the matched components to the
            `max` function, and choose only the result — though with regexes, there may
            be distinct partial matches, each of which will pass only its corresponding
            components to `max`, resulting in multiple matches.

        Returns
        -------
        matched_paths : Set[str]

        See Also
        --------
        sxs.utilities.select_by_path_component

        Examples
        --------
        First, we can choose the `h` waveform with `n=2` extrapolation in the
        highest-resolution (Lev) run from the simulation SXS:BBH:0002 with

            >>> catalog.select("SXS:BBH:0002/Lev/h_ex.*_n2")
            {"SXS:BBH:0002v7/Lev6/h_extrapolated_n2.h5"}

        Because the "Lev" component of the input string only matched the "Lev" portion
        of the path components "Lev4", "Lev5", and "Lev6", all of those components were
        passed to `max` together, and only "Lev6" was returned.  Similarly, the highest
        version of SXS:BBH:0002 (currently, v7) was chosen automatically.

        We could, instead, use a regular expression to include the numbers in the
        matches:

            >>> catalog.select("SXS:BBH:0002/Lev./h_ex.*_n2")
            {"SXS:BBH:0002v7/Lev4/h_extrapolated_n2.h5",
             "SXS:BBH:0002v7/Lev5/h_extrapolated_n2.h5",
             "SXS:BBH:0002v7/Lev6/h_extrapolated_n2.h5"}

        In this case, the pattern "Lev." matches each component entirely, so they are
        each retained.  Similarly, "SXS:BBH:0002v." would match all versions, and
        "h_ex.*_n." would match all extrapolation orders — and "h.*" would match all
        extrapolation orders, as well as the outermost extraction.

        Another helpful pattern is to use alternation in the regular expression, to
        explicitly list acceptable matches:

            >>> catalog.select("SXS:BBH:0002/Lev(4|6)/h_ex.*_n2")
            {"SXS:BBH:0002v7/Lev4/h_extrapolated_n2.h5",
             "SXS:BBH:0002v7/Lev6/h_extrapolated_n2.h5"}

        """
        from ..utilities import select_by_path_component
        files = self.files
        selections = select_by_path_component(location, set(self.files))

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
