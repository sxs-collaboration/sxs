from __future__ import absolute_import, division, print_function

def read_catalog(catalog_root_directory='.', exclude_patterns=['Attic', '.*Links'],
                 ignore_invalid_lines=False, suppress_errors=False,
                 cache_results=False, error_on_cache_failure=True, indent=4, separators=(',', ': ')):
    """Parse the catalog's metadata into a single ordered dictionary

    Parameters
    ----------
    catalog_root_directory: str, default='.'
        Relative or absolute path to the root of the directory to be searched for metadata files.
    exclude_patterns: list of str, default=['Attic', '.*Links']
        List of regex patterns which, if matched, will exclude a directory and all subdirectories
        from the search.
    ignore_invalid_lines: bool, default=False
        This option is passed to `Metadata.from_txt_file`, when it is used.  If True, individual
        lines that are not well constructed are simply skipped, allowing remaining lines to be
        processed.  Otherwise, an Exception will be raised for that file, though it may be caught by
        the next argument.
    suppress_errors: bool, default=False
        If True, any errors will be raised immediately, causing this function to halt; if False,
        warnings will be issued, but the function will continue.  Note that warnings may be filtered
        with the `warnings` python module.  This argument does not apply to errors related to
        writing the cache; see below.
    cache_results: bool, default=False
        If True, cache the results to `catalog_metadata.json` in the `catalog_root_directory`.
    error_on_cache_failure: bool, default=True
        If True, raise an error if the cache file cannot be made.
    indent: int, default=4
        Indentation level of the json file, if `cache_results` is True.  If zero or negative, only
        newlines will be added.  If `None`, no newlines will be used.
    separators: tuple of str, default=(',', ': ')
        Tuple of `(item_separator, key_separator)` to use in the json file, if `cache_results` is
        True.

    Returns
    -------
    catalog_metadata: OrderedDict
        Collection of Metadata objects arranged by `simulation_group` and `lev` numbers of the
        Metadata objects.  For example, to access the simulation `d16_q1_s0_s0/Lev6`, use
        `catalog_metadata['d16_q1_s0_s0']['6']['metadata']`.  The lowest-level dictionary also
        stores the directory in which that metadata was found:
        `catalog_metadata['d16_q1_s0_s0']['6']['directory']`.

    """
    import os
    import re
    import warnings
    import json
    from collections import OrderedDict
    from . import Metadata

    exclude_patterns = [re.compile(pattern) for pattern in exclude_patterns]

    def exclude(directory):
        for pattern in exclude_patterns:
            if pattern.match(directory):
                return True
        return False

    catalog_metadata = OrderedDict()

    for root, dirs, files in os.walk(catalog_root_directory, topdown=True):
        try:
            # From help(os.walk): "When topdown is true, the caller can modify the dirnames list
            # in-place (e.g., via del or slice assignment), and walk will only recurse into the
            # subdirectories whose names remain in dirnames; this can be used to prune the search"
            dirs[:] = [d for d in dirs if not exclude(d)]

            if 'metadata.txt' in files or 'metadata.json' in files:
                # Either metadata.txt or metadata.json (or both) is acceptable; the
                # `Metadata.from_file` function will decide which one to use based on modification
                # times.
                metadata = Metadata.from_file(os.path.join(root, 'metadata'), ignore_invalid_lines=ignore_invalid_lines)
                simulation_group = metadata.simulation_group
                lev = metadata.lev
                catalog_metadata.setdefault(simulation_group, OrderedDict())
                catalog_metadata[simulation_group][str(lev)] = {'directory': root, 'metadata': metadata}
        except Exception as e:
            if suppress_errors:
                message = "\nWarning: While trying to read metadata in {0}, an exception of type {1} occurred.\n".format(root, type(e).__name__)
                message += "\tException arguments: {0!r}\n".format(e.args)
                message += "Continuing to read other metadata files."
                warnings.warn(message)
            else:
                raise e

    if cache_results:
        try:
            with open(os.path.join(catalog_root_directory, 'catalog_metadata.json'), 'w') as f:
                json.dump(catalog_metadata, f, indent=indent, separators=separators)
        except Exception as e:
            if error_on_cache_failure:
                raise e

    return catalog_metadata


def drop_all_but_highest_levs(catalog, in_place=False):
    new_catalog = OrderedDict()

    for key in catalog:
        levs = list(catalog[key])
        highest_lev = sorted(levs)[-1]
        if in_place:
            for lev in levs:
                if lev != highest_lev:
                    del catalog[key][lev]
        else:
            new_catalog[key] = {
                highest_lev: {
                    'directory': catalog[key][highest_lev]['directory'],
                    'metadata': catalog[key][highest_lev]['metadata']
                }
            }

    if in_place:
        return catalog
    else:
        return new_catalog


def key_by_alternative_name(catalog_metadata, alternative_name_patterns=[r"""^SXS:""",]):
    import re

    new_catalog = {}

    alternative_name_patterns = [re.compile(pattern) for pattern in alternative_name_patterns]

    def first_match(names):
        for pattern in alternative_name_patterns:
            for name in names:
                if pattern.match(name):
                    return name
        return ''

    for key in catalog:
        levs = list(catalog[key])
        highest_lev = sorted(levs)[-1]
        alternative_names = catalog[key][highest_lev]['metadata'].alternative_names
        if not isinstance(alternative_names, list):
            alternative_names = [alternative_names,]
        new_key = first_match(alternative_names)
        if new_key:
            new_catalog[new_key] = {
                'directory': catalog[key][highest_lev]['directory'],
                'metadata': catalog[key][highest_lev]['metadata']
            }

    return collections.OrderedDict([(key, new_catalog[key]) for key in sorted(new_catalog)])
