from __future__ import absolute_import, division, print_function

def read_catalog(catalog_root_directory='.', exclude_patterns=[r'^\.', r'^Attic', r'.*Links$'],
                 ignore_invalid_lines=False, suppress_errors=False,
                 cache_results=False, error_on_cache_failure=True,
                 indent=4, separators=(',', ': '), verbosity=1):
    """Parse the catalog's metadata into a single ordered dictionary

    Parameters
    ----------
    catalog_root_directory: str (default: '.')
        Relative or absolute path to the root of the directory to be searched for metadata files.
    exclude_patterns: list of str (default: [r'^\.', '^Attic', '.*Links$'])
        List of regex patterns which, if matched, will exclude a directory and all subdirectories
        from the search.
    ignore_invalid_lines: bool (default: False)
        This option is passed to `Metadata.from_txt_file`, when it is used.  If True, individual
        lines that are not well constructed are simply skipped, allowing remaining lines to be
        processed.  Otherwise, an Exception will be raised for that file, though it may be caught by
        the next argument.
    suppress_errors: bool (default: False)
        If True, any errors will be raised immediately, causing this function to halt; if False,
        warnings will be issued, but the function will continue.  Note that warnings may be filtered
        with the `warnings` python module.  This argument does not apply to errors related to
        writing the cache; see below.
    cache_results: bool (default: False)
        If True, cache the results to `catalog_metadata.json` in the `catalog_root_directory`.
    error_on_cache_failure: bool (default: True)
        If True, raise an error if the cache file cannot be made.
    indent: int (default: 4)
        Indentation level of the json file, if `cache_results` is True.  If zero or negative, only
        newlines will be added.  If `None`, no newlines will be used.
    separators: tuple of str (default: (',', ': '))
        Tuple of `(item_separator, key_separator)` to use in the json file, if `cache_results` is
        True.
    verbosity: int
        If greater than 0, print each simulation's directory as it is processed.

    Returns
    -------
    catalog_metadata: OrderedDict
        Collection of Metadata objects, with keys given by the directory (relative to the
        `catalog_root_directory`) in which the metadata file is found.

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
                # times if both are present, or will simply use the only one present.
                metadata = Metadata.from_file(os.path.join(root, 'metadata'), ignore_invalid_lines=ignore_invalid_lines)
                key = os.path.relpath(root, catalog_root_directory)
                if verbosity>0:
                    print(key)
                catalog_metadata[key] = metadata
        except Exception as e:
            if suppress_errors:
                message = "\nWhile trying to read metadata in {0}, an exception of type {1} occurred.\n".format(root, type(e).__name__)
                message += "\tException arguments: {0!r}\n".format(e.args)
                message += "Continuing to read other metadata files.\n"
                warnings.warn(message)
            else:
                message = "\nWhile trying to read metadata in {0}, an exception occurred:\n".format(root)
                warnings.warn(message)
                raise e

    if cache_results:
        try:
            with open(os.path.join(catalog_root_directory, 'catalog_metadata.json'), 'w') as f:
                json.dump(catalog_metadata, f, indent=indent, separators=separators)
        except Exception as e:
            if error_on_cache_failure:
                raise e

    return catalog_metadata


def drop_all_but_highest_levs(catalog):
    """Return a new catalog with all but the highest Levs removed

    Sorts the catalog into groups, and keeps only the simulation with the highest Lev.

    See also
    --------
    drop_all_but_selected_resolutions

    """
    import itertools
    highest_levs = [max(g) for k, g in itertools.groupby(sorted(catalog), key=lambda k: k.split('/Lev')[0])]
    return {k: catalog[k] for k in highest_levs}


def drop_all_but_selected_resolutions(catalog,
                                      group_and_resolution_parser=lambda key, val: (val.simulation_group, val.lev),
                                      resolution_selector=lambda resolution_list: max(resolution_list)):
    """Return a new catalog with all but the selected resolutions removed

    Parameters
    ----------
    catalog: dict-like
        Any mapping, but presumably an OrderedDict with keys given by some simulation name and
        values of Metadata objects.
    group_and_resolution_parser: function (default: lambda key, val: (val.simulation_group, val.lev))
        Function taking arguments of the key and value of each item in the `catalog`.  The return
        value should be a tuple with first element being the group into which the given item should
        sorted, and the second being the resolution assigned to the item.  This group name will be
        the key in the output catalog, and the resolution will be one of the inputs to the following
        function.  The default function assumes the input dict has values of Metadata objects that
        follow the standard SpEC naming conventions.
    resolution_selector: function (default: lambda resolution_list: max(resolution_list))
        Function taking a list of resolutions for a given group.  The output should be the single
        resolution selected to represent this group.  The default function simply sorts the list
        using python's built-in `sorted` function (which is alphabetical or numerical), and takes
        the last element.

    Returns
    -------
    catalog_metadata: OrderedDict
        Collection of objects that were the values of the input dict, with keys given by the first
        element returned by the `group_and_resolution_parser` input function.  If the input
        dict-like object was unordered, the order in this output will be meaningless.

    """
    import collections
    hierarchical_dict = collections.OrderedDict()
    for key in catalog:
        group, resolution = group_and_resolution_parser(key, catalog[key])
        hierarchical_dict.setdefault(group, dict())
        hierarchical_dict[group][resolution] = key

    keys_with_highest_resolutions = [hierarchical_dict[group][resolution_selector(list(hierarchical_dict[group]))]
                                     for group in hierarchical_dict]

    new_catalog = collections.OrderedDict([(key_with_highest_resolution, catalog[key_with_highest_resolution])
                                           for key_with_highest_resolution in keys_with_highest_resolutions])

    return new_catalog


def key_by_alternative_name(catalog, alternative_name_patterns=[r"""^SXS:""",],
                            error_on_duplicate_keys=True, error_on_missing_key=False,
                            warn_on_missing_key=False, allow_ugly_keys=False):
    """Return a new catalog, with keys replaced by first-matched alternative name found in metadata

    NOTE: You almost certainly want to run `drop_all_but_highest_levs` before running this function.
    If you don't, it is likely that multiple entries with different Levs will have the same
    alternative name, which will raise an exception by default.

    Parameters
    ----------
    catalog: OrderedDict
        Output from `drop_all_but_highest_levs`, for example, mapping a directory to a Metadata
        object.
    alternative_name_patterns: list of strings (default: ['^SXS:'])
        List of regular expressions to match in `alternative-names`.  If a match is found, the first
        matching name is used as the key in the output OrderedDict.
    error_on_duplicate_keys: bool (default: True)
        If False, only issue a warning when duplicate alternative-name keys are found; otherwise
        raise an Exception.  The last metadata found with that alternative name will be the only one
        that appears in the output.
    error_on_missing_key: bool (default: False)
        If True, raise an exception whenever no match is found in the list of alternative names.
    warn_on_missing_key: bool (default: False)
        If True, raise a warning whenever no match is found in the list of alternative names.  The
        code will not arrive at the point of the warning unless the previous argument was False.
    allow_ugly_keys: bool (default: False)
        If no match is found in the list of alternative names, simply use the ugly key.

    Returns
    -------
    catalog_metadata: OrderedDict
        Collection of Metadata objects, with keys given by the first match to the
        `alternative_name_patterns` found in the metadata's `alternative_names` field.  The values
        are copies of the input dictionary's values, and are given a new key named
        '_original_catalog_key'.

    """
    import re
    import warnings
    import copy
    import collections

    new_catalog = collections.OrderedDict()

    alternative_name_patterns = [re.compile(pattern) for pattern in alternative_name_patterns]

    def first_match(names, key):
        for pattern in alternative_name_patterns:
            for name in names:
                if pattern.match(name):
                    return name
        if allow_ugly_keys:
            return key
        return ''

    for key in catalog:
        alternative_names = catalog[key].alternative_names
        if not isinstance(alternative_names, list):
            alternative_names = [alternative_names,]
        new_key = first_match(alternative_names, key)
        if not new_key:
            message = ("\nCould not find a matching pattern among the alternative names for {0}:\n".format(key)
                       + "    alternative_names = {0!r}\n".format(alternative_names)
                       + "    alternative_name_patterns = {0!r}\n".format([p.pattern for p in alternative_name_patterns]))
            if error_on_missing_key:
                raise ValueError(message)
            elif warn_on_missing_key:
                warnings.warn(message)
        else:
            if new_key in new_catalog:
                message = "Name {0} is duplicated in input directories\n    {1}\nand\n    {2}"
                message = message.format(new_key, new_catalog[new_key]['_original_catalog_key'], key)
                if error_on_duplicate_keys:
                    raise ValueError(message)
                else:
                    warnings.warn(message)
            new_catalog[new_key] = copy.copy(catalog[key])
            try:
                new_catalog[new_key]['_original_catalog_key'] = key
            except Exception as e:
                pass  # Don't bother

    return new_catalog


def symlink_runs(source_directory='Catalog', target_directory='CatalogLinks',
                 remove_old_target_dir=False, alternative_name_patterns=[r"""^SXS:""",],
                 exclude_patterns=[r'^\.', r'^Attic', r'.*Links$'],
                 use_relative_links=False, relative_directory_path=None, verbosity=1):
    """Make nicely named symbolic links for entries in the SpEC waveform catalog

    Use this function to write a directory of symlinks with "nice" names like `SXS:BBH:0001` for
    entries in the `Catalog` directory, or names like `PRIVATE:BBH:0001` for entries in the
    `Incoming` directory.  These "nice" names are taken from any `metadata.json` or `metadata.txt`
    file.

    Parameters
    ----------
    source_directory: str (default: 'Catalog')
        Name of directory (either absolute or relative to working directory) to traverse looking for
        metadata files containing the `alternative-names` key with a matching a pattern (given below).
    target_directory: str (default: 'CatalogLinks')
        Name of directory (either absolute or relative to working directory) in which to store the
        resulting links.  If this does not exist, it is created.
    remove_old_target_dir: bool (default: False)
        If True, remove and recreate the `target_directory`, so that old links are gone.
    alternative_name_patterns: list of strings (default: ['^SXS:'])
        List of regular expressions to match in `alternative-names`.  If a match is found, the first
        matching name is used as the name of the "nice" link.
    exclude_patterns: list of str (default: [r'^\.', '^Attic', '.*Links$'])
        List of regex patterns which, if matched, will exclude a directory and all subdirectories
        from the search.
    use_relative_links: bool (default: False)
        If True, use relative paths in the links instead of absolute paths.  This can be helpful,
        for example, when the catalog is mounted as a volume in a docker container, so that the
        names can change.
    relative_directory_path: str (default: None)
        If this path is not None, use relative links and replace the relative part of the path
        between `source_directory` and `target_directory` with this string.  This forces
        `use_relative_links` to be True.  Again, this can be helpful when the source and target
        directories are mounted as volumes in a docker container, but under different names.
    verbosity: int 0, 1, or 2 (default: 1)
        Output only exceptions if verbosity is 0, only directories missing a matching
        alternative_name if verbosity is 1, and all directories if verbosity is 2.

    """
    import os
    import sys
    import shutil
    import re
    from . import _mkdir_recursively

    def symlink_force(target, link_name):
        import errno
        try:
            os.symlink(target, link_name)
        except OSError as e:
            if e.errno == errno.EEXIST:
                os.remove(link_name)
                os.symlink(target, link_name)
            else:
                raise e

    if not os.path.isdir(source_directory):
        raise ValueError("source_directory '{0}' not found".format(source_directory))
    if remove_old_target_dir and os.path.isdir(target_directory):
        shutil.rmtree(target_directory)
    _mkdir_recursively(target_directory)
    if relative_directory_path is not None:
        use_relative_links = True
        relative_directory_replaced = os.path.relpath(source_directory, target_directory)
    if verbosity < 1:
        verbosity = 0
    elif verbosity > 1:
        verbosity = 2

    catalog = read_catalog(source_directory, ignore_invalid_lines=True, suppress_errors=True,
                           exclude_patterns=exclude_patterns, verbosity=verbosity)
    catalog = drop_all_but_highest_levs(catalog)
    catalog = key_by_alternative_name(catalog, alternative_name_patterns=alternative_name_patterns,
                                      error_on_duplicate_keys=True, warn_on_missing_key=(verbosity>0))

    for name in catalog:
        source = os.path.join(source_directory, catalog[name]['_original_catalog_key'])
        target = os.path.join(target_directory, name)
        lev_index = source.rindex(os.sep + 'Lev')
        if lev_index > 0:
            source = source[:lev_index]

        if use_relative_links:
            source_link = os.path.relpath(source, os.path.dirname(target))
            if relative_directory_path is not None:
                source_link = source_link.replace(relative_directory_replaced, relative_directory_path)
        else:
            source_link = os.path.abspath(source)

        if verbosity > 1:
            print(target, '->', source_link)

        symlink_force(source_link, target)

    return catalog
