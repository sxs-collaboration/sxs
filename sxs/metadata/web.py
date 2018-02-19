from __future__ import absolute_import, division, print_function


def create_web_files(catalog_root_directory='.',
                     relative_directory_path=None,
                     public_json_directory='web/public',
                     public_file_listings_directory='web/public_file_listings',
                     public_links_directory='web/links',
                     private_json_directory='web/private',
                     public_directory_patterns=[r'Catalog'],
                     excluded_directory_patterns=[r'^\.', r'^Attic', r'.*Links$'],
                     public_altname_patterns=[r"""^SXS:""", r"""^BHNS:"""],
                     private_altname_patterns=[r"""^SXS:""", r"""^BHNS:""", r"""^PRIVATE:""", r""".+"""]):
    """Function to create the files needed by the website's catalog pages

    1) Create a directory symlinking public runs to data files
    2) Generate JSON for public runs
    3) Generate JSON for private runs
    4) Generate JSON describing metadata fields

    This function totally ignores hidden files, the `Attic` directory, and any directory ending in
    `Links`.  Only directories found within top-level directories matching any of the regexes given
    in `public_directory_patterns` will be put into the public dataset; all others will be put into
    the private dataset.

    """
    import os
    from sys import platform
    import subprocess
    import math
    import re
    import collections
    import json
    from . import (read_catalog, drop_all_but_highest_levs,
                   key_by_alternative_name, symlink_runs, metadata_fields, _mkdir_recursively)

    if relative_directory_path is None:
        relative_directory_path = os.path.join('..', '..')

    # Compile regex patterns and create functions for matching the various inputs
    public_directory_patterns = [re.compile(pattern) for pattern in public_directory_patterns]
    def is_public(directory):
        for pattern in public_directory_patterns:
            if pattern.match(directory):
                return True
        return False
    excluded_directory_patterns = [re.compile(pattern) for pattern in excluded_directory_patterns]
    def exclude(directory):
        for pattern in excluded_directory_patterns:
            if pattern.match(directory):
                return True
        return False
    public_altname_patterns = [re.compile(pattern) for pattern in public_altname_patterns]
    def first_public_match(names):
        for pattern in public_altname_patterns:
            for name in names:
                if pattern.match(name):
                    return name
        return ''
    private_altname_patterns = [re.compile(pattern) for pattern in private_altname_patterns]
    def first_private_match(names):
        for pattern in private_altname_patterns:
            for name in names:
                if pattern.match(name):
                    return name
        return ''

    # Figure out which directories are public and which are private
    catalog_root_directory = os.path.expanduser(catalog_root_directory)
    public_dirs = [os.path.normpath(os.path.join(catalog_root_directory, d))
                   for d in os.listdir(catalog_root_directory)
                   if os.path.isdir(os.path.join(catalog_root_directory, d)) and is_public(d) and not exclude(d)]
    private_dirs = [os.path.normpath(os.path.join(catalog_root_directory, d))
                    for d in os.listdir(catalog_root_directory)
                    if os.path.isdir(os.path.join(catalog_root_directory, d)) and not is_public(d) and not exclude(d)]

    # Assemble the public parts of the catalog, symlinking for each directory, and getting the catalogs
    public_catalog = collections.OrderedDict()
    for directory in public_dirs:
        replacement = os.path.join(relative_directory_path, os.path.relpath(directory, catalog_root_directory))
        sub_catalog = symlink_runs(source_directory=directory,
                                   target_directory=public_links_directory,
                                   remove_old_target_dir=False, alternative_name_patterns=public_altname_patterns,
                                   exclude_patterns=excluded_directory_patterns,
                                   use_relative_links=True, relative_directory_path=replacement, verbosity=1)
        public_catalog.update(sub_catalog)

    # Get the private catalogs
    private_catalog = collections.OrderedDict()
    for directory in private_dirs:
        sub_catalog = read_catalog(directory, ignore_invalid_lines=True, suppress_errors=True, exclude_patterns=excluded_directory_patterns)
        sub_catalog = drop_all_but_highest_levs(sub_catalog)
        sub_catalog = key_by_alternative_name(sub_catalog, alternative_name_patterns=private_altname_patterns,
                                              allow_ugly_keys=True)
        private_catalog.update(sub_catalog)

    # Rearrange the catalogs to be lists of OrderedDicts for JSON
    def modify_metadata(key, metadata):
        """Add 'name', 'object_types', and '*_mass_ratio' fields, and expand three-vectors to four separate fields"""
        m = [('name', key)]
        for k in metadata:
            v = metadata[k]
            if isinstance(v, list) and len(v)==3 and (
                    isinstance(v[0], float) and isinstance(v[1], float) and isinstance(v[2], float)):
                m += [(k+'_mag', math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)),
                      (k+'_x', v[0]), (k+'_y', v[1]), (k+'_z', v[2])]
            else:
                if k == 'object1':  # Insert combined type first
                    m += [('object_types', ''.join(sorted([metadata[k].upper(), metadata['object2'].upper()])))]
                if k == 'initial_mass1':  # Insert mass ratio first
                    m += [('initial_mass_ratio', metadata[k]/metadata['initial_mass2'])]
                if k == 'relaxed_mass1':  # Insert mass ratio first
                    m += [('relaxed_mass_ratio', metadata[k]/metadata['relaxed_mass2'])]
                m += [(k, v)]
        return collections.OrderedDict(m)
    public_catalog = [modify_metadata(key, metadata)
                      for key in public_catalog for metadata in [public_catalog[key]]]
    private_catalog = [modify_metadata(key, metadata)
                       for key in private_catalog for metadata in [private_catalog[key]]]

    # public_catalog = [collections.OrderedDict([('name', key)] + [(k, val[k]) for k in val])
    #                   for key in public_catalog for val in [public_catalog[key]]]
    # private_catalog = [collections.OrderedDict([('name', key)] + [(k, val[k]) for k in val])
    #                   for key in private_catalog for val in [private_catalog[key]]]

    # Get date of last git change
    try:
        on_windows = ('win' in platform.lower() and not 'darwin' in platform.lower())
        use_shell = not on_windows
        git_revision = subprocess.check_output("""git show -s --format="%ci" HEAD""",
                                               cwd=catalog_root_directory,
                                               shell=use_shell).decode('ascii').rstrip()
        date, time, utc_offset = git_revision.split(' ')
        last_changed = date + ' ' + time
    except Exception as e:
        last_changed = ''

    # Make sure the directories we want are present
    _mkdir_recursively(public_json_directory)
    _mkdir_recursively(private_json_directory)

    # Finally, output the JSON files
    with open(os.path.join(public_json_directory, 'catalog_info.json'), 'w') as f:
        catalog_info = {
            'title': 'SXS Gravitational Waveform Database',
            'subtitle': 'Completed Simulations',
            'documentation': True,
            'news': True,
            'last_changed': last_changed,
            'fields': metadata_fields
        }
        json.dump(catalog_info, f, indent=4)
    with open(os.path.join(public_json_directory, 'catalog.json'), 'w') as f:
        json.dump(public_catalog, f, indent=4)
    write_file_listings(public_catalog, public_file_listings_directory, public_links_directory)
    with open(os.path.join(private_json_directory, 'catalog.json'), 'w') as f:
        json.dump(private_catalog, f, indent=4)


def write_file_listings(catalog, file_listings_directory, catalog_root_directory,
                        file_extension_whitelist = ['.h5', '.txt', '.out', '.perl', '.tgz']):
    import os
    import os.path
    import json
    from . import _mkdir_recursively

    def lsdir(path):
        path = os.path.realpath(path)
        for entry in os.listdir(path):
            if (not entry.startswith('.')
                and (os.path.isdir(os.path.join(path, entry))
                     or os.path.splitext(entry)[1].lower() in file_extension_whitelist)):
                yield entry

    def path_to_dict(annex_root, path, depth=-1):
        d = {
            "basename": os.path.basename(path),
            "depth": depth,
            "path": path,
        }
        try:
            if os.path.isdir(os.path.join(annex_root, path)):
                children = [y
                            for x in lsdir(os.path.join(annex_root, path))
                            for y in [path_to_dict(annex_root, os.path.join(path, x), depth+1)]
                            if y is not None]
                if len(children) == 0:
                    return None
                d['type'] = "directory"
                d['children'] = sorted(children, key=lambda d: d['basename'].lower())
            else:
                split_extension = os.path.splitext(d['basename'])
                if len(split_extension) < 1:
                    return None
                d['type'] = split_extension[1].lower()
                d['size'] = os.stat(os.path.join(annex_root, path)).st_size
            return d
        except:
            return None

    _mkdir_recursively(file_listings_directory)

    for group in catalog:
        name = group['name']
        d = path_to_dict(catalog_root_directory, name)
        if d is None:
            continue
        d['basename'] = name
        #print(d['basename'])
        with open(os.path.join(file_listings_directory, name + '.json'), 'w') as f:
            json.dump(d, f)

