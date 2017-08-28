from __future__ import absolute_import, division, print_function


def create_web_files(catalog_root_directory='.', output_directory='web', relative_directory_path=None,
                     public_directory_patterns=[r'Catalog'], excluded_directory_patterns=[r'^\.', r'^Attic', r'.*Links$'],
                     public_altname_patterns=[r"""^SXS:""",], private_altname_patterns=[r"""^SXS:""", r"""^PRIVATE:""", r""".+"""]):
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
    import re
    import collections
    import json
    from . import (read_catalog, drop_all_but_highest_levs,
                   key_by_alternative_name, symlink_runs, metadata_field_mapping, _mkdir_recursively)

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
    public_dirs = [os.path.join(catalog_root_directory, d) for d in os.listdir(catalog_root_directory)
                   if os.path.isdir(os.path.join(catalog_root_directory, d)) and is_public(d) and not exclude(d)]
    private_dirs = [os.path.join(catalog_root_directory, d) for d in os.listdir(catalog_root_directory)
                    if os.path.isdir(os.path.join(catalog_root_directory, d)) and not is_public(d) and not exclude(d)]

    # Assemble the public parts of the catalog, symlinking for each directory
    public_catalog = collections.OrderedDict()
    for directory in public_dirs:
        sub_catalog = symlink_runs(source_directory=directory,
                                   target_directory=os.path.join(output_directory, 'public'),
                                   remove_old_target_dir=False, alternative_name_patterns=public_altname_patterns,
                                   exclude_patterns=excluded_directory_patterns,
                                   use_relative_links=True, relative_directory_path=os.path.join('..', '..', directory), verbosity=1)
        public_catalog.update(sub_catalog)
    
    private_catalog = collections.OrderedDict()
    for directory in private_dirs:
        sub_catalog = read_catalog(directory, ignore_invalid_lines=True, suppress_errors=True, exclude_patterns=excluded_directory_patterns)
        sub_catalog = drop_all_but_highest_levs(sub_catalog)
        sub_catalog = key_by_alternative_name(sub_catalog, alternative_name_patterns=private_altname_patterns,
                                              allow_ugly_keys=True)
        private_catalog.update(sub_catalog)

    # Get date of last git change
    try:
        on_windows = ('win' in platform.lower() and not 'darwin' in platform.lower())
        use_shell = not on_windows
        git_revision = subprocess.check_output("""git show -s --format="%ci" HEAD""", shell=use_shell).decode('ascii').rstrip()
        date, time, utc_offset = git_revision.split(' ')
        last_changed = date + ' ' + time
    except Exception as e:
        last_changed = ''

    # Make sure the directories we want are present
    _mkdir_recursively(os.path.join(output_directory, 'public'))
    _mkdir_recursively(os.path.join(output_directory, 'private'))

    # Finally, output the JSON files
    with open(os.path.join(output_directory, 'public', 'catalog_info.json'), 'w') as f:
        catalog_info = {
            'last_changed': last_changed,
            'metadata_fields': metadata_field_mapping
        }
        json.dump(catalog_info, f, indent=4)
    with open(os.path.join(output_directory, 'public', 'catalog.json'), 'w') as f:
        json.dump(public_catalog, f, indent=4)
    with open(os.path.join(output_directory, 'private', 'catalog.json'), 'w') as f:
        json.dump(private_catalog, f, indent=4)
