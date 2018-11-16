
catalog_file_description = """
        This JSON file has the following format.  Comments are, of course, not present (since JSON does not support
        comments).  Single quotes here are, of course, double quotes in the rest of the file (since JSON encloses
        strings in double quotes).  Anything inside <angle brackets> is just standing in for the relevant value.  An
        ellipsis ... indicates that the preceding structure can be repeated.  Also note that the metadata entries for
        simulations may not be present if the record on zenodo is closed-access; see catalog_private_metadata.json if
        you have access to those simulations, which will contain the missing information.  That file should be read
        and written automatically by functions in this module, so that the catalog dict returned will contain all
        available information.

        {
            'catalog_file_description': '<this description>',
            'modified': '<YYYY-MM-DDThh:mm:ss.ssssss>',  # UTC time of last-modified record in this file
            'records': {  # Includes *all* versions of *all* records published in the 'sxs' community, not just simulations
                '<id>': {  # This Zenodo ID key is a *string* containing the 'id' value inside this object (JSON requires keys to be strings)
                    # More details about this 'representation' object at http://developers.zenodo.org/#depositions
                    'conceptdoi': '10.5281/zenodo.<conceptrecid>',  # Permanent DOI for all versions of this record
                    'conceptrecid': '<conceptrecid>',  # ~7-digit integer (as string) collectively identifying all versions of this record
                    'created': '<YYYY-MM-DDThh:mm:ss.ssssss>',  # UTC time of creation of this record on Zenodo
                    'doi': '10.5281/zenodo.<id>',  # Permanent DOI for this record
                    'doi_url': 'https://doi.org/10.5281/zenodo.<id>',  # URL for permanent DOI of this record
                    'id': <id>,  # ~7-digit integer uniquely identifying this record
                    'links': {
                         'badge': 'https://zenodo.org/badge/doi/10.5281/zenodo.<id>.svg',
                         'bucket': 'https://zenodo.org/api/files/<uuid>',  # Base URL for file uploads and downloads
                         'conceptbadge': 'https://zenodo.org/badge/doi/10.5281/zenodo.<conceptrecid>.svg',
                         'conceptdoi': 'https://doi.org/10.5281/zenodo.<conceptrecid>',  # Permanent link to webpage for most-recent version
                         'discard': 'https://zenodo.org/api/deposit/depositions/<id>/actions/discard',  # API action to discard a draft
                         'doi': 'https://doi.org/10.5281/zenodo.<id>',  # Permanent URL for this version
                         'edit': 'https://zenodo.org/api/deposit/depositions/<id>/actions/edit',  # API action to edit this record
                         'files': 'https://zenodo.org/api/deposit/depositions/<id>/files',  # Only present for author
                         'html': 'https://zenodo.org/deposit/<id>',  # Webpage for this version
                         'latest': 'https://zenodo.org/api/records/<id>',  # API endpoint for most-recent version
                         'latest_html': 'https://zenodo.org/record/<id>',  # Webpage for most-recent version
                         'publish': 'https://zenodo.org/api/deposit/depositions/<id>/actions/publish',  # Only present for author
                         'record': 'https://zenodo.org/api/records/<id>',  # Only present for author
                         'record_html': 'https://zenodo.org/record/<id>',  # Webpage for this particular version; only present for author
                         'self': 'https://zenodo.org/api/deposit/depositions/<id>'
                    },
                    'metadata': {  # Note that this is Zenodo metadata, and is different from the SXS metadata found below
                        'access_right': '<access>',  # Can be 'open', 'closed', 'embargoed', or 'restricted'
                        'communities': [
                            {'identifier': '<community_name>'},  # Names may include 'sxs' and 'zenodo'
                            ...
                        ],
                        'creators': [
                            {
                                'name': '<name>',  # Name of this creator in the format Family name, Given names
                                'affiliation': '<affiliation>',  # (Optional) Affiliation of this creator
                                'orcid': '<orcid>',  # (Optional) ORCID identifier of this creator
                                'gnd': '<gnd>'  # (Optional) GND identifier of this creator
                            },
                            ...
                        ],
                        'description': '<description>',  # Text description of this record
                        'doi': '10.5281/zenodo.<id>',  # Permanent DOI of this record
                        'keywords': [
                            '<keyword>',  # Optional; this array may be empty
                            ...
                        ],
                        'license': '<license_type>',  # Usually 'CC-BY-4.0' for SXS
                        'prereserve_doi': {'doi': '10.5281/zenodo.<id>', 'recid': <id>},
                        'publication_date': '<YYYY-MM-DD>',  # Possibly meaningless date (UTC)
                        'title': '<title>',
                        'upload_type': 'dataset'
                    },
                    'modified': '<YYYY-MM-DDThh:mm:ss.ssssss>',  # (UTC) Last modification of this record (possibly just Zenodo metadata modified)
                    'owner': <user_id>,  # ~5-digit integer identifying the user who owns this record
                    'record_id': <id>,  # Same as 'id'
                    'state': '<state>',  # Can be 'done', 'inprogress', 'error', 'unsubmitted', possibly others
                    'submitted': <submitted>,  # True or false (always true for published records)
                    'title': '<title>'  # Same as ['metadata']['title'],
                    'files': [  # May not be present if this simulation is closed-access; see catalog_private_metadata.json as noted above
                        # See https://data.black-holes.org/waveforms/documentation.html for
                        # detailed descriptions of the *contents* of the files in each record.
                        {
                            'checksum': '<checksum>',  # MD5 checksum of file on Zenodo
                            'filename': '<filename>',  # Name of file; may contain slashes denoting directories
                            'filesize': <filesize>,  # Number of bytes in the file
                            'id': '<fileid>',  # A standard UUID (hexadecimal with characters in the pattern 8-4-4-4-12)
                            'links': {
                                'download': 'https://zenodo.org/api/files/<bucket>/<filename>',  # The URL to use to download this file
                                'self': 'https://zenodo.org/api/deposit/depositions/<deposition_id>/files/<fileid>'  # Ignore this
                            }
                        },
                        ...  # Other file descriptions in the order in which they were uploaded (not necessarily a meaningful order)
                    ]
                },
                ...
            },
            'simulations': {
                '<sxs_id>': {  # The SXS ID is a string like SXS:BHNS:0001 or SXS:BBH:1234
                    'conceptrecid': '<conceptrecid>',  # The Zenodo ID of the 'concept' record, which *resolves to* the most-recent version
                    'versions': [  # Zenodo IDs of each version.  There may only be one; index this list with [-1] to always get the most recent
                        '<id1>',  # Oldest version first
                        ...
                    ],
                    'metadata': {  # May not be present if this simulation is closed-access; see catalog_private_metadata.json as noted above
                        # Variable content describing (mostly) physical parameters of the system.  It's basically a
                        # python-compatible version of the information contained in 'metadata.txt' from the
                        # highest-resolution run in the most-recent version of this simulation.  That file is meant to
                        # be more-or-less as suggested in <https://arxiv.org/abs/0709.0093>.  The conversion to a
                        # python-compatible format means that keys like 'simulation-name' have had hyphens replaced by
                        # underscores so that they can be used as variable names in python and any other sane language
                        # (with apologies to Lisp).  As far as possible, values that are just strings in that file
                        # have been converted into the relevant types -- like numbers, integers, and arrays.  Note
                        # that some keys like eccentricity are sometimes numbers and sometimes the string '<number'
                        # (meaning that the eccentricity is less than the number), which is necessarily a string.
                        #
                        # Below are just the first few keys that *may* be present.  Note that closed-access
                        # simulations will have empty dictionaries here.
                        #
                        'simulation_name': '<directory_name>',  # This may be distinctly uninformative
                        'alternative_names': '<sxs_id>',  # This may be a list of strings
                        'initial_data_type': '<type>',  # Something like 'BBH_CFMS'
                        'number_of_orbits': <number>,  # This is a float
                        'relaxed_mass1': <m2>,
                        'relaxed_mass2': <m1>,
                        'relaxed_dimensionless_spin1': [
                            <chi1_x>,
                            <chi1_y>,
                            <chi1_z>
                        ],
                        'relaxed_dimensionless_spin2': [
                            <chi2_x>,
                            <chi2_y>,
                            <chi2_z>
                        ],
                        'relaxed_eccentricity': <eccentricity>,  # Or maybe a string...
                        'relaxed_orbital_frequency': [
                            <omega_x>,
                            <omega_y>,
                            <omega_z>
                        ],
                        'relaxed_measurement_time': <time>,
                        ...
                    },
                },
                ...
            }
        }
"""


def split_to_public_and_private(catalog):
    from collections import OrderedDict
    from copy import deepcopy
    public = deepcopy(catalog)
    private = {'records': OrderedDict(), 'simulations': OrderedDict()}
    for record_id in catalog['records']:
        is_public = (catalog['records'][record_id]['metadata']['access_right'] == 'open')
        if not is_public and 'files' in catalog['records'][record_id]:
            private['records'][record_id] = {'files': deepcopy(catalog['records'][record_id]['files'])}
            public['records'][record_id].pop('files', None)
    for sxs_id in catalog['simulations']:
        version = str(catalog['simulations'][sxs_id]['versions'][-1])
        record = catalog['records'][version]
        is_public = (record['metadata']['access_right'] == 'open')
        if not is_public and 'metadata' in catalog['simulations'][sxs_id]:
            private['simulations'][sxs_id] = {'metadata': deepcopy(catalog['simulations'][sxs_id]['metadata'])}
            public['simulations'][sxs_id].pop('metadata', None)
    return public, private


def join_public_and_private(public, private):
    from copy import deepcopy
    catalog = deepcopy(public)
    for record_id in private['records']:
        if record_id in public['records']:
            catalog['records'][record_id]['files'] = deepcopy(private['records'][record_id]['files'])
    for sxs_id in private['simulations']:
        if sxs_id in public['simulations']:
            catalog['simulations'][sxs_id]['metadata'] = deepcopy(private['simulations'][sxs_id]['metadata'])
    return catalog


def update(path='~/.sxs/catalog/catalog.json', verbosity=1):
    """Update a local copy of the SXS catalog

    Because git has better handling of revision history with incremental updates, and because most
    users will have set up their credentials for github, we prefer git to simply re-downloading the
    catalog from black-holes.  Specifically, we try to update the catalog in the following order.

    1) Private copy via github (git scheme)
    2) Private copy via github (https scheme)
    3) Public copy via direct download (https://data.black-holes.org/catalog.json)

    Parameters
    ==========
    path: str, defaults to '~/.sxs/catalog/catalog.json'
        Absolute or relative path to JSON file containing the catalog.  If the path does not end
        with '.json', it is assumed to be a directory containing a 'catalog.json' file.
    verbosity: int, defaults to 1

        Amount of information to output.  Less than 1 corresponds to no output; 1 to only print a
        notice if the private file cannot be retrieved; greater than 1 to print a notice about
        wherever the file is retrieved; greater than 2 shows the stdout/stderr from external calls;
        and greater than 3 also asks git to be verbose.

    """
    from os.path import expanduser, isdir, join, dirname, basename, exists
    from os import makedirs, chdir, getcwd, remove
    from shutil import copyfile
    from subprocess import call, check_call
    from warnings import warn
    from .api.utilities import download
    path = expanduser(path)
    if not path.endswith('.json'):
        path = join(path, 'catalog.json')
    directory = dirname(path)
    if not exists(directory):
        makedirs(directory)
    chdir(directory)
    if verbosity > 2:
        stdout = None
        stderr = None
    else:
        stdout = subprocess.DEVNULL
        stderr = subprocess.DEVNULL
    if verbosity > 3:
        git_verbosity = '-v'
    else:
        git_verbosity = ''
    if exists(path):
        copyfile(path, path+'.bak')
    try:
        git_success = False
        try:
            if call("git status {0} .".format(git_verbosity), shell=True, stdout=stdout, stderr=stderr):
                check_call("git init {0} .".format(git_verbosity), shell=True, stdout=stdout, stderr=stderr)
            call("git remote {0} add origin_git git@github.com:sxs-collaboration/zenodo_catalog.git".format(git_verbosity),
                 shell=True, stdout=stdout, stderr=stderr)
            call("git remote {0} add origin_https https://github.com/sxs-collaboration/zenodo_catalog.git".format(git_verbosity),
                 shell=True, stdout=stdout, stderr=stderr)
            for remote in ["origin_git", "origin_https"]:
                if not call("git pull {0} {1} master".format(git_verbosity, remote), shell=True, stdout=stdout, stderr=stderr):
                    call("git reset {0} --hard HEAD".format(git_verbosity), shell=True, stdout=stdout, stderr=stderr)
                    git_success = True
                    if verbosity>1:
                        print('Retrieved catalog from {0}.'.format(remote))
                    break
        except:  # If for *any* reason git failed...
            pass
        if not git_success:  # ...fall back to direct download
            print("Failed to pull private copy of catalog; downloading public version.")
            if verbosity >= 1:
                verbose = 1
            download('https://data.black-holes.org/catalog.json', basename(path), verbose)
    except:  # If for *any* reason that failed...
        if exists(path+'.bak'):  # ... move the original file (if it exists) back into place
            rename(path+'.bak', path)
        raise
    else:  # If everything went well...
        if exists(path+'.bak'):  # ... remove the backup
            remove(path+'.bak')


def read(path=None):
    from os.path import expanduser, exists, join, dirname
    from json import load
    import sxs
    if path is None:
        if exists('catalog.json'):
            path = 'catalog.json'
        elif exists(expanduser('~/.sxs/catalog/catalog.json')):
            path = expanduser('~/.sxs/catalog/catalog.json')
        else:
            raise ValueError("Cannot find 'catalog.json' file in current directory or ~/.sxs/catalog/catalog.json.")
    else:
        path = expanduser(path)
    with open(path, 'r') as f:
        catalog = load(f)
    return catalog


def write(catalog, catalog_file_name=None, private_metadata_file_name=None):
    """Write catalog dictionary to file
    
    This function separates the catalog into a public part and any private SXS metadata.  If the
    latter exists, it gets written to the file given as the third parameter, or simply the file
    'catalog_private_metadata.json' in the same directory as the public part.

    Parameters
    ==========
    catalog: dict
        The catalog information in the format described by the string
        `sxs.zenodo.catalog.catalog_file_description`.
    catalog_file_name: str or None
        Path to the output public JSON file describing this catalog.  If None, the file is
        'catalog.json' in the working directory.  If the string is precisely
        'sxs/zenodo/catalog.json', the file will be placed in the sxs module's path, which is
        typically in some directory like .../lib/python3.x/site-packages/sxs/zenodo.
    private_metadata_file_name: str or None
        Path to the output private JSON file describing any private metadata.  If None, the file
        will be placed alongside the 'catalog.json' file, and named 'catalog_private_metadata.json'.
        Note that this file will not be written at all if there are no private metadata sets.

    """
    from os.path import join, dirname
    from json import dump
    import sxs
    if catalog_file_name == 'sxs/zenodo/catalog.json':
        catalog_file_name = join(dirname(sxs.__file__), 'zenodo', 'catalog.json')
    elif catalog_file_name is None:
        catalog_file_name = 'catalog.json'
    if private_metadata_file_name is None:
        private_metadata_file_name = join(dirname(catalog_file_name), 'catalog_private_metadata.json')
    public, private = split_to_public_and_private(catalog)
    with open(catalog_file_name, 'w') as f:
        dump(public, f, indent=4, separators=(',', ': '), ensure_ascii=True)
    if private['simulations']:
        with open(private_metadata_file_name, 'w') as f:
            dump(private, f, indent=4, separators=(',', ': '), ensure_ascii=True)


def modification_time(representation_list):
    return sorted([r['modified'] for r in representation_list])[-1]


def sxs_metadata_file_description(representation):
    """Find metadata file from highest Lev for this simulation"""
    from os.path import basename
    files = representation['files']
    metadata_files = [f for f in files if basename(f['filename'])=='metadata.json']
    metadata_files = sorted(metadata_files, key=lambda f: f['filename'])
    if not metadata_files:
        return None
    return metadata_files[-1]


def fetch_metadata(url, login, *args, **kwargs):
    """Get the json file from zenodo"""
    from .api import Login
    login = login or Login(*args, **kwargs)
    r = login.session.get(url)
    if r.status_code != 200:
        return {}
    try:
        return r.json()
    except:
        return {}


def order_version_list(representation_dict, versions):
    return sorted([str(v) for v in versions], key=lambda v: representation_dict[v]['created'])


def update_simulations(catalog, representation_list, login=None, *args, **kwargs):
    """Update list of simulations (and SXS metadata) from list of representations

    This function can be used to refresh the information about simulations found in the catalog from
    a list of "representation" objects as returned by zenodo.  Thus, if you search for new records
    on zenodo, and get that list back, you can simply run this function on the catalog and that
    list, and it will update all the information about simulations.  Note that this function returns
    the updated 'simulations' dictionary, rather than modifying it in place.

    Parameters
    ==========
    catalog: dict
        The catalog information in the format described by the string
        `sxs.zenodo.catalog.catalog_file_description`.
    representation_list: list
        List of "representation" dictionaries as returned by Zenodo.

    """
    from copy import deepcopy
    import re
    from collections import OrderedDict
    # from .. import sxs_identifier_regex
    from sxs import sxs_identifier_regex
    sxs_identifier_regex = re.compile(sxs_identifier_regex)
    simulations = deepcopy(catalog['simulations'])
    verbosity = kwargs.pop('verbosity', 2)
    for i, r in enumerate(representation_list, 1):
        print('{0:6} of {1}: {2}'.format(i, len(representation_list), r['id']))
        sxs_id_match = sxs_identifier_regex.search(r['title'])
        if sxs_id_match:
            sxs_id = sxs_id_match['sxs_identifier']
            zenodo_id = r['id']
            conceptrecid = r['conceptrecid']
            if sxs_id in simulations:
                if zenodo_id not in simulations[sxs_id]['versions']:
                    # First, get the information for the current most-recent metadata
                    last_record_id = str(simulations[sxs_id]['versions'][-1])
                    last_representation = catalog['records'][last_record_id]
                    last_metadata_file_description = sxs_metadata_file_description(last_representation)
                    # Now, create the new sorted version list
                    simulations[sxs_id]['versions'] = order_version_list(catalog['records'], simulations[sxs_id]['versions'] + [zenodo_id])
                    new_last_record_id = str(simulations[sxs_id]['versions'][-1])
                    if last_record_id != new_last_record_id:
                        # Only if the most-recent Zenodo ID has changed do we need to check any more
                        new_last_representation = catalog['records'][new_last_record_id]
                        new_metadata_file_description = sxs_metadata_file_description(new_last_representation)
                        if last_metadata_file_description['checksum'] != new_last_metadata_file_description['checksum']:
                            # Only if the metadata checksum has changed do we need to change the metadata
                            url = new_last_metadata_file_description['links']['download']
                            simulations[sxs_id]['metadata'] = fetch_metadata(url, login, *args, **kwargs)
                elif not simulations[sxs_id]['metadata']:
                    # Try to download the most-recent metadata; nothing else has changed.
                    last_record_id = str(simulations[sxs_id]['versions'][-1])
                    last_representation = catalog['records'][last_record_id]
                    last_metadata_file = sxs_metadata_file_description(last_representation)
                    url = last_metadata_file['links']['download']
                    simulations[sxs_id]['metadata'] = fetch_metadata(url, login, *args, **kwargs)
                else:
                    pass  # Nothing for 'simulations' changed.  Something might have changed for 'records', though.
            else:
                metadata_file_description = sxs_metadata_file_description(r)
                url = metadata_file_description['links']['download']
                simulations[sxs_id] = {
                    'conceptrecid': conceptrecid,
                    'versions': [zenodo_id],
                    'metadata': fetch_metadata(url, login, *args, **kwargs)
                }
    return OrderedDict([(s, simulations[s]) for s in sorted(simulations)])


def catalog_from_representation_list(representation_list, simulation_dict={}, login=None, *args, **kwargs):
    """Convert list of representations from Zenodo into catalog dictionary

    Given a list of "representation" dictionaries as returned by Zenodo, this function returns a "catalog" dictionary,
    as described by `catalog_file_description`.  In brief, this dictionary contains the `catalog_file_description`
    itself under that key, the UTC timestamp of the last modification of any item in the list, and then a series of keys
    given by the SXS identifiers of all the simulations in the input list, values for which are just the representations
    themselves.

    """
    from copy import deepcopy
    from collections import OrderedDict
    from textwrap import dedent
    catalog = OrderedDict()
    catalog['catalog_file_description'] = dedent(catalog_file_description).split('\n')[1:-1]
    catalog['modified'] = modification_time(representation_list)
    catalog['records'] = OrderedDict(sorted([(str(r['id']), r) for r in representation_list], key=lambda kv:kv[0]))
    catalog['simulations'] = deepcopy(simulation_dict)
    catalog['simulations'] = update_simulations(catalog, representation_list=representation_list, login=login, *args, **kwargs)
    return catalog



# def compare_catalogs(c1, c2):
#     from copy import deepcopy
#     # First, copy the catalogs so we can modify them without screwing other things up
#     c1 = deepcopy(c1)
#     c2 = deepcopy(c2)
#     # Remove keys that are allowed to differ
#     for c in [c1, c2]:
#         del c['catalog_file_description']
#         del c['modified']
#         for sxs_id in c:
#             if 'versions' in c[sxs_id]:
#                 for i in range(len(c[sxs_id]['versions'])):
#                     files = c[sxs_id]['versions'][i]['representation']['files']
#                     for j in range(len(files)):
#                         del files[j]['links']['self']
#     # Now we're ready to compare
#     return c1 == c2


# def catalog_from_representation_list(representation_list):
#     """Convert list of representations from Zenodo into catalog dictionary

#     Given a list of "representation" dictionaries as returned by Zenodo, this function returns a "catalog" dictionary,
#     as described by `catalog_file_description`.  In brief, this dictionary contains the `catalog_file_description`
#     itself under that key, the UTC timestamp of the last modification of any item in the list, and then a series of keys
#     given by the SXS identifiers of all the simulations in the input list, values for which are just the representations
#     themselves.

#     """
#     import re
#     from collections import OrderedDict, defaultdict
#     from sxs import sxs_identifier_regex
#     sxs_identifier_regex = re.compile(sxs_identifier_regex)
#     repr_by_sxs_id = defaultdict(list)
#     for r in representation_list:
#         m = sxs_identifier_regex.search(r['title'])
#         if m:
#             sxs_id = m['sxs_identifier']
#             repr_by_sxs_id[sxs_id].append(r)
#     catalog = OrderedDict()
#     catalog['catalog_file_description'] = catalog_file_description
#     catalog['modified'] = sorted([r['modified'] for sxs_id in repr_by_sxs_id for r in repr_by_sxs_id[sxs_id]])[-1]
#     for sxs_id in sorted(repr_by_sxs_id):
#         catalog[sxs_id] = {}
#         catalog[sxs_id]['conceptrecid'] = repr_by_sxs_id[sxs_id][0]['conceptrecid']
#         catalog[sxs_id]['versions'] = [
#             {'representation': r, 'sxs_metadata': {}}
#             for r in sorted(repr_by_sxs_id[sxs_id], key=lambda r: r['created'])
#         ]
#     return catalog


# def update(catalog_file_name='catalog.json'):
#     import json
#     from . import Login
#     with open(catalog_file_name, 'r') as f:
#         catalog = json.load(f)
#     modified = catalog['modified']
#     l = Login()
#     new_records = l.list_deposits(q='updated:["{0}" TO *]'.format(modified), status='published', all_versions=True)


def map(catalog_file_name='complete_catalog.json', map_file_name='sxs_to_zenodo.map'):
    """Create a mapping from SXS identifiers to Zenodo record numbers for nginx

    The input must be a catalog file.  The output is formatted for inclusion into an nginx
    configuration.  Note that this map file includes a `map_hash_max_size` directive, and thus must
    precede any `map` directives, or you will get an error that this directive is a duplicate (even
    if you never explicitly gave it previously).

    The output map matches both the plain SXS identifier (with nothing following it) and the
    identifier followed by an arbitrary file path.

    Parameters
    ==========
    catalog_file_name: string [defaults to 'complete_catalog.json']
        Relative or absolute path to catalog JSON file.  This is expected to have been created by
        the `update` function.
    map_file_name: string [defaults to 'sxs_to_zenodo.map']
        Relative or absolute path to output file.

    """
    import json
    import math
    with open(catalog_file_name, 'r') as f:
        catalog = json.load(f)
    size = 256 * 2**math.ceil(math.log2(len(catalog)+1))
    def file_prefix(sxs_id):
        prefix = sxs_id + '/'
        files = catalog[sxs_id]['files']
        for f in files:
            if not f['filename'].startswith(prefix):
                return ''
        return prefix
    record_string = "    /waveforms/data/{0} record/{1};\n"
    file_string = "    ~/waveforms/data/{0}/(.*) record/{1}/files/{2}$1;\n"
    with open(map_file_name, 'w') as f:
        f.write("map_hash_max_size {0};\n".format(size))
        f.write("map $uri $zenodo_identifier {\n")
        f.write("    default communities/sxs;\n")
        for sxs_identifier in sorted(catalog):
            f.write(record_string.format(sxs_identifier, catalog[sxs_identifier]['id']))
        for sxs_identifier in sorted(catalog):
            f.write(file_string.format(sxs_identifier, catalog[sxs_identifier]['id'], file_prefix(sxs_identifier)))
        f.write("}\n")


# def catalog(complete_catalog_file_name='complete_catalog.json', public_catalog_file_name='public_catalog.json', *args, **kwargs):
#     """Update (or construct) a catalog of SXS simulations hosted by Zenodo

#     NOTE: This function returns two catalogs: one "complete" catalog containing all records
#     available to the user which may include closed-access datasets (including their files and SXS
#     metadata); and a second catalog containing all information only for datasets marked 'open',
#     along with basic information for all others.

#     This function creates or updates a JSON file containing information for each SXS simulation on
#     Zenodo.  It first looks through all records in the 'sxs' community on Zenodo.  Each one whose
#     title includes an SXS identifier like SXS:BBH:nnnn, etc., is included in the catalog.  The
#     catalog itself is a dictionary mapping the SXS identifier to a "representation".  This is mostly
#     the same as the Zenodo representation described at <http://developers.zenodo.org/#depositions>,
#     but also includes a 'files' field of <http://developers.zenodo.org/#deposition-files>, as well
#     as a 'sxs_metadata' field, containing the metadata from the highest Lev in the dataset.

#     """
#     import re
#     import json
#     import os.path
#     from .. import sxs_identifier_regex
#     from . import records, Login

#     verbosity = kwargs.pop('verbosity', 2)

#     sxs_identifier_regex = re.compile(sxs_identifier_regex)
#     # Download all records in the 'sxs' community from Zenodo
#     kwargs['sxs'] = True
#     r = records(*args, **kwargs)
#     del kwargs['sxs']
#     # Now start a Login object for better interaction with the website
#     l = Login(*args, **kwargs)
#     # Just to make sure the input file contains at least an empty dictionary
#     if not os.path.isfile(complete_catalog_file_name) or not os.path.getsize(complete_catalog_file_name) > 1:
#         with open(complete_catalog_file_name, 'w') as f:
#             f.write('{}')
#     # Load the catalog from the input file
#     with open(complete_catalog_file_name, 'r') as f:
#         complete_catalog = json.load(f)
#     public_catalog = {}
#     # Initialize the catalogs
#     complete_catalog['catalog_file_description'] = catalog_file_description
#     public_catalog['catalog_file_description'] = catalog_file_description
    
#     # Step through the records, making sure we've got everything
#     for i, record in enumerate(r, 1):
#         complete_catalog['records'][record['id']] = record
#         # Get the SXS identifier
#         sxs_identifier_match = sxs_identifier_regex.search(record['title'])
#         if not sxs_identifier_match:
#             print('No SXS identifier found in {0}; skipping.'.format(record['title']))
#             continue
#         sxs_identifier = sxs_identifier_match['sxs_identifier']
#         simulation_type = sxs_identifier_match['simulation_type']
#         if verbosity > 0:
#             print('Working on {0} ({1}/{2})'.format(sxs_identifier, i, len(r)))
#         if sxs_identifier not in complete_catalog:
#             update = True
#         else:
#             # Check to make sure that the local copy is the latest one
#             local_latest = complete_catalog['simulations'][sxs_identifier]['links']['latest_html']
#             zenodo_latest = record['links']['latest_html']
#             if local_latest != zenodo_latest:
#                 update = True
#             else:
#                 complete_catalog['simulations'][sxs_identifier].update(record)  # Just in case there've been metadata changes
#                 update = False
#         if update or 'files' not in complete_catalog[sxs_identifier]:
#             if verbosity > 1:
#                 print('\tUpdating {0}'.format(sxs_identifier))
#             # Get the most recent Deposit object for this record
#             d = l.deposit(record['id'], ignore_deletion=True)
#             if not d.published:
#                 print('Record titled "{0}" has not yet been published; skipping.'.format(record['title']))
#                 continue
#             if not d.is_latest:
#                 d = d.get_latest()
#             # Add the Zenodo representation
#             complete_catalog[sxs_identifier] = d.representation
#             # Add the list of files, along with their MD5 checksums and list of links
#             complete_catalog[sxs_identifier]['files'] = d.files
#         if update or 'sxs_metadata' not in complete_catalog[sxs_identifier]:
#             if verbosity > 1 and not update:
#                 print('\tGetting SXS metadata for {0}'.format(sxs_identifier))
#             # Download and add the SXS metadata
#             metadata_url = sorted([f['links']['download']
#                                    for f in complete_catalog[sxs_identifier]['files']
#                                    if 'metadata.json' in f['filename']])[-1]  # Use highest Lev, for no good reason
#             complete_catalog[sxs_identifier]['sxs_metadata'] = l.session.get(metadata_url).json()
#         if complete_catalog[sxs_identifier]['metadata']['access_right'] == 'open':
#             keys_to_exclude = []
#         else:
#             keys_to_exclude = ['sxs_metadata', 'files']
#         public_catalog[sxs_identifier] = {key:complete_catalog[sxs_identifier][key]
#                                           for key in complete_catalog[sxs_identifier]
#                                           if key not in keys_to_exclude}
#         # # Write a temporary copy of the JSON file, to ensure that the work isn't lost
#         # with open(complete_catalog_file_name+'.tmp', 'w') as f:
#         #     json.dump(complete_catalog, f, indent=4, separators=(',', ': '))
#     with open(complete_catalog_file_name, 'w') as f:
#         json.dump(complete_catalog, f, indent=4, separators=(',', ': '))
#     with open(public_catalog_file_name, 'w') as f:
#         json.dump(public_catalog, f, indent=4, separators=(',', ': '))
#     if verbosity > 2:
#         return complete_catalog
#     else:
#         return
