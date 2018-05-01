

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


def update(complete_catalog_file_name='complete_catalog.json', public_catalog_file_name='public_catalog.json', *args, **kwargs):
    """Update (or construct) a catalog of SXS simulations hosted by Zenodo

    NOTE: This function returns two catalogs: one "complete" catalog containing all records
    available to the user which may include closed-access datasets (including their files and SXS
    metadata); and a second catalog containing all information only for datasets marked 'open',
    along with basic information for all others.

    This function creates or updates a JSON file containing information for each SXS simulation on
    Zenodo.  It first looks through all records in the 'sxs' community on Zenodo.  Each one whose
    title includes an SXS identifier like SXS:BBH:nnnn, etc., is included in the catalog.  The
    catalog itself is a dictionary mapping the SXS identifier to a "representation".  This is mostly
    the same as the Zenodo representation described at <http://developers.zenodo.org/#depositions>,
    but also includes a 'files' field of <http://developers.zenodo.org/#deposition-files>, as well
    as a 'sxs_metadata' field, containing the metadata from the highest Lev in the dataset.

    """
    import re
    import json
    import os.path
    from .. import sxs_identifier_regex
    from . import records, Login

    verbosity = kwargs.pop('verbosity', 2)

    sxs_identifier_regex = re.compile(sxs_identifier_regex)
    # Download all records in the 'sxs' community from Zenodo
    kwargs['sxs'] = True
    r = records(*args, **kwargs)
    del kwargs['sxs']
    # Now start a Login object for better interaction with the website
    l = Login(*args, **kwargs)
    # Just to make sure the input file contains at least an empty dictionary
    if not os.path.isfile(complete_catalog_file_name) or not os.path.getsize(complete_catalog_file_name) > 1:
        with open(complete_catalog_file_name, 'w') as f:
            f.write('{}')
    # Load the catalog from the input file
    with open(complete_catalog_file_name, 'r') as f:
        complete_catalog = json.load(f)
    public_catalog = {}
    # Step through the records, making sure we've got everything
    for i, record in enumerate(r, 1):
        # Get the SXS identifier
        sxs_identifier_match = sxs_identifier_regex.search(record['title'])
        if not sxs_identifier_match:
            print('No SXS identifier found in {0}; skipping.'.format(record['title']))
            continue
        sxs_identifier = sxs_identifier_match['sxs_identifier']
        simulation_type = sxs_identifier_match['simulation_type']
        if verbosity > 0:
            print('Working on {0} ({1}/{2})'.format(sxs_identifier, i, len(r)))
        if sxs_identifier not in complete_catalog:
            update = True
        else:
            # Check to make sure that the local copy is the latest one
            local_latest = complete_catalog[sxs_identifier]['links']['latest_html']
            zenodo_latest = record['links']['latest_html']
            if local_latest != zenodo_latest:
                update = True
            else:
                complete_catalog[sxs_identifier].update(record)  # Just in case there've been metadata changes
                update = False
        if update or 'files' not in complete_catalog[sxs_identifier]:
            if verbosity > 1:
                print('\tUpdating {0}'.format(sxs_identifier))
            # Get the most recent Deposit object for this record
            d = l.deposit(record['id'], ignore_deletion=True)
            if not d.published:
                print('Record titled "{0}" has not yet been published; skipping.'.format(record['title']))
                continue
            if not d.is_latest:
                d = d.get_latest()
            # Add the Zenodo representation
            complete_catalog[sxs_identifier] = d.representation
            # Add the list of files, along with their MD5 checksums and list of links
            complete_catalog[sxs_identifier]['files'] = d.files
        if update or 'sxs_metadata' not in complete_catalog[sxs_identifier]:
            if verbosity > 1 and not update:
                print('\tGetting SXS metadata for {0}'.format(sxs_identifier))
            # Download and add the SXS metadata
            metadata_url = sorted([f['links']['download']
                                   for f in complete_catalog[sxs_identifier]['files']
                                   if 'metadata.json' in f['filename']])[-1]  # Use highest Lev, for no good reason
            complete_catalog[sxs_identifier]['sxs_metadata'] = l.session.get(metadata_url).json()
        if complete_catalog[sxs_identifier]['metadata']['access_right'] == 'open':
            keys_to_exclude = []
        else:
            keys_to_exclude = ['sxs_metadata', 'files']
        public_catalog[sxs_identifier] = {key:complete_catalog[sxs_identifier][key]
                                          for key in complete_catalog[sxs_identifier]
                                          if key not in keys_to_exclude}
        # # Write a temporary copy of the JSON file, to ensure that the work isn't lost
        # with open(complete_catalog_file_name+'.tmp', 'w') as f:
        #     json.dump(complete_catalog, f, indent=4, separators=(',', ': '))
    with open(complete_catalog_file_name, 'w') as f:
        json.dump(complete_catalog, f, indent=4, separators=(',', ': '))
    with open(public_catalog_file_name, 'w') as f:
        json.dump(public_catalog, f, indent=4, separators=(',', ': '))
    if verbosity > 2:
        return complete_catalog
    else:
        return
