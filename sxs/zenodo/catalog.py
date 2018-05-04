
catalog_file_description = """
        This JSON file has the following format.  Comments are, of course, not present (since JSON does not support
        comments).  Single quotes here are, of course, double quotes in the rest of the file (since JSON encloses
        strings in double quotes).  Anything inside <angle brackets> is just standing in for the relevant value.  An
        ellipsis ... indicates that the preceding structure can be repeated.

          {
              'catalog_file_description': '<this description>',
              'last_file_modification': '<YYYY-MM-DDThh:mm:ss.ssssss>',  # UTC time of last modification made to this file
              '<sxs_id>': {  # The SXS ID is a string like SXS:BHNS:0001 or SXS:BBH:1234
                  'conceptrecid': '<conceptrecid>',  # The Zenodo ID of the 'concept' record, which resolves to the current version
                  'versions': [
                      {  # Oldest version first
                          'representation': {  # More details about this object at http://developers.zenodo.org/#depositions
                              'conceptdoi': '10.5281/zenodo.<conceptrecid>',  # Permanent DOI for all versions of this record
                              'conceptrecid': '<conceptrecid>',  # ~7-digit integer identifying collectively all versions of this record
                              'created': '<YYYY-MM-DDThh:mm:ss.ssssss>',  # UTC time of creation of this record on Zenodo
                              'doi': '10.5281/zenodo.<id>',  # Permanent DOI for this record
                              'doi_url': 'https://doi.org/10.5281/zenodo.<id>',  # URL for permanent DOI of this record
                              'id': <id>,  # ~7-digit integer uniquely identifying this record
                              'links': {
                                   'badge': 'https://zenodo.org/badge/doi/10.5281/zenodo.<id>.svg',
                                   'bucket': 'https://zenodo.org/api/files/<uuid>',  # Base URL for file uploads and downloads
                                   'conceptbadge': 'https://zenodo.org/badge/doi/10.5281/zenodo.<conceptrecid>.svg',
                                   'conceptdoi': 'https://doi.org/10.5281/zenodo.<conceptrecid>',
                                   'discard': 'https://zenodo.org/api/deposit/depositions/<id>/actions/discard',
                                   'doi': 'https://doi.org/10.5281/zenodo.<id>',  # Permanent URL for this version
                                   'edit': 'https://zenodo.org/api/deposit/depositions/<id>/actions/edit',
                                   'files': 'https://zenodo.org/api/deposit/depositions/<id>/files',
                                   'html': 'https://zenodo.org/deposit/<id>',  # Webpage for this version
                                   'latest': 'https://zenodo.org/api/records/<id>',  # 
                                   'latest_html': 'https://zenodo.org/record/<id>',  # Webpage for most recent version
                                   'publish': 'https://zenodo.org/api/deposit/depositions/<id>/actions/publish',
                                   'record': 'https://zenodo.org/api/records/<id>',
                                   'record_html': 'https://zenodo.org/record/<id>',
                                   'self': 'https://zenodo.org/api/deposit/depositions/<id>'
                              },
                              'metadata': {  # Note that this is different from the SXS metadata
                                  'access_right': '<access>',  # Can be 'open', 'closed', 'embargoed', or 'restricted'
                                  'communities': [
                                      {'identifier': '<community_name>'},  # Names may include 'sxs' and 'zenodo'
                                      ...
                                  ],
                                  'creators': [
                                      {
                                          'name': '<name>',  # Name of this creator in the format Family name, Given names
                                          'affiliation': '<affiliation>',  # Affiliation of this creator (optional)
                                          'orcid': '<orcid>',  # ORCID identifier of this creator (optional)
                                          'gnd': '<gnd>'  # GND identifier of this creator (optional)
                                      },
                                      ...
                                  ],
                                  'description': '<description>',  # Text description of this record
                                  'doi': '10.5281/zenodo.<id>',  # Permanent DOI of this record
                                  'keywords': [
                                      '<keyword>',  # Optional; this array may be empty
                                      ...
                                  ],
                                  'prereserve_doi': {'doi': '10.5281/zenodo.<id>', 'recid': <id>},
                                  'publication_date': '<YYYY-MM-DD>',  # Possibly meaningless date (UTC)
                                  'title': '<title>',
                                  'upload_type': 'dataset'
                              },
                              'modified': '<YYYY-MM-DDThh:mm:ss.ssssss>',  # (UTC) Last modification of this record (possibly just Zenodo metadata)
                              'owner': <user_id>,  # ~5-digit integer identifying the user who owns this record
                              'record_id': <id>,  # Same as 'id'
                              'state': '<state>',  # Can be 'done', 'inprogress', 'error', 'unsubmitted', possibly others
                              'submitted': <submitted>,  # True or false (always true for published records)
                              'title': '<title>'  # Same as ['metadata']['title']
                          }
                          'files': [
                              # See https://data.black-holes.org/waveforms/documentation.html for
                              # detailed descriptions of the *contents* of the files in each record.
                              {
                                  'checksum': '<checksum>',  # MD5 checksum of file on Zenodo
                                  'filename': '<filename>',  # Name of file; may contain slashes denoting directories
                                  'filesize': <filesize>,  # Number of bytes in the file
                                  'id': '<id>',  # A standard UUID (hexadecimal with characters in the pattern 8-4-4-4-12)
                                  'links': {
                                      'download': 'https://zenodo.org/api/files/<bucket>/<filename>',  # The URL to use to download this file
                                      'self': 'https://zenodo.org/api/deposit/depositions/<deposition_id>/files/<id>'  # Ignore this
                                  }
                              },
                              ...  # Other file descriptions in the order in which they were uploaded (not necessarily a meaningful order)
                          ]
                          'sxs_metadata': {
                              # Variable content describing (mostly) physical parameters of the system.  It's basically
                              # a python-compatible version of the information contained in 'metadata.txt' from the
                              # highest-resolution run in this record.  That file is meant to be more-or-less as
                              # suggested in <https://arxiv.org/abs/0709.0093>.  The conversion to a python-compatible
                              # format means that keys like 'simulation-name' have had hyphens replaced by underscores
                              # so that they can be used as variable names in python and any other sane language (with
                              # apologies to Lisp).  As far as possible, values that are just strings in that file have
                              # been converted into the relevant types -- like numbers, integers, and arrays.  Note that
                              # some keys like eccentricity are sometimes numbers and sometimes the string '<number'
                              # (meaning that the eccentricity is less than the number), which is necessarily a string.
                              #
                              # Below are just the first few keys that *may* be present.
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
                          }
                      },
                      ...  # Progressively newer versions.  There may only be one.  To get the newest in every case, just use ['versions'][-1]
                  ]
              },
              ...  # More SXS IDs
          }
          """

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
