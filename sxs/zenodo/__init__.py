from __future__ import print_function  # Because apparently print(..., file=f) is the only thing stopping python 2.7
from .api import Login, Deposit, Records

# See https://github.com/moble/nb-marine-science for other examples using the Zenodo API
# The other python API interface I found is here: https://github.com/moble/zenodo-python

def map(catalog_file_name, map_file_name):
    """Create a mapping from SXS identifiers to Zenodo record numbers

    Parameters
    ==========
    catalog_file_name: string
        Relative or absolute path to catalog JSON file.  This is expected to have been created by
        the `catalog` function.
    map_file_name: string
        Relative or absolute path to output file.

    """
    import json
    import math
    with open(catalog_file_name, 'r') as f:
        catalog = json.load(f)
    size = 2**math.ceil(math.log2(len(catalog)))
    with open(map_file_name, 'w') as f:
        print("map_hash_max_size {0};".format(size), file=f)
        print("map $sxs_identifier $zenodo_identifier {", file=f)
        # print("    default https://zenodo.org/communities/sxs/search?page=1&size=20;", file=f)
        for sxs_identifier in sorted(catalog):
            print("    {0} {1};".format(sxs_identifier, catalog[sxs_identifier]['id']), file=f)
        print("}", file=f)


def catalog(catalog_file_name, public_catalog_file_name, *args, **kwargs):
    """Construct a catalog of SXS simulations hosted by Zenodo

    NOTE: This function returns all records available to the user.  This may include closed-access
    datasets.  If you want to create a public catalog, just loop through the output from this
    function, and filter by each entry's ["metadata"]["access_right"] value.

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

    verbosity = kwargs.pop('verbosity', 2)

    sxs_identifier_regex = re.compile(sxs_identifier_regex)
    # Download all records in the 'sxs' community from Zenodo
    kwargs['sxs'] = True
    r = records(*args, **kwargs)
    del kwargs['sxs']
    # Now start a Login object for better interaction with the 
    l = Login(*args, **kwargs)
    # Just to make sure the input file contains at least an empty dictionary
    if not os.path.isfile(catalog_file_name) or not os.path.getsize(catalog_file_name) > 1:
        with open(catalog_file_name, 'w') as f:
            f.write('{}')
    # Load the catalog from the input file
    with open(catalog_file_name, 'r') as f:
        local_catalog = json.load(f)
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
        if sxs_identifier not in local_catalog:
            update = True
        else:
            # Check to make sure that the local copy is the latest one
            local_latest = local_catalog[sxs_identifier]['links']['latest_html']
            zenodo_latest = record['links']['latest_html']
            if local_latest != zenodo_latest:
                update = True
            else:
                local_catalog[sxs_identifier].update(record)  # Just in case there've been metadata changes
                update = False
        if update:
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
            local_catalog[sxs_identifier] = d.representation
            # Add the list of files, along with their MD5 checksums and list of links
            local_catalog[sxs_identifier]['files'] = d.files
        if update or 'sxs_metadata' not in local_catalog[sxs_identifier]:
            if verbosity > 1 and not update:
                print('\Getting SXS metadata for {0}'.format(sxs_identifier))
            # Download and add the SXS metadata
            metadata_url = sorted([f['links']['download']
                                   for f in local_catalog[sxs_identifier]['files']
                                   if 'metadata.json' in f['filename']])[-1]  # Use highest Lev, for no good reason
            local_catalog[sxs_identifier]['sxs_metadata'] = l.session.get(metadata_url).json()
        if local_catalog[sxs_identifier]['metadata']['access_right'] == 'open':
            public_catalog[sxs_identifier] = local_catalog[sxs_identifier]
        # # Write a temporary copy of the JSON file, to ensure that the work isn't lost
        # with open(catalog_file_name+'.tmp', 'w') as f:
        #     json.dump(local_catalog, f, indent=4, separators=(',', ': '))
    with open(catalog_file_name, 'w') as f:
        json.dump(local_catalog, f, indent=4, separators=(',', ': '))
    with open(public_catalog_file_name, 'w') as f:
        json.dump(public_catalog, f, indent=4, separators=(',', ': '))
    if verbosity > 2:
        return local_catalog
    else:
        return

def records(*args, **kwargs):
    """List all deposits"""
    import sys
    import json
    json_output = kwargs.pop('json', False)
    sxs = kwargs.pop('sxs', False)
    q = kwargs.pop('q', '')
    if sxs:
        q += ' communities: "sxs" '
    status = kwargs.pop('status', None)
    sort = kwargs.pop('sort', None)
    page = kwargs.pop('page', None)
    size = kwargs.pop('size', 9999)
    l = Login(*args, **kwargs)
    r = l.list_deposits(q, status, sort, page, size)
    if page is None and len(r) == size:
        # Continue until we've gotten all the records
        last_response = r
        last_page = 1
        while len(last_response) == size:
            last_page = 1
            last_response = l.list_deposits(q, status, sort, last_page+1, size)
            if last_response:
                last_page += 1
                r += last_response
    if json_output:
        return json.dumps(r, indent=4, separators=(',', ': '))
    else:
        return r


def upload(directory, exclude=[],
           sandbox=False, access_token_path=None,
           skip_existing=True, deposition_id=None, ignore_deletion=False,
           access_right='closed', license='CC-BY-4.0',
           creators=[], description='', keywords=[]):
    """Publish or edit a Zenodo entry for an SXS simulation

    This is essentially a wrapper around many of the Zenodo API's functions, specialized for SXS
    systems and intended to account for various possible errors or special conditions

    This function should be able to safely handle
      1) new deposits that Zenodo has not seen previously;
      2) drafts that were started previously but failed for some reason, like a spurious Zenodo
         server failure, or some problem in the data that has now been fixed; or
      3) systems that have been published on Zenodo previously and have not changed at all, but you
         want to verify that the local copy and the version on Zenodo are in sync.

    Most of the parameters to this function are simply passed to other functions.  For more
    explanation of these parameters, see the relevant function's documentation.  Most commonly, the
    only parameter you really need to pass is the first.  You may also wish to pass the last
    parameter if you want the deposit to be published automatically.

    This function returns a Deposit object, which may be used to examine the deposit, change it, or
    publish it if the final parameter is not given as True.


    Parameters to `.api.utilities.find_files`
    =========================================
    directory: string
        Absolute or relative path to a directory containing 'common-metadata.txt' listing an SXS
        identifier starting with 'SXS:BBH:', 'SXS:BHNS:', or 'SXS:NSNS:' and containing at least one
        'metadata.txt' file somewhere in its file hierarchy.
    exclude: list of strings [defaults to an empty list]

    Parameters to `.api.login.Login`
    ================================
    sandbox: bool [defaults to False]
    access_token_path: string or None [defaults to None]

    skip_existing: bool [defaults to True]
        If a record with this name exists already, skip this upload.

    Parameters to `.api.deposit.Deposit`
    ====================================
    deposition_id: string, int, or None [defaults to None]
    ignore_deletion: bool [defaults to False]
        If True and this function call does not succeed in publishing the record, the returned
        Deposit object will issue a warning that it has not been published when it is deleted (which
        may happen after the function returns).
   
    Parameters to `.api.deposit.Deposit.update_metadata`
    ====================================================
    access_right: string [defaults to 'open']
    license: string [defaults to 'cc-by']
    creators: string [defaults to empty list]
    description: string [defaults to '']
    keywords: string [defaults to empty list]
        Note that the last three parameters, if not passed to this function, will be derived
        automatically from the 'metadata.txt' files found in the SXS system directory; they will be
        the union of the parameters found in each file if there are multiple such files.

    """
    import re
    import os
    from .api.utilities import md5checksum, find_files
    from .creators import known_creators, creators_emails, default_creators
    from ..metadata import Metadata
    from .. import sxs_identifier_regex
    default_creators = [{'name': 'SXS Collaboration'}]

    if not os.path.isdir(directory):
        print('The input directory name "{0}" does not appear to be a directory.'.format(directory))
        raise ValueError(directory)
    if not os.path.isfile(os.path.join(directory, 'common-metadata.txt')):
        print('Could not find common-metadata.txt in this directory')
        raise ValueError(directory)
    sxs_system_re = re.compile(sxs_identifier_regex)
    sxs_system_type = None
    sxs_system_number = None
    with open(os.path.join(directory, 'common-metadata.txt')) as f:
        for line in f.readlines():
            line = line.strip()
            if 'alternative-names' in line:
                m = sxs_system_re.search(line)
                if m:
                    sxs_system_type = m['simulation_type']
                    sxs_system_number = m['sxs_number']
                    break
    if not sxs_system_type or not sxs_system_number:
        raise ValueError("No SXS identifier found in common-metadata.txt")
    sxs_system = 'SXS:{0}:{1}'.format(sxs_system_type, sxs_system_number)
    if sxs_system_type == 'BBH':
        title = 'Binary black-hole simulation {0}'.format(sxs_system)
        default_description = """Simulation of a black-hole binary system evolved by the <a href="{0}">SpEC code</a>."""
    elif sxs_system_type == 'BHNS':
        title = 'Black-hole neutron-star binary simulation {0}'.format(sxs_system)
        default_description = """Simulation of a black-hole neutron-star binary system evolved by the <a href="{0}">SpEC code</a>."""
    elif sxs_system_type == 'NSNS':
        title = 'Binary neutron-star simulation {0}'.format(sxs_system)
        default_description = """Simulation of a neutron-star binary system evolved by the <a href="{0}">SpEC code</a>."""
    else:
        raise ValueError('Did not recognize SXS system type "{0}"; should be BBH, BHNS, or NSNS.'.format(sxs_system_type))
    print("Beginning work on {0}".format(sxs_system))

    # Log in to zenodo
    l = Login(sandbox=sandbox, access_token_path=access_token_path)

    # Get this deposit and the title
    if deposition_id is not None:
        if skip_existing:
            print('Asking to skip existing uploads, but also asked for deposition id {0}.'.format(deposition_id))
            print('These are contradictory.  Raise an exception.')
            raise ValueError(deposition_id)
        d = l.deposit(deposition_id, ignore_deletion=ignore_deletion)
        title = d.title
    else:
        # Check to see if this simulation exists in the list of the user's deposits or in the sxs community
        matching_deposits = l.list_deposits(q='title: "{0}"'.format(title))
        if len(matching_deposits) == 1:
            deposition_id = matching_deposits[0]['id']
            print('A deposit with title "{0}"'.format(title))
            if skip_existing:
                print('has already been started with id {0}.'.format(deposition_id))
                raise ValueError(title)
            print('has already been started with this login.  Opening it for editing.')
            d = l.deposit(deposition_id, ignore_deletion=ignore_deletion)
        elif len(matching_deposits) > 1:
            print('Multiple deposits titled "{0}" have been found.'.format(title))
            raise ValueError(title)
        elif len(matching_deposits) == 0:
            # Check to see if this simulation is already in sxs but not owned by this login
            records = Records.search('title: "{0}"'.format(title))
            records = [r for r in records if r['title'] == title]  # Ensure *exact* match
            communities = [community['identifier']
                           for representation in records
                           for community in representation['metadata']['communities']]
            if 'sxs' in communities:
                print('A record exists on Zenodo with exactly the name "{0}",'.format(title))
                print('but is apparently owned by someone else in the "sxs" community.')
                print('Please contact the owner of that record if you wish to update it.')
                print('Web link: {0}'.format(records[0]['links']['html']))
                raise ValueError(title)
            elif len(communities) > 0:
                print('A record exists on Zenodo with exactly the name "{0}",'.format(title))
                print('but is apparently not owned by you, nor is it in the "sxs" community.')
                print('Web link: {0}'.format(records[0]['links']['html']))
                raise ValueError(title)
            d = l.deposit(deposition_id=None, ignore_deletion=ignore_deletion)
    print('Working on deposit "{0}"'.format(title))

    # Convert each metadata.txt file to a metadata.json file sorted with interesting stuff at the
    # top of the file, so it appears prominently on Zenodo's preview without scrolling.  Do this
    # before checking for new files in case these are new or get changed in the process.
    paths_and_names = find_files(directory, exclude=exclude, include_top_directory_in_name=False)
    authors_emails = set()
    point_of_contact_email = ''
    keywords = set(keywords)
    for path,_ in paths_and_names:
        if os.path.basename(path) == 'metadata.txt':
            json_path = os.path.join(os.path.dirname(path), 'metadata.json')
            print('Converting metadata.txt to JSON in {0}'.format(json_path))
            m = Metadata.from_txt_file(path, cache_json=False).reorder_keys()
            del m['metadata_path']  # Don't bother recording the local path to the metadata file
            m.to_json_file(json_path)
            authors_emails |= set(m.get('authors_emails', []))
            point_of_contact_email = m.get('point_of_contact_email', point_of_contact_email)
            keywords |= set(m.get('keywords', []))
                
    # Get list of creators, keywords, and description
    print('Constructing metadata')
    if not creators:
        creators = d.metadata.get('creators', [])
        if not creators:
            authors_emails = list(authors_emails)
            if not authors_emails:
                print("No creators found in input arguments, on Zenodo, or in any metadata.txt file.")
                if point_of_contact_email in creators_emails:
                    creators.append(creators_emails[point_of_contact_email])
                    print('Using point-of-contact email to add', creators_emails[point_of_contact_email]['name'])
                elif point_of_contact_email:
                    print('Unknown point-of-contact email: {0}'.format(point_of_contact_email))
                creators.extend(default_creators)
                print('Adding default creators')
            else:
                for author_email in authors_emails:
                    name = ' '.join(author_email.split()[:-1])
                    if name in known_creators:
                        creators.append(known_creators[name])
                    else:
                        # We tried our best; let's just get this over with.  Sorry Dr. van Whatever.
                        name_parts = name.split()
                        first_name = ' '.join(name_parts[:-1])
                        last_name = name_parts[-1]
                        if first_name:
                            creators.append({'name': '{0}, {1}'.format(last_name, first_name)})
                        else:
                            creators.append({'name': last_name})
    # print('Creators: {0}'.format(creators))
    keywords = list(set(keywords) | set(d.metadata.get('keywords', [])))
    # print('Keywords: {0}'.format(keywords))
    if not description:
        description = d.metadata.get('description', '')
        if not description:
            spec_url = "https://www.black-holes.org/code/SpEC.html"
            description = default_description.format(spec_url)
    # print('Description: {0}'.format(description))
    communities = d.metadata.get('communities', [])
    if 'sxs' not in [c['identifier'] for c in communities]:
        communities.append({'identifier': 'sxs'})
    print('Finished constructing metadata')

    # Send Zenodo the metadata before messing with files, in case this deposit is interrupted (e.g.,
    # due to long upload times)
    new_metadata = {
        'title': title,
        'upload_type': 'dataset',
        'access_right': access_right,
        'license': license,
        'communities': communities,
        'description': description,
        'keywords': keywords,
        'creators': creators,
    }
    # print('New metadata: {0}'.format(new_metadata))
    metadata = d.metadata.copy()
    metadata.update(new_metadata)  # Ensure that fields we haven't changed are still present
    unchanged_metadata = (metadata == d.metadata)
    if unchanged_metadata:
        print('No metadata changed.  Updating it on Zenodo would produce an error, so skipping that.')
    else:
        try:
            d.edit()
        except:
            pass
        d.update_metadata(metadata)
        print('Uploaded metadata')

    # Get the list of files we'll be uploading and compare to files already in the deposit to see if
    # any have changed.  If so, we need to create a new version.  Otherwise, we can just edit this
    # version.
    zenodo_filenames = d.file_names
    local_paths_and_names = find_files(directory, exclude=exclude, include_top_directory_in_name=False)
    if len(local_paths_and_names) == 0:
        print('Zenodo requires that there be at least one file.  None found in {0}.'.format(directory))
        raise ValueError('No files found')
    local_filenames = [name for path, name in local_paths_and_names]
    zenodo_filenames_to_delete = [zf for zf in zenodo_filenames if not zf in local_filenames]
    file_checksums = d.file_checksums  # {filename: checksum}
    print('Comparing MD5 checksums')
    for path, name in local_paths_and_names.copy():
        if name in file_checksums:
            zenodo_checksum = file_checksums[name]
            local_checksum = md5checksum(path)
            if zenodo_checksum == local_checksum:
                local_paths_and_names.remove([path, name])

    # Now, if needed do the file deletions and/or uploads, and publish
    if not zenodo_filenames_to_delete and not local_paths_and_names and unchanged_metadata:
        print('Nothing will change in this deposit.  Just checking that it is published.')
    elif not zenodo_filenames_to_delete and not local_paths_and_names:
        print('Only metadata will change in this deposit.  Proceeding to publication.')
    else:
        if d.published:
            # If this deposit already has files, we need to create a new deposit to change the files
            print('Changing files that are already present in a published deposit.  Getting a new version.')
            d = d.get_new_version()
        else:
            # Otherwise this is presumably a new deposit, so we can add files directly with no trouble
            print('Uploading files to an unpublished deposit.')
        print('Files to delete: {0}'.format(zenodo_filenames_to_delete))
        print('Files to upload: {0}\n'.format([name for path, name in local_paths_and_names]))
        for file_name in zenodo_filenames_to_delete:
            print('\tDeleting {0}'.format(file_name))
            d.delete_file(file_name)
        for path, name in local_paths_and_names:
            print('\tUploading {0}'.format(name))
            d.upload_file(path, name=name, skip_checksum=True)

    # Publish this version
    d.publish()
    print('Finished publishing "{0}" to {1}.'.format(title, d.website))

    return d
