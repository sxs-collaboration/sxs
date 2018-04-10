from .api import Login, Deposit, Records

# See https://github.com/moble/nb-marine-science for other examples using the Zenodo API
# The other python API interface I found is here: https://github.com/moble/zenodo-python


def deposit_sxs_simulation(sxs_system_directory_name, exclude=[],
                           sandbox=False, access_token_path=None,
                           deposition_id=None, ignore_deletion=False,
                           access_right='open', license='CC-BY-4.0',
                           creators=[], description='', keywords=[],
                           publish=False):
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
    sxs_system_directory_name: string
        Absolute or relative path to a directory starting with 'SXS:BBH:', 'SXS:BHNS:', or
        'SXS:NSNS:' and containing at least one 'metadata.txt' file somewhere in its file hierarchy.
    exclude: list of strings [defaults to an empty list]

    Parameters to `.api.login.Login`
    ================================
    sandbox: bool [defaults to False]
    access_token_path: string or None [defaults to None]

    Parameters to `.api.deposit.Deposit`
    ====================================
    deposition_id: string, int, or None [defaults to None]
    ignore_deletion: bool [defaults to False]
   
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

    Final parameter
    ==================
    publish: bool [defaults to False]
        If True, and everything before it succeeds, publish the deposit on Zenodo.

    """
    import re
    import os
    from .api import md5checksum, find_files
    from .creators import known_creators, creators_emails, default_creators
    from ..metadata import Metadata

    if not os.path.isdir(sxs_system_directory_name):
        print('The input directory name "{0}" does not appear to be a directory.'.format(sxs_system_directory_name))
        raise ValueError(sxs_system_directory_name)
    sxs_system = os.path.basename(os.path.normpath(sxs_system_directory_name))
    if sxs_system.startswith('BHNS:') or sxs_system.startswith('NSNS:'):
        sxs_system = 'SXS:' + sxs_system
    if not ((sxs_system.startswith('SXS:BBH:') and len(sxs_system) > 8)
            or (sxs_system.startswith('SXS:BHNS:') and len(sxs_system) > 9)
            or (sxs_system.startswith('SXS:NSNS:') and len(sxs_system) > 9)):
        print('Input directory "{0}" does not start with "SXS:BBH:", "SXS:BHNS:", or "SXS:NSNS:" followed by some identifier.'.format(sxs_system))
        print('This function is only designed to handle SXS directories.  Since the input directory does not conform to this')
        print('specification, this function cannot construct the data with any robustness.  Feel free to copy this function\'s')
        print('source to publish to Zenodo as you wish, and submit a Pull Request if you feel that it would make a useful')
        print('addition to this module.')
        raise ValueError(sxs_system_directory_name)
    print("Beginning work on {0}".format(sxs_system))

    # Log in to zenodo
    l = Login(sandbox=sandbox, access_token_path=access_token_path)

    # Get this deposit and the title
    new_deposit = False
    if deposition_id is not None:
        d = l.deposit(deposition_id, ignore_deletion=ignore_deletion)
        title = d.title
        new_deposit = False
    else:
        # Check to see if this simulation exists in the list of the user's deposits or in the sxs community
        title = 'Binary black-hole simulation {0}'.format(sxs_system)
        matching_deposits = l.list_deposits(q='title: "{0}"'.format(title))
        if len(matching_deposits) == 1:
            deposition_id = matching_deposits[0]['id']
            print('A deposit with title "{0}"'.format(title))
            print('has already been started with this login.  Opening it for editing.')
            d = l.deposit(deposition_id, ignore_deletion=ignore_deletion)
            new_deposit = False
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
            new_deposit = True
    print('Working on deposit "{0}"'.format(title))

    # Convert each metadata.txt file to a metadata.json file sorted with interesting stuff at the
    # top of the file, so it appears prominently on Zenodo's preview without scrolling.  Do this
    # before checking for new files in case these are new or get changed in the process.
    paths_and_names = find_files(sxs_system_directory_name, exclude=exclude)
    authors_emails = set()
    point_of_contact_email = ''
    keywords = set(keywords)
    for path,_ in paths_and_names:
        if os.path.basename(path) == 'metadata.txt':
            json_path = os.path.join(os.path.dirname(path), 'metadata.json')
            print('Converting metadata.txt to JSON in {0}'.format(json_path))
            m = Metadata.from_txt_file(path, cache_json=False).reorder_keys()
            del m['metadata_path']
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
            description = """Simulation of a black-hole binary system evolved by the <a href="{0}">SpEC code</a>."""
            description = description.format(spec_url)
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
    local_paths_and_names = find_files(sxs_system_directory_name, exclude=exclude)
    if len(local_paths_and_names) == 0:
        print('Zenodo requires that there be at least one file.  None found in {0}.'.format(sxs_system_directory_name))
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
    if publish:
        if d.published:
            print('Nothing has changed in this deposit, and it has already been published.')
        else:
            d.publish()
            print('Finished publishing "{0}" to {1}.'.format(title, d.website))
    else:
        print('As requested, {0} has not yet been published.'.format(title))
        print('If you want to publish it, you can simply take the object returned from')
        print('this function and run its `.publish()` method.  Alternatively, you can')
        print('publish this deposit via the web interface at\n')
        print('    {0}\n'.format(d.website))

    return d
