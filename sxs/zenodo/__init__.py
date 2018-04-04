from .api import Login, Deposit, Records


def deposit_sxs_bbh_simulation(sxs_bbh_directory_name, exclude=[],
                               sandbox=False, deposition_id=None, access_token_path=None,
                               access_right='open', license='cc-by',
                               creators=[], description='', keywords=[]):
    """Publish or edit a Zenodo entry for an SXS:BBH simulation


    """
    import re
    import os
    from .api import md5checksum, find_files
    from .creators import known_creators
    from ..metadata import Metadata

    if not os.path.isdir(sxs_bbh_directory_name):
        print('The input directory name "{0}" does not appear to be a directory.'.format(sxs_bbh_directory_name))
        raise ValueError(sxs_bbh_directory_name)
    sxs_bbh = os.path.basename(os.path.normpath(sxs_bbh_directory_name))
    if not sxs_bbh.startswith('SXS:BBH:') or len(sxs_bbh) <= 8:
        print('Input directory "{0}" does not start with "SXS:BBH:" followed by some identifier.'.format(sxs_bbh))
        print('This function is only designed to handle SXS:BBH directories.')
        print('Since the input directory does not conform to this specification,')
        print('this function cannot construct the data with any robustness.')
        print('Feel free to copy this function\'s source to publish to Zenodo')
        print('as you wish, and submit a Pull Request if you feel that it would')
        print('make a useful addition to this module.')
        raise ValueError(sxs_bbh_directory_name)

    # Log in to zenodo
    l = Login(sandbox=sandbox, access_token_path=access_token_path)

    # Get this deposit and the title
    if deposition_id is not None:
        d = l.deposit(deposition_id)
        title = d.title
    else:
        # Check to see if this simulation exists in the list of the user's deposits or in the sxs community
        title = 'Binary black-hole simulation {0}'.format(sxs_bbh)
        matching_deposits = l.list_deposits(q='title: "{0}"'.format(title))
        if len(matching_deposits) == 1:
            deposition_id = matching_deposits[0]['id']
            print('A deposit with title "{0}"'.format(title))
            print('has already been started with this login.  Opening it for editing.')
            d = l.deposit(deposition_id)
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
            d = l.new_deposit
    print('Working on deposit "{0}"'.format(title))

    # Convert each metadata.txt file to a metadata.json file sorted with interesting stuff at the
    # top of the file, so it appears prominently on Zenodo's preview without scrolling.  Do this
    # before checking for new files in case these are new or get changed in the process.
    paths_and_names = find_files(sxs_bbh_directory_name, exclude=exclude)
    authors_emails = set()
    keywords = set(keywords)
    for path,_ in paths_and_names:
        if os.path.basename(path) == 'metadata.txt':
            json_path = os.path.join(os.path.dirname(path), 'metadata.json')
            print('Converting metadata.txt to JSON in {0}'.format(json_path))
            m = Metadata.from_txt_file(path, cache_json=False).reorder_keys()
            m.to_json_file(json_path)
            authors_emails |= set(m.get('authors_emails', []))
            keywords |= set(m.get('keywords', []))

    # Get the list of files we'll be uploading and compare to files already in the deposit to see if
    # any have changed.  If so, we need to create a new version.  Otherwise, we can just edit this
    # version.
    zenodo_filenames = d.file_names
    local_paths_and_names = find_files(sxs_bbh_directory_name, exclude=exclude)
    if len(local_paths_and_names) == 0:
        print('Zenodo requires that there be at least one file.  None found in {0}.'.format(sxs_bbh_directory_name))
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
    if not zenodo_filenames_to_delete and not local_paths_and_names:
        # No files will change, so we just want to edit this Deposit
        print('No files will change; just editing this deposit.')
        d.edit()
    else:
        # We need to create a new deposit to change the files
        d = d.get_new_version()
        print('Changing files.')
        print('Deleting {0}'.format(zenodo_filenames_to_delete))
        print('Uploading {0}'.format([name for path, name in local_paths_and_names]))
        for file_name in zenodo_filenames_to_delete:
            print('\tDeleting {0}'.format(file_name))
            d.delete_file(file_name)
        for path, name in local_paths_and_names:
            print('\tUploading {0}'.format(name))
            d.upload_file(path, name=name, skip_checksum=True)

    # Get list of creators, keywords, and description
    if not creators:
        creators = d.metadata.get('creators', [])
        if not creators:
            authors_emails = list(authors_emails)
            if not authors_emails:
                print("No creators found in input arguments, on Zenodo, or in any metadata.txt file.")
                raise ValueError("Missing creators")
            creators = []
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
    print('Creators: {0}'.format(creators))
    keywords = list(set(keywords) | set(d.metadata.get('keywords', [])))
    print('Keywords: {0}'.format(keywords))
    if not description:
        description = d.metadata.get('description', '')
        if not description:
            spec_url = "https://www.black-holes.org/code/SpEC.html"
            description = """Simulation of a black-hole binary system evolved by the <a href="{0}">SpEC code</a>."""
            description = description.format(spec_url)
    print('Description: {0}'.format(description))

    # Construct the Zenodo metadata
    new_metadata = {
        'title': title,
        'description': description,
        'keywords': keywords,
        'creators': creators,
        'upload_type': 'dataset',
        'access_right': access_right,
        'license': 'cc-by',
        'communities': [{'identifier': 'sxs'}],
    }
    print('New metadata: {0}'.format(new_metadata))
    metadata = d.metadata
    metadata.update(new_metadata)  # Ensure that fields we haven't changed are still present
    d.update_metadata(metadata)

    # Publish this version
    d.publish()
    print('Finished publishing {0} to {1}.'.format(title, d.website))

    return d
