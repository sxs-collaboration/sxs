"""Functions for syncing the SurrogateModeling data to Zenodo

"""

def sync(annex_dir='./', exclude=None, publish='if_pending',
         sandbox=False, access_token_path=None, skip_checksums='if_file_is_older'):
    """Sync the SurrogateModeling annex to Zenodo

    This function first searches for any changes in the data files; if there are
    any, a new version is created on Zenodo, and the relevant files are uploaded,
    deleted, etc.  It then collects the Zenodo metadata (author list, related
    identifiers, etc.) from zenodo_metadata.json, along with the description to be
    placed on the Zenodo record from zenodo_description.html.  These changes are
    made to the Zenodo record, whether the existing record (if no files changed) or
    the new record (if files changed).

    Parameters
    ----------
    annex_dir : str, optional [defaults to './']
        Absolute or relative path to base directory of SurrogateModeling annex with
        data files.

    """
    import os
    import json
    import datetime
    import pytz
    from ..utilities import md5checksum, find_files, fit_to_console
    from .api import Login

    if exclude is None:
        exclude = []

    l = Login(sandbox=sandbox, access_token_path=access_token_path)
    surrogate_record = l.search(q="title:'Binary black-hole surrogate waveform catalog'", size=1, max_pages=1)[0]
    surrogate_id = surrogate_record['id']
    d = l.deposit(surrogate_id, ignore_deletion=True)

    # Get the time at which this deposit was created, for possible comparison to file times below
    deposit_created = datetime.datetime.strptime(d.representation['created'], '%Y-%m-%dT%H:%M:%S.%f%z')  # in UTC

    # Read the metadata from the top-level file 'zenodo_metadata.json'
    stored_keys = ['creators', 'related_identifiers', 'title']
    metadata = d.metadata.copy()
    with open(os.path.join(annex_dir, 'zenodo_metadata.json'), 'r') as f:
        local_metadata = json.load(f)
    unchanged_metadata = True
    for key in stored_keys:
        if d.metadata.get(key, '') != local_metadata.get(key, ''):
            unchanged_metadata = False
            metadata[key] = local_metadata.get(key, '')

    # Read the description from the top-level file 'zenodo_description.html'
    zenodo_description = d.metadata['description']
    with open(os.path.join(annex_dir, 'zenodo_description.html'), 'r') as f:
        local_description = f.read()
    unchanged_description = (zenodo_description == local_description)
    if not unchanged_description:
        metadata.update({'description': local_description})

    # Make any changes
    if unchanged_description and unchanged_metadata:
        print('No zenodo metadata changed.')
    else:
        try:
            d.edit()
        except Exception:
            print("Failure to unlock deposit is probably acceptable")
        d.update_metadata(metadata)
        print('Uploaded new metadata')

    # Get the list of files we'll be uploading and compare to files already in the deposit to see if
    # any have changed.  If so, we need to create a new version.  Otherwise, we can just edit this
    # version.
    zenodo_filenames = [fn for fn in d.file_names if not fn.endswith('.html')]
    local_paths_and_names = find_files(os.path.join(annex_dir, 'data'), exclude=exclude,
                                       include_top_directory_in_name=False)
    if len(local_paths_and_names) == 0:
        print('Zenodo requires that there be at least one file.  None found in {0}.'.format(annex_dir))
        raise ValueError('No files found')
    names_to_delete = sorted(set(zenodo_filenames) - set(name for path, name in local_paths_and_names))
    zenodo_file_sizes = d.file_sizes  # formatted as {filename: size_in_bytes}
    zenodo_file_checksums = d.file_checksums  # formatted as {filename: md5checksum}
    print('Comparing sizes and MD5 checksums')
    local_timezone = pytz.timezone("America/Los_Angeles")  # assuming this is the server's timezone
    for path, name in local_paths_and_names.copy():
        if name in zenodo_file_sizes:
            zenodo_filesize = zenodo_file_sizes[name]
            local_filesize = os.path.getsize(path)
            if zenodo_filesize != local_filesize:
                continue  # short-circuit the md5 check because we know they'll be different
            elif skip_checksums == 'if_file_is_older':
                local_mtime = datetime.datetime.fromtimestamp(os.path.getmtime(path))
                if local_timezone.localize(local_mtime) < deposit_created:
                    local_paths_and_names.remove([path, name])
                    continue  # Just assume these files are the same and skip the checksum
            elif skip_checksums:
                local_paths_and_names.remove([path, name])
                continue  # Just assume these files are the same and skip the checksum
        if name in zenodo_file_checksums:
            zenodo_checksum = zenodo_file_checksums[name]
            local_checksum = md5checksum(path)
            if zenodo_checksum == local_checksum:
                local_paths_and_names.remove([path, name])
    name_to_path_map = {name: path for path, name in local_paths_and_names}
    names_to_upload = sorted(set(name_to_path_map) - set(zenodo_filenames))
    names_to_replace = sorted(set(name_to_path_map) & set(zenodo_filenames))
    paths_and_names_to_upload = [[name_to_path_map[name], name] for name in names_to_upload]
    paths_and_names_to_replace = [[name_to_path_map[name], name] for name in names_to_replace]

    # Now, if needed do the file deletions and/or uploads, and publish
    if not names_to_delete and not names_to_upload and not names_to_replace and unchanged_description:
        if publish == 'if_pending':
            publish = True
        print('Nothing will change in this deposit.  Just checking that it is published.')
    elif not names_to_delete and not names_to_upload and not names_to_replace:
        if publish == 'if_pending':
            publish = False
        print('Only metadata will change in this deposit.  Proceeding to publication.')
    else:
        if publish == 'if_pending':
            publish = False
        if d.published:
            # If this deposit already has files, we need to create a new deposit to change the files
            print('Changing files that are already present in a published deposit.  Getting a new version.')
            d = d.get_new_version()
        else:
            # Otherwise this is presumably a new deposit, so we can add files directly with no trouble
            print('Uploading files to an unpublished deposit.')
        # Set the new publication_date
        publication_date = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d")
        metadata = d.metadata.copy()
        metadata.update({'publication_date': publication_date})
        d.update_metadata(metadata, refresh_information=False)
        # Make the changes to the files
        print(fit_to_console(names_to_delete, 'Files to delete: ', subsequent_indent='    '))
        print(fit_to_console(names_to_upload, 'Files to upload: ', subsequent_indent='    '))
        print(fit_to_console(names_to_replace, 'Files to replace: ', subsequent_indent='    '))
        for name in names_to_delete:
            print('    Deleting {0}'.format(name))
            d.delete_file(name, refresh_information=False)
        for path, name in paths_and_names_to_upload:
            print('    Uploading {0}'.format(name))
            d.upload_file(path, name=name, skip_checksum=True, refresh_information=False)
        for path, name in paths_and_names_to_replace:
            print('    Replacing {0}'.format(name))
            d.delete_file(name, refresh_information=False)
            d.upload_file(path, name=name, skip_checksum=True, refresh_information=False)

    # Before returning or checking for published state, make sure our copy of the deposit is up to date with zenodo
    # noinspection PyStatementEffect
    d.refresh_information

    # Publish this version
    if publish:
        if not d.published:
            d.publish()
            print('Finished publishing "{0}" to {1}.'.format(d.title, d.website))
        else:
            print('Already published "{0}" to {1}.'.format(d.title, d.website))
    else:
        print('Skipping publication, as requested.')

    return d
