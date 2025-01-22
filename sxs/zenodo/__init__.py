"""Interface to Zenodo or similar Invenio-backed repositories

"""

import os.path
from pathlib import Path

from .api import Login, Deposit, Records
from . import catalog, simannex, surrogatemodeling

# See https://github.com/moble/nb-marine-science for other examples using the Zenodo API
# The other python API interface I found is here: https://github.com/moble/zenodo-python


def path_to_invenio(file_path):
    """Convert a file path to an invenio-compatible name"""
    return str(file_path).replace(os.path.sep, ":")


def invenio_to_path(file_name):
    """Convert an invenio-compatible name to a file path"""
    return Path(file_name.replace(":", os.path.sep))


def translate(sxs_identifier, url=False):
    """Query data.black-holes.org to get the current Zenodo equivalent of the given SXS ID

    Parameters
    ----------
    sxs_identifier: string
        The SXS identifier of the simulation (something like SXS:BBH:1234, or
        SXS:BHNS:0001).
    url: bool [defaults to False]
        If True, return the complete URL to the Zenodo record.  Otherwise, only the
        Zenodo identifier of the record (the last part of the URL and DOI) will be
        returned.

    """
    import requests
    r = requests.get('https://data.black-holes.org/waveforms/data/{0}'.format(sxs_identifier), allow_redirects=False)
    if url:
        return r.headers['Location']
    else:
        return r.headers['Location'][len('https://zenodo.org/record/'):]


def records(*args, **kwargs):
    """List all published records

    By default, this function lists all deposits by the current user, logging in by
    default as with the sxs.zenodo.Login class.  Optional parameters allow for
    different searches.

    A list of dicts is returned, each of which contains the "representation" of the
    record, as described in the API documentation:
    http://developers.zenodo.org/#depositions.

    Parameters
    ----------
    json_output: bool [defaults to False]
        If True, this function returns a JSON string; otherwise, it returns a
        python dict.
    sxs: bool [defaults to False]
        If True, this function looks for all records in the 'sxs' community on
        Zenodo.
    all_versions: bool [defaults to False]
        If True, return all versions, not just the most recent version of each
        record.

    Parameters for .api.Login.search
    --------------------------------
    q: string
        Search query, using Elasticsearch query string syntax.  See
        https://help.zenodo.org/guides/search/ for details.  If the above argument
        `sxs` is True, the string ' communities: "sxs" ' is added to the query.
    status: string
        Filter result based on deposit status (either 'draft' or 'published')
    sort: string
        Sort order ('bestmatch' or 'mostrecent').  Prefix with minus to change form
        ascending to descending (e.g., '-mostrecent').
    page: int
        Page number for pagination
    size: int
        Number of results to return per page.  Note that Zenodo (as of this
        writing) seems to place a hard limit of 9999 responses.  Anything more will
        result in an error.

    All remaining parameters are passed to .api.Login.

    """
    import json
    json_output = kwargs.pop('json_output', kwargs.pop('json_output', False))
    sxs = kwargs.pop('sxs', False)
    all_versions = kwargs.pop('all_versions', False)
    q = kwargs.pop('q', '')
    if sxs:
        q += ' communities: "sxs" '
    status = kwargs.pop('status', None)
    sort = kwargs.pop('sort', None)
    page = kwargs.pop('page', None)
    size = kwargs.pop('size', 9999)
    l = Login(*args, **kwargs)
    r = l.search(q, status, sort, page, size, all_versions=all_versions)
    if json_output:
        return json.dumps(r, indent=4, separators=(',', ': '))
    else:
        return r


def related_identifier_formatter(identifier, relation='isSupplementTo', scheme='url'):
    """Convert an identifier to a representation acceptable to Zenodo (INCOMPLETE)

    NOTE: Because of missing documentation, this function is likely incomplete.  In
    particular, the scheme is assumed to be 'url' for everything except simple DOI
    URLs and arXiv links, which are converted.  It is unclear what Zenodo else
    actually converts, and how.

    """
    if 'https://dx.doi.org/' in identifier:
        identifier = identifier.replace('https://dx.doi.org/', '')
        scheme = 'doi'
    elif identifier.lower().startswith('arxiv:'):
        identifier = 'arXiv:' + identifier[6:]
        scheme = 'arxiv'
    elif identifier.startswith('https://arxiv.org/'):
        identifier = 'arXiv:'+identifier.replace('https://arxiv.org/abs/', '').replace('https://arxiv.org/pdf/', '').replace('https://arxiv.org/', '')
        scheme = 'arxiv'
    elif identifier.startswith('http://arxiv.org/'):
        identifier = 'arXiv:'+identifier.replace('http://arxiv.org/abs/', '').replace('http://arxiv.org/pdf/', '').replace('http://arxiv.org/', '')
        scheme = 'arxiv'
    return {
        'identifier': identifier,
        'relation': relation,
        'scheme': scheme,
    }


def upload(directory,
           exclude=[
               'HorizonsDump.h5', 'RedshiftQuantities.h5', 'SpEC.out',
               'rh_FiniteRadii_CodeUnits.h5', 'rPsi4_FiniteRadii_CodeUnits.h5',
               'rhOverM_Asymptotic_GeometricUnits.h5', 'rMPsi4_Asymptotic_GeometricUnits.h5',
           ],
           sandbox=False, access_token_path=None, skip_checksums='if_file_is_older',
           skip_existing=True, deposition_id=None, ignore_deletion=False,
           access_right='closed', license='CC-BY-4.0',
           creators=[], description='', keywords=[], related_identifiers=[],
           error_on_existing=True, publish=True):
    """Publish or edit a Zenodo entry for an SXS simulation

    This is essentially a wrapper around many of the Zenodo API's functions,
    specialized for SXS systems and intended to account for various possible errors
    or special conditions.

    This function should be able to safely handle
      1) new deposits that Zenodo has not seen previously;
      2) drafts that were started previously but failed for some reason, like a
         spurious Zenodo server failure, or some problem in the data that has now
         been fixed; or
      3) systems that have been published on Zenodo previously but have changed in
         some way, so you want to ensure that the local copy and the version on
         Zenodo are in sync.
      4) systems that have been published on Zenodo previously and have not changed
         at all, but you want to verify that the local copy and the version on
         Zenodo are in sync.

    Most of the parameters to this function are simply passed to other functions.
    For more explanation of these parameters, see the relevant function's
    documentation.  Most commonly, the only parameter you really need to pass is
    the first.  You may also wish to pass the last parameter if you want the
    deposit to be published automatically.

    This function returns a Deposit object, which may be used to examine the
    deposit, change it, or publish it if the final parameter is not given as True.

    Parameters only used in this function
    -------------------------------------
    skip_checksums : bool or 'if_file_is_older' [defaults to 'if_file_is_older']
        If False, an MD5 checksum is run for any file that exists on zenodo and
        locally (unless the file sizes are different).  If 'if_file_is_older', the
        files are assumed to match if the local modification time is earlier than
        the creation date of the deposit on zenodo; otherwise, the checksum is run
        (unless the file sizes are different).  If the file sizes or checksums are
        different, the local file is uploaded to zenodo.
    skip_existing : bool [defaults to True]
        If a record with this name exists already, skip this upload.
    error_on_existing : bool [defaults to False]
        If True, and a record already exists for this system, raise an error.
    publish : bool or 'if_pending' [defaults to True]
        If True and the current version on zenodo is not published, publish it.  If
        the input value is the string 'if_pending', it will be changed to True if
        no other changes are made during this run of the function, or to False if
        other changes are made.  This allows for a delay, to ensure that the files
        have stopped changing before the system is actually published.


    Parameters to `.utilities.find_files`
    -------------------------------------
    directory : string
        Absolute or relative path to a directory containing 'common-metadata.txt'
        listing an SXS identifier starting with 'SXS:BBH:', 'SXS:BHNS:', or
        'SXS:NSNS:' and containing at least one 'metadata.txt' file somewhere in
        its file hierarchy.
    exclude : list of strings [defaults to an empty list]

    Parameters to `.api.login.Login`
    --------------------------------
    sandbox : bool [defaults to False]
    access_token_path : string or None [defaults to None]

    Parameters to `.api.deposit.Deposit`
    ------------------------------------
    deposition_id : string, int, or None [defaults to None]
    ignore_deletion : bool [defaults to False]
        If True and this function call does not succeed in publishing the record,
        the returned Deposit object will issue a warning that it has not been
        published when it is deleted (which may happen after the function returns).
   
    Parameters to `.api.deposit.Deposit.update_metadata`
    ----------------------------------------------------
    access_right : string [defaults to 'open']
    license : string [defaults to 'cc-by']
    creators : string [defaults to empty list]
    description : string [defaults to '']
    keywords : string [defaults to empty list]
        Note that the last three parameters, if not passed to this function, will
        be derived automatically from the 'metadata.txt' files found in the SXS
        system directory; they will be the union of the parameters found in each
        file if there are multiple such files.

    """
    import re
    import os
    import datetime
    import pytz
    from ..utilities import sxs_identifier_regex, inspire, md5checksum, find_files, fit_to_console
    from ..metadata import Metadata
    from .creators import known_creators, creators_emails, default_creators
    default_creators = [{'name': 'SXS Collaboration'}]

    directory = os.path.normpath(directory)
    if not os.path.isdir(directory) and os.path.isfile(directory):
        directory = os.path.dirname(directory)
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
        raise ValueError("No SXS identifier found in {0}".format(os.path.join(directory, 'common-metadata.txt')))
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
        raise ValueError('Did not recognize SXS system type "{0}"'.format(sxs_system_type)
                         + 'in directory "{0}"; should be BBH, BHNS, or NSNS.'.format(directory))
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
        matching_deposits = l.search(q='title: "{0}"'.format(title))
        if len(matching_deposits) == 1:
            deposition_id = matching_deposits[0]['id']
            print('A deposit with title "{0}"'.format(title))
            if skip_existing:
                print('has already been started with id {0}.'.format(deposition_id))
                if error_on_existing:
                    raise ValueError(title)
                else:
                    return
            print('has already been started with this login.  Opening it for editing.')
            d = l.deposit(deposition_id, ignore_deletion=ignore_deletion)
        elif len(matching_deposits) > 1:
            print('Multiple deposits titled "{0}" have been found.'.format(title))
            raise ValueError(title)
        elif len(matching_deposits) == 0:
            # Check to see if this simulation is already in sxs but not owned by this login
            records = Records.search('title: "{0}"'.format(title))
            records = [r for r in records if r.get('title', '') == title]  # Ensure *exact* match
            communities = [community.get('identifier', '')
                           for representation in records
                           for community in representation.get('metadata', {}).get('communities', {})]
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

    # Get the time at which this original deposit was created, for possible comparison to file times below
    deposit_created = datetime.datetime.strptime(d.representation['created'], '%Y-%m-%dT%H:%M:%S.%f%z')  # in UTC

    # Convert each metadata.txt file to a metadata.json file sorted with interesting stuff at the
    # top of the file, so it appears prominently on Zenodo's preview without scrolling.  Do this
    # before checking for new files in case these are new or get changed in the process.
    paths_and_names = find_files(directory, exclude=exclude, include_top_directory_in_name=False)
    authors_emails = set()
    point_of_contact_email = ''
    keywords = set(keywords)
    related_identifiers = set(related_identifiers)
    for path,_ in paths_and_names:
        if os.path.basename(path) == 'metadata.txt':
            json_path = os.path.join(os.path.dirname(path), 'metadata.json')
            print('Converting metadata.txt to JSON in {0}'.format(json_path))
            m = Metadata.from_txt_file(path, cache_json=False).add_extras().reorder_keys()
            del m['metadata_path']  # Don't bother recording the local path to the metadata file
            m.to_json_file(json_path)
            authors_emails |= set(m.get('authors_emails', []))
            point_of_contact_email = m.get('point_of_contact_email', point_of_contact_email)
            keywords |= set(m.get('keywords', []))
            bibtex_keys = m.get('bibtex_keys', m.get('simulation_bibtex_keys', []))
            bib_to_id = inspire.map_bibtex_keys_to_identifiers(bibtex_keys)
            related_identifiers |= set(bib_to_id[bib] for bib in bib_to_id)

    # Get list of creators, keywords, and description
    print('Constructing zenodo metadata')
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
    keywords = sorted(set(keywords) | set(d.metadata.get('keywords', [])))
    if not description:
        description = d.metadata.get('description', '')
        if not description:
            spec_url = "https://www.black-holes.org/code/SpEC.html"
            description = default_description.format(spec_url)
    communities = d.metadata.get('communities', [])
    if 'sxs' not in [c['identifier'] for c in communities]:
        communities.append({'identifier': 'sxs'})
    related_identifiers = [related_identifier_formatter(rid) for rid in sorted(related_identifiers)]
    print('Finished constructing zenodo metadata')

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
    if related_identifiers:
        new_metadata['related_identifiers'] = related_identifiers
    metadata = d.metadata.copy()
    metadata.update(new_metadata)  # Ensure that fields we haven't changed are still present
    unchanged_metadata = (metadata == d.metadata)
    if unchanged_metadata:
        print('No zenodo metadata changed.')
    else:
        try:
            d.edit()
        except:
            print("Failure to unlock deposit is probably acceptable")
        d.update_metadata(metadata)
        print('Uploaded metadata')

    # Get the list of files we'll be uploading and compare to files already in the deposit to see if
    # any have changed.  If so, we need to create a new version.  Otherwise, we can just edit this
    # version.
    zenodo_filenames = d.file_names
    zenodo_name_prefix = ''
    for file_name in zenodo_filenames:
        if file_name.startswith(sxs_system + '/'):
            zenodo_name_prefix = sxs_system + '/'
            break
    local_paths_and_names = find_files(directory, exclude=exclude, include_top_directory_in_name=False)
    if zenodo_name_prefix:
        local_paths_and_names = [[path, zenodo_name_prefix+name] for path, name in local_paths_and_names]
    if len(local_paths_and_names) == 0:
        print('Zenodo requires that there be at least one file.  None found in {0}.'.format(directory))
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
    zenodo_names_to_upload_or_replace = [zf for zf in zenodo_filenames if zf in name_to_path_map]
    names_to_upload = sorted(set(name_to_path_map) - set(zenodo_filenames))
    names_to_replace = sorted(set(name_to_path_map) & set(zenodo_filenames))
    paths_and_names_to_upload = [[name_to_path_map[name], name] for name in names_to_upload]
    paths_and_names_to_replace = [[name_to_path_map[name], name] for name in names_to_replace]

    # Now, if needed do the file deletions and/or uploads, and publish
    if not names_to_delete and not names_to_upload and not names_to_replace and unchanged_metadata:
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
    d.refresh_information

    # Publish this version
    if publish:
        if not d.published:
            d.publish()
            print('Finished publishing "{0}" to {1}.'.format(title, d.website))
        else:
            print('Already published "{0}" to {1}.'.format(title, d.website))
    else:
        print('Skipping publication, as requested.')

    return d
