"""Class encapsulating "deposit" objects on Zenodo

Also called "depositions", these are only accessible to the user who creates
them.  Thus, you need the correct login information to access a deposit.

"""

class Deposit(object):

    def __init__(self, login, deposition_id=None, ignore_deletion=False):
        """Initialize a Deposit object for creating a new zenodo entry

        This object encapsulates all the actions you might want to take when creating, publishing,
        updating, or replacing an entry in zenodo.

        Parameters
        ----------
        login: Login object
            Object containing a `requests.Session` object that can successfully interact with the
            Zenodo API.  See the help string for `Login` in this module.

        deposition_id: string, int, or None [default: None]
            If present, this is used as the id of the deposit to edit.  If `None`, a new
            deposit is created.

        ignore_deletion: bool [default: False]
            If True and this object is deleted before the record is published, issue a warning
            explaining that and suggesting how it might be published by the user.

        """
        self.login = login
        self.ignore_deletion = ignore_deletion
        self._representation = {}  # Just during construction, to prevent errors if construction fails

        # Now, create or reacquire the specific deposit we're looking for
        if deposition_id is not None:
            # If the deposit id was given, check that we can access it
            url = "{0}api/deposit/depositions/{1}".format(self.base_url, deposition_id)
            r = self._get(url)
            if r.status_code != 200:
                print('The input deposition id "{0}" could not be accessed on {1}.'.format(deposition_id, url))
                try:
                    print(r.json())
                except:
                    pass
                r.raise_for_status()
                raise RuntimeError()  # Will only happen if the response was not strictly an error
        else:
            url = "{0}api/deposit/depositions".format(self.base_url)
            r = self._post(url, data="{}")
            if r.status_code != 201:
                print('Unable to create a new deposit on {0}.'.format(url))
                try:
                    print(r.json())
                except:
                    pass
                r.raise_for_status()
                raise RuntimeError()  # Will only happen if the response was not strictly an error

        # Now, using the response generated above, set some data describing this deposit
        self._representation = r.json()

    @property
    def _get(self):
        return self.login.session.get

    @property
    def _post(self):
        return self.login.session.post

    @property
    def _put(self):
        return self.login.session.put

    @property
    def _delete(self):
        return self.login.session.delete

    @property
    def refresh_information(self):
        """Retrieve current information about this Deposit from Zenodo

        This function updates this object's `representation` data, which contains information like
        the id, metadata, submission status, and various links for this deposit.  That
        information is also used by many of this object's other member functions.

        """
        url = "{0}api/deposit/depositions/{1}".format(self.base_url, self.deposition_id)
        r = self._get(url)
        if r.status_code != 200:
            print('This deposit (id "{0}") could not be accessed on {1}.'.format(self.deposition_id, url))
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        self._representation = r.json()
        return self.representation

    @property
    def base_url(self):
        return self.login.base_url

    @property
    def representation(self):
        return self._representation

    @property
    def metadata(self):
        return self.representation['metadata']

    @property
    def links(self):
        return self.representation['links']
    
    @property
    def deposition_id(self):
        """Return id number of this deposit"""
        return self.id

    @property
    def id(self):
        """Return id number of this deposit"""
        return self.representation['id']

    @property
    def id_latest(self):
        """Return id number of the most recent version of this deposit"""
        return self.links['latest'].split('/')[-1]

    @property
    def id_latest_draft(self):
        """Return id number of the most recent draft (unpublished) version of this deposit

        Note: There may be no draft version, in which case this function will raise a KeyError.

        """
        return self.links['latest_draft'].split('/')[-1]

    @property
    def website(self):
        """URL of web page on which this deposit is found"""
        return self.links['html']

    @property
    def website_latest(self):
        """URL of web page on which the current version of this deposit is found"""
        return self.links['latest_html']

    @property
    def website_latest_draft(self):
        """URL of web page on which the current draft (unpublished) version of this deposit is found

        Note: There may be no draft version, in which case this function will raise a KeyError.

        """
        return self.links['latest_draft_html']

    @property
    def record(self):
        """URL of API endpoint for this deposit"""
        return self.links['record']

    @property
    def record_latest(self):
        """URL of API endpoint for most recent version of this deposit"""
        return self.links['latest']

    @property
    def record_latest_draft(self):
        """URL of API endpoint for most recent draft (unpublished) version of this deposit

        Note: There may be no draft version, in which case this function will raise a KeyError.

        """
        return self.links['latest_draft']

    @property
    def state(self):
        """Current status of this deposit

        The documentation states that this may be one of:
            * inprogress: Deposit metadata can be updated. If deposit is also unsubmitted (see
              submitted) files can be updated as well.
            * done: Deposit has been published. 
            * error: Deposit is in an error state - contact Zenodo support.
        However, 'inprogress' seems to have been replaced by 'unsubmitted'.

        """
        return self.representation.get('state', 'error')
    
    @property
    def submitted(self):
        """Return True if this deposit has been submitted/published"""
        return bool(self.representation.get('submitted', False))
    
    @property
    def published(self):
        """Return True if this deposit has been submitted/published"""
        return self.submitted and self.state=='done'

    @property
    def is_latest(self):
        """Return True if this deposit is the most recent version"""
        return (self.links.get('latest', 1) == self.links.get('record', 2))

    def get_latest(self):
        """Return a new Deposit object pointing to the latest version of this deposit
        
        Note: This deposit may already be the latest version, in which case a new object pointing
        to this deposit is returned.

        """
        return self.login.deposit(self.id_latest)

    def get_latest_draft(self):
        """Return a new Deposit object pointing to the latest draft (unpublished) version of this deposit

        Note: There may be no draft version, in which case this function will raise a KeyError.

        """
        return self.login.deposit(self.id_latest_draft)

    @property
    def versions(self):
        """Get information for all versions of this deposit

        The versions are returned as a list of record data from oldest to newest.  The only change
        that this function makes to the data is to add a version number to the metadata if it is not
        present already.

        """
        url = "{0}api/deposit/depositions".format(self.base_url)
        conceptrecid = self.representation['conceptrecid']
        params = {
            'q': 'conceptrecid:{0}'.format(conceptrecid),
            'all_versions': '',
            'sort': 'version',
        }
        r = self._get(url, params=params)
        if r.status_code != 200:
            print('The list of versions for this deposit (id "{0}") could not be accessed on {1}.'.format(self.deposition_id, url))
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        version_list = r.json()
        for v, version in enumerate(version_list, 1):
            if 'version' not in version['metadata']:
                version['metadata']['version'] = v
        return version_list

    @property
    def title(self):
        return self.representation['title']

    @property
    def files(self):
        """Return list of dictionaries describing uploaded files

        Each file is described by a dictionary containing these keys:

            * checksum (MD5 fingerprint)
            * filename
            * filesize (in bytes)
            * id (Zenodo-generated hex id of this file)
            * links
                * download
                * self

        The last item, 'links' is another dictionary giving a couple URLs; it is not clear how to
        use these URLs.  For example, trying to get the 'download' link results in a 404 message.

        """
        return self.representation['files']

    @property
    def file_checksums(self):
        return {d['filename']: d['checksum'] for d in self.files}

    @property
    def file_sizes(self):
        return {d['filename']: d['filesize'] for d in self.files}

    @property
    def file_ids(self):
        return {d['filename']: d['id'] for d in self.files}

    @property
    def file_names(self):
        return [d['filename'] for d in self.files]

    def update_metadata(self, metadata, refresh_information=True):
        """Update this deposit with the given metadata

        The `metadata` argument should be a dictionary representing the metadata.

        Note that a deposit cannot be published unless all required metadata fields are present,
        including things like 'doi', which are automatically assigned by the system.  So if you are
        changing the metadata of a deposit, it's usually simpler to copy the metadata, modify it
        as desired, and then use this function to update it on Zenodo.


        Metadata keys and allowed values
        --------------------------------
        title: string, required

        creators: list of dictionaries, required
            Each entry should be a dictionary of these entries:
                * name: Name of creator in the format Family name, Given names
                * affiliation: Affiliation of creator (optional).
                * orcid: ORCID identifier of creator (optional).
                * gnd: GND identifier of creator (optional).

        upload_type: string, required
            * 'publication'
            * 'poster'
            * 'presentation'
            * 'dataset'
            * 'image'
            * 'video'
            * 'software'
            * 'lesson'
            * 'other'

        publication_type: string, required if 'upload_type' is 'publication'
            * 'book'
            * 'section'
            * 'conferencepaper'
            * 'article'
            * 'patent'
            * 'preprint'
            * 'deliverable'
            * 'milestone'
            * 'proposal'
            * 'report'
            * 'softwaredocumentation'
            * 'thesis'
            * 'technicalnote'
            * 'workingpaper'
            * 'other'

        image_type: string, required if if 'upload_type' is 'image'
            * 'figure'
            * 'plot'
            * 'drawing'
            * 'diagram'
            * 'photo'
            * 'other'
        
        publication_date: string [defaults to current date]
            Date of publication in ISO8601 format (YYYY-MM-DD).

        description: string, required
            Abstract or description for deposit. The following HTML tags are allowed: a, p, br,
            blockquote, strong, b, u, i, em, ul, ol, li, sub, sup, div, strike.

        access_right: string [defaults to 'open']
            * 'open'
            * 'embargoed'
            * 'restricted'
            * 'closed'

        license: string [defaults to 'cc-by' for non-datasets and 'cc-zero' for datasets]

        embargo_date: date [defaults to current date]
            When the deposited files will be made automatically made publicly available by the
            system.

        access_conditions: string, required if 'access_right' is 'restricted'
            Specify the conditions under which you grant users access to the files in your
            upload. User requesting access will be asked to justify how they fulfil the
            conditions. Based on the justification, you decide who to grant/deny access. You are not
            allowed to charge users for granting access to data hosted on Zenodo. Following HTML
            tags are allowed: a, p, br, blockquote, strong, b, u, i, em, ul, ol, li, sub, sup, div,
            strike.

        keywords: list of strings
            Free form keywords for this deposit.

        notes: string
            Additional notes. No HTML allowed.

        communities: list of dictionaries
            List of communities in which you wish the deposit to appear. The owner of the
            community will be notified, and can either accept or reject your request. Each array
            element is a dictionary with the key 'identifier' followed by the name of a Zenodo
            community.

        doi: string
            Digital Object Identifier. Did a publisher already assign a DOI to your deposited files?
            If not, leave the field empty and we will register a new DOI for you when you publish. A
            DOI allow others to easily and unambiguously cite your deposit.

        prereserve_doi: bool
            Set to `True`, to reserve a Digital Object Identifier (DOI). The DOI is automatically
            generated by our system and cannot be changed. Also, The DOI is not registered with
            DataCite until you publish your deposit, and thus cannot be used before
            then. Reserving a DOI is useful, if you need to include it in the files you upload, or
            if you need to provide a dataset DOI to your publisher but not yet publish your
            dataset. The response from the REST API will include the reserved DOI.

        related_identifiers: list of dictionaries
            Persistent identifiers of related publications and datasets. Supported identifiers
            include: DOI, Handle, ARK, PURL, ISSN, ISBN, PubMed ID, PubMed Central ID, ADS
            Bibliographic Code, arXiv, Life Science Identifiers (LSID), EAN-13, ISTC, URNs and
            URLs. Each array element is an object with the attributes:
                * identifier: The persistent identifier
                * relation: Relationship. Controlled vocabulary (isCitedBy, cites, isSupplementTo,
                  isSupplementedBy, isNewVersionOf, isPreviousVersionOf, isPartOf, hasPart,
                  compiles, isCompiledBy, isIdenticalTo, isAlternateIdentifier).

        contributors: list of dictionaries
            The contributors of the deposit (e.g., editors, data curators, etc.). Each array
            element is a dictionary with the keys:
                * name: Name of creator in the format Family name, Given names
                * type: Contributor type. Controlled vocabulary (ContactPerson, DataCollector,
                  DataCurator, DataManager, Editor, Researcher, RightsHolder, Sponsor, Other)
                * affiliation: Affiliation of creator (optional).
                * orcid: ORCID identifier of creator (optional).
                * gnd: GND identifier of creator (optional).

        references: list of strings
            List of references.

        grants: list of dictionaries
            List of OpenAIRE-supported grants, which have funded the research for this
            deposit. Each array element is an object with the key 'id' giving the grant ID.

        journal_title: string
            Journal title, if deposit is a published article.

        journal_volume: string
            Journal volume, if deposit is a published article.

        journal_issue: string
            Journal issue, if deposit is a published article.

        journal_pages: string
            Journal pages, if deposit is a published article.

        conference_title: string
            Title of conference (e.g., 20th International Conference on Computing in High Energy and
            Nuclear Physics).

        conference_acronym: string
            Acronym of conference (e.g., CHEP'13).

        conference_dates: string
            Dates of conference (e.g., 14-18 October 2013). Conference title or acronym must also be
            specified if this field is specified.

        conference_place: string
            Place of conference in the format city, country (e.g., Amsterdam, The
            Netherlands). Conference title or acronym must also be specified if this field is
            specified.

        conference_url: string
            URL of conference (e.g., http://www.chep2013.org/).

        conference_session: string
            Number of session within the conference (e.g., VI).

        conference_session_part: string
            Number of part within a session (e.g., 1).

        imprint_publisher: string
            Publisher of a book/report/chapter

        imprint_isbn: string
            ISBN of a book/report

        imprint_place: string
            Place of publication of a book/report/chapter in the format city, country.

        partof_title: string
            Title of book for chapters

        partof_pages: string
            Pages numbers of book

        thesis_supervisors: array of objects
            Supervisors of the thesis. Same format as for creators.

        thesis_university: string
            Awarding university of thesis.

        subjects: array of objects
            Specify subjects from a taxonomy or controlled vocabulary. Each term must be uniquely
            identified (e.g., a URL). For free form text, use the keywords field. Each array element
            is an object with the attributes:
                * term: Term from taxonomy or controlled vocabulary.
                * identifier: Unique identifier for term.
                * scheme: Persistent identifier scheme for id (automatically detected).

        version: string
            Version of the resource. Any string will be accepted, however the suggested format is a
            semantically versioned tag (see more details on semantic versioning at semver.org)

        language: string
            Specify the main language of the record as ISO 639-2 or 639-3 code, see Library of
            Congress ISO 639 codes list.

        """
        import json
        url = "{0}api/deposit/depositions/{1}".format(self.base_url, self.deposition_id)
        r = self._put(url, data=json.dumps({'metadata': metadata}))
        if r.status_code != 200:
            print('Updating deposit {0} failed.'.format(self.deposition_id))
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        if refresh_information:
            self.refresh_information
        return r

    def edit(self, refresh_information=True):
        """Unlock a previously submitted deposit for editing."""
        url = "{0}api/deposit/depositions/{1}/actions/edit".format(self.base_url, self.deposition_id)
        r = self._post(url)
        if r.status_code != 201:
            print('Unlocking deposit {0} for editing failed.'.format(self.deposition_id))
            if r.status_code == 400:
                print('Deposit state does not allow for editing (e.g., deposits in state `inprogress`).')
            if r.status_code == 409:
                print('Deposit is in the process of being integrated.  Please wait 5 minutes before trying again.')
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        if refresh_information:
            self.refresh_information
        return r
        
    def discard(self, refresh_information=True):
        """Discard changes in the current editing session."""
        url = "{0}api/deposit/depositions/{1}/actions/discard".format(self.base_url, self.deposition_id)
        r = self._post(url)
        if r.status_code != 201:
            print('Discarding changes from the current editing session to deposit {0} failed.'.format(self.deposition_id))
            if r.status_code == 400:
                print('Deposit is not being edited.')
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        if refresh_information:
            self.refresh_information
        return r

    def get_new_version(self, ignore_deletion=True, refresh_information=True):
        """Create a new Deposit object describing a new version of this deposit.

        To publish the new version, its files must differ from all previous versions.

        This action will create a new deposit, which will be a snapshot of the current resource,
        inheriting the metadata as well as snapshot of files. The new version deposit will have a
        state similar to a new, unpublished deposit, most importantly its files will be modifiable
        as for a new deposit.

        Only one unpublished new version deposit can be available at any moment, i.e.: calling new
        version action multiple times will have no effect, as long as the resulting new version
        deposit from the first call is not published or deleted.

        """
        r = self.register_new_version(refresh_information=refresh_information)
        id_latest_draft = r.json()['links']['latest_draft'].split('/')[-1]
        return self.login.deposit(id_latest_draft, ignore_deletion=ignore_deletion)
        
    def register_new_version(self, refresh_information=True):
        """Create a new version of a deposit.

        To publish the new version, its files must differ from all previous versions.

        This action will create a new deposit, which will be a snapshot of the current resource,
        inheriting the metadata as well as snapshot of files. The new version deposit will have a
        state similar to a new, unpublished deposit, most importantly its files will be modifiable
        as for a new deposit.

        Only one unpublished new version deposit can be available at any moment, i.e.: calling new
        version action multiple times will have no effect, as long as the resulting new version
        deposit from the first call is not published or deleted.

        NOTE: The response body of this action is NOT the new version deposit, but the original
        resource. The new version deposit can be accessed through the "latest_draft" under
        "links" in the response body.

        """
        url = "{0}api/deposit/depositions/{1}/actions/newversion".format(self.base_url, self.deposition_id)
        r = self._post(url)
        if r.status_code != 201:
            print('Failed to register new version of deposit {0}.'.format(self.deposition_id))
            if r.status_code == 403:
                print('This deposit may not have been published, in which case new versions are not allowed.')
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        if refresh_information:
            self.refresh_information
        return r

    def delete_file(self, file_name, refresh_information=True):
        matching_files = [(i, f) for i, f in enumerate(self.files) if f['filename'] == file_name]
        if not matching_files:
            print('File name "{0}" not found on Zenodo.'.format(file_name))
            raise ValueError(file_name)
        file_index, file = matching_files[0]
        file_name = file['filename']
        file_id = file['id']
        url = '{0}api/deposit/depositions/{1}/files/{2}'.format(self.base_url, self.id, file_id)
        r = self._delete(url)
        if r.status_code != 204:
            print('Deleting file "{0}" from deposit "{1}" failed.'.format(file_name, self.id))
            if r.status_code == 403:
                print('Server replied with "Forbidden: Deleting an already published deposition file."')
                print('Try get_new_version, and then delete the file from that version.')
            if r.status_code == 404:
                print('File id "{0}" does not exist in deposition "{1}".'.format(file_id, self.id))
                print('Try refreshing the local information in this Deposit object, and trying again as relevant.')
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        if refresh_information:
            self.refresh_information
        else:
            self.representation['files'].pop(file_index)
        return r

    def upload_file(self, path, name=None, relpath_start=None, skip_checksum=False, refresh_information=True):
        """Upload a single file to the deposit

        The current list of files uploaded to Zenodo is checked.  If `name` is in that list, the MD5
        checksum of `path` is evaluated and compared to the MD5 checksum of the file on Zenodo.  If
        they match, the upload is skipped and `None` is returned.

        Parameters
        ----------
        path: string
            Relative or absolute path to the file to upload
        name: string or None [default: None]
            Name of the file as it should appear in the deposit.  This can be the same as or
            different from the path, and can contain directories.  If this is None, the name will be
            the relative path from `relpath_start` to the file itself.  Note that if the absolute
            path to `relpath_start` is not contained within the absolute path to `path`, then that
            name would start with '../' or something.  All such parts are removed from the `name`.
        relpath_start: string [default: current directory]
            Relative or absolute path at which to start the relative path to the `name` if required.
        skip_checksum: bool [default: False]
            If True, ignore the checksum and always upload the file.  If False, see if any file with
            exactly this name and checksum has been uploaded to Zenodo; if so, don't bother
            uploading it again.

        """
        import os
        from ...utilities import md5checksum
        if relpath_start is None:
            relpath_start = os.curdir
        if name is None:
            abspath = os.path.abspath(path)
            absstart = os.path.abspath(relpath_start)
            name = os.path.relpath(abspath, absstart)
            pardir_sep = os.pardir + os.sep
            while name.startswith(pardir_sep):
                name = name[len(pardir_sep):]
        if not skip_checksum:
            file_checksums = self.file_checksums
            if name in file_checksums:
                if md5checksum(path) == file_checksums[name]:
                    print('File {0} has already been uploaded.  Skipping this upload.'.format(name))
                    return None
        if name in self.file_ids:
            self.delete_file(name, refresh_information=False)
        url = '{0}/{1}'.format(self.links['bucket'], name)
        r = self._put(url, data=open(path, 'rb'), headers={"Content-Type": "application/octet-stream"})
        if r.status_code != 201:
            print('Uploading {0} to deposit {1} failed.'.format(path, self.deposition_id))
            print('Upload url was {0}.'.format(url))
            if r.status_code == 400:
                if os.stat(path).st_size == 0:
                    print('This file has size zero, which leads to an error response.')
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        if refresh_information:
            self.refresh_information
        return r

    def rename_file(self, old_file_name, new_file_name, refresh_information=True):
        import json
        file_ids = self.file_ids
        if old_file_name not in file_ids:
            print('File name "{0}" not found on Zenodo.'.format(old_file_name))
            raise ValueError(old_file_name)
        file_id = file_ids[old_file_name]
        url = '{0}api/deposit/depositions/{1}/files/{2}'.format(self.base_url, self.id, file_id)
        rename_data = {'name': new_file_name}
        r = self._put(url, data=json.dumps(rename_data))
        if r.status_code != 200:
            print('Renaming file "{0}" to "{1}" in deposit "{2}" failed.'.format(old_file_name, new_file_name, self.id))
            if r.status_code == 403:
                print('Server replied with "Forbidden: Renaming an already published deposition file."')
                print('Try get_new_version, and then delete the file from that version.')
            if r.status_code == 404:
                print('File id "{0}" does not exist in deposition "{1}".'.format(file_id, self.id))
                print('Try refreshing the local information in this Deposit object, and trying again as relevant.')
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        if refresh_information:
            self.refresh_information
        return r
        
    def upload_all_files(self, top_directory, exclude=[], skip_checksum=False, refresh_information=True):
        """Recursively upload all files found in `top_directory`

        This function simply calls the `sxs.utilities.find_files` function, and then calls
        `upload_file` for each one.  See those functions for explanations of the parameters to this
        function.

        """
        from ...utilities import find_files
        paths_and_names = find_files(top_directory, exclude=exclude)
        for path, name in paths_and_names:
            print("Uploading\n    {0}\nas\n    {1}".format(path, name))
            self.upload_file(path, name=name, skip_checksum=skip_checksum, refresh_information=False)
            print("Upload succeeded\n")
        if refresh_information:
            self.refresh_information

    def prereserve_doi(self, refresh_information=True):
        """Pre-reserve the DOI for an unpublished deposit, so that it does not change upon publication

        Note that this is undocumented behavior of the zenodo API, so it may change without notice;
        this function was reverse-engineered by analyzing requests made by the web interface.  The
        API itself is not explicit about how this could work, but the approach that is vaguely
        suggested does not work.

        """
        from copy import deepcopy
        metadata = deepcopy(self.metadata)
        metadata.update({'doi': metadata['prereserve_doi']['doi']})
        return self.update_metadata(metadata, refresh_information=refresh_information)
    
    def publish(self, refresh_information=True):
        """Publish this deposit on Zenodo.

        In order to be published successfully, a deposit must contain at least one file, and all
        of the required metadata fields must be present.  In particular if you have changed the
        metadata, be sure that you have not removed one of the fields that is automatically
        assigned, like 'doi' or 'prereserve_doi', etc.

        Note that you will not be able to change the files after publishing, unless you create a new
        version of this deposit, which will result in a new DOI -- though anyone looking for this
        deposit will see a notice that there is a newer version.  You will still be able to edit
        the metadata (including the description) without changing the DOI.

        """
        url = '{0}api/deposit/depositions/{1}/actions/publish'.format(self.base_url, self.deposition_id)
        r = self._post(url)
        if r.status_code != 202:
            print('Publishing deposit {0} failed.'.format(self.deposition_id))
            if r.status_code == 500:
                print('Server returned the code "500 Internal Server Error".')
                print('This can have any number of causes, but frequently it is because not')
                print('all of the required information is present.  In particular, every')
                print('deposit must contain at least one file, and all of its metadata')
                print('must be present.  See the warning in `Deposit.update_metadata`.')
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        if refresh_information:
            self.refresh_information
        return r
    
    def delete_deposit(self, confirmed=False):
        """Permanently delete this deposit from Zenodo

        If you are sure that you want to delete this deposit, you can pass `confirmed=True`.
        Otherwise, you will be prompted to input 'yes' in order to complete the deletion.  This
        prompt will only last for 60 seconds.  If no input is given, we assume that the answer is
        'no', a warning is printed, and the program continues.

        """
        if not confirmed:
            import sys, select
            timeout = 60
            print("Please confirm that you want to delete the deposit {0}.".format(self.deposition_id))
            print("You have {0} seconds to confirm by entering 'yes'.".format(timeout))
            i, o, e = select.select([sys.stdin], [], [], timeout)
            if not i or sys.stdin.readline().strip().lower() != 'yes':
                print('No confirmation received.  Aborting deletion of zenodo deposit {0}.'.format(self.deposition_id))
                raise RuntimeError('No confirmation')
        url = '{0}api/deposit/depositions/{1}'.format(self.base_url, self.deposition_id)
        r = self._delete(url)
        if r.status_code != 204:
            print('Deleting deposit {0} failed.'.format(self.deposition_id))
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        return r

    def __del__(self):
        if not self.published and not self.ignore_deletion:
            from textwrap import dedent
            from warnings import warn

            warning = r"""
            The Zenodo Deposit object has not been published.  The deposit id is '{deposition_id}'.
            If you want to publish this deposit, you can do it manually from the website by
            going to

                {base_url}deposit/{deposition_id}

            filling in any remaining information that is needed, and clicking "Publish".  If you do not
            want to publish it at any point, you can go to that page and click "Delete".

            You can also perform these actions using this interface by running

                >>> from sxs.zenodo import Login
                >>> l = Login(<YourLoginInfo>)
                >>> d = l.deposit({deposition_id})
                >>> # d.upload_file(...), d.update_metadata(...), etc.
                >>> d.publish()

            or

                >>> d.delete()

            """.format(deposition_id=self.deposition_id, base_url=self.base_url)
            warn(dedent(warning))
