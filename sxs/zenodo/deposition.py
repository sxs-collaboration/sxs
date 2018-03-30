class Deposition(object):
    from http.client import responses

    def __init__(self, login, deposition_id=None):
        """Initialize a Deposition object for creating a new zenodo entry

        This object encapsulates all the actions you might want to take when creating, publishing,
        updating, or replacing an entry in zenodo.

        Parameters
        ==========
        login: Login object
            Object containing a `requests.Session` object that can successfully interact with the
            Zenodo API.  See the help string for `Login` in this module.

        deposition_id: string, int, or None [default: None]
            If present, this is used as the id of the deposition to edit.  If `None`, a new
            deposition is created.

        """
        self.login = login

        # Now, create or reacquire the specific deposition we're looking for
        if deposition_id is not None:
            # If the deposition id was given, check that we can access it
            url = "{0}api/deposit/depositions/{1}".format(self.base_url, deposition_id)
            r = self.login.session.get(url)
            if r.status_code != 200:
                print('The input deposition id "{0}" could not be accessed on {1}.'.format(self.deposition_id, url))
                print('The returned HTTP status code was "{0} {1}".'.format(r.status_code, responses[r.status_code]))
                r.raise_for_status()
                raise RuntimeError()  # Will only happen if the response was not strictly an error
        else:
            url = "{0}api/deposit/depositions".format(self.base_url, self.deposition_id), data="{}"
            r = self.login.session.post(url)
            if r.status_code != 201:
                print('Unable to create a new deposition on {0}.'.format(url))
                print('The returned HTTP status code was "{0} {1}".'.format(r.status_code, responses[r.status_code]))
                r.raise_for_status()
                raise RuntimeError()  # Will only happen if the response was not strictly an error

        # Now, using the response generated above, set some data describing this deposition
        r_json = r.json()
        self._deposition_id = r_json['id']
        self._links = r_json['links']
        self._state = r_json['state']
        self._submitted = bool(r_json['submitted'])

    @property
    def id(self):
        return self._deposition_id

    @property
    def deposition_id(self):
        return self._deposition_id
    
    @property
    def links(self):
        return self._links
    
    @property
    def state(self):
        return self._state
    
    @property
    def published(self):
        return self._submitted
    
    def publish(self):
        """Publish this deposition on Zenodo.

        Note that you will not be able to change the files after publishing, unless you create a new
        version of this deposition, which will result in a new DOI -- though anyone looking for this
        deposition will see a notice that there is a newer version.  You will still be able to edit
        the metadata (including the description) without changing the DOI.

        """
        url = '{0}api/deposit/depositions/{1}/actions/publish'.format(self.login.base_url, deposition_id)
        r = self.login.session.post(url)
        if r.status_code != 202:
            print('Publishing deposition {0} failed.'.format(self.deposition_id))
            print('The returned HTTP status code was "{0} {1}".'.format(r.status_code, responses[r.status_code]))
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        r_json = r.json()
        self._state = r_json['state']
        self._submitted = bool(r_json['submitted'])
        return r

    def delete(self, confirmed=False):
        """Permanently delete this deposition from Zenodo

        If you are sure that you want to delete this deposition, you can pass `confirmed=True`.
        Otherwise, you will be prompted to input 'yes' in order to complete the deletion.  This
        prompt will only last for 10 seconds.  If no input is given, we assume that the answer is
        'no', a warning is printed, and the program continues.

        """
        raise NotImplementedError()

        url = '{0}api/deposit/depositions/{1}'.format(self.login.base_url, deposition_id)
        r = self.login.session.delete(url)
        # TODO: Check if this really should be 204; the documentation is contradictory, and says '201 Created' somewhere else
        if r.status_code != 204:
            print('Deleting deposition {0} failed.'.format(self.deposition_id))
            print('The returned HTTP status code was "{0} {1}".'.format(r.status_code, responses[r.status_code]))
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        r_json = r.json()
        self._state = r_json['state']
        self._submitted = bool(r_json['submitted'])
        return r

    def update(self, data):
        """Update this deposition with the given data

        The `data` argument should be a dictionary representing the metadata

        Metadata keys and allowed values
        ================================
        upload_type: string, required
            * 'publication'
            * 'poster'
            * 'presentation'
            * 'dataset'
            * 'image'
            * 'video'
            * 'software'

        publication_type: string, required if 'upload_type' is 'publication'
            * 'book'
            * 'section'
            * 'conferencepaper'
            * 'article'
            * 'patent'
            * 'preprint'
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

        title: string, required

        creators: list of dictionaries, required
            Each entry should be a dictionary of these entries:
                * name: Name of creator in the format Family name, Given names
                * affiliation: Affiliation of creator (optional).
                * orcid: ORCID identifier of creator (optional).
                * gnd: GND identifier of creator (optional).

        description: string, required
            Abstract or description for deposition. The following HTML tags are allowed: a, p, br,
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
            Free form keywords for this deposition.

        notes: string
            Additional notes. No HTML allowed.

        communities: list of dictionaries
            List of communities in which you wish the deposition to appear. The owner of the
            community will be notified, and can either accept or reject your request. Each array
            element is a dictionary with the key 'identifier' followed by the name of a Zenodo
            community.

        doi: string
            Digital Object Identifier. Did a publisher already assign a DOI to your deposited files?
            If not, leave the field empty and we will register a new DOI for you when you publish. A
            DOI allow others to easily and unambiguously cite your deposition.

        prereserve_doi: bool
            Set to `True`, to reserve a Digital Object Identifier (DOI). The DOI is automatically
            generated by our system and cannot be changed. Also, The DOI is not registered with
            DataCite until you publish your deposition, and thus cannot be used before
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
            The contributors of the deposition (e.g., editors, data curators, etc.). Each array
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
            deposition. Each array element is an object with the key 'id' giving the grant ID.

        journal_title: string
            Journal title, if deposition is a published article.

        journal_volume: string
            Journal volume, if deposition is a published article.

        journal_issue: string
            Journal issue, if deposition is a published article.

        journal_pages: string
            Journal pages, if deposition is a published article.

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
        url = "{0}api/deposit/depositions/{1}".format(self.login.base_url, self.deposition_id)
        r = self.login.session.put(url, data=data)
        if r.status_code != 200:
            print('Updating deposition {0} failed.'.format(self.deposition_id))
            print('The returned HTTP status code was "{0} {1}".'.format(r.status_code, responses[r.status_code]))
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        r_json = r.json()
        self._state = r_json['state']
        self._submitted = bool(r_json['submitted'])
        return r

    def upload_file(self, path, name=None, relpath_start=os.curdir):
        """Upload a single file to the deposition

        Parameters
        ==========
        path: string
            Relative or absolute path to the file to upload
        name: string or None [default: None]
            Name of the file as it should appear in the deposition.  This can be the same as or
            different from the path, and can contain directories.  If this is None, the name will be
            the relative path from `relpath_start` to the file itself.  Note that if the absolute
            path to `relpath_start` is not contained within the absolute path to `path`, then that
            name would start with '../' or something.  All such parts are removed from the `name`.
        relpath_start: string [default: current directory]
            Relative or absolute path at which to start the relative path to the `name` if required.

        """
        if name is None:
            abspath = os.path.abspath(path)
            absstart = os.path.abspath(relpath_start)
            name = os.path.relpath(abspath, absstart)
            pardir_sep = os.pardir + os.sep
            while name.startswith(pardir_sep):
                name = name[len(pardir_sep):]
        url = '{0}/{1}'.format(self.links['bucket_url'], name)
        r = requests.put(url, data=open(path, 'rb'),  headers={"Content-Type":"application/octet-stream"})
        if r.status_code != 200:
            print('Uploading {0} to deposition {1} failed.'.format(path, self.deposition_id))
            print('Upload url was {0}.'.format(url))
            print('The returned HTTP status code was "{0} {1}".'.format(r.status_code, responses[r.status_code]))
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        r_json = r.json()
        self._state = r_json['state']
        self._submitted = bool(r_json['submitted'])
        return r

    def __del__(self):
        if not self.published:
            from textwrap import dedent
            from warnings import warn

            self.sandbox, self.access_argument
            
            warning = r"""\
            The Zenodo Deposition object has not been published.  The deposition id is '{deposition_id}'.
            If you want to publish this deposition, you can do it manually from the website by
            going to

                {base_url}deposit/{deposition_id}

            filling in any remaining information that is needed, and clicking "Publish".  If you do not
            want to publish it at any point, you can go to that page and click "Delete".

            You can also perform these actions using this interface by running

                >>> from sxs.zenodo import Login
                >>> l = Login(<YourLoginInfo>)
                >>> d = l.deposition({deposition_id})
                >>> # d.upload_file(...), d.update(...), etc.
                >>> d.publish()

            or

                >>> d.delete()

            """.format(deposition_id=self.deposition_id, base_url=self.l.base_url)
            warn(dedent(warning))
