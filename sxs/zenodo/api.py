url_standard = 'https://zenodo.org/'
url_sandbox = 'https://sandbox.zenodo.org/'

class Login(object):

    def __init__(self, sandbox=False, access_token=None, access_token_path=None, session=None):
        """Initialize a Login object for interacting with zenodo

        This object encapsulates the credentials needed to interact with the Zenodo API.  It can be
        used for generic requests, but note that other objects in this module make certain tasks
        easier -- such as creating or modifying a "deposition", which is Zenodo's name for a new
        upload.

        Parameters
        ==========
        sandbox: bool [default: False]
            If True, use the zenodo sandbox site, which is intended solely for testing purposes.
            This site is cleaned out regularly, so you cannot expect any entry to be here for very
            long.

        access_token: string or None [default: None]
            If present this is used as the zenodo API access token.  These are obtained from the
            website -- either
                https://zenodo.org/account/settings/applications/tokens/new/
            or
                https://sandbox.zenodo.org/account/settings/applications/tokens/new/
            Note that zenodo.org and sandbox.zenodo.org use separate login systems and separate
            access tokens.

        access_token_path: string or None [default: None]
            If `access_token` is not given, this file is read and the first line is used as the
            access token.  If this argument is None, it defaults to either
            '~/.credentials/zenodo/access_token' for the regular website or
            '~/.credentials/zenodo/access_token_sandbox' for the sandbox website.  As a basic
            security measure, please ensure that these files are not readable by anyone but the
            user.

        session: requests.Session or None [default: None]
            This is the object that handles all of the requests made to the API.  If `None`, a
            Session is created for you, and sensible default headers (including the access token)
            are created.  If you need to adjust some of the Session parameters like proxies or SSL
            verification, you can simply create your own and pass it in here.  Remember to set the
            access token using a header like
                {"Authorization": "Bearer <YourAccessTokenHere>"}

        """
        import requests
        import os

        self.sandbox = sandbox
        if self.sandbox:
            self.base_url = url_sandbox
        else:
            self.base_url = url_standard

        # The Session object will handle all requests we make.
        if session is None:
            self.session = requests.Session()
        else:
            self.session = session

        # If the input session object succeeds at a `get`, we can skip a lot of the following
        r = self.session.get("{0}api/deposit/depositions".format(self.base_url))
        if r.status_code != 200:  # That's okay, we mostly expected this

            # Set the Zenodo API access token
            if access_token is not None:
                self.access_token = access_token
                self.access_argument = 'access_token=<YourTokenHere>'
            else:
                if access_token_path is None:
                    if self.sandbox:
                        access_token_path = '~/.credentials/zenodo/access_token_sandbox'
                    else:
                        access_token_path = '~/.credentials/zenodo/access_token'
                path = os.path.expanduser(access_token_path)
                try:
                    with open(path, 'r') as f:
                        self.access_token = f.readline().strip()
                    self.access_argument = "access_token_path='{0}'".format(access_token_path)
                except IOError:
                    print('Unable to find the Zenodo access token needed to make a Deposition.')
                    print('Failed to open file "{0}" for reading.'.format(path))
                    raise
                if not self.access_token:
                    print('The file "{0}" did not contain any text on the first line.'.format(path))
                    print('This is should be a Zenodo access token, which is need to make a Deposition.')
                    raise ValueError('Deposition requires a Zenodo access token')

            # Ensure that, by default, this session sends the Authorization header
            self.session.headers.update({"Authorization": "Bearer {0}".format(self.access_token)})

        # Note that some requests require different choices for 'Accept' and 'Content-Type'; these
        # are altered in the corresponding methods below.
        default_headers= {
            "Accept": "application/json",
            "Content-Type": "application/json",
        }
        self.session.headers.update(default_headers)

        # Test to see if we can use the given access token
        url = "{0}api/deposit/depositions".format(self.base_url)
        r = self.session.get(url)
        if r.status_code == 401:
            print('The given Zenodo access token was not accepted by {0}.  Please ensure that it is still valid.'.format(self.base_url))
            print('Also note that zenodo.org and sandbox.zenodo.org use separate logins and separate access tokens.')
            print(r.json())
            r.raise_for_status()
        elif r.status_code != 200:
            print('An unknown error occurred when trying to access {0}.'.format(self.base_url))
            print(r.json())
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error

    @property
    def new_deposition(self):
        """Create a new Deposition object using this login"""
        return self.deposition()

    def deposition(self, deposition_id=None):
        """Retrieve a deposition created with this login"""
        return Deposition(self, deposition_id)

    @property
    def list_depositions(self, q=None, status=None, sort=None, page=None, size=None):
        """Return list of dictionaries describing each deposition created with this login

        It is possible to filter the results use the optional parameters

        Optional parameters
        ===================
        q: string
            Search query (using Elasticsearch query string syntax)

        status: string
            Filter result based on deposit status (either 'draft' or 'published')

        sort: string
            Sort order ('bestmatch' or 'mostrecent').  Prefix with minus to change form ascending to
            descending (e.g., '-mostrecent').

        page: int
            Page number for pagination

        size: int
            Number of results to return per page

        """
        params={}
        if q is not None:
            params['q'] = q
        if status is not None:
            params['status'] = status
        if sort is not None:
            params['sort'] = sort
        if page is not None:
            params['page'] = page
        if size is not None:
            params['size'] = size
        
        url = "{0}api/deposit/depositions".format(self.base_url)
        r = self.session.get(url, params=params)
        if r.status_code != 200:
            print('An unknown error occurred when trying to access {0}.'.format(url))
            print('The search parameters were "{0}"'.format(params))
            print(r.json())
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        return r.json()


class Deposition(object):

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
            r = self.get(url)
            if r.status_code != 200:
                print('The input deposition id "{0}" could not be accessed on {1}.'.format(self.deposition_id, url))
                print(r.json())
                r.raise_for_status()
                raise RuntimeError()  # Will only happen if the response was not strictly an error
        else:
            url = "{0}api/deposit/depositions".format(self.base_url)
            r = self.post(url, data="{}")
            if r.status_code != 201:
                print('Unable to create a new deposition on {0}.'.format(url))
                print(r.json())
                r.raise_for_status()
                raise RuntimeError()  # Will only happen if the response was not strictly an error

        # Now, using the response generated above, set some data describing this deposition
        r_json = r.json()
        self._deposition_id = r_json['id']
        self._links = r_json['links']
        self._state = r_json['state']
        self._submitted = bool(r_json['submitted'])

    @property
    def base_url(self):
        return self.login.base_url

    @property
    def get(self):
        return self.login.session.get

    @property
    def post(self):
        return self.login.session.post

    @property
    def put(self):
        return self.login.session.put

    @property
    def delete(self):
        return self.login.session.delete

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
    def submitted(self):
        return self._submitted
    
    @property
    def published(self):
        return self._submitted

    @property
    def files(self):
        url = '{0}api/deposit/depositions/{1}/files'.format(self.base_url, self.deposition_id)
        r = self.get(url)
        if r.status_code != 200:
            print('Getting file listing for deposition {0} failed.'.format(self.deposition_id))
            print(r.json())
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        return r.json()

    def update(self, metadata):
        """Update this deposition with the given metadata

        The `metadata` argument should be a dictionary representing the metadata

        Metadata keys and allowed values
        ================================
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
        import json
        url = "{0}api/deposit/depositions/{1}".format(self.base_url, self.deposition_id)
        r = self.put(url, data=json.dumps({'metadata': metadata}))
        if r.status_code != 200:
            print('Updating deposition {0} failed.'.format(self.deposition_id))
            print(r.json())
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        return r

    def upload_file(self, path, name=None, relpath_start=None):
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
        import os

        if relpath_start is None:
            relpath_start = os.curdir

        if name is None:
            abspath = os.path.abspath(path)
            absstart = os.path.abspath(relpath_start)
            name = os.path.relpath(abspath, absstart)
            pardir_sep = os.pardir + os.sep
            while name.startswith(pardir_sep):
                name = name[len(pardir_sep):]
        url = '{0}/{1}'.format(self.links['bucket'], name)
        r = self.put(url, data=open(path, 'rb'),  headers={"Content-Type":"application/octet-stream"})
        if r.status_code != 200:
            print('Uploading {0} to deposition {1} failed.'.format(path, self.deposition_id))
            print('Upload url was {0}.'.format(url))
            print(r.json())
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error

        return r

    def upload_all_files(self, top_directory, exclude=[]):
        """Recursively upload all files found in `top_directory`

        Each file is named with its path relative to the parent directory of `top_directory`.

        Parameters
        ==========
        top_directory: string
            Absolute or relative path to the top directory from which the recursive search for files
            begins.
        exclude: list of strings
            Each string is compiled as a regular expression.  The path to each directory and file
            relative to `top_directory` is searched for a match, and if found that item is excluded.
            In particular, if a directory matches, no files from that directory will be uploaded.

        """
        import os
        import re
        exclude = [re.compile(exclusion) for exclusion in exclude]
        top_directory = os.path.abspath(top_directory)
        parent_directory = os.path.dirname(top_directory)
        for root, dirs, files in os.walk(top_directory, topdown=True):
            dirs.sort(key=str.lower)  # Go in case-insensitive alphabetical order
            files.sort(key=str.lower)  # Go in case-insensitive alphabetical order
            for exclusion in exclude:
                for d in dirs:
                    if exclusion.search(os.path.relpath(d, top_directory)):
                        dirs.remove(d)
                for f in files:
                    if exclusion.search(os.path.relpath(f, top_directory)):
                        files.remove(f)
            for f in files:
                path = os.path.join(root, f)
                name = os.path.relpath(path, parent_directory)
                print("Uploading\n    {0}\nas\n    {1}\n".format(path, name))
                self.upload_file(path, name=name)
                print("Upload succeeded\n")

    def publish(self):
        """Publish this deposition on Zenodo.

        Note that you will not be able to change the files after publishing, unless you create a new
        version of this deposition, which will result in a new DOI -- though anyone looking for this
        deposition will see a notice that there is a newer version.  You will still be able to edit
        the metadata (including the description) without changing the DOI.

        """
        url = '{0}api/deposit/depositions/{1}/actions/publish'.format(self.base_url, deposition_id)
        r = self.post(url)
        if r.status_code != 202:
            print('Publishing deposition {0} failed.'.format(self.deposition_id))
            print(r.json())
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
        prompt will only last for 60 seconds.  If no input is given, we assume that the answer is
        'no', a warning is printed, and the program continues.

        """
        if not confirmed:
            import sys, select
            timeout = 60
            print("Please confirm that you want to delete the deposition {0}.".format(self.deposition_id))
            print("You have {0} seconds to confirm by entering 'yes'.".format(timeout))
            i, o, e = select.select([sys.stdin], [], [], timeout)
            if not i or sys.stdin.readline().strip().lower() != 'yes':
                print('No confirmation received.  Aborting deletion of zenodo deposition {0}.'.format(self.deposition_id))
                raise RuntimeError('No confirmation')
        url = '{0}api/deposit/depositions/{1}'.format(self.base_url, deposition_id)
        r = self.delete(url)
        # TODO: Check if this really should be 204; the documentation is contradictory, and says '201 Created' somewhere else
        if r.status_code != 204:
            print('Deleting deposition {0} failed.'.format(self.deposition_id))
            print(r.json())
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
