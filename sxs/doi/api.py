class DOI(object):
    resourceTypeGeneral_options = [  
        "Audiovisual", "Collection", "Data", "Paper", "Dataset", "Event", "Image", "InteractiveResource",
        "Model", "PhysicalObject", "Service", "Software", "Sound", "Text", "Workflow", "Other"
    ]

    def _validate_resourceTypeGeneral(resourceTypeGeneral_string):
        import textwrap
        if resourceTypeGeneral_string in resourceTypeGeneral_options:
            return
        message = textwrap.fill(
            f"The input 'resourceTypeGeneral' string, '{resourceTypeGeneral_string}' is not allowed.  "
            f"These strings must come from the following list: {str(resourceTypeGeneral_options)}.  "
            "This list is determined by the DataCite Metadata Schema "
            "<https://schema.datacite.org/meta/kernel-4.3/>."
        )
        raise ValueError(message)

    def __init__(self, information={}):
        import copy
        self._info = copy.deepcopy(information)

    @property
    def data(self):
        return self._info["data"]

    @property
    def id(self):
        return self.data.get("id", "")

    @property
    def attributes(self):
        return self.data.get("attributes", {})

    @property
    def doi(self):
        return self.attributes.get("doi", "")

    @property
    def prefix(self):
        return self.attributes.get("prefix", "")

    @property
    def suffix(self):
        return self.attributes.get("suffix", "")

    @property
    def url(self):
        return self.attributes.get("url", "")

    @url.setter
    def url(self, url_string):
        self.attributes["url"] = str(url_string)

    @property
    def creators(self):
        return self.attributes.get("creators", [])

    @creators.setter
    def creators(self, creator_list):
        self.attributes["creators"] = creator_list

    @property
    def title(self):
        return self.attributes.get("titles", [{}])[0].get("title", "")

    @title.setter
    def title(self, title_string):
        self.attributes["titles"] = [{"title": str(title_string)}]

    @property
    def publisher(self):
        return self.attributes.get("publisher", "")

    @publisher.setter
    def publisher(self, publisher_string):
        self.attributes["publisher"] = str(publisher_string)

    @property
    def publicationYear(self):
        return self.attributes.get("publicationYear", 0)

    @publicationYear.setter
    def publicationYear(self, year):
        self.attributes["publicationYear"] = int(year)

    @property
    def resourceTypeGeneral(self):
        return self.attributes.get("types", {}).get("resourceTypeGeneral", "")

    @resourceTypeGeneral.setter
    def resourceTypeGeneral(self, resourceTypeGeneral_string):
        resourceTypeGeneral_string = str(resourceTypeGeneral_string)
        _validate_resourceTypeGeneral(resourceTypeGeneral_string)
        self.attributes["types"]["resourceTypeGeneral"] = resourceTypeGeneral_string

    @property
    def state(self):
        # Note that this cannot be set directly; use the `draft`, `register`, `publish`,
        # or `hide` methods, along with a Login object, to change state.
        return self.attributes.get("state", "")

    def draft(self, login):
        return login.draft(self)

    def register(self, login):
        return login.register(self)

    def publish(self, login):
        return login.publish(self)

    def hide(self, login):
        return login.hide(self)


class Login(object):

    def __init__(self, username, password=None, password_path='~/.credentials/datacite/password',
                 test=False, total_retry_count=50, backoff_factor=0.1, session=None):
        """Initialize a Login object for interacting with the DataCite REST API"""
        import os
        import requests
        from requests.adapters import HTTPAdapter
        from requests.packages.urllib3.util.retry import Retry

        self.base_url = (
            "https://api.test.datacite.org/"
            if test else
            "https://api.datacite.org/"
        )

        # The Session object will handle all requests we make.
        self.session = session or requests.Session()

        # If the input session object succeeds at a `get`, we can skip a lot of the following
        r = self.session.get(f"{self.base_url}dois", headers={'accept': 'application/vnd.api+json'},
                             params={"page[size]":"1"})
        if r.status_code != 200:  # That's okay, we mostly expected this

            # Set the API access token
            if password is None:
                path = os.path.expanduser(password_path)
                try:
                    with open(path, 'r') as f:
                        password = f.readline().strip()
                except IOError:
                    print('Unable to find the API access token needed to create DOIs.')
                    print(f'Failed to open file "{path}" for reading.')
                    raise
                if not password:
                    print(f'The file "{path}" did not contain any text on the first line.')
                    print('This is should be an API access token, which is need to make a Deposit.')
                    raise ValueError('Deposit requires an API access token')

            # Ensure that this session sends the Authorization header with every request to the base_url
            class DataCiteAuth(requests.auth.AuthBase):
                # This is like the HTTPBasicAuth class, but checks to make sure we're only sending
                # auth to the base url.
                def __init__(self, base_url, username, password):
                    from requests.auth import _basic_auth_str
                    self.base_url = base_url
                    self.auth_str = _basic_auth_str(username, password)
                    super(InvenioRDMAuth, self).__init__()
                def __call__(self, r):
                    if r.url.startswith(self.base_url):
                        r.headers["Authorization"] = self.auth_str
                    return r
            self.session.auth = DataCiteAuth(self.base_url, username, password)

        # Note that some requests require different choices for 'Accept' and 'Content-Type'; these
        # are altered in the corresponding methods below.
        default_headers= {
            "accept": "application/vnd.api+json",
            "content-type": "application/vnd.api+json",
        }
        self.session.headers.update(default_headers)

        ## Retry automatically on certain types of errors
        retry = Retry(
            total=total_retry_count,
            backoff_factor=backoff_factor,
            status_forcelist=[500, 502, 503, 504,],
        )
        adapter = HTTPAdapter(max_retries=retry)
        self.session.mount(self.base_url, adapter)

        # Test to see if we can use the given access token
        r = self.session.get(f"{self.base_url}dois", params={"page[size]":"1"})
        if r.status_code != 200:
            if r.status_code == 401:
                print(f'The given Datacite access token was not accepted by {self.base_url}.'
                      'Please ensure that it is still valid.')
            else:
                print(f'An unknown error occurred when trying to access {self.base_url}.')
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error

    @property
    def get(self):
        return self.session.get

    @property
    def post(self):
        return self.session.post

    @property
    def put(self):
        return self.session.put

    @property
    def delete(self):
        return self.session.delete

    def list(self, **kwargs):
        """Return a list of dois

        The following query params can be given as keyword arguments to this function:

        query: string
            See <https://support.datacite.org/docs/api-queries>
        created: float
            metadata where year of DOI creation is {created}
        registered: float
            metadata where year of DOI registration is {year}
        provider-id: string
            metadata associated with a specific DataCite provider
        client-id: string
            metadata associated with a specific DataCite client
        person-id: string
            metadata associated with a specific person's ORCID iD
        resource-type-id: string
            metadata for a specific resourceTypeGeneral
        subject: string
            subject, keyword, classification code, or key phrase describing the resource
        schema-version: string
            metadata where schema version of the deposited metadata is {schema-version}
        random: boolean
            if "true", return random subset of results
        sample-size: float
            number of random samples to return for each sample-group
        sample-group: string
            return random samples grouped by "client", "provider", or "resource-type"
        page[number]: float
            page number
        page[size]: float
            results per per page
        page[cursor]: float
            see <https://support.datacite.org/docs/pagination>
        include: string
            side-load associations
        sort: string
            sort results by a certain field

        """
        r = self.get(f"{self.base_url}dois")
        if r.status_code != 200:
            print(f"The input query params resulted in an error")
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        return r.json()

    def update(self, doi):
        doi = DOI(doi)
        r = self.put(f"{self.base_url}dois/{doi.id}", data=doi)
        if r.status_code != 200:
            print(f"The input query params resulted in an error")
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        return r

    def draft(self, doi):
        self.update(doi)

    def register(self, doi):
        doi = DOI(doi._info)
        doi.attributes["event"] = "register"
        self.update(doi)

    def publish(self, doi):
        doi = DOI(doi._info)
        doi.attributes["event"] = "publish"
        self.update(doi)

    def hide(self, doi):
        doi = DOI(doi._info)
        doi.attributes["event"] = "hide"
        self.update(doi)
