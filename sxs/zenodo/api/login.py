class Login(object):

    def __init__(self, sandbox=False, access_token=None, access_token_path=None, session=None):
        """Initialize a Login object for interacting with zenodo

        This object encapsulates the credentials needed to interact with the Zenodo API, and exposes
        a Session object that can be used to make requests which automatically include the
        credentials.  It can be used for generic requests, but note that other objects in this
        module make certain tasks easier -- such as creating or modifying a "deposit", which is
        Zenodo's name for a new upload.  The Deposit object should be created from this object.

        These actions require a Zenodo API access token.  These are obtained from the website --
        either
            https://zenodo.org/account/settings/applications/tokens/new/
        or
            https://sandbox.zenodo.org/account/settings/applications/tokens/new/
        Note that zenodo.org and sandbox.zenodo.org use separate login systems and separate access
        tokens.  The access token may either be passed as a string to this function (though that
        means it will probably be found in some script file somewhere, which is probably not a good
        idea for security), or can be read from a file.  By default the file from which the token is
        read is '~/.credentials/zenodo/access_token' or the same name with '_sandbox' appended.
        Thus, it is probably easiest to simply place your access tokens in those files, so that no
        arguments need to be passed to this function.  As a basic security measure, please ensure
        that those files are not readable by anyone but the user.


        Parameters
        ==========
        sandbox: bool [default: False]
            If True, use the zenodo sandbox site, which is intended solely for testing purposes.
            This site is cleaned out regularly, so you cannot expect any entry to be here for very
            long.

        access_token: string or None [default: None]
            If present this is used as the Zenodo API access token.

        access_token_path: string or None [default: None]
            If `access_token` is not given, this file is read and the first line is used as the
            access token.  If this argument is None, it defaults to either
            '~/.credentials/zenodo/access_token' for the regular website or
            '~/.credentials/zenodo/access_token_sandbox' for the sandbox website.

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
        from . import url_sandbox, url_standard

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
                    print('Unable to find the Zenodo access token needed to make a deposit.')
                    print('Failed to open file "{0}" for reading.'.format(path))
                    raise
                if not self.access_token:
                    print('The file "{0}" did not contain any text on the first line.'.format(path))
                    print('This is should be a Zenodo access token, which is need to make a Deposit.')
                    raise ValueError('Deposit requires a Zenodo access token')

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
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
        elif r.status_code != 200:
            print('An unknown error occurred when trying to access {0}.'.format(self.base_url))
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error

    @property
    def new_deposit(self):
        """Create a new Deposit object using this login"""
        return self.deposit()

    def deposit(self, deposition_id=None):
        """Retrieve a deposit created with this login"""
        from .deposit import Deposit
        return Deposit(self, deposition_id)

    def list_deposits(self, q=None, status=None, sort=None, page=None, size=None):
        """Return list of dictionaries describing each deposit created with this login

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
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        return r.json()

    def community_curate_accept(community_id, record_id):
        url = "https://zenodo.org/communities/{0}/curate/".format(community_id)
        data = {"action": "accept", "recid": record_id}
        r = self.session.post(url, data=data)
        if r.status_code != 201:
            print('Unable to accept record id {0} into community {1}.'.format(record_id, community_id))
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        return r.json()
