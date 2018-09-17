class Login(object):

    def __init__(self, sandbox=False, access_token=None, access_token_path=None,
                 total_retry_count=50, backoff_factor=0.1, backoff_max=20.0, session=None):
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

        total_retry_count: int [default: 50]
            Total number of times to retry requests that fail for retry-able reasons.

        backoff_factor: float [default: 0.1]
            A delay factor to apply between requests after the second try (most errors are resolved
            immediately by a second try without a delay).  After a certain number of total retries,
            the request Session will sleep for:

                {backoff factor} * (2 ^ ({number of total retries} - 1))

            seconds before trying again. For example, if the `backoff_factor` is 0.1, then the
            session will sleep for [0.0s, 0.2s, 0.4s, ...] between retries.  It will never be longer
            than `backoff_max`.

        backoff_max: float [default: 20.0]
            Longest time (in seconds) to wait between retries.

        session: requests.Session or None [default: None]
            This is the object that handles all of the requests made to the API.  If `None`, a
            Session is created for you, and sensible default headers (including the access token)
            are created.  If you need to adjust some of the Session parameters like proxies or SSL
            verification, you can simply create your own and pass it in here.  Note that any `auth`
            property on the passed object will be replaced by one that adds this to the header of
            each request to the chosen Zenodo domain:
                {"Authorization": "Bearer <YourAccessTokenHere>"}

        """
        import os
        import requests
        from requests.adapters import HTTPAdapter
        from requests.packages.urllib3.util.retry import Retry
        from . import url_sandbox, url_standard

        self.sandbox = sandbox
        if self.sandbox:
            self.base_url = url_sandbox
        else:
            self.base_url = url_standard

        # The Session object will handle all requests we make.
        self.session = session or requests.Session()

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

            # Ensure that this session sends the Authorization header with every request to the base_url
            class ZenodoAuth(requests.auth.AuthBase):
                def __init__(self, base_url, access_token):
                    self.base_url = base_url
                    self.access_token = access_token
                    super(ZenodoAuth, self).__init__()
                def __call__(self, r):
                    if r.url.startswith(self.base_url):
                        r.headers.update({"Authorization": "Bearer {0}".format(self.access_token)})
                    return r
            self.session.auth = ZenodoAuth(self.base_url, self.access_token)

        # Note that some requests require different choices for 'Accept' and 'Content-Type'; these
        # are altered in the corresponding methods below.
        default_headers= {
            "Accept": "application/json",
            "Content-Type": "application/json",
        }
        self.session.headers.update(default_headers)

        ## Retry automatically on certain types of errors
        Retry.BACKOFF_MAX = backoff_max  # Must be set on the class, not the instance
        retry = Retry(
            total=total_retry_count,
            backoff_factor=backoff_factor,
            status_forcelist=[500, 502, 503, 504,],
        )
        adapter = HTTPAdapter(max_retries=retry)
        self.session.mount(self.base_url, adapter)

        # Test to see if we can use the given access token
        url = "{0}api/deposit/depositions".format(self.base_url)
        r = self.session.get(url)
        if r.status_code != 200:
            if r.status_code == 401:
                print('The given Zenodo access token was not accepted by {0}.  Please ensure that it is still valid.'.format(self.base_url))
                print('Also note that zenodo.org and sandbox.zenodo.org use separate logins and separate access tokens.')
            else:
                print('An unknown error occurred when trying to access {0}.'.format(self.base_url))
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error

    def download(self, url, path):
        """Download large file efficiently

        Parameters
        ==========
        url: string
            The URL to download from.  Redirects are followed.
        path: string
            Relative or absolute path to the file in which the download will be stored.  If this is
            an existing directory or ends in a path separator, the "path" component of the URL will
            be used as the file name, and the full directory path will be created.

        """
        from shutil import copyfileobj
        from os import makedirs
        from os.path import split, exists, join, isdir
        from functools import partial
        try:
            from urllib.parse import urlparse
        except ImportError:
            from urlparse import urlparse
        url_path = urlparse(url).path
        if isdir(path):
            path = join(path, url_path[1:])
            directory, filename = split(path)
            if not exists(directory):
                makedirs(directory)
            local_filename = join(directory, filename)
        else:
            directory, filename = split(path)
            if not exists(directory):
                makedirs(directory)
            if not filename:
                filename = url_path
            local_filename = join(directory, filename)
        r = self.session.get(url, stream=True, allow_redirects=True)
        if r.status_code != 200:
            print('An error occurred when trying to access <{0}>.'.format(url))
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        r.raw.read = partial(r.raw.read, decode_content=True)
        # r.raw.decode_content = True
        with open(local_filename, 'wb') as f:
            copyfileobj(r.raw, f)
        return local_filename

    @property
    def new_deposit(self):
        """Create a new Deposit object using this login"""
        return self.deposit()

    def deposit(self, deposition_id=None, ignore_deletion=False):
        """Retrieve a deposit created with this login"""
        from .deposit import Deposit
        return Deposit(self, deposition_id, ignore_deletion)

    def list_deposits(self, q=None, status=None, sort=None, page=None, size=9999, all_versions=False):
        """Return list of dictionaries describing each deposit created with this login

        It is possible to filter the results using the optional parameters.

        Optional parameters
        ===================
        q: string
            Search query, using Elasticsearch query string syntax.  See
            https://help.zenodo.org/guides/search/ for details.
        status: string [default: None]
            Filter result based on deposit status (either 'draft' or 'published').  If None, don't
            filter.
        sort: string
            Sort order ('bestmatch' or 'mostrecent').  Prefix with minus to change form ascending to
            descending (e.g., '-mostrecent').
        page: int
            Page number for pagination
        size: int
            Number of results to return per page.  Note that Zenodo (as of this writing) seems to
            place a hard limit of 9999 responses.  Anything more will result in an error.  Use
            multiple pages to get more results.
        all_versions: bool [defaults to False]
            If True return all records, including older versions of published records.

        """
        params={}
        if q is not None and q:
            params['q'] = q
        if status is not None:
            params['status'] = status
        if sort is not None:
            params['sort'] = sort
        if page is not None:
            params['page'] = page
        if size is not None:
            params['size'] = size
        if all_versions:
            params['all_versions'] = ''
        
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

    def delete_untitled_empty_deposits(self):
        deleted_deposits = 0
        deposits = self.list_deposits(size=9999)
        for d in deposits:
            try:
                if d['title'] == '':
                    d = l.deposit(d['id'], ignore_deletion=True)
                    if not d.files:
                        d.delete_deposit(confirmed=True)
                        deleted_deposits += 1
            except:
                pass
        print('Deleted {0} deposits'.format(deleted_deposits))

    def discard_all_drafts(self):
        discarded_drafts = 0
        deposits = self.list_deposits(size=9999)
        for d in deposits:
            try:
                if d['state'] == 'inprogress':
                    d = l.deposit(d['id'], ignore_deletion=True)
                    d.discard()
                    discarded_drafts += 1
            except:
                pass
        print('Discarded {0} drafts'.format(discarded_drafts))

    def awaiting_approval(self, community_id):
        """List all records awaiting approval for the given community"""
        url = '{0}/api/records/?q=provisional_communities:{1}'.format(self.base_url, community_id)
        r = self.session.get(url)
        if r.status_code != 200:
            print('Unable to find any records for community {0}.'.format(community_id))
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        return r.json()

    def community_curate_accept(self, community_id, record_id):
        """Accept a record into the community"""
        import json
        import requests
        url = "{0}/communities/{1}/curaterecord/".format(self.base_url, community_id)
        data = {"recid": int(record_id), "action": "accept"}
        r = self.session.post(url, json=data)
        if r.status_code != 200:
            print('Unable to accept record id {0} into community {1}; status code={2}.'.format(record_id, community_id, r.status_code))
            try:
                r_json = r.json()
                print('Response JSON:')
                print(r_json)
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        return r.json()

    def total_deposit_size(self, deposition_id=None, human_readable=True):
        import math

        def convert_size(size_bytes):
            if size_bytes == 0:
                return "0B"
            size_name = ("B  ", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB", "YiB")
            i = int(math.floor(math.log(size_bytes, 1024)))
            p = math.pow(1024, i)
            s = round(size_bytes / p, 3)
            return "{0:>8.3f} {1}".format(s, size_name[i])

        if deposition_id is None:
            depositions = self.list_deposits(page=1, size=9999)
        else:
            depositions = [self.deposit(deposition_id, ignore_deletion=True).representation]
        total_size = 0
        for deposition in depositions:
            id = deposition['id']
            d = self.deposit(id, ignore_deletion=True)
            d_total_size = sum([f['filesize'] for f in d.files])
            print('{1} in "{2}" (Zenodo ID {0})'.format(id, convert_size(d_total_size), d.title))
            total_size += d_total_size
        print('{0} in {1} deposits'.format(convert_size(total_size), len(depositions)))
        if human_readable:
            return convert_size(total_size)  # Note: the return type will be str
        else:
            return total_size  # Note: the return type will be int
