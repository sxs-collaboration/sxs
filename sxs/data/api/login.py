class Login(object):

    def __init__(self, sandbox=False, access_token=None, access_token_path=None,
                 total_retry_count=50, backoff_factor=0.1, backoff_max=20.0, session=None):
        """Initialize a Login object for interacting with InvenioRDM

        This object encapsulates the credentials needed to interact with the InvenioRDM API, and
        exposes a Session object that can be used to make requests which automatically include the
        credentials.  It can be used for generic requests, but note that other objects in this
        module make certain tasks easier -- such as creating or modifying a "deposit", which is
        InvenioRDM's name for a new upload.  The Deposit object should be created from this object.

        These actions require an InvenioRDM API access token.  These are obtained from the website
        -- typically under the account/settings/applications/tokens menu.  Note that live sites and
        sandbox sites usually use separate login systems and separate access tokens.  The access
        token may either be passed as a string to this function (though that means it will probably
        be found in some script file somewhere, which is probably not a good idea for security), or
        can be read from a file.  By default the file from which the token is read is
        '~/.credentials/sxs/data_access_token' or the same name with '_sandbox' appended.  Thus, it
        is probably easiest to simply place your access tokens in those files, so that no arguments
        need to be passed to this function.  As a basic security measure, please ensure that those
        files are not readable by anyone but the user.


        Parameters
        ==========
        sandbox: bool [default: False]
            If True, use the sandbox site, which is intended solely for testing purposes.  This site
            is cleaned out regularly, so you cannot expect any entry to be here for very long.

        access_token: string or None [default: None]
            If present, this is used as the API access token.

        access_token_path: string or None [default: None]
            If `access_token` is not given, this file is read and the first line is used as the
            access token.  If this argument is None, it defaults to either
            '~/.credentials/sxs/data_access_token' for the regular website or
            '~/.credentials/sxs/data_access_token_sandbox' for the sandbox website.

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
        r = self.session.get(f"{self.base_url}api/deposit/depositions")
        if r.status_code != 200:  # That's okay, we mostly expected this

            # Set the API access token
            if access_token is not None:
                self.access_token = access_token
            else:
                if access_token_path is None:
                    if self.sandbox:
                        access_token_path = '~/.credentials/sxs/data_access_token_sandbox'
                    else:
                        access_token_path = '~/.credentials/sxs/data_access_token'
                path = os.path.expanduser(access_token_path)
                try:
                    with open(path, 'r') as f:
                        self.access_token = f.readline().strip()
                except IOError:
                    print('Unable to find the API access token needed to make a deposit.')
                    print(f'Failed to open file "{path}" for reading.')
                    raise
                if not self.access_token:
                    print(f'The file "{path}" did not contain any text on the first line.')
                    print('This is should be an API access token, which is need to make a Deposit.')
                    raise ValueError('Deposit requires an API access token')

            # Ensure that this session sends the Authorization header with every request to the base_url
            class InvenioRDMAuth(requests.auth.AuthBase):
                def __init__(self, base_url, access_token):
                    self.base_url = base_url
                    self.access_token = access_token
                    super(InvenioRDMAuth, self).__init__()
                def __call__(self, r):
                    if r.url.startswith(self.base_url):
                        r.headers.update({"Authorization": f"Bearer {self.access_token}"})
                    return r
            self.session.auth = InvenioRDMAuth(self.base_url, self.access_token)

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
        url = f"{self.base_url}api/deposit/depositions"
        r = self.session.get(url)
        if r.status_code != 200:
            if r.status_code == 401:
                print(f'The given InvenioRDM access token was not accepted by {url}.'
                      'Please ensure that it is still valid.')
                print('Also note that live sites and sandbox sites typically use separate logins '
                      'and separate access tokens.')
            else:
                print(f'An unknown error occurred when trying to access {url}.')
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
            print(f'An error occurred when trying to access <{url}>.')
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

    def search(self, q=None, status=None, sort=None, page=1, size=1000, all_versions=False, max_pages=10):
        """Return list of dictionaries describing each deposit created with this login

        It is possible to filter the results use the optional parameters.  Note that the web interface
        can sometimes be used to find search parameters by looking in the `search-hidden-params`
        parameter of the `invenio-search` tag.

        Example queries
        ===============
        'title:"SXS:BBH:0003"'  # Finds titles with given string; use quotes for robustness
        'owners: 38418'  # Find records by id number of owner

        Optional parameters
        ===================
        q: string [optional]
            Search query, using Elasticsearch query string syntax.  See
            https://help.zenodo.org/guides/search/ for details.
        status: string, either 'draft' or 'published' [optional]
            Filter result based on deposit status.
        sort: string [optional]
            Sort order ('bestmatch' or 'mostrecent').  Prefix with minus to change from ascending to
            descending (e.g., '-mostrecent').
        page: int [optional, defaults to 1]
            Page number for pagination
        size: int [optional, defaults to 1000]
            Number of results to return per page.  Note that Zenodo (as of this writing) seems to
            place a hard limit of 9999 responses.  Anything more will result in an error.  Use
            multiple pages to get more results.
        all_versions: bool [optional, defaults to False]
            If True return all records, including older versions of published records.
        max_pages: int [optional, defaults to 10]
            If the query returns a number of records equal to `size`, it is evidently incomplete.
            This function will attempt to retrieve successive pages until the number of records is
            less than `size`.  If the query is still incomplete after this many pages, just return
            what we've got.

        """
        params={}
        if q is not None:
            params['q'] = q
        if status is not None and status in ['draft', 'published']:
            params['status'] = status
        if sort is not None and sort in ['bestmatch', 'mostrecent', '-bestmatch', '-mostrecent']:
            params['sort'] = sort
        params['page'] = page
        params['size'] = size
        if all_versions:
            params['all_versions'] = ''

        url = f"{self.base_url}api/deposit/depositions"
        r = self.session.get(url, params=params)
        if r.status_code != 200:
            print(f'An unknown error occurred when trying to access {url}.')
            print(f'The search parameters were "{params}"')
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error

        json = r.json()
        if len(json) == size:
            page += 1
            if page > max_pages:
                print(f'Search is not yet complete after {max_pages} pages; returning with what we have.')
                return json  # This will percolate back up the recursion to include other results
            return json + self.search(q=q, sort=sort, page=page, size=size,
                                      all_versions=all_versions, max_pages=max_pages)

        return json

    def delete_untitled_empty_deposits(self):
        deleted_deposits = 0
        deposits = self.search()
        for d in deposits:
            try:
                if d['title'] == '':
                    d = l.deposit(d['id'], ignore_deletion=True)
                    if not d.files:
                        d.delete_deposit(confirmed=True)
                        deleted_deposits += 1
            except:
                pass
        print(f'Deleted {deleted_deposits} deposits')

    def discard_all_drafts(self):
        discarded_drafts = 0
        deposits = self.search()
        for d in deposits:
            try:
                if d['state'] == 'inprogress':
                    d = l.deposit(d['id'], ignore_deletion=True)
                    d.discard()
                    discarded_drafts += 1
            except:
                pass
        print(f'Discarded {discarded_drafts} drafts')

    def awaiting_approval(self, community_id):
        """List all records awaiting approval for the given community"""
        url = f'{self.base_url}/api/records/?q=provisional_communities:{community_id}'
        r = self.session.get(url)
        if r.status_code != 200:
            print(f'Unable to find any records for community {community_id}.')
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
        url = f"{self.base_url}/communities/{community_id}/curaterecord/"
        data = {"recid": int(record_id), "action": "accept"}
        r = self.session.post(url, json=data)
        if r.status_code != 200:
            print(f'Unable to accept record id {record_id} into community {community_id};'
                  f'status code={r.status_code}.')
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
            depositions = self.search()
        else:
            depositions = [self.deposit(deposition_id, ignore_deletion=True).representation]
        total_size = 0
        for deposition in depositions:
            deposition_id = deposition['id']
            d = self.deposit(deposition_id, ignore_deletion=True)
            d_total_size = sum([f['filesize'] for f in d.files])
            print(f'{convert_size(d_total_size)} in "{d.title}" (deposition ID {deposition_id})')
            total_size += d_total_size
        print(f'{convert_size(total_size)} in {len(depositions)} deposits')
        if human_readable:
            return convert_size(total_size)  # Note: the return type will be str
        else:
            return total_size  # Note: the return type will be int
