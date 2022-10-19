"""Class encapsulating interactions with CaltechDATA

This class handles requests through the CaltechDATA web API.  In particular, it manages authorization,
and retries any failed requests automatically.

"""

import requests

class CaltechDATAAuth(requests.auth.AuthBase):
    def __init__(self, base_url, access_token):
        self.base_url = base_url
        self.access_token = access_token
        #super(CaltechDATAAuth, self).__init__()

    def __call__(self, request):
        if request.url.startswith(self.base_url):
            request.headers["Authorization"] = f"Bearer {self.access_token}"
        return request

    def __eq__(self, other):
        return all(
            [
                self.base_url == getattr(other, "base_url", None),
                self.access_token == getattr(other, "access_token", None),
            ]
        )

    def __ne__(self, other):
        return not self == other


class Login(object):
    def __init__(
        self,
        url="https://cd-sandbox.tind.io/",
        access_token=None,
        access_token_path=None,
        total_retry_count=50,
        backoff_factor=0.1,
        backoff_max=20.0,
        session=None,
    ):
        """Initialize a Login object for interacting with CaltechDATA

        This object encapsulates the credentials needed to interact with the CaltechDATA API, and
        exposes a Session object that can be used to make requests which automatically include the
        credentials.  It can be used for generic requests, but note that other objects in this
        module make certain tasks easier -- such as creating or modifying a "deposit", which is
        CaltechDATA's name for a new upload.  The Deposit object should be created from this object.

        These actions require a CaltechDATA API personal access token.  These are obtained from the
        website -- either
            https://data.caltech.edu/account/settings/applications/tokens/new/
        or
            https://cd-sandbox.tind.io/account/settings/applications/tokens/new/
        Note that these two options use separate login systems and separate access tokens.  The
        access token may either be passed as a string to this function (though that means it will
        probably be found in some script file somewhere, which is probably not a good idea for
        security), or can be read from a file.  By default the file from which the token is read is
        '~/.credentials/caltechdata/access_token' or the same name with '_sandbox' appended.  Thus,
        it is probably easiest to simply place your access tokens in those files, so that no
        arguments need to be passed to this function.  As a basic security measure, please ensure
        that those files are not readable by anyone but the user.

        Parameters
        ----------
        url : str [default: "https://cd-sandbox.tind.io/"]
            The base URL of the archive.  Note that the default URL is the "sandbox" version for
            CaltechDATA, which is just for testing purposes, and will likely be cleaned out
            regularly.  To upload to the archival site, pass its URL: "https://data.caltech.edu/".

        access_token: string or None [default: None]
            If present, this is used as the CaltechDATA API access token.

        access_token_path: string or None [default: None]
            If `access_token` is not given, this file is read and the first line is used as the
            access token.  If this argument is None, it defaults to either
            '~/.credentials/caltechdata/access_token' for the regular website or
            '~/.credentials/caltechdata/access_token_sandbox' for the sandbox website.

        total_retry_count: int [default: 50]
            Total number of times to retry requests that fail for retry-able reasons.

        backoff_factor: float [default: 0.1]
            A delay factor to apply between requests after the second try (most errors are resolved
            immediately by a second try without a delay).  After a certain number of total retries,
            the request Session will sleep for:

                {backoff factor} * (2 ^ ({number of total retries} - 1))

            seconds before trying again. For example, if the `backoff_factor` is 0.1, then the
            session will sleep for [0.0s, 0.2s, 0.4s, 0.8s, ...] between retries.  It will never be
            longer than `backoff_max`.

        backoff_max: float [default: 20.0]
            Longest time (in seconds) to wait between retries.

        session: requests.Session or None [default: None]
            This is the object that handles all of the requests made to the API.  If `None`, a
            Session is created for you, and sensible default headers (including the access token)
            are created.  If you need to adjust some of the Session parameters like proxies or SSL
            verification, you can simply create your own and pass it in here.  Note that any `auth`
            property on the passed object will be replaced by one that adds this to the header of
            each request to the chosen CaltechDATA domain:
                {"Authorization": "Bearer <YourAccessTokenHere>"}

        """
        import os
        import requests
        from requests.adapters import HTTPAdapter
        from urllib3.util.retry import Retry
        from datacite import DataCiteRESTClient

        self.base_url = url

        # The `session` object will handle all requests we make.
        self.session = session or requests.Session()

        # Set the CaltechDATA API access token
        if "sandbox" in url:
            access_token_path = "~/.credentials/caltechdata/access_token_sandbox"
            doi_auth_path = "~/.credentials/caltechdata/doi_auth_sandbox"
            self.doi_prefix = "10.80269"
            test_mode = True
        else:
            access_token_path = "~/.credentials/caltechdata/access_token"
            doi_auth_path = "~/.credentials/caltechdata/doi_auth"
            self.doi_prefix = "10.22002"
            test_mode = False
        with open(os.path.expanduser(doi_auth_path), "r") as f:
            user, password = f.readline().strip().split(":", 1)
        self.datacite = DataCiteRESTClient(
            username=user,
            password=password,
            prefix=self.doi_prefix,
            test_mode=test_mode,
        )
        path = os.path.expanduser(access_token_path)
        try:
            with open(path, "r") as f:
                self.access_token = f.readline().strip()
        except IOError:
            print("Unable to find the CaltechDATA access token needed to change a record.")
            print(f"Failed to open file '{path}' for reading.")
            raise
        if not self.access_token:
            print(f"The file '{path}' did not contain any text on the first line.")
            print("This is should be a CaltechDATA access token, which is needed to change a record.")
            raise ValueError("Deposit requires a CaltechDATA access token")

        # Ensure that this session sends the Authorization header with every request to the base_url
        self.session.auth = CaltechDATAAuth(self.base_url, self.access_token)

        # Note that some requests require different choices for 'Accept' and 'Content-Type'; these
        # are altered in the corresponding methods below.
        default_headers = {
            "Accept": "application/json",
            "Content-Type": "application/json",
        }
        self.session.headers.update(default_headers)

        ## Retry automatically on certain types of errors
        Retry.BACKOFF_MAX = backoff_max  # Must be set on the class, not the instance
        retry = Retry(
            total=total_retry_count,
            backoff_factor=backoff_factor,
            status_forcelist=[
                500,
                502,
                503,
                504,
            ],
        )
        adapter = HTTPAdapter(max_retries=retry)
        self.session.mount(self.base_url, adapter)

        # Test to see if we can use the given access token
        url = "{0}api/records".format(self.base_url)
        r = self.session.get(url)
        if r.status_code != 200:
            if r.status_code == 401:
                print(
                    f"The given CaltechDATA access token was not accepted by {self.base_url}.  Please ensure that it is still valid."
                )
                print(
                    "Also note that the standard site and the sandbox site use separate logins and separate access tokens."
                )
            else:
                print(f"An unknown error occurred when trying to access {self.base_url}.")
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an HTTP error

    def send_s3(self, path, name=None, verbose=False):
        import sys
        import pathlib
        #import tqdm

        if name is None:
            name = str(path)
        path = pathlib.Path(path).expanduser().resolve()
        size = path.stat().st_size

        if verbose:
            print(f"  Uploading {name} ({size:_} B) ", end="", flush=True)

        s3url = f"{self.base_url}tindfiles/sign_s3/"  # Note trailing slash
        chkurl = f"{self.base_url}tindfiles/md5_s3"

        r = self.session.get(s3url)
        if r.status_code != 200:
            if r.status_code == 401:
                print(
                    f"The given CaltechDATA access token was not accepted by {self.base_url}.  Please ensure that it is still valid."
                )
                print(
                    "Also note that the standard site and the sandbox site use separate logins and separate access tokens."
                )
                print(f"Used headers {r.request.headers}")
            else:
                print(f"An unknown error occurred when trying to access {self.base_url}.")
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an HTTP error
        s3 = r.json()
        data = s3["data"]
        bucket = s3["bucket"]

        key = data["fields"]["key"]
        policy = data["fields"]["policy"]
        aid = data["fields"]["AWSAccessKeyId"]
        signature = data["fields"]["signature"]
        url = data["url"]

        s3headers = {
            "Host": f"{bucket}.s3.amazonaws.com",
            "Date": "date",
            "x-amz-acl": "public-read",
            "Access-Control-Allow-Origin": "*",
        }

        with path.open("rb") as f:
            #with tqdm.tqdm.wrapattr(f, "read", total=size, desc="    ") as fw:  # To monitor upload progress
            form = (
                ("key", key),
                ("acl", "public-read"),
                ("AWSAccessKeyID", aid),
                ("policy", policy),
                ("signature", signature),
                ("file", f),
            )
            response = requests.session().post(url, files=form, headers=s3headers)
        if response.status_code != 204:
            if response.status_code == 400:
                print(f"Bad request: Probably caused by incorrectly formed input")
                print(f"Used headers {response.request.headers}")
            else:
                print(f"An unknown error occurred when trying to access {self.base_url}.")
            try:
                print(response.json())
            except:
                pass
            try:
                print(response.text)
            except:
                pass
            response.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an HTTP error

        response = self.session.get(f"{chkurl}/{bucket}/{key}/")
        md5 = response.json()["md5"]

        fileinfo = {"url": key, "filename": name, "md5": md5, "size": size}

        if verbose:
            print(f"✓")

        return fileinfo

    def download(self, url, path):
        """Download large file efficiently

        Parameters
        ----------
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
        from urllib.parse import urlparse

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
            print("An error occurred when trying to access <{0}>.".format(url))
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        r.raw.read = partial(r.raw.read, decode_content=True)
        # r.raw.decode_content = True
        with open(local_filename, "wb") as f:
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

    def search(self, q=None, sort=None, size=1000, page=1, allversions=False, max_pages=10):
        """Return list of dictionaries describing each deposit created with this login

        It is possible to filter the results use the optional parameters.  Note that the web interface
        can sometimes be used to find search parameters by looking in the `search-hidden-params`
        parameter of the `invenio-search` tag.

        Example queries
        ---------------
        'title:"SXS:BBH:0003"'  # Finds titles with given string; use quotes for robustness
        'communities:sxs'  # Records in the 'sxs' CaltechDATA community
        'provisional_communities:sxs'  # Records awaiting approval by the community curator
        'owners: 38418'  # Find records by id number of owner

        Optional parameters
        -------------------
        q: string [optional]
            Search query, using Elasticsearch query string syntax.  See
            https://help.zenodo.org/guides/search/ for details.
        sort: string [optional]
            Sort order ('bestmatch' or 'mostrecent').  Prefix with minus to change from ascending to
            descending (e.g., '-mostrecent').
        size: int [optional, defaults to 1000]
            Number of results to return per page.  Note that CaltechDATA (as of this writing) seems to
            place a hard limit of 9999 responses.  Anything more will result in an error.  Use
            multiple pages to get more results.
        page: int [optional, defaults to 1]
            Page number for pagination
        allversions: bool [optional, defaults to False]
            If True return all records, including older versions of published records.
        max_pages: int [optional, defaults to 10]
            If the query returns a number of records equal to `size`, it is evidently incomplete.
            This function will attempt to retrieve successive pages until the number of records is
            less than `size`.  If the query is still incomplete after this many pages, just return
            what we've got.

        """
        params = {}
        if q is not None:
            params["q"] = q
        if sort is not None and sort in ["bestmatch", "mostrecent", "-bestmatch", "-mostrecent"]:
            params["sort"] = sort
        params["page"] = page
        params["size"] = size
        if allversions:
            params["allversions"] = ""

        url = "{0}api/records".format(self.base_url)
        r = self.session.get(url, params=params)
        if r.status_code != 200:
            print("An unknown error occurred when trying to access {0}.".format(url))
            print('The search parameters were "{0}"'.format(params))
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
                print("Search is not yet complete after {0} pages; returning with what we have.".format(max_pages))
                return json  # Note: This will percolate back up the recursion to include other results
            return json + self.search(
                q=q, sort=sort, page=page, size=size, allversions=allversions, max_pages=max_pages
            )

        return json

    def delete_untitled_empty_deposits(self):
        deleted_deposits = 0
        deposits = self.search()
        for d in deposits:
            try:
                if d["title"] == "":
                    d = self.deposit(d["id"], ignore_deletion=True)
                    if not d.files:
                        d.delete_deposit(confirmed=True)
                        deleted_deposits += 1
            except:
                pass
        print("Deleted {0} deposits".format(deleted_deposits))

    def discard_all_drafts(self):
        discarded_drafts = 0
        deposits = self.search()
        for d in deposits:
            try:
                if d["state"] == "inprogress":
                    d = self.deposit(d["id"], ignore_deletion=True)
                    d.discard()
                    discarded_drafts += 1
            except:
                pass
        print("Discarded {0} drafts".format(discarded_drafts))

    def awaiting_approval(self, community_id):
        """List all records awaiting approval for the given community"""
        url = "{0}/api/records/?q=provisional_communities:{1}".format(self.base_url, community_id)
        r = self.session.get(url)
        if r.status_code != 200:
            print("Unable to find any records for community {0}.".format(community_id))
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        return r.json()

    def community_curate_accept(self, community_id, record_id):
        """Accept a record into the community"""
        url = "{0}/communities/{1}/curaterecord/".format(self.base_url, community_id)
        data = {"recid": int(record_id), "action": "accept"}
        r = self.session.post(url, json=data)
        if r.status_code != 200:
            print(
                "Unable to accept record id {0} into community {1}; status code={2}.".format(
                    record_id, community_id, r.status_code
                )
            )
            try:
                r_json = r.json()
                print("Response JSON:")
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
            id = deposition["id"]
            d = self.deposit(id, ignore_deletion=True)
            d_total_size = sum([f["filesize"] for f in d.files])
            print('{1} in "{2}" (CaltechDATA ID {0})'.format(id, convert_size(d_total_size), d.title))
            total_size += d_total_size
        print("{0} in {1} deposits".format(convert_size(total_size), len(depositions)))
        if human_readable:
            return convert_size(total_size)  # Note: the return type will be str
        else:
            return total_size  # Note: the return type will be int
