

class Login(object):
    import re
    import requests
    import json
    import os
    from http.client import responses

    url_standard = 'https://zenodo.org/'
    url_sandbox = 'https://sandbox.zenodo.org/'

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
        r = self.session.get(self.base_url + "api/deposit/depositions")
        if r.status_code != 200:

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

            # Note that some requests require different choices for 'Accept' and 'Content-Type'; these
            # are altered in the corresponding methods below.
            default_headers= {
                "Accept": "application/json",
                "Content-Type": "application/json",
                "Authorization": "Bearer {0}".format(self.access_token),
            }
            self.session.headers.update(default_headers)

        # Test to see if we can use the given access token
        r = self.session.get(self.base_url + "api/deposit/depositions")
        if r.status_code == 401:
            print('The given Zenodo access token was not accepted.  Please ensure that it is still valid.')
            print('Also note that zenodo.org and sandbox.zenodo.org use separate logins and separate access tokens.')
            r.raise_for_status()
        else if r.status_code != 200:
            print('An unknown error occurred when trying to access {0}.'.format(self.base_url))
            print('The returned HTTP status code was "{0} {1}".'.format(r.status_code, responses[r.status_code]))
            r.raise_for_status()



class Deposition(object):
    import re
    import requests
    import json
    import os
    from http.client import responses

    url_standard = 'https://zenodo.org/'
    url_sandbox = 'https://sandbox.zenodo.org/'

    def __init__(self, deposition_id=None, sandbox=False, access_token=None, access_token_path=None, **kwargs):
        """Initialize a Deposition object for creating a new zenodo entry

        This object encapsulates all the actions you might want to take when creating, publishing,
        updating, or replacing an entry in zenodo.

        Parameters
        ==========
        deposition_id: string, int, or None [default: None]
            If present, this is used as the id of the deposition to edit.  If `None`, a new
            deposition is created.

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
        r = self.session.get(self.base_url + "api/deposit/depositions")
        if r.status_code != 200:

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

            # Note that some requests require different choices for 'Accept' and 'Content-Type'; these
            # are altered in the corresponding methods below.
            default_headers= {
                "Accept": "application/json",
                "Content-Type": "application/json",
                "Authorization": "Bearer {0}".format(self.access_token),
            }
            self.session.headers.update(default_headers)

        # Test to see if we can use the given access token
        r = self.session.get(self.base_url + "api/deposit/depositions")
        if r.status_code == 401:
            print('The given Zenodo access token was not accepted.  Please ensure that it is still valid.')
            print('Also note that zenodo.org and sandbox.zenodo.org use separate logins and separate access tokens.')
            r.raise_for_status()
        else if r.status_code != 200:
            print('An unknown error occurred when trying to access {0}.'.format(self.base_url))
            print('The returned HTTP status code was "{0} {1}".'.format(r.status_code, responses[r.status_code]))
            r.raise_for_status()

        # Now, create or reacquire the specific deposition we're looking for
        if deposition_id is not None:
            # If the deposition id was given, check that we can access it
            r = self.session.get("{0}api/deposit/depositions/{1}".format(self.base_url, deposition_id))
            if r.status_code != 200:
                print('The input deposition id "{0}" could not be accessed on {1}.'.format(self.deposition_id, self.base_url))
                print('The returned HTTP status code was "{0} {1}".'.format(r.status_code, responses[r.status_code]))
                r.raise_for_status()
                raise RuntimeError()  # Will only happen if the response was not strictly an error
        else:
            r = self.session.post("{0}api/deposit/depositions".format(self.base_url, self.deposition_id), data="{}")
            if r.status_code != 201:
                print('Unable to create a new deposition on {0}.'.format(self.base_url))
                print('The returned HTTP status code was "{0} {1}".'.format(r.status_code, responses[r.status_code]))
                r.raise_for_status()
                raise RuntimeError()  # Will only happen if the response was not strictly an error

        # Now, using the response generated above, set some data for this deposition
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
        raise NotImplementedError

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
        r = requests.put('{0}/{1}'.format(self.links['bucket_url'], name),
                         data=open(path, 'rb'),  headers={"Content-Type":"application/octet-stream"})
        raise NotImplementedError


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

                >>> from sxs.zenodo import Deposition
                >>> d = Deposition(id='{deposition_id}', {arguments})

            """.format(deposition_id=self.deposition_id, base_url=self.base_url, arguments=self.arguments)
            warn(dedent(warning))
