"""Search NASA's ADS system

Documentation can be found here: <http://adsabs.github.io/help/>

"""


api_url = 'https://api.adsabs.harvard.edu/v1/'


class Login(object):

    def __init__(self, api_token=None, api_token_path=None,
                 total_retry_count=50, backoff_factor=0.1, backoff_max=20.0, session=None):
        """Initialize a Login object for interacting with NASA ADS

        This object encapsulates the credentials needed to interact with the NASA/ADS API, and
        exposes a Session object that can be used to make requests which automatically include the
        credentials.

        These actions require a NASA/ADS API access token.  These are obtained from the website:
            https://ui.adsabs.harvard.edu/user/settings/token
        The access token may either be passed as a string to this function (though that means it
        will probably be found in some script file somewhere, which is probably not a good idea for
        security), or can be read from a file.  By default the file from which the token is read is
        '~/.credentials/ads/api_token'.  Thus, it is probably easiest to simply place your access
        tokens in those files, so that no arguments need to be passed to this function.  As a basic
        security measure, please ensure that those files are not readable by anyone but the user.


        Parameters
        ==========
        api_token: string or None [default: None]
            If present, this is used as the NASA/ADS API access token.

        api_token_path: string or None [default: None]
            If `api_token` is not given, this file is read and the first line is used as the access
            token.  If this argument is None, it defaults to '~/.credentials/ads/api_token'.

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
            each request to the NASA/ADS domain:
                {"Authorization": "Bearer <YourAccessTokenHere>"}

        """
        import os
        import requests
        from requests.adapters import HTTPAdapter
        from requests.packages.urllib3.util.retry import Retry

        self.base_url = api_url

        # The Session object will handle all requests we make.
        self.session = session or requests.Session()

        # If the input session object succeeds at a `get`, we can skip a lot of the following
        r = self.session.get(self.base_url+'search/query', params={'q': 'star'})
        if r.status_code != 200:  # That's okay, we mostly expected this

            # Set the NASA/ADS API access token
            if api_token is not None:
                self.api_token = api_token
            else:
                if api_token_path is None:
                    api_token_path = '~/.credentials/ads/api_token'
                path = os.path.expanduser(api_token_path)
                try:
                    with open(path, 'r') as f:
                        self.api_token = f.readline().strip()
                except IOError:
                    print('Unable to find the NASA/ADS access token needed to query the system.')
                    print('Failed to open file "{0}" for reading.'.format(path))
                    raise
                if not self.api_token:
                    print('The file "{0}" did not contain any text on the first line.'.format(path))
                    print('This is should be a NASA/ADS access token, which is need to query the system.')
                    raise ValueError('Query requires a NASA/ADS access token')

            # Ensure that this session sends the Authorization header with every request to the base_url
            class NASAADSAuth(requests.auth.AuthBase):
                def __init__(self, base_url, api_token):
                    self.base_url = base_url
                    self.api_token = api_token
                    super(NASAADSAuth, self).__init__()
                def __call__(self, r):
                    if r.url.startswith(self.base_url):
                        r.headers.update({"Authorization": "Bearer {0}".format(self.api_token)})
                    return r
            self.session.auth = NASAADSAuth(self.base_url, self.api_token)

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
        r = self.session.get(self.base_url+'search/query', params={'q': 'star'})
        if r.status_code != 200:
            if r.status_code == 401:
                print('The given NASA/ADS access token was not accepted by {0}.  Please ensure that it is still valid.'.format(self.base_url))
            else:
                print('An unknown error occurred when trying to access {0}.'.format(self.base_url))
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error

    @property
    def export_formats(self):
        """List of formats in which ADS can export records"""
        return ['bibtex', 'bibtexabs', 'ads', 'endnote', 'procite', 'ris', 'refworks', 'rss', 'medlars',
                'dcxml', 'refxml', 'refabsxml', 'aastex', 'icarus', 'mnras', 'soph', 'votable', 'custom']

    def export(self, format='bibtex', ):
        pass

    def query(self, query_dict):
        query_url = api_url + 'search/query'
        r = self.session.get(query_url, params={'q': query_dict})
        if r.status_code != 200:
            print('An error occurred when trying to access <{0}>.'.format(query_url), file=sys.stderr)
            try:
                print(r.json(), file=sys.stderr)
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        return r.json()


def ads_bearer_token(error_if_not_found=True):
    import os.environ
    if 'ADS_TOKEN' not in os.environ:
        message = 'Environment variable "ADS_TOKEN" not found'
        if error_if_not_found:
            raise EnvironmentError(message)
        else:
            from warnings import warn
            warn(message)
            return ''
    return os.environ('ADS_TOKEN')


def authorization_header(error_if_not_found=True):
    token = ads_bearer_token(error_if_not_found)
    return 'Authorization: Bearer:{0}'.format(token)



#curl -H 'Authorization: Bearer:Co3woTdgh00kgY46RrP51sVHZBFEiml6LeQq2xNK' 'http://api.adsabs.harvard.edu/v1/search/query?q=doi:10.1111/j.1365-2966.2010.18127.x&fl=pub+volume+issue+year+page'
