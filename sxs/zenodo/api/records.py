
class Records(object):

    @classmethod
    def search(cls, q=None, sort=None, page=None, size=None, sandbox=False):
        """Search public records

        It is possible to filter the results use the optional parameters.  Note that the web interface
        can sometimes be used to find search parameters by looking in the `search-hidden-params`
        parameter of the `invenio-search` tag.

        Example queries
        ===============
        'title "SXS:BBH:0003"'  # Finds titles with given string; use quotes for robustness
        'communities:sxs'  # Records in the 'sxs' Zenodo community
        'provisional_communities:sxs'  # Records awaiting approval by the community curator
        'owners: 38418'  # Find records by id number of owner

        Optional parameters
        ===================
        q: string
            Search query, using Elasticsearch query string syntax.  See
            https://help.zenodo.org/guides/search/ for details.
        sort: string
            Sort order ('bestmatch' or 'mostrecent').  Prefix with minus to change form ascending to
            descending (e.g., '-mostrecent').
        page: int
            Page number for pagination
        size: int
            Number of results to return per page
        sandbox: bool [defaults to False]
            If True use sandbox.zenodo.org instead of the standard site.

        """
        import requests
        from . import url_sandbox, url_standard

        if sandbox:
            base_url = url_sandbox
        else:
            base_url = url_standard
        url = base_url + "api/records/"

        params={}
        if q is not None:
            params['q'] = q
        if sort is not None:
            params['sort'] = sort
        if page is not None:
            params['page'] = page
        if size is not None:
            params['size'] = size

        r = requests.get(url, params=params, headers={'Accept': 'application/json'})
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
