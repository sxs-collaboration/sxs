"""Class encapsulating "record" objects on Zenodo

These are the published versions of "deposits".  In particular, if a record is publicly available,
no Login object is required.

"""


class Records(object):

    @classmethod
    def search(cls, q=None, sort=None, page=1, size=1000, sandbox=False, all_versions=False, max_pages=10):
        """Search public records

        It is possible to filter the results use the optional parameters.  Note that the web interface
        can sometimes be used to find search parameters by looking in the `search-hidden-params`
        parameter of the `invenio-search` tag.

        Example queries
        ---------------
        'title:"SXS:BBH:0003"'  # Finds titles with given string; use quotes for robustness
        'communities:sxs'  # Records in the 'sxs' Zenodo community
        'provisional_communities:sxs'  # Records awaiting approval by the community curator
        'owners: 38418'  # Find records by id number of owner

        Optional parameters
        -------------------
        q: string
            Search query, using Elasticsearch query string syntax.  See
            https://help.zenodo.org/guides/search/ for details.
        sort: string
            Sort order ('bestmatch' or 'mostrecent').  Prefix with minus to change from ascending to
            descending (e.g., '-mostrecent').
        page: int
            Page number for pagination
        size: int
            Number of results to return per page.  Note that Zenodo (as of this writing) seems to
            place a hard limit of 9999 responses.  Anything more will result in an error.  Use
            multiple pages to get more results.
        sandbox: bool [defaults to False]
            If True use sandbox.zenodo.org instead of the standard site.
        all_versions: bool [defaults to False]
            If True return all records, including older versions of published records.
        max_pages: int [defaults to 10]
            If the query returns a number of records equal to `size`, it is evidently incomplete.
            This function will attempt to retrieve successive pages until the number of records is
            less than `size`.  If the query is still incomplete after this many pages, just return
            what we've got.

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
        params['page'] = page
        params['size'] = size
        if all_versions:
            params['all_versions'] = ''

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

        json = r.json()
        if len(json) == size:
            page += 1
            if page > max_pages:
                print('Search is not yet complete after {0} pages; returning with what we have.'.format(max_pages))
                return r.json()
            return json + cls.search(q=q, sort=sort, page=page, size=size, sandbox=sandbox, all_versions=all_versions, max_pages=max_pages)

        return json
