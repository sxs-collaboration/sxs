"""Search the INSPIRE database

Documentation can be found here: <https://github.com/inspirehep/rest-api-doc>

"""

import warnings

api_url = "https://inspirehep.net/api/{record_type}"


def inspire2doi(inspire_bibtex_key, raise_exceptions=False):
    try:
        result = query(inspire_bibtex_key, fields="dois.value,arxiv_eprints")
    except KeyboardInterrupt:
        raise
    except Exception as e:
        warnings.warn(f"Error querying for INSPIRE key {inspire_bibtex_key}")
        if raise_exceptions:
            raise e
    else:
        if not result:
            warnings.warn(f"No entry found for INSPIRE key {inspire_bibtex_key}")
            if raise_exceptions:
                raise ValueError(f"No entry found for INSPIRE key {inspire_bibtex_key}")
        else:
            try:
                return result[0]["metadata"]["dois"][0]["value"]
            except Exception as e:
                # DOI not present; INSPIRE only has an arXiv number, so we'll convert to an arXiv DOI
                try:
                    eprint_value = result[0]["metadata"]["arxiv_eprints"][0]["value"]
                    # https://info-arxiv-org.proxy.library.cornell.edu/help/arxiv_identifier.html
                    # says that identifiers starting with 9107 to 0703 need their category included
                    if eprint_value.startswith("9") or int(eprint_value[:4])<704:
                        category = result[0]["metadata"]["arxiv_eprints"][0]["categories"][0]
                        return f"10.48550/arXiv.{category}/{eprint_value}"
                    else:
                        return f"10.48550/arXiv.{eprint_value}"
                except Exception as e:
                    warnings.warn(f"Unexpected result format for INSPIRE key {inspire_bibtex_key}: \n{result}")
                    if raise_exceptions:
                        raise e


def query(
    query,
    sort="mostrecent",
    page=1,
    size=1000,
    fields=None,
    page_limit=10,
    record_type="literature"
):
    """Fetch records from the INSPIRE API.

    Parameters
    ----------
    query : str
        The search query.  For details, see
        https://github.com/inspirehep/rest-api-doc/tree/master#search-query
    fields : str, optional
        The fields to retrieve from INSPIRE.  This can be a
        comma-separated list of field names.  For example, to
        retrieve the list of DOIs: "dois.value".  Other options at
        https://inspire-schemas.readthedocs.io/en/latest/schemas/elements/reference.html
        Note that other fields will also be returned for each record,
        even when you just ask for one — including "id", "links",
        "metadata", "created", and "updated".  The default value of
        None will return all fields.
    sort : str, optional
        The sorting order. Default is "mostrecent".
    size : int, optional
        The number of records per page. Default is 1000.
    page : int, optional
        The starting page number. Default is 1.
    page_limit : int, optional
        Optional limit to the number of pages to collect.
    record_type : str, optional
        The type of record to fetch (e.g., "literature"). Default is
        "literature".  Other options can be found here:
        https://github.com/inspirehep/rest-api-doc/tree/master#obtaining-a-record
    
    Returns
    -------
    list
        The JSON response from the API.

    Raises
    ------
    requests.exceptions.HTTPError :
        If the HTTP request to INSPIRE failed

    """
    import sys
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    session = requests.Session()
    collected_results = []

    ## Retry automatically on certain types of errors
    retry = Retry(
        total=10,
        backoff_factor=0.1,
        status_forcelist=[
            413,  # Request Entity Too Large
            429,  # Too Many Requests
            500,  # Internal Server Error
            502,  # Bad Gateway
            503,  # Service Unavailable
            504,  # Gateway Timeout
        ],
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)

    url = api_url.format(record_type=record_type)
    params = {
        "q": query,
        "sort": sort,
        "size": size,
        "page": page,
    }
    if fields:
        params["fields"] = fields

    while params["page"] - page <= page_limit:
        response = session.get(url, params=params)
        if response.status_code != 200:
            print(f"An error occurred when trying to access <{api_url}>.", file=sys.stderr)
            try:
                print(response.json(), file=sys.stderr)
            except:
                pass
            response.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error

        try:
            json_response = response.json()
        except ValueError:
            print(f"Response from {url} does not contain valid JSON; returning early.")
            return collected_results

        new_hits = json_response.get("hits", {}).get("hits", [])
        if not new_hits:
            break
        collected_results.extend(new_hits)

        params["page"] += 1

    return collected_results
