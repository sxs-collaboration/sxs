"""Search the INSPIRE database

Documentation can be found here: <https://inspirehep.net/info/hep/api>

"""

api_url = "https://inspirehep.net/search"


def query(pattern, output_format='recjson', output_tags=None, records_in_groups_of=None, jump_to_records=None):
    """Search the INSPIRE database

    API documentation can be found here: <https://inspirehep.net/info/hep/api>
    More search documentation can be found here: <https://inspirehep.net/info/hep/search-tips>

    Parameters
    ==========
    pattern: str
        This is the query in the Inspire search syntax. All search features and operators familiar
        from the Inspire web interface and documented in the online help are supported and complex
        queries are possible.  Documentation here: <https://inspirehep.net/info/hep/search-tips>.

    output_format: str, optional [defaults to 'recjson']
        The format of the response sent back to the client. There are two choices, 'xm' for
        (MARC-)XML or 'recjson' for JSON. The XML response format is MARCXML or fragments thereof
        when individual fields are selected via the output_tags parameter

    output_tags: str, optional [defaults to None]
        If present, this selects (filter) specific tags from the response. If output_format is 'xm',
        this option takes a comma separated list of MARC tags; valid MARC tags for Inspire records
        can be found here <https://twiki.cern.ch/twiki/bin/view/Inspire/DevelopmentRecordMarkup>.
        If output_format is 'recjson', this is similar to selecting MARC tags, however by name
        instead of numerical value.  In addition the JSON response can contain derived or
        dynamically calculated values which are not available in MARC.  If this argument is None,
        the default output will be returned, which generally includes all available information for
        each record.

    records_in_groups_of: int, optional [defaults to None]
        If present, this parameter specifies the number of records per chunk for long
        responses. Note that the default setting is 25 records per chunk. The maximum number is 250.

    jump_to_records: int, optional [defaults to None]
        Long responses are split into several chunks. To access subsequent chunks specify the record
        offset with this parameter.

    Returns
    =======
    json: list or dict
        Usually, this will be a list of the results from the query, even if there is only one.  Each
        result will be a dictionary mapping the requested 'output_tags' to their corresponding
        values.

    Raises
    ======
    requests.exceptions.HTTPError
        If the HTTP request to INSPIRE failed for any reason.

    """
    import sys
    import requests

    params = {
        'p': pattern,
        'of': output_format,
    }
    if output_tags is not None:
        params['ot'] = output_tags
    if records_in_groups_of is not None:
        params['rg'] = records_in_groups_of
    if jump_to_records is not None:
        params['jrec'] = jump_to_records
    
    r = requests.get(api_url, params=params)
    if r.status_code != 200:
        print('An error occurred when trying to access <{0}>.'.format(api_url), file=sys.stderr)
        try:
            print(r.json(), file=sys.stderr)
        except:
            pass
        r.raise_for_status()
        raise RuntimeError()  # Will only happen if the response was not strictly an error
    return r.json()


def extract_bibtex_key(system_control_numbers):
    """Get bibtex key from 'system_control_numbers' field
    
    Unfortunately, this seems to be the only way to get the bibtex key.  I have seen suggestions
    around the github issues for inspirehep/invenio that this should always be present
    
    """
    bibtex_keys = [number.get('value', '') for number in system_control_numbers
                   if number['institute']=='INSPIRETeX' or number['institute']=='SPIRESTeX']
    bibtex_keys = [key for key in bibtex_keys if key]
    if not bibtex_keys:
        return ''
    return bibtex_keys[0]


def extract_doi(doi):
    """Ensure that 'doi' is a single string
    
    Occasionally, INSPIRE returns a list of identical DOIs.  This just extracts the first
    element if it is such a list, or returns the input otherwise.
    
    """
    if isinstance(doi, list) and len(doi)>0:
        return doi[0]
    return doi


def map_bibtex_keys_to_doi(bibtex_keys):
    """Map a list of INSPIRE bibtex keys to DOIs

    This function queries the INSPIRE database, searching for each of the input bibtex keys, and
    extracting the corresponding DOIs (if present on INSPIRE).

    Parameter
    =========
    bibtex_keys: str, or list of str
        Each string should be precisely one bibtex key from INSPIRE.  Note that any bibtex keys that
        are not found by INSPIRE are simply ignored.

    Returns
    =======
    key_to_doi: dict
        The output is a dictionary mapping each input bibtex key itself to the corresponding DOI.
        These are raw DOIs; to get a URL, just append the DOI to 'https://dx.doi.org/', which will
        resolve to the appropriate web page for that DOI.  Note that any record that is found by its
        bibtex key but does not have a DOI will simply be absent from this dictionary.

    Raises
    ======
    requests.exceptions.HTTPError
        If the HTTP request to INSPIRE failed for any reason.

    """
    if not isinstance(bibtex_keys, list):
        bibtex_keys = [bibtex_keys,]
    pattern = 'find texkey ' + ' or texkey '.join(bibtex_keys)
    output_tags = 'system_control_number,doi'
    results = query(pattern, output_tags=output_tags)
    mapping = {
        bibtex_key: doi
        for result in results
        for bibtex_key in [extract_bibtex_key(result['system_control_number']),]
        for doi in [extract_doi(result['doi']),]
        if bibtex_key and doi
    }
    return mapping
