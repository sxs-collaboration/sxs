"""Search the INSPIRE database

Documentation can be found here: <https://inspirehep.net/info/hep/api>

"""

api_url = "https://inspirehep.net/search"


def query(pattern, output_format='recjson', output_tags=None, records_in_groups_of=None, jump_to_records=None):
    """Search the INSPIRE database

    API documentation can be found here: <https://inspirehep.net/info/hep/api>

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
        instead of numerical value. In addition the JSON response can contain derived or dynamically
        calculated values which are not available in MARC.

    records_in_groups_of: int, optional [defaults to None]
        If present, this parameter specifies the number of records per chunk for long
        responses. Note that the default setting is 25 records per chunk. The maximum number is 250.

    jump_to_records: int, optional [defaults to None]
        Long responses are split into several chunks. To access subsequent chunks specify the record
        offset with this parameter.

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


def bibtex_key_to_doi(bibtex_key):
    """Translate INSPIRE bibtex key to DOI"""
    import sys
    pattern = 'find texkey ' + bibtex_key
    output_tags = 'doi'
    r_json = query(pattern, output_tags=output_tags)
    if len(r_json) < 1:
        print('No records found on INSPIRE with bibtex key "{0}".'.format(bibtex_key), file=sys.stderr)
        raise ValueError()
    return r_json[0]['doi']
