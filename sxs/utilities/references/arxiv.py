arxiv_api_url = 'https://export.arxiv.org/api/query'



def request_data(search_query='', id_list='', start=0, max_results=10):
    import requests
    data = {
        'search_query': search_query,
        'id_list': id_list,
        'start': start,
        'max_results': max_results,
    }
    return requests.post(arxiv_api_url, data=data).text


def get_entry_by_arxiv_id(id):
    import feedparser
    response = request_data(id_list=id)
    parsed = feedparser.parse(response)
    return parsed['entries'][0]


def get_all_entries(search_query='', id_list='', start=0, max_results=100):
    import requests
    data = {
        'search_query': search_query,
        'id_list': id_list,
        'start': start,
        'max_results': max_results,
    }
    result = requests.post(arxiv_api_url, data=data)
    entries = result['entries']
    tr = result['opensearch_totalresults']
    si = result['opensearch_startindex']
    ip = result['opensearch_itemsperpage']
    while si+ip < tr:
        data['start'] = si+ip
        result = requests.post(arxiv_api_url, data=data)
        entries += result['entries']
        tr = result['opensearch_totalresults']
        si = result['opensearch_startindex']
        ip = result['opensearch_itemsperpage']
    return entries


def get_journal_reference(entry, string=''):
    if 'journal_ref' in entry:
        return entry['journal_ref']
    if 'arxiv_journal_ref' in entry:
        return entry['arxiv_journal_ref']
    if 'arxiv_doi' in entry:
        try:
            return get_journal_reference_from_doi(entry['arxiv_doi'], string)
        except:
            pass  # Didn't work for some reason; oh well...
    return ''


def get_submission_comment(entry):
    import re
    if 'arxiv_comment' in entry:
        publication_comment_regex = re.compile(
            r"""(?P<comment>[Aa]ccepted (?:for publication )?(?:by|in|to) |[Ss]ubmitted to |[Ii]n press (?:with )?)(?P<publication>[^,;]*)""")
        comment = entry['arxiv_comment']
        search = publication_comment_regex.search(comment)
        if search:
            comment = search['comment']
            publication = search['publication']
            publication_comment = ', ' + comment[0].lower() + comment[1:] + publication
            return publication_comment
    return ''


def get_journal_reference_from_doi(doi, string=''):
    if 'PhysRev' in doi:
        jr = get_journal_reference_from_phys_rev(doi)
    else:
        jr = get_journal_reference_from_ads(doi)
    if string:
        return string.format(jr)
    else:
        return jr


def get_journal_reference_from_phys_rev(doi):
    import re
    import urllib.request
    reference = {}
    pr = re.search('PhysRev([A-Z])', doi)
    fields = re.compile(' (' + '|'.join(['title', 'pub', 'volume', 'issue', 'year', 'pages'])
                        + ') = {(.*)}')
    if pr:
        reference['pub'] = 'Phys. Rev. {0}'.format(pr.groups()[0])
        url = 'https://journals.aps.org/pr{0}/export/{1}'.format(pr.groups()[0].lower(), doi)
        with urllib.request.urlopen(url) as response:
            html = response.read().decode()
        for line in html.split('\n'):
            f = fields.search(line)
            if f:
                field, value = f.groups()[0], f.groups()[1]
                reference[field] = value
    return reference


def get_journal_reference_from_ads(doi):
    import sys
    import difflib
    import ads
    from .journal_abbreviations import journal_abbreviation_pairs
    reference = {'title':'', 'pub':'', 'volume':'', 'issue':'', 'year':'', 'page':''}
    try:
        article = ads.SearchQuery(q="doi:{0}".format(doi), fl=list(reference)).next()
        reference['title'] = article.title
        reference['pub'] = article.pub
        reference['volume'] = article.volume
        reference['issue'] = article.issue
        reference['year'] = article.year
        reference['page'] = article.page
    except Exception as e:
        print('Failed to get ADS entry for doi "{0}"'.format(doi), file=sys.stderr)
    closest_match = difflib.get_close_matches(reference['pub'], journal_abbreviation_pairs, n=1)
    if closest_match:
        reference['pub'] = journal_abbreviation_pairs[closest_match[0]]
    return reference
