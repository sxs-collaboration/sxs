from .references import limit_author_list
from .arxiv import get_entry_by_arxiv_id, get_journal_reference, get_submission_comment


def format_paper_list_entry(arxiv_id, submission_comment='', desired_authors=None, author_list_length_limit=10):
    import re
    from pylatexenc.latexencode import utf8tolatex
    import textwrap
    wrapper = textwrap.TextWrapper(initial_indent='  ', subsequent_indent='    ', width=70,
                                   break_long_words=False, break_on_hyphens=False)
    entry = get_entry_by_arxiv_id(arxiv_id)
    authors = [author['name'] for author in entry['authors']]
    authors = limit_author_list(authors, desired_authors=desired_authors, author_list_length_limit=author_list_length_limit)
    authors = [utf8tolatex(a).replace(r'{\textasciitilde}', '~') for a in authors]
    author_list = ', '.join(authors)
    title = r'\textit{{{0}}}'.format(entry['title'])
    submission_comment = submission_comment.strip()
    if not submission_comment:
        journal_ref = get_journal_reference(entry, r'{pub} \textbf{{{volume}}}, {issue} ({year})')
        if journal_ref:
            submission_comment = ', published as ' + journal_ref.strip()
        else:
            submission_comment = get_submission_comment(entry)
    if submission_comment and not submission_comment.startswith(','):
        submission_comment = ', ' + submission_comment
    arxiv_url = re.sub('v[0-9]*$', '', entry['id'].replace('http://', 'https://'))
    entry = r'\item {0}, {1}{2}; available at \url{{{3}}}'.format(author_list, title, submission_comment, arxiv_url)
    entry = wrapper.fill(' '.join(entry.split()))
    return entry
