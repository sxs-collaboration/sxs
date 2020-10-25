from .journal_abbreviations import (
    journal_macros, journal_abbreviation_pairs, journal_full_names, journal_abbreviations
)

from .arxiv import (
    request_data, get_entry_by_arxiv_id, get_all_entries, get_journal_reference,
    get_submission_comment, get_journal_reference_from_doi
)

from .fairchild_report import format_paper_list_entry

from .references import limit_author_list
