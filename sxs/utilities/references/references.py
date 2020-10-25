def limit_author_list(all_authors, desired_authors=[], author_list_length_limit=10):
    if len(all_authors) <= author_list_length_limit:
        return all_authors
    author_list = [a for a in all_authors if a in desired_authors]
    if not author_list:
        author_list = all_authors[:author_list_length_limit]
    if len(author_list) < len(all_authors):
        author_list += ['et~al.']
    return author_list
