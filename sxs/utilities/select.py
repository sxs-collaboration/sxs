"""Select items from a list"""


def select_by_path_component(path_pattern, possible_matches, recursion_index=0):
    """Select from a list of paths by progressively matching path components

    Parameters
    ----------
    path_pattern : str
        A pattern to search for in the sequence of possible matches.  If
        `recursion_index` is 0, this is first searched as an exact partial match,
        meaning paths that start with exactly this string are chosen.  If that does
        not produce any matches, the pattern is split into components, separated by
        "/".  This is then treated as a series of patterns to match corresponding
        path components in the `possible_matches`.  Each component pattern will be
        interpreted first as a literal match, then — if no match is found — as a
        regex, compiled by python's `re` module.
    possible_matches : Iterable[str]
        A sequence of paths to select matches from.
    recursion_index : int, optional
        This is presumably mostly for use from inside this function, and denotes
        the number of path components to try to match at once.  Defaults to 0.

    Returns
    -------
    matched_paths : Set[str]
        
    Notes
    -----
    Here, we consider a path as a sequence of components separated by "/".  Each
    path component is just a string without any "/" character.

    This function progressively matches paths, starting with the first path
    component, then searching the remaining possibilities for the second path
    component, and so on.  In each search, the pattern is first treated as a
    literal string to match to *the beginning of* a path.  If nothing is found, the
    pattern is compiled as a regex (by python's standard `re` module), and the
    search is run again.  In each case, if the current component was only a partial
    match (meaning there are characters in this path component *after* whatever was
    matched by the pattern), then the "largest" result is used — meaning literally
    that the matching path components are fed into the `max` function.  The idea is
    that if there are, for example, multiple versions or multiple resolutions of a
    simulation, then the largest is usually the "best", so we can default to that.

    An example may help to clarify.  Suppose we have this list of possible matches:

        possible_matches = [
            "abc/lm/xyz",
            "abc/ln/xyz",
            "abd/lm/xyz",
            "abd/ln/xyz",
        ]

    Then, we select from them with this:

        >>> select_by_path_component("ab/ln/xyz", possible_matches)
        {"abd/ln/xyz"}

    Compare this with the result from

        >>> select_by_path_component("ab./ln/xyz", possible_matches)
        {"abc/ln/xyz", "abd/ln/xyz"}

    In the first case, only "ab" matched to the first component in any of
    possibilities, so the function narrowed down the list to just the largest of
    these: those paths where the first element is "abd".  In the second case, the
    "." did not literally match anything, so "ab." was compiled to a regex, which
    matched both "abc" and "abd" — which leave no remaining characters, and so
    taking the `max` element in each case does nothing.  On the other hand, if we
    had used a regex that only matched part of the first component, the `max` stage
    would still have eliminated some possibilities:

        >>> select_by_path_component("a./ln/xyz", possible_matches)
        {"abd/ln/xyz"}

    Only the "ab" portion was matched by this regex, so we continue searching only
    the paths starting with `max("abc", "abd")` — which is "abd".

    This process is repeated successively for each component, and any type of match
    may be used in any component.  For example, we can select with

        >>> select_by_path_component("ab./l./x", possible_matches)
        {"abc/lm/xyz", "abc/ln/xyz", "abd/lm/xyz", "abd/ln/xyz"}

    Here, both the first and second components use regexes.  Or we could use
    partial matching in the second component:

        >>> select_by_path_component("ab./l/x", possible_matches)
        {"abc/ln/xyz", "abd/ln/xyz"}

    Note that both "abc" and "abd" match for the first component, and their matches
    are evaluated separately when selecting the second component.  In both cases
    "ln" is the max matching component.

    """
    import collections
    import re

    # Look for exact matches
    if recursion_index == 0:
        matches = {p for p in possible_matches if p.startswith(path_pattern)}
        if matches:
            return matches

    # Separate the pattern into components
    split_path_pattern = path_pattern.split("/")

    # If we've searched all components of the pattern, return all possible matches
    if recursion_index > len(split_path_pattern):
        return possible_matches

    # Split the pattern into the part we're currently working on, and the rest
    path_working = "/".join(split_path_pattern[:recursion_index+1])
    path_suffix = "/".join(split_path_pattern[recursion_index+1:])
    if path_suffix:
        path_suffix = "/" + path_suffix

    # First, look for exact matches to beginning of string
    matches = collections.defaultdict(list)
    for f in possible_matches:
        if f.startswith(path_working):
            next_slash = f.find("/", f.index(path_working) + len(path_working) - 1)
            if next_slash < 0:
                next_slash = len(f)
            matched_component = f[:next_slash]
            matches[path_working].append((matched_component, f))

    # If that didn't give us anything, compile to a regex
    if not matches:
        path_regex = re.compile(path_working)
        component_regex = re.compile(path_working + r"[^/]*")
        for f in possible_matches:
            m = path_regex.match(f)
            if m:
                match = m.group()
                matched_component = component_regex.match(f).group()
                matches[match].append((matched_component, f))

    # Go through, finding the `max` for any group of partial matches
    best_matches = {}
    for match in matches:
        best_component = max(set(m[0] for m in matches[match]))
        best_possible_matches = [m[1] for m in matches[match] if m[0] == best_component]
        best_matches[best_component] = best_possible_matches

    # Recurse, if necessary
    selected = {
        path
        for best_component in best_matches
        for path in select_by_path_component(
            best_component + path_suffix,  # Search for found components place rest of pattern
            possible_matches=best_matches[best_component],
            recursion_index=recursion_index+1
        )
    }

    return selected
