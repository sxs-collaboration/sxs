def cite(*sxs_ids, bibtex=True):
    """Cite this package and/or data

    Prints out a citation for this package, the most recent paper
    describing the catalog, the catalog data itself, and optionally
    individual simulations.

    Note that this function makes web requests to obtain the DOIs and
    corresponding BibTeX entries.

    If you want to cite a specific version of the SXS catalog data, be
    sure to load it with `sxs.load("dataframe", tag="v3.0.0")` or
    similar *before* calling this function.  Otherwise, the most recent
    version will be used.
    
    Parameters
    ----------
    sxs_ids : str
        Any number of SXS IDs to include in the citation.
    bibtex : bool, optional
        If True — the default — returns full BibTeX entries.
        Otherwise, just returns a list of DOIs to cite.
    
    Returns
    -------
    str or list[str]
        A string containing BibTeX entries, or a list of DOI strings.
    
    """
    from . import doi_prefix, load, __version__
    from . import sxs_id as sxs_identifier

    simulations = load("simulations")

    # Get the DOI for this version of this package
    package_version = f"v{__version__}"
    package_doi = get_zenodo_doi("The sxs package", package_version)

    # Get the DOI for the paper
    paper_doi = f"{doi_prefix}/SXSCatalog3"

    # Get the DOI for this version of the catalog
    catalog_tag = getattr(simulations, "tag", "v3.0.0")
    catalog_doi = get_zenodo_doi(
        "The SXS catalog of simulations",
        catalog_tag
    )

    sxs_id_references = {
        doi
        for sxs_id in sxs_ids
        for doi in simulations[sxs_identifier(sxs_id)].get("citation_dois", [])
        if doi.startswith("10.")
    } - {package_doi, paper_doi, catalog_doi}

    if bibtex:
        package_bibtex = doi2bibtex(package_doi, key=f"SXSPackage_{package_version}")
        paper_bibtex = doi2bibtex(paper_doi, key="SXSCatalogPaper_3")
        catalog_bibtex = doi2bibtex(
            catalog_doi,
            key=f"SXSCatalogData_{catalog_tag[1:]}",
            title_suffix=f" {catalog_tag}"
        )
        sxs_id_references_bibtex = [
            doi2bibtex(doi) for doi in sxs_id_references
        ]
        sxs_id_bibtex = [
            doi2bibtex(f"{doi_prefix}/{sxs_id}", key=sxs_id)
            for sxs_id in sxs_ids
        ]
        return "\n".join(
            [package_bibtex, paper_bibtex, catalog_bibtex]
            + sxs_id_references_bibtex
            + sxs_id_bibtex
        )
    else:
        sxs_id_references_dois = list(sxs_id_references)
        sxs_id_dois = [f"{doi_prefix}/{sxs_id}" for sxs_id in sxs_ids]
        return [package_doi, paper_doi, catalog_doi] + sxs_id_references_dois + sxs_id_dois


def doi2bibtex(doi, key="", title_suffix=""):
    """Convert a DOI to a BibTeX entry

    This function queries doi.org — and possibly the service it
    redirects to — for the BibTeX entry corresponding to the DOI.

    Note that there are some services that do not provide good BibTeX.
    In particular, MNRAS does not support this service very well.

    Parameters
    ----------
    doi : str
        The DOI to convert.
    key : str, optional
        The BibTeX key to use.  If not provided, the key will be
        generated from the DOI.  If multiple entries are found, this
        function will raise an exception.
    
    Returns
    -------
    str
        The BibTeX entry as a string.
    
    """
    import requests
    import bibtexparser
    from .utilities import sxs_identifier_re

    url = f"https://doi.org/{doi}"
    headers = {"Accept": "application/x-bibtex"}
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    bibtex = response.text

    bibtex_format = bibtexparser.BibtexFormat()
    bibtex_format.indent = "  "
    bibtex_format.block_separator = "\n"

    try:
        library = bibtexparser.parse_string(bibtex)

        for i,entry in enumerate(library.entries):
            # Replace the key, if requested
            if key:
                if i>1:
                    raise Exception("Multiple entries found in BibTeX, but one key was given.")
                entry.key = key

            # Ensure that the "author" field does not end in ","
            if (author := entry.fields_dict.get("author", "")):
                if author.value.rstrip().endswith(","):
                    author.value = author.value.rstrip()[:-1]
                if author.value == "SXS Collaboration":
                    author.value = "{SXS Collaboration}"
                entry.set_field(author)

            # Now, look for an SXS ID in the `title` field; if present
            # surround it with curly braces to prevent BibTeX from
            # downcasing it.
            title = entry.fields_dict.get("title", "")
            if (m:=sxs_identifier_re.search(title.value)):
                # Surround the SXS ID with curly braces
                replacement = f"{{{m.group('sxs_identifier')}}}"

                # The title will not (currently) include the version,
                # but if the requested DOI contains a *versioned* SXS
                # ID, we add the version (e.g., "v3.0") to the
                # replacement string
                if (
                    (match:=sxs_identifier_re.search(doi))
                    and match.group("sxs_identifier") == m.group("sxs_identifier")
                    and match.group("version")
                ):
                    replacement += "v" + match.group("version")

                title.value = title.value.replace(m.group("sxs_identifier"), replacement)
                entry.set_field(title)
            elif title_suffix:
                title.value += title_suffix
                entry.set_field(title)

        bibtex = bibtexparser.write_string(library, bibtex_format=bibtex_format)

    except Exception as e:
        print("Failed to replace key and/or escape titles in BibTeX entry:\n")
        for line in bibtex.split("\n"):
            print("    " + line)
        print("")
        raise

    return bibtex


def get_zenodo_doi(title, version):
    import requests

    params = {
        "q": f'metadata.title:"{title}" AND metadata.version:"{version}"',
        "size": 1,
        "page": 1,
        "sort": "mostrecent",
        "all_versions": "",  # Need to include this to get older versions
    }
    url = "https://zenodo.org/api/records"
    r = requests.get(url, params=params)

    if r.status_code != 200:
        print(f"An unknown error occurred when trying to access '{url}'.")
        print(f"The search parameters were '{params}'")
        try:
            print(r.json())
        except:
            pass
        r.raise_for_status()
        raise RuntimeError()  # Will only happen if the response was not strictly an error

    json = r.json()

    if json.get("hits", {}).get("total", 0) == 0:
        print(f"No results found for '{title}' version '{version}'.")
        return None
    else:
        hit = json["hits"]["hits"][0]
        return hit["doi"]
