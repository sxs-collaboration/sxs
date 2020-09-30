"""Create a new catalog of SXS data"""


def _create(login=None):
    """Create a new catalog of SXS data

    WARNING: This function is private; you almost certainly don't want to use it.
    It completely reconstructs the catalog from scratch, which can take a long
    time, and requires a Zenodo login key.  And unless your login key gives you
    access to records restricted to SXS members, the resulting catalog will
    probably have incorrect version numbers.

    Unless you are sure you really want to reconstruct the catalog, you probably
    just want to use the existing version of the catalog, which can be downloaded
    and loaded automatically using `sxs.load("catalog")`.

    Parameters
    ----------
    login : {None, sxs.zenodo.Login}, optional
        If not present, this uses the default constructor of sxs.zenodo.Login.

    See Also
    --------
    sxs.load : Call this with "catalog" to load the existing version.
    sxs.zenodo.Login

    """
    import collections
    from ..zenodo import Login
    from .description import catalog_file_description
    from .catalog import Catalog

    def create_simulations(records, login):
        """Create a dictionary of simulations with information for downloading SXS metadata"""
        from collections import defaultdict
        from tqdm.auto import tqdm
        from .. import sxs_id, load, Metadata

        def get_highest_metadata(sxs_sim_id, record, login):
            if record.get("metadata", {}).get("access_right", "") == "open":
                metadata = load(f"{sxs_sim_id}/Lev/metadata.json")
            else:
                highest_lev_metadata_json = max(
                    [f for f in record.get("files", []) if "/metadata.json" in f["filename"]],
                    default={}, key=lambda f: f["filename"]
                )
                download_url = highest_lev_metadata_json.get("links", {}).get("download", "")
                if not download_url:
                    return {}
                metadata = Metadata(login.session.get(download_url).json())
            url = record.get("links", {}).get("conceptdoi", record["doi_url"])
            if url:
                metadata["url"] = url
            return metadata.reorder_keys()

        version_map = defaultdict(list)
        for r in records.values():
            sxs_sim_id = sxs_id(r["title"])
            if sxs_sim_id:  # and r.get("metadata", {}).get("access_right", "closed") == "open":
                version_map[sxs_sim_id].append(r)
        simulations = {
            sxs_sim_id: m
            for sxs_sim_id, versions in tqdm(version_map.items(), total=len(version_map), dynamic_ncols=True)
            for m in [dict(get_highest_metadata(sxs_sim_id, max(versions, key=lambda r: r["id"]), login))] if m
        }
        return simulations

    # If login is None, this creates a Login object to use
    l = login or Login()

    # Search for *all* versions — even unpublished ones, to get the versions right.  Note that we
    # have to run these queries separately because there are more than 10,000 if combined, which
    # exceeds zenodo's limit.  Currently, this hack works to get them all — though it will fail if
    # more drafts are published.
    print("Searching for all versions of all records...", flush=True)
    published = l.search(q="communities:sxs", all_versions=True, status="published")
    unpublished = l.search(q="communities:sxs", all_versions=True, status="draft")
    all_versions = published + unpublished

    # Figure out the version numbers
    concept_to_versions = collections.defaultdict(list)
    for record in all_versions:
        concept_to_versions[record["conceptrecid"]].append(record["doi_url"])
    doi_url_to_version = {}
    for conceptrecid, doi_urls in concept_to_versions.items():
        for version, doi_url in enumerate(sorted(doi_urls), start=1):
            doi_url_to_version[doi_url] = version
    for record in all_versions:
        record["version"] = doi_url_to_version[record["doi_url"]]

    # Make it into a dictionary sorted by title and version, dropping unpublished
    records = {
        r["doi_url"]: r
        for r in sorted(all_versions, key=lambda rec: (rec["title"], rec["version"]))
        if r.get("state", "error") == "done" and bool(r.get("submitted", False)) == True
    }

    # Get the latest modification time of
    modified = max(r.get("modified", "") for r in records.values())

    print("Downloading metadata files:", flush=True)
    simulations = create_simulations(records, l)

    return Catalog({
        "catalog_file_description": catalog_file_description,
        "modified": modified,
        "records": records,
        "simulations": simulations,
    })
