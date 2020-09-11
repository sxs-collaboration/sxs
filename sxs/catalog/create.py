"""Create a new catalog of SXS data"""


def create(login=None):
    """Create a new catalog of SXS data

    Note that this completely reconstructs the catalog from scratch, which can take
    a long time, and requires a Zenodo login key.  Unless you are sure you really
    want to reconstruct the catalog, you probably just want to use the existing
    version, which can be downloaded and loaded automatically using
    `sxs.load("catalog")`.

    Parameters
    ----------
    login : sxs.zenodo.Login

    """
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

    print("Searching for all published records", flush=True)
    all_versions = l.search(q="communities:sxs", all_versions=True, status="published", size=9999)

    modified = max(r.get("modified", "") for r in all_versions)

    records = {
        r["doi_url"]: r
        for r in sorted(all_versions, key=lambda rec: rec["title"])
    }

    simulations = create_simulations(records, l)

    return Catalog({
        "catalog_file_description": catalog_file_description,
        "modified": modified,
        "records": records,
        "simulations": simulations,
    })
