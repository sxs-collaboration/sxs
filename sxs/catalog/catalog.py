import functools

class Catalog(object):
    url = "https://data.black-holes.org/catalog.json"

    def __init__(self, catalog=None, **kwargs):
        import collections

        self._dict = catalog or type(self).load(**kwargs)

        # Add version numbers to all records
        concept_to_versions = collections.defaultdict(list)
        for doi_url, record in self.records.items():
            concept_to_versions[record["conceptrecid"]].append(doi_url)
        for conceptrecid, doi_urls in concept_to_versions.items():
            for v, doi_url in enumerate(sorted(doi_urls), start=1):
                self.records[doi_url]["version"] = v

    @classmethod
    @functools.lru_cache()
    def load(cls, download=None, progress=None, **kwargs):
        """Load the SXS catalog

        Note that — unlike most SXS data files — the catalog file is updated
        frequently.  As a result, this function — unlike the loading functions for most
        SXS data files — will download the catalog by default each time it is called.
        However, also note that this function is itself cached, meaning that the same
        dict will be returned on each call in a given python session.  If you want to
        avoid that behavior, use `Catalog.reload`.

        Parameters
        ----------
        download : {None, bool}, optional
            If False, this function will look for the catalog in the sxs cache and
            raise an error if it is not found.  If True, this function will download
            the catalog and raise an error if the download fails.  If None (the
            default), it will try to download the file, warn but fall back to the cache
            if that fails, and only raise an error if the catalog is not found in the
            cache.  Note that this ignores the sxs configuration file entirely.
        progress : {None, bool}, optional
            If True, and a nonzero Content-Length header is returned, a progress bar
            will be shown during any downloads.  Default is False unless
            `read_config("download_progress")` is True.

        See Also
        --------
        sxs.sxs_directory : Locate cache directory
        Catalog.reload : Avoid caching the result of this function

        """
        import json
        import pathlib
        import tempfile
        import zipfile
        from datetime import datetime, timezone
        from .. import sxs_directory
        from ..utilities import download_file

        if progress is None:
            progress = read_config("download_progress")

        cache_path = sxs_directory("cache") / "catalog.zip"

        if cache_path.exists():
            if_newer = datetime.utcfromtimestamp(
                cache_path.stat().st_mtime
            ).replace(tzinfo=timezone.utc)
        else:
            if_newer = False

        if download or download is None:
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_path = pathlib.Path(temp_dir) / "catalog.json"
                zip_path = pathlib.Path(temp_dir) / "catalog.zip"
                try:
                    downloaded_path = download_file(cls.url, temp_path, progress=progress, if_newer=if_newer)
                except Exception as e:
                    if download:
                        raise RuntimeError(f"Failed to download '{cls.url}'") from e
                    download_failed = e  # We'll try the cache
                else:
                    download_failed = False
                    if downloaded_path == temp_path:
                        with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_BZIP2) as catalog_zip:
                            catalog_zip.write(temp_path, arcname="catalog.json")
                        zip_path.replace(cache_path)

        if not cache_path.exists():
            if download_failed:
                raise ValueError(f"Catalog not found in '{cache_path}' and download failed") from download_failed
            elif download is False:  # Test if it literally *is* False, rather than just casts to False
                raise ValueError(f"The catalog was not found in '{cache_path}', and downloading was turned off")
            else:
                raise ValueError(f"Catalog not found in '{cache_path}' for unknown reasons")

        try:
            with zipfile.ZipFile(cache_path, "r") as catalog_zip:
                try:
                    with catalog_zip.open("catalog.json") as catalog_json:
                        try:
                            catalog = json.load(catalog_json)
                        except Exception as e:
                            raise ValueError(f"Failed to parse 'catalog.json' in '{cach_path}'") from e
                except Exception as e:
                    raise ValueError(f"Failed to open 'catalog.json' in '{cache_path}'") from e
        except Exception as e:
            raise ValueError(f"Failed to open '{cache_path}' as a ZIP file") from e

        return cls(catalog)

    @classmethod
    def reload(cls, download=True, progress=None, **kwargs):
        """Reload the SXS catalog, without caching

        Clears the cache of `Catalog.load` and returns the result of calling it again.
        Note that in this function, the default value of `download` is `True`, rather
        than `None` as in `Catalog.load` — though both behaviors are available.

        Parameters
        ----------
        download : {None, bool}, optional
            If False, this function will look for the catalog in the sxs cache and
            raise an error if it is not found.  If True (the default), this function
            will download the catalog and raise an error if the download fails.  If
            None (the default), it will try to download the file, warn but fall back to
            the cache if that fails, and only raise an error if the catalog is not
            found in the cache.  Note that this ignores the sxs configuration file
            entirely.
        progress : {None, bool}, optional
            If True, and a nonzero Content-Length header is returned, a progress bar
            will be shown during any downloads.  Default is False unless
            `read_config("download_progress")` is True.

        See Also
        --------
        sxs.sxs_directory : Locate cache directory
        Catalog.load : Caching version of this function

        """
        cls.load.cache_clear()
        return cls.load(download=download, progress=progress, **kwargs)

    def save(self, **kwargs):
        raise NotImplementedError()

    def select_files(path_pattern):
        """Return paths and file information from all files available in the catalog

        This function is essentially the same as the `select` function, except that
        — whereas that function returns only the matching paths — this function
        returns matching paths and file information about those paths, including
        things like checksum, filesize, and download links.

        Parameters
        ----------
        path_pattern : str
            A pattern to search for among the catalog files.  See the docstring of
            `select` for details about how this pattern is used.

        Returns
        -------
        selections : dict
            This is a dictionary with keys given by the paths to selected files in
            the catalog, and values given by `file_info` dicts described below.

        See Also
        --------
        select : Return only the paths, not the file info

        Notes
        -----
        Because multiple versions of the same simulation may exist with some
        identical files, and because we don't want to reproduce those identical
        files in each version, many files in the catalog should just be links to
        those files in previous versions of their simulations.  This relationship
        is encapsulated by the `truepath` element of the `file_info` dict.

        The `file_info` dict returned for each selection contains the following
        key:value pairs:

          * checksum: MD5 sum of the file
          * filename: Name of this file in this version of the simulation.
          * filesize: Size of the file in bytes
          * download: The URL to use to download this file
          * truepath: The full SXS path to the original file

        """
        files = self.files
        selection = self.select(path_pattern, files=files)
        return {s: files[s] for s in selection}

    def select(self, path_pattern, files=None):
        """Select from all catalog files by progressively matching path components

        Parameters
        ----------
        path_pattern : str
            A pattern to search for among the catalog files.  This is first searched as
            a literal substring, which will return a set of one or more matches.  If no
            matches were found, the pattern is split into components by "/", which are
            used progressively to match corresponding path components — first as
            literal substrings, then as (python-style) regular expressions.  Each
            partial match in this step will pass all the matched components to the
            `max` function, and choose only the result — though with regexes, there may
            be distinct partial matches, each of which will pass only its corresponding
            components to `max`, resulting in multiple matches.
        files : Iterable[str], optional
            A list of files to select from.  Default is all files in the catalog.

        Returns
        -------
        matched_paths : Set[str]

        See Also
        --------
        sxs.utilities.select_by_path_component

        Examples
        --------
        First, we can choose the `h` waveform with `n=2` extrapolation in the
        highest-resolution (Lev) run from the simulation SXS:BBH:0002 with

            >>> catalog.select("SXS:BBH:0002/Lev/h_ex.*_n2")
            {"SXS:BBH:0002v7/Lev6/h_extrapolated_n2.h5"}

        Because the "Lev" component of the input string only matched the "Lev" portion
        of the path components "Lev4", "Lev5", and "Lev6", all of those components were
        passed to `max` together, and only "Lev6" was returned.  Similarly, the highest
        version of SXS:BBH:0002 (currently, v7) was chosen automatically.

        We could, instead, use a regular expression to include the numbers in the
        matches:

            >>> catalog.select("SXS:BBH:0002/Lev./h_ex.*_n2")
            {"SXS:BBH:0002v7/Lev4/h_extrapolated_n2.h5",
             "SXS:BBH:0002v7/Lev5/h_extrapolated_n2.h5",
             "SXS:BBH:0002v7/Lev6/h_extrapolated_n2.h5"}

        In this case, the pattern "Lev." matches each component entirely, so they are
        each retained.  Similarly, "SXS:BBH:0002v." would match all versions, and
        "h_ex.*_n." would match all extrapolation orders — and "h.*" would match all
        extrapolation orders, as well as the outermost extraction.

        Another helpful pattern is to use alternation in the regular expression, to
        explicitly list acceptable matches:

            >>> catalog.select("SXS:BBH:0002/Lev(4|6)/h_ex.*_n2")
            {"SXS:BBH:0002v7/Lev4/h_extrapolated_n2.h5",
             "SXS:BBH:0002v7/Lev6/h_extrapolated_n2.h5"}

        """
        from ..utilities import select_by_path_component
        files = files or self.files
        selections = select_by_path_component(location, set(files))

    @property
    @functools.lru_cache()
    def files(self):
        import re
        import pathlib
        from ..utilities import sxs_identifier_regex, sxs_id

        sxs_id_regex = re.compile(sxs_identifier_regex + r"/?")

        file_infos = {}
        for record in self.records.values():
            sxs_sim_id = sxs_id(record.get("title", ""))
            if not sxs_sim_id:
                continue
            version = record["version"]
            prefix = pathlib.Path(f"{sxs_sim_id}v{version}")
            files = record["files"]
            for file in files:
                path_str = str(prefix / sxs_id_regex.sub("", file["filename"], count=1))
                file_info = {
                    "checksum": file["checksum"],
                    "filename": file["filename"],
                    "filesize": int(file["filesize"]),
                    "download": file["links"]["download"],
                }
                if path_str in file_infos:
                    raise ValueError(
                        f"\nFile collision involving '{path_str}':\n"
                        f"    file_info 1: {file_infos[path_str]}\n"
                        f"    file_info 2: {file_info}"
                    )
                file_infos[path_str] = file_info

        unique_files = collections.defaultdict(list)
        for k, v in file_infos.items():
            unique_files[f"{v['checksum']}{v['filesize']}"].append(k)

        original_paths = {k: min(v) for k, v in unique_files.items()}

        for v in file_infos.values():
            v["truepath"] = original_paths[f"{v['checksum']}{v['filesize']}"]

        return file_infos

    @property
    def description(self):
        return self._dict['description']

    @property
    def modified(self):
        return self._dict['modified']

    @property
    def records(self):
        return self._dict['records']

    @property
    def simulations(self):
        return self._dict['simulations']
