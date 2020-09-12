import functools


# noinspection SpellCheckingInspection
class Catalog(object):
    url = "https://data.black-holes.org/catalog.json"

    def __init__(self, catalog=None, **kwargs):
        self._dict = catalog or type(self).load(**kwargs)

    @classmethod
    @functools.lru_cache()
    def load(cls, download=None):
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
        from .. import sxs_directory, read_config
        from ..utilities import download_file

        progress = read_config("download_progress") or False

        cache_path = sxs_directory("cache") / "catalog.zip"

        if cache_path.exists():
            if_newer = cache_path
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
                    if temp_path.exists():
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
                            raise ValueError(f"Failed to parse 'catalog.json' in '{cache_path}'") from e
                except Exception as e:
                    raise ValueError(f"Failed to open 'catalog.json' in '{cache_path}'") from e
        except Exception as e:
            raise ValueError(f"Failed to open '{cache_path}' as a ZIP file") from e

        return cls(catalog)

    @classmethod
    def reload(cls, download=True):
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

        See Also
        --------
        sxs.sxs_directory : Locate cache directory
        Catalog.load : Caching version of this function

        """
        cls.load.cache_clear()
        return cls.load(download=download)

    def save(self, **kwargs):
        raise NotImplementedError()

    def select_files(self, path_pattern):
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

            >>> catalog = Catalog.load()
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
        selections = select_by_path_component(path_pattern, set(files))
        return sorted(selections)

    @property
    @functools.lru_cache()
    def files(self):
        import collections
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
            files = sorted(record["files"], key=lambda f: f["filename"])
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
        return self._dict['catalog_file_description']

    @property
    def modified(self):
        return self._dict['modified']

    @property
    def records(self):
        return self._dict['records']

    @property
    def simulations(self):
        return self._dict['simulations']

    @property
    def open_access(self):
        """Return a catalog containing only open-access information

        Note that the returned object contains references to the original; meaning that
        changing the returned object could change the original.  Use `copy.deepcopy` to
        construct a copy with separate data.

        """
        from .. import sxs_id
        from . import catalog_file_description

        open_records = {
            k: v for k, v in self.records.items()
            if v.get("metadata", {}).get("access_right", "closed") == "open"
        }

        open_sxs_ids = set(sxs_id(r["title"]) for r in open_records.values())

        open_simulations = {
            k: v for k, v in self.simulations.items()
            if k in open_sxs_ids
        }

        open_modified = max(r.get("modified", "") for r in open_records.values())

        open_catalog = type(self)({
            "catalog_file_description": catalog_file_description,
            "modified": open_modified,
            "records": open_records,
            "simulations": open_simulations,
        })

        return open_catalog

    def split_and_write(self, directory, set_atime_and_mtime=True):
        """Write public and private JSON files

        This function writes four JSON files to the given directory: complete and
        open-access-only versions of the full catalog (including records) and
        simulations-only files.

        Parameters
        ----------
        directory : {str, pathlib.Path}
            Directory in which to place all four files.  If it doesn't exist it will be
            created.
        set_atime_and_mtime : bool, optional
            If True, we try to set the access and modification times of the output
            files to match the last modification time of the data itself.  If this
            fails, an exception will be raised.

        """
        import json
        import pathlib
        import os
        from datetime import datetime

        dir_path = pathlib.Path(directory).expanduser().resolve()
        public_catalog_path = dir_path / "public_catalog.json"
        public_simulations_path = dir_path / "public_simulations.json"
        private_catalog_path = dir_path / "private_catalog.json"
        private_simulations_path = dir_path / "private_simulations.json"

        dir_path.mkdir(parents=True, exist_ok=True)
        if not os.access(str(dir_path), os.W_OK) or not dir_path.is_dir():
            raise OSError(f"Could not create writable directory '{dir_path}'")

        public_catalog = self.open_access
        private_catalog = self

        with public_catalog_path.open("w") as f:
            json.dump(public_catalog._dict, f, indent=2, separators=(",", ": "), ensure_ascii=True)
        with public_simulations_path.open("w") as f:
            json.dump(public_catalog._dict["simulations"], f, indent=None, separators=(",", ":"), ensure_ascii=True)
        with private_catalog_path.open("w") as f:
            json.dump(private_catalog._dict, f, indent=2, separators=(",", ": "), ensure_ascii=True)
        with private_simulations_path.open("w") as f:
            json.dump(private_catalog._dict["simulations"], f, indent=None, separators=(",", ":"), ensure_ascii=True)

        if set_atime_and_mtime:
            public_modified = public_catalog.modified
            private_modified = private_catalog.modified
            if public_modified:
                public_modified = int(datetime.timestamp(datetime.strptime(public_modified, "%Y-%m-%dT%H:%M:%S.%f")))
                times = (public_modified, public_modified)
                os.utime(str(public_catalog_path), times)
                os.utime(str(public_simulations_path), times)
            if private_modified:
                private_modified = int(datetime.timestamp(datetime.strptime(private_modified, "%Y-%m-%dT%H:%M:%S.%f")))
                times = (private_modified, private_modified)
                os.utime(str(private_catalog_path), times)
                os.utime(str(private_simulations_path), times)
