# Copyright (c) 2020, Michael Boyle
# See LICENSE file for details:
# <https://github.com/sxs-collaboration/sxs/blob/master/LICENSE>

try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:  # pragma: no cover
    import importlib_metadata

__version__ = importlib_metadata.version(__name__)

from . import catalog, data, utilities


def load(location, /, download=False, cache=None, **kwargs):
    """Load an SXS-format dataset, optionally downloading and caching

    Parameters
    ----------
    location : local file path, SXS path, or URL
        This can be a relative or absolute path to a local file or an SXS ID
        followed by a path supplied with that SXS data -- for example,
        'SXS:BBH:1234/Lev5/h_Extrapolated_N2.h5'.  In the former case, all
        following parameters are ignored.  In the latter case, this file is first
        sought in the cache directory, or optionally downloaded from CaltechDATA.
    download : bool, optional
        If this is True and the data is recognized as starting with an SXS ID
        but cannot be found in the cache, the data will be downloaded from
        CaltechDATA.  Note that if this is True but `cache` is None, it will
        automatically be switched to True.
    cache: {None, str, bool}, optional
        If `location` is not found directly, the cache is used, with `location`
        appended to the cache path.  If this is a string, it is interpreted as a
        path to the cache directory.  If it is True, the cache directory is
        obtained from `sxs.utilities.get_sxs_directory`.

    Keyword Parameters
    ------------------
    All remaining parameters are passed to the `load` function responsible for the
    requested data.

    Notes
    -----
    This function can load data in various ways.

      1) Given an absolute or relative path to a local file, it just loads the
         data directly.

      2) Given an SXS path -- like 'SXS:BBH:1234/Lev5/h_Extrapolated_N2.h5' --
         it first looks in a local cache directory (see details below).

      3) If the SXS path is not found in the cache directory and `download` is
         set to `True` this function attempts to download the data from
         CaltechDATA.  Note that `download` must be explicitly set in this
         case, or a ValueError will be raised.

      4) An arbitrary URL to download.  Note that the scheme (https://, etc.)
         must be included

    Note that downloading is switched off by default, but if it is switched on
    (set to True), the cache is also switched on by default.

    """
    raise NotImplementedError()

    import contextlib
    import warnings
    import re
    import pathlib
    import tempfile
    import json
    import h5py

    sxs_formats = {
        # 'format name': (format_load_function, requires_json)
        'corotating_paired_xor': (rotating_paired_xor.load, False),
        'rotating_paired_xor': (rotating_paired_xor.load, False),
    }

    if download and cache is None:
        cache = True

    warnings.warn('I promised to parse URLs here and switch download to True if needed')
    path = pathlib.Path(location).expanduser().resolve()
    h5_path = path.with_suffix('.h5')
    json_path = path.with_suffix('.json')

    # See if `location` *starts with* something like SXS:BBH:1234
    sxs_matches = re.match(sxs.sxs_identifier_regex, location)

    with contextlib.ExitStack() as exit_stack:
        if not cache:
            temporary_directory = exit_stack.enter_context(tempfile.TemporaryDirectory())
            cache_path = pathlib.Path(temporary_directory)
        else:
            if cache is True:  # We want it to literally *be* True, not just equal True
                warnings.warn('Logic is missing here!!!')
                cache_path = pathlib.Path('~/shouldntbehere/../dir/name').expanduser().resolve()
            else:  # It's apparently a non-empty string
                cache_path = pathlib.Path(cache).expanduser().resolve()
            if not cache_path.exists():
                warnings.warn(f'Creating directory "{path}" to serve as cache for SXS files')
                cache_path.mkdir(parents=True, exist_ok=True)
                cached_path = cache_path / location

        if not h5_path.exists():
            if not sxs_matches:
                raise ValueError(f'File "{h5_path}" does not exist and is not recognized as an SXS dataset')
            h5_cached_path = cached_path.with_suffix('.h5')
            if not h5_cached_path.exists():
                if not cache:
                    warnings.warn('\nDownloading "{h5_path}" but not caching it; set `cache=True` if possible.')
                download(doi_prefix + h5_path, h5_cached_path)
            h5_file = exit_stack.enter_context(h5py.File(h5_cached_path, 'r'))
        else:
            h5_file = exit_stack.enter_context(h5py.File(h5_path, 'r'))

        sxs_format = h5_file.get('sxs_format', h5_file.get('format', None))
        if sxs_format is None:
            raise ValueError(f'No format information found in "{h5_path}"')
        load_function, requires_json = sxs_formats[sxs_format]

        if requires_json:
            if not json_path.exists():
                if not sxs_matches:
                    raise ValueError(f'File "{json_path}" does not exist and is not recognized as an SXS dataset')
                json_cached_path = cached_path.with_suffix('.json')
                if not json_cached_path.exists():
                    if not cache:
                        warnings.warn('\nDownloading "{json_path}" but not caching it; set `cache=True` if possible.')
                    download(doi_prefix + json_path, json_cached_path)
                with open(json_cached_path, 'r') as f:
                    json_file = json.load(f)
            else:
                with open(json_path, 'r') as f:
                    json_file = json.load(f)
            return load_function(h5_file, json_file)
        else:
            return load_function(h5_file)
