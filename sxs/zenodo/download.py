def matching(*args, **kwargs):
    """Download all open-access SXS files matching the input arguments

    Files will be output into subdirectories named for the SXS ID of each matching record found on
    Zenodo.  For example, if the Zenodo record contains 'SXS:BBH:1234', the files will be placed
    into a subdirectory named 'SXS_BBH_1234'.  If only the highest Lev is desired (the default), the
    files will be placed directly into that subdirectory; otherwise, additional levels of
    subdirectories will be created as needed.

    Parameters
    ==========
    file: string or multiple strings as non-keyword arguments
        Zenodo file name to match.  This will be compiled as a regex, so it can include python-style
        regex matches.
    sxs_ids: list of strings (defaults to ['SXS:BBH:'])
        Keyword argument only.  SXS IDs to be searched for in the Zenodo deposit's title.  Each will
        be compiled as a regex, so it can include python-style regex matches.
    highest_lev_only: bool (defaults to True)
        Keyword argument only.  If True, only download only files from the highest-numbered Lev.

    """
    import re
    import requests
    from .. import sxs_id as sxs_id_finder
    from .api.records import Records
    from tqdm import tqdm

    print('args:', args)
    print('kwargs:', kwargs)

    file_name_matches = [re.compile(f) for f in args]
    sxs_ids = [re.compile(i) for i in kwargs.pop('sxs_ids', ['SXS:BBH:'])]
    highest_lev_only = kwargs.pop('highest_lev_only', True)
    lev_path_re = re.compile(r'/Lev[-0-9]*')

    def local_path(sxs_id, filename):
        """Return the local filename where you want to save this file"""
        if not filename.startswith(sxs_id):
            filename = sxs_id + '/' + filename
        filename = filename.replace(':', '_')
        if download_only_highest_lev:
            filename = re.sub(lev_path_re, '', filename)
        return filename

    def title_matches(title):
        for sxs_id in sxs_ids:
            if sxs_id.search(title):
                return True
        return False

    def file_matches(file_name):
        for file_name_match in file_name_matches:
            if file_name_match.search(file_name):
                return True
        return False

    all_records = Records.search(q='communities:sxs AND access_right:open')
    catalog = [record for record in all_records if title_matches(record.get('metadata', {}).get('title', {}))]

    for simulation in tqdm(catalog):
        try:  # We probably don't want this entire script to abort if something goes wrong with one simulation
            doi_url = simulation['doi_url']
            resolver = requests.get(doi_url)
            if resolver.status_code != 200:
                resolver.raise_for_status()
                raise RuntimeError()
            response = requests.get(resolver.url.replace('/record/', '/api/records/'), headers={"Accept": "application/json"})
            if response.status_code != 200:
                response.raise_for_status()
                raise RuntimeError()

            title = response.json()['metadata']['title']
            # if 'SXS:BBH:' not in title:
            #     continue  # Skip non-BBH systems
            sxs_id = sxs_id_finder(title)

            all_files = response.json()['files']

            horizons_files = [f for f in all_files if f['filename'].endswith('Horizons.h5')]
            highest_lev_horizons_file = sorted(horizons_files, key=lambda f: f['filename'])[-1]

            h_files = [f for f in all_files if f['filename'].endswith('rhOverM_Asymptotic_GeometricUnits_CoM.h5')]
            highest_lev_h_file = sorted(h_files, key=lambda f: f['filename'])[-1]

            print('Working on "{0}"'.format(sxs_id))

            if download_only_highest_lev:
                files_to_download = [highest_lev_h_file, highest_lev_horizons_file]
            else:
                files_to_download = h_files+horizons_files

            for file_description in files_to_download:
                url = file_description['links']['download']
                filename = file_description['filename']
                path = decide_where_to_put_this(sxs_id, filename)
                print('\tDownloading to "{0}"'.format(path))
                #download(url, path)

        except KeyboardInterrupt:  # Don't catch Ctrl-C, so you can actually interrupt this loop if you want
            raise

        except Exception as e:  # For anything else, just print the error, and continue
            traceback.print_tb(e.__traceback__)
