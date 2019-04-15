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
        regex matches or partial completions.
    sxs_ids: list of strings (defaults to ['SXS:BBH:'])
        Keyword argument only.  SXS IDs to be searched for in the Zenodo deposit's title.  Each will
        be compiled as a regex, so it can include python-style regex matches or partial completions.
    highest_lev_only: bool (defaults to True)
        Keyword argument only.  If True, only download only files from the highest-numbered Lev.

    """
    import traceback
    import re
    import requests
    from .. import sxs_id as sxs_id_finder
    from .api.records import Records
    from tqdm import tqdm

    file_name_matches = [re.compile(f) for f in args]
    sxs_ids = [re.compile(i) for i in kwargs.pop('sxs_ids', ['SXS:BBH:'])]
    highest_lev_only = kwargs.pop('highest_lev_only', True)
    lev_path_re = re.compile(r'Lev(?P<lev>[-0-9]*)/')

    def local_path(sxs_id, filename):
        """Return the local filename where you want to save this file"""
        if not filename.startswith(sxs_id):
            filename = sxs_id + '/' + filename
        filename = filename.replace(':', '_')
        if highest_lev_only:
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
            title = simulation['metadata']['title']
            sxs_id = sxs_id_finder(title)
            print('Working on "{0}"'.format(sxs_id))

            all_files = simulation['files']

            if highest_lev_only:
                files_to_download = {}
                for file_description in all_files:
                    filename = file_description['filename']
                    if file_matches(filename):
                        search = lev_path_re.search(filename)
                        if search:
                            lev = search['lev']
                            generic_filename = filename.replace('Lev{0}/'.format(lev), 'Lev{0}/')
                            files_to_download[generic_filename] = files_to_download.get(generic_filename, []) + [lev,]
                        else:
                            files_to_download[filename] = ['']
                files_to_download = [key.format(sorted(files_to_download[key])[-1]) for key in files_to_download]
                files_to_download = [f for f in all_files if f['filename'] in files_to_download]
            else:
                files_to_download = [f for f in all_files if file_matches(f['filename'])]

            for file_description in files_to_download:
                url = file_description['links']['download']
                filename = file_description['filename']
                path = local_path(sxs_id, filename)
                print('\tDownloading "{0}" to "{1}"'.format(filename, path))
                #download(url, path)

        except KeyboardInterrupt:  # Don't catch Ctrl-C, so you can actually interrupt this loop if you want
            raise

        except:  # For anything else, just print the error, and continue
            traceback.print_exc()
