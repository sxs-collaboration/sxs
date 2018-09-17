

def backup(catalog_file_name='complete_catalog.json', checksum_map_file_name='checksums.json',
           backup_directory='.', top_directory=None, dry_run=False):
    """Back up all Zenodo records, including files

    * catalog file - complete catalog JSON file
    * checksums file - JSON file mapping checksums to file sizes and full paths relative to this file
    * records directory - complete copies of each record in the 'sxs' Zenodo community
    * current directory - links named by SXS identifier pointing to current record directory
    * versions directory - links named by SXS identifier and version number pointing to corresponding record directory

    Parameters
    ==========
    catalog_file_name: str [defaults to 'complete_catalog.json']
        Relative or absolute path to catalog JSON file.  This is expected to have been created by
        the `catalog` function.
    checksum_map_file_name: str [defaults to 'checksums.json']
        Absolute or relative path to JSON file containing the checksums of all the files within
        `top_directory`.  This will be updated as part of this function.
    backup_directory: str [defaults to '.']
        Absolute or relative path to directory in which the backup will be stored.
    top_directory: str [defaults to value of backup_directory]
        Absolute or relative path to directory that will be searched for files to link to in the
        backup.  This must enclose the backup directory itself and may be used, for example, if the
        SimulationAnnex contains the original files, in which case the backup will not require
        downloads.
    dry_run: bool [defaults to False]
        If True, don't actually create any directories or file, or download anything.  Instead, just
        print all actions that will be taken.

    """
    import os.path

    if top_directory is None:
        top_directory = backup_directory
    else:
        if not os.path.realpath(backup_directory).startswith(os.path.realpath(top_directory)):
            error_message = 'Input backup directory "{0}" is not a subdirectory of input top directory "{1}"'
            raise ValueError(error_message.format(backup_directory, top_directory))


def main(*args, **kwargs):
    # __doc__ = backup.__doc__
    return backup(*args, **kwargs)
