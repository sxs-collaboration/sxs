#! /usr/bin/env python

"""Make nicely named symbolic links for entries in the SpEC waveform catalog

This script is a wrapper around the `symlink_runs` function.  See that function's documentation for
details.

NOTE: This is a modified version of the older SymlinkRuns script found in the SimulationAnnex
itself.  The modifications:
  1) allow this to run under python 3
  2) add an option to use relative links rather than absolute (so that the catalog can be
     mounted as a volume in docker, for example)
  3) add an option to place links in another directory (useful for docker)
  4) require the "nice" name's prefix to be given as a regex to match

If the symlinks are desired for the `Catalog` directory *and* the `Incoming` directory (probably
with prefixes `^SXS:` and `^PRIVATE:`, respectively, either run the command twice or allow for
either possibility in the regex.

"""

def symlink_runs(source_directory='Catalog', target_directory='CatalogLinks',
                 remove_old_target_dir=False, prefix='^SXS:',
                 use_relative_links=False, relative_directory_path=None, verbosity=2):
    """Make nicely named symbolic links for entries in the SpEC waveform catalog

    Use this function to write a directory of symlinks with "nice" names like `SXS:BBH:0001` for
    entries in the `Catalog` directory, or names like `PRIVATE:BBH:0001` for entries in the
    `Incoming` directory.  These "nice" names are taken from any `metadata.json` or `metadata.txt`
    file.

    Parameters
    ----------
    source_directory: str (default: 'Catalog')
        Name of directory (either absolute or relative to working directory) to traverse looking for
        metadata files containing the `alternative-names` key with a matching prefix (given below).
    target_directory: str (default: 'CatalogLinks')
        Name of directory (either absolute or relative to working directory) in which to store the
        resulting links.  If this does not exist, it is created.
    remove_old_target_dir: bool (default: False)
        If True, remove and recreate the `target_directory`, so that old links are gone.
    prefix: str (default: '^SXS:')
        Name of prefix to be searched for in `alternative-names`.  If a match is found, that name is
        used as the name of the "nice" link.
    use_relative_links: bool (default: False)
        If True, use relative paths in the links instead of absolute paths.  This can be helpful,
        for example, when the catalog is mounted as a volume in a docker container, so that the
        names can change.
    relative_directory_path: str (default: None)
        If this path is not None, use relative links and replace the relative part of the path
        between `source_directory` and `target_directory` with this string.  This forces
        `use_relative_links` to be True.  Again, this can be helpful when the source and target
        directories are mounted as volumes in a docker container, but under different names.
    verbosity: int 0, 1, or 2 (default: 2)
        Output nothing if verbosity is 0, only duplicate or missing directories if verbosity is 1,
        and all directories if verbosity is 2.

    """
    import os
    import sys
    import shutil
    import re
    from parse_metadata import ParseMetadataTree, error, warning

    if not os.path.isdir(source_directory):
        raise ValueError("source_directory '{0}' not found".format(source_directory))
    if remove_old_target_dir and os.path.isdir(target_directory):
        shutil.rmtree(target_directory)
    if not os.path.isdir(target_directory):
        os.mkdir(target_directory)
    prefix_regex = re.compile(prefix)
    if relative_directory_path is not None:
        use_relative_links = True
    if verbosity < 1:
        verbosity = 0
    else if verbosity > 1:
        verbosity = 2
    

    # DirsToLink = {"Catalog": "SXS"}
    # # DirsToLink = {"Catalog": "SXS", "Incoming": "PRIVATE"}
    # LinkDirs = [x+'Links' for x in DirsToLink]
    
    # Dups = {}
    # NoID = []
    # for Dir in DirsToLink:
    #     LinkDir = Dir+'Links'

    #     if store_links_in is not None:
    #         target_directory = os.path.join(store_links_in, LinkDir)
    #     else:
    #         target_directory = LinkDir

        # # Recreate Link directories
        # if os.path.isdir(target_directory):
        #     shutil.rmtree(target_directory)
        # os.mkdir(target_directory)

        # Parse the metadata
        print("Parsing the root for {}".format(Dir))
        reader = ParseMetadataTree(Dir, array=True)

        for run in reader:
            RunMetadata = reader.HighestLev(run)
            RunPath = RunMetadata['full-path']
            src = os.path.realpath(os.path.join(RunPath,".."))

            try:
                prefix = DirsToLink[Dir]
                RunName = next(ID for ID in RunMetadata['alternative-names']
                               if re.match(prefix, ID) or re.match("SXS", ID))
            except StopIteration:
                # If no ID was found, use the relative path to Dir since it is unique
                NoID.append(RunPath)
                RunName = os.path.relpath(src, os.path.realpath(Dir)).replace('/', '_')

            # Set the symlink destination (accounting for duplicate RunNames)
            dest = os.path.join(LinkDir, RunName)
            while os.path.exists(dest):
                duplicates = [src, os.path.realpath(dest)]
                try:
                    Dups[RunName].extend(duplicates)
                except KeyError:
                    Dups[RunName] = duplicates
                dest += "A"

            if use_relative_links:
                if use_relative_links_from:
                    src = os.path.join(use_relative_links_from, os.path.relpath(src, Dir))
                else:
                    src = os.path.relpath(src, os.path.dirname(dest))
                dest = os.path.join(target_directory, os.path.relpath(dest, LinkDir))
            os.symlink(src,dest)
            print(dest + " <- " + src)

    for ID in Dups:
        warning("Duplicate runs named {}:\n{}".format(ID,"\n".join(set(Dups[ID]))))

    if NoID:
        warning("Missing case numbers in:\n{}".format("\n".join(NoID)))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument("--use_relative_links",
                        help="Use relative links rather than absolute (default: absolute links)",
                        action="store_true")
    parser.add_argument("--use_relative_links_from",
                        help="Allow specification of the relative link prefix",
                        default=None)
    parser.add_argument("--store_links_in",
                        help="Place the links in the given directory, rather than the current directory",
                        default=None)
    args = parser.parse_args()

    main(use_relative_links=args.use_relative_links,
         use_relative_links_from=args.use_relative_links_from,
         store_links_in=args.store_links_in)
