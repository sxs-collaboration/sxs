from pathlib import Path
from .. import sxs_id, Metadata, sxs_directory
from ..utilities import sxs_identifier_re
from ..zenodo import path_to_invenio as p2i

def file_upload_allowed(file, directory_listing):
    """Return True if the file should be uploaded
    
    A file should be uploaded if
        * it is named "metadata.json" or "Horizons.h5"
        * it is named "Strain_*.json" or "ExtraWaveforms.json" and the corresponding
          ".h5" file is in the directory listing
        * it is named "Strain_*.h5" or "ExtraWaveforms.h5" and the corresponding
          ".json" file is in the directory listing
    
    """
    # Check `file.name` to ignore the directory
    if file.name in ["metadata.json", "Horizons.h5"]:
        return True
    if file.name.startswith("Strain_") or file.name.startswith("ExtraWaveforms"):
        # Ensure that both `.h5` and `.json` exist for all such files
        if file.suffix == ".json":
            return file.with_suffix(".h5") in directory_listing
        elif file.suffix == ".h5":
            return file.with_suffix(".json") in directory_listing
        else:
            return False
    return False


def files_to_upload(directory, annex_dir="."):
    """Return a list of files to upload

    The files to upload are those that are in the directory listing
    and pass the `file_upload_allowed` function.

    """
    full_directory = annex_dir / Path(directory)
    files = []
    for lev in full_directory.resolve().glob("Lev*"):
        directory_listing = list(lev.iterdir())
        files.extend([
            file for file in directory_listing
            if file_upload_allowed(file, directory_listing)
        ])
    return sorted(files, key=lambda x: str(x).lower())


def extract_id_from_common_metadata(file, annex_dir):
    """Extract the SXS ID from a common-metadata.txt file
    
    If the ID doesn't exist, return the directory path, relative to
    the `annex_dir`.
    """
    file = Path(file)
    annex_dir = Path(annex_dir)
    key = str(file.resolve().parent.relative_to(annex_dir.resolve()))
    with file.open("r") as f:
        for line in f.readlines():
            line = line.strip()
            if "alternative-names" in line:
                if (m := sxs_identifier_re.search(line)):
                    key = m["sxs_identifier"]
                    break
    return key


def local_simulations(annex_dir):
    """
    Walk the annex directory to find and process all simulations

    For each `common-metadata.txt` file found:
        - Ensures that at least one directory starting with "Lev"
          exists; if not, the process is skipped.
        - Defines a key for the metadata, which is either:
            - The SXS ID contained in that file's "alternative-names"
              field, if present.
            - The directory path relative to `annex_dir`.
        - Chooses the highest "Lev" directory and extracts the
          metadata.
        - Finds all files to upload in the directory; if none are
          found, the process is skipped.
        - Adds the "files" dictionary to the metadata, pointing to
          each file that would be uploaded if the simulation were
          published.

    Args:
        annex_dir (str or Path): The path to the annex directory to be
        processed.

    Returns:
        dict: A dictionary containing the processed simulations
        metadata.
    """
    from os import walk

    simulations = {}
    annex_dir = Path(annex_dir).resolve()

    # The `walk` method can be made *much* faster than the `glob` method
    for dirpath, dirnames, filenames in walk(annex_dir, topdown=True):
        dirpath = Path(dirpath)

        # Ignore hidden directories
        if dirpath.name.startswith("."):
            dirnames[:] = []
            continue

        if "common-metadata.txt" in filenames:
            if not any(d.startswith("Lev") for d in dirnames):
                continue

            key = extract_id_from_common_metadata(dirpath / "common-metadata.txt", annex_dir)

            # Find the highest Lev directory and extract the metadata
            highest_lev = sorted(
                [d for d in dirnames if d.startswith("Lev")]
            )[-1]
            metadata = Metadata.load(dirpath / highest_lev / "metadata")

            metadata["files"] = {
                p2i(file.relative_to(dirpath)): {"link": str(file)}
                for file in files_to_upload(dirpath, annex_dir)
            }

            simulations[key] = metadata

            dirnames[:] = []  # Don't keep looking for common-metadata.txt files under this directory

    return simulations


def write_local_simulations(annex_dir):
    """Write the local simulations to a file for use when loading `Simulations`

    Args:
        annex_dir (str or Path): The path to the annex directory to be
        processed.

    Returns:
        None
    """
    from json import dump

    simulations = local_simulations(annex_dir)
    with open(sxs_directory("cache") / "local_simulations.json", "w") as f:
        dump(simulations, f, indent=2, separators=(",", ": "), ensure_ascii=True)
