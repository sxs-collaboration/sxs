import copy
import time
import warnings
import re
import datetime
import json
import pathlib

from .login import Login
from .. import sxs_id, Metadata
from ..utilities import sxs_identifier_regex, SimpleVersion


# To do:
# - Finish up the code to create a new version from an existing one


def mtime(f):
    """Look for git or filesystem modification time

    If the git call fails for any reason, we fall back to the filesystem's `mtime`.

    The input `f` is expected to be a path to the file.  The return value is a DateTime object

    """
    from subprocess import check_output
    try:
        timestamp = check_output(f"""git log -1 --format="%ad" --date=unix -- {f}""", shell=True).decode().strip()
    except:
        timestamp = f.stat().st_mtime
    finally:
        mt = datetime.datetime.fromtimestamp(timestamp, tz=datetime.timezone.utc)
    return mt


def filter_files(files, utime):
    return [f for f in files if mtime(f) > utime]


def upload_simannex_dir(login, directory, date=None):
    """Update/upload data from the SimAnnex to CaltechDATA

    Parameters
    ----------
    login : caltechdata.Login
        This is used to communicate with CaltechDATA.
    directory : str or pathlib.Path
        This must contain either "Public" or "Private", followed by the specific
        directory in SimAnnex to upload.  It can also contain higher path elements
        — for example, you can give the absolute path.
    date : str, optional
        If given, this is used as the `publicationDate` field.  It must be
        formatted as "%Y-%m-%d".  The default is just the current UTC time.

    """
    raise Exception("This interface is for the old TIND API; use zenodo for the new API")
    from caltechdata_api import customize_schema, decustomize_schema

    allowed_files = [
        "metadata.json",
        "Strain_Outer.h5",
        "Strain_N2.h5",
        "Strain_N3.h5",
        "Strain_N4.h5",
        "Horizons.h5",
        "Matter.h5",
    ]

    # Check that `directory` is a valid SXS dataset
    # - Name must contain either "Public" or "Private" followed by at least one subdirectory
    # - Must contain common-metadata.txt (though this will not be uploaded)
    # - Must contain at least one Lev* subdirectory
    # - Each Lev* subdirectory must contain these files (which will be uploaded)
    #   - metadata.json
    #   - Strain_Outer.h5
    #   - Strain_N2.h5
    #   - Strain_N3.h5
    #   - Strain_N4.h5
    #   - Horizons.h5 and/or Matter.h5
    directory = pathlib.Path(directory).expanduser().resolve()
    if "/Public/" not in str(directory) and "/Private/" not in str(directory):
        raise ValueError(f"""Either "Public" or "Private" must be in the input directory "{directory}".""")
    if not (directory / "common-metadata.txt").exists():
        raise ValueError(f"Missing common-metadata.txt in {directory}")
    if not directory.glob("Lev*"):
        raise ValueError(f"Missing Lev* in {directory}")
    for lev in directory.glob("Lev*"):
        warnings.warn("Temporarily skipping file checks")
        # for f in ["Strain_Outer.h5", "Strain_N2.h5", "Strain_N3.h5", "Strain_N4.h5"]:
        #     if not (lev / f).exists():
        #         raise ValueError(f"""Missing file "{f}" in "{lev}".""")
        # if not (lev / "metadata.json").exists():
        #     mtxt = lev / "metadata.txt"
        #     if not mtxt.exists():
        #         raise ValueError(f"""Missing both "metadata.txt" and "metadata.json" in "{lev}".""")
        #     else:
        #         m = Metadata.from_txt_file(str(mtxt), cache_json=False).add_extras().reorder_keys()
        #         m.to_json_file(str(mtxt.with_suffix(".json")))
        # if not (lev / "Horizons.h5").exists() and not (lev / "Matter.h5").exists():
        #     raise ValueError(f"""Missing both "Horizons.h5" and "Matter.h5" in "{lev}".""")

    # Figure out the "original name" of this simulation — something like "q4_7d/022" or "HighSpinMethods/BBH_SKS_d15.4_q1_..."
    root_dir = [par for par in directory.parents if par.match("Public") or par.match("Private")][0]
    original_name = str(directory.relative_to(root_dir))

    # Get the SXS ID from common-metadata.txt
    sxs_system_re = re.compile(sxs_identifier_regex)
    sxs_system_type = None
    sxs_system_number = None
    with (directory / "common-metadata.txt").open("r") as f:
        for line in f.readlines():
            line = line.strip()
            if "alternative-names" in line:
                m = sxs_system_re.search(line)
                if m:
                    sxs_system_type = m["simulation_type"]
                    sxs_system_number = m["sxs_number"]
                    break
    if not sxs_system_type or not sxs_system_number:
        raise ValueError(f"No SXS identifier found in {directory / 'common-metadata.txt'}")
    sxs_system = f"SXS:{sxs_system_type}:{sxs_system_number}"
    spec_url = "https://www.black-holes.org/code/SpEC.html"
    if sxs_system_type == "BBH":
        title = f"Binary black-hole simulation {sxs_system}"
        description = f"""Simulation of a black-hole binary system evolved by the <a href="{spec_url}">SpEC code</a>."""
    elif sxs_system_type == "BHNS":
        title = f"Black-hole neutron-star binary simulation {sxs_system}"
        description = f"""Simulation of a black-hole neutron-star binary system evolved by the <a href="{spec_url}">SpEC code</a>."""
    elif sxs_system_type == "NSNS":
        title = f"Binary neutron-star simulation {sxs_system}"
        description = f"""Simulation of a neutron-star binary system evolved by the <a href="{spec_url}">SpEC code</a>."""
    else:
        raise ValueError(
            f"""Did not recognize SXS system type "{sxs_system_type}" in directory "{directory}"; should be BBH, BHNS, or NSNS."""
        )
    print(f"Beginning work on {sxs_system}")

    # Find our list of files
    files = sorted(
        [f for f in directory.glob("Lev*/*") if f.name in allowed_files],
        key=lambda f: f"{f.parent}/{allowed_files.index(f.name)}"
    )

    # Search for an existing record with this name
    existing = login.search(f'"{sxs_system}"')
    exists = existing.get("hits", {}).get("total", 0) > 0

    if exists:
        # Find the latest version
        latest = max(
            existing["hits"]["hits"],
            key=lambda h: SimpleVersion(h["metadata"]["version"])
        )
        latest_link = latest["links"]["self"].replace("/api/", "/")
        latest_doi = latest["metadata"]["doi"]

        # Check to see if any files have changed since the latest version
        utime = datetime.datetime.fromisoformat(latest["updated"])
        files_to_upload = filter_files(files, utime)
        if not files:
            print(f"No changes to {latest_link}.")
            return latest_link

        # Update version, relatedIdentifiers, and updated date
        metadata = latest["metadata"]
        old_version = SimpleVersion(metadata["version"])
        metadata["version"] = str(old_version.increment())
        metadata["relatedIdentifiers"] =[
            {
                "relatedIdentifierRelation": "IsNewVersionOf",
                "relatedIdentifier": latest_doi,
                "relatedIdentifierScheme": "DOI"
            }
        ]
        metadata["relevantDates"] = [
            {
                "relevantDateType": "Updated",
                "relevantDateValue": datetime.datetime.utcnow().strftime("%Y-%m-%d"),
            }
        ]

        # Upload files
        raise NotImplementedError(
            "Adding new files to an existing record is not yet supported.  We will need some logic like this:\n" +
            "  https://github.com/caltechlibrary/caltechdata_api/blob/840a359a017d1257ac344b6af80038f591fd8a97/caltechdata_api/caltechdata_edit.py#L189-L233\n" +
            "More specifically, if we want to create a new version, we will probably need to essentially copy\n" +
            "the existing record, publish it, and then follow steps like the above to upload new files and\n" +
            "delete the old ones.  It is not clear to me if the deletion would destroy files being used by the\n" +
            "older version, or just delete that entry in the new record."
        )
        existing_names = {
            f["electronic_name"][0]: f["uniform_resource_identifier"].split("/")[-2]
            for f in metadata["electronic_location_and_access"]
        }
        for f in files_to_upload:
            fileinfo = login.send_s3(f, str(f.relative_to(directory)), verbose=True)
            name = fileinfo["filename"]
            if name in existing_names:
                metadata["files"][existing_names.index(name)] = fileinfo
            else:
                metadata["files"].append(fileinfo)

    else:
        # Set up the metadata.  Note that CaltechDATA internally uses an undocumented
        # customized schema similar to — but different from — the standard Invenio
        # schema.  Fortunately, `caltechdata_api` includes a function to convert from
        # DataCite 4.0 or 4.3 to their schema.  So, we set ours up as DataCite 4.0 and
        # convert it: <https://schema.datacite.org/meta/kernel-4.0/>

        # Generate a little table of info to append to the description
        highest_lev = sorted(set(f.relative_to(directory).parent for f in files))[-1]
        sim_metadata = Metadata.from_file(directory / highest_lev / "metadata")
        object_types = sim_metadata.get("object_types", "")
        mass_ratio = sim_metadata.get("reference_mass_ratio", "Unknown")
        if "reference_dimensionless_spin1" in sim_metadata:
            chi1 = sim_metadata.reference_dimensionless_spin1
            chi1 = f"[{chi1[0]:0.4f}, {chi1[1]:0.4f}, {chi1[2]:0.4f}]"
        else:
            chi1 = "N/A"
        if "reference_dimensionless_spin2" in sim_metadata:
            chi2 = sim_metadata.reference_dimensionless_spin2
            chi2 = f"[{chi2[0]:0.4f}, {chi2[1]:0.4f}, {chi2[2]:0.4f}]"
        else:
            chi2 = "N/A"
        n_orbits = sim_metadata.get("number_of_orbits", "Unknown")
        reference_eccentricity = sim_metadata.get("reference_eccentricity", "Unknown")
        description = f"""{description}
        <br/>
        <table>
          <tbody>
            <tr><td style="padding:0 15px;">Mass ratio</td><td>{mass_ratio:.4f}</td></tr>
            <tr><td style="padding:0 15px;">Spin 1</td><td>{chi1}</td></tr>
            <tr><td style="padding:0 15px;">Spin 2</td><td>{chi2}</td></tr>
            <tr><td style="padding:0 15px;">Number of orbits</td><td>{n_orbits:.4f}</td></tr>
            <tr><td style="padding:0 15px;">Eccentricity</td><td>{reference_eccentricity}</td></tr>
          </tbody>
        </table>
        <br/>
        Originally named "{original_name}"
        """
        files_to_upload = files

        metadata = {
            "resourceType": {
                "resourceTypeGeneral": "Dataset"
            },
            "titles": [
                {"title": title},
            ],
            "version": "2.0",
            "creators": [{"creatorName": "SXS Collaboration"}],
            "descriptions": [
                {
                    "descriptionType": "Abstract",
                    "description": description,
                }
            ],
            "subjects": [
                {"subject": "Astronomy"},
                {"subject": "Astrophysics"},
                {"subject": "Gravitational Waves"},
                {"subject": "Numerical Relativity"},
            ],
            "dates": [
                {
                    "date": datetime.datetime.utcnow().strftime("%Y-%m-%d"),
                    "dateType": "Updated"
                },
            ],
            "rightsList": [
                {
                    "rights": "Creative Commons Attribution 4.0 International License",
                    "rightsURI": "https://creativecommons.org/licenses/by/4.0/",
                }
            ],
        }
        if "bh" in object_types.lower():
            metadata["subjects"].append({"subject": "Black Holes"})
        if "ns" in object_types.lower():
            metadata["subjects"].append({"subject": "Neutron Stars"})
        metadata = customize_schema(metadata, schema="40")
        metadata["publicationDate"] = date or datetime.datetime.utcnow().strftime("%Y-%m-%d")
        metadata["files"] = [login.send_s3(f, str(f.relative_to(directory)), verbose=True) for f in files]

    metadata["titles"] = [
        {"title": title},
        {"title": f"{original_name}", "titleType": "AlternativeTitle"},
    ]

    # if doi is None:
    #     # We want tind to generate the identifier
    #     metadata["final_actions"] = [
    #         {
    #             "type": "create_doi",
    #             "parameters": {"type": "records", "field": "doi"},
    #         }
    #     ]
    # else:
    metadata["doi"] = f"""{login.doi_prefix}/{sxs_system}v{metadata["version"]}"""

    # Now tell CaltechDATA about it
    recordurl = f"{login.base_url}submit/api/create/"
    response = login.session.post(recordurl, data=json.dumps({"record": metadata}))
    if response.status_code != 200:
        print(f"An error occurred when trying to create a new record for '{directory}'.")
        try:
            print(response.text)
        except:
            pass
        try:
            print(response.json())
        except:
            pass
        response.raise_for_status()
        raise RuntimeError()  # Will only happen if the response was not strictly an error
    else:
        print(f"""  Publishing as "{title}".""")
        print(f"""  {response.text}""")

    url = response.text[response.text.find(login.base_url):].rstrip().rstrip(".")
    pid = url.rsplit("/", 1)[1]

    # Create the DOI
    print("  Creating/updating DOIs")
    api_url = f"{login.base_url}api/record/{pid}"
    for retry in range(20):
        r = login.session.get(api_url)
        if r.status_code == 200:
            break
        time.sleep(1.1)  # Let things settle down at CaltechDATA, so we can get the metadata
    if r.status_code != 200:
        print(f"""An error occurred when trying to access "{api_url}".""")
        try:
            print(r.json())
        except:
            pass
        r.raise_for_status()
        raise RuntimeError()  # Will only happen if the response was not strictly an error
    output_metadata = r.json()
    doi_metadata = decustomize_schema(output_metadata, schema="43")
    if not exists:
        doi = f"""{login.doi_prefix}/{sxs_system}"""
        login.datacite.public_doi(doi_metadata, url, doi)
    doi = f"""{login.doi_prefix}/{sxs_system}v{metadata["version"]}"""
    login.datacite.public_doi(doi_metadata, url, doi)

    return url
