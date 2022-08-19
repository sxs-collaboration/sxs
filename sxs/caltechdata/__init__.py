from .api import Login, Deposit, Records

# from . import catalog, simannex, surrogatemodeling


import datetime
import json
import pathlib
from caltechdata_api import customize_schema
from .. import sxs_id, Metadata


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

    """

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
        for f in ["Strain_Outer.h5", "Strain_N2.h5", "Strain_N3.h5", "Strain_N4.h5"]:
            if not (lev / f).exists():
                raise ValueError(f"""Missing file "{f}" in "{lev}".""")
        if not (lev / "metadata.json").exists():
            mtxt = lev / "metadata.txt"
            if not mtxt.exists():
                raise ValueError(f"""Missing both "metadata.txt" and "metadata.json" in "{lev}".""")
            else:
                m = Metadata.from_txt_file(str(mtxt), cache_json=False).add_extras().reorder_keys()
                m.to_json_file(str(mtxt.with_suffix(".json")))
        if not (lev / "Horizons.h5").exists() and not (lev / "Matter.h5").exists():
            raise ValueError(f"""Missing both "Horizons.h5" and "Matter.h5" in "{lev}".""")

    # Figure out the "original name" of this simulation — something like "q4_7d/022" or "HighSpinMethods/BBH_SKS_d15.4_q1_..."
    root_dir = [par for par in directory.parents if par.match("Public") or par.match("Private")][0]
    original_name = str(directory.relative_to(root_dir))

    # Get the SXS ID from common-metadata.txt
    with (directory / "common-metadata.txt").open("r") as f:
        common_metadata = f.read()
    sxsid = sxs_id(common_metadata)
    if not sxsid:
        raise ValueError(f"""Could not find SXS ID in {directory / "common-metadata.txt"}""")

    # Find our list of files
    files = [f for f in directory.glob("Lev*/*") if f.name in allowed_files]
    highest_lev = sorted(set(f.relative_to(directory).parent for f in files))[-1]

    # Set up the metadata.  Note that CaltechDATA internally uses an undocumented
    # customized schema similar to — but different from — the standard Invenio
    # schema.  Fortunately, `caltechdata_api` includes a function to convert from
    # DataCite 4.0 or 4.3 to their schema.  So, we set ours up as DataCite 4.0 and
    # convert it: <https://schema.datacite.org/meta/kernel-4.0/>
    metadata = customize_schema({
        "resourceType": {
            "resourceTypeGeneral": "Dataset"
        },
        "titles": [
            {"title": f"Binary black-hole simulation {sxsid}"},
            {"title": f"{original_name}", "titleType": "AlternativeTitle"},
        ],
        "version": "10",
        "creators": [{"creatorName": "SXS Collaboration"}],
        "descriptions": [
            {
                "description": 'Simulation of a black-hole binary system evolved by the <a href="https://www.black-holes.org/code/SpEC.html">SpEC code</a>.',
                "descriptionType": "Abstract"
            }
        ],
        "subjects": [
            {"subject": "Astronomy"},
            {"subject": "Astrophysics"},
            {"subject": "Gravitational Waves"},
            {"subject": "Numerical Relativity"},
            {"subject": "Black Holes"},
            {"subject": "Neutron Stars"},
        ],
        "dates": [
            {
                "date": date or datetime.datetime.utcnow().strftime("%Y-%m-%d"),
                "dateType": "Created"
            },
            {
                "date": datetime.datetime.utcnow().strftime("%Y-%m-%d"),
                "dateType": "Updated"
            },
        ],
        "rightsList": [
            {
                "rights": "Creative Commons Attribution 4.0 International License",
                "rightsURI": "https://creativecommons.org/licenses/by/4.0/",
                "rightsIdentifier": "CC-BY-4.0",
            }
        ],
    }, schema="40")

    recordurl = f"{login.base_url}submit/api/create/"
    s3surl = f"{login.base_url}tindfiles/sign_s3/"
    chkurl = f"{login.base_url}tindfiles/md5_s3"

    # Search for an existing record with this name
    existing = login.search(f'title:"{sxsid}"')
    exists = existing.get("hits", {}).get("total", 0) > 0

    if exists:
        raise NotImplementedError()
        # find the latest version
        # update `metadata` with its content, and bump `publication_date` and `version`
        # check to see if any files have changed
        # if files have changed, edit the existing record and proceed
        # otherwise, exit early
    else:
        pass
        # create a new record and proceed
        # d = login.deposit(deposition_id=None)

        # url = f"{login.base_url}api/records"
        # r = login.session.post(url, data=json.dumps(metadata))
        # if r.status_code != 201:
        #     print(f"Unable to create a new deposit on {url}.")
        #     try:
        #         print(r.json())
        #     except:
        #         pass
        #     r.raise_for_status()
        #     raise RuntimeError()  # Will only happen if the response was not strictly an error
        # links = r.json()["links"]
        # url = "{0}api/records".format(login.base_url)
        # r = login.session.get(url, data=json.dumps(metadata))
        # assert r.status_code == 201, \
        #     f"Failed to create record (code: {r.status_code})"

    # d.update_metadata(metadata)
    # links = r.json()["links"]

    # Upload files
    metadata["files"] = [login.send_s3(f, f.relative_to(directory), verbose=True) for f in files]

    # Now tell CaltechDATA about it
    response = login.session._post(recordurl, data=json.dumps({"record": metadata}))

    return response.text
