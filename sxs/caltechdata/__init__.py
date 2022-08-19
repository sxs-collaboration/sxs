from .api import Login, Deposit, Records

# from . import catalog, simannex, surrogatemodeling


import datetime
import json
import pathlib
from .. import sxs_id, Metadata


"""Update/upload data


"""


def upload(directory, login=None):
    if login is None:
        login = Login()

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
    files = [f.relative_to(directory) for f in directory.glob("Lev*/*") if f.name in allowed_files]
    highest_lev = sorted(set(f.parent for f in files))[-1]

    # Set up the metadata; see https://inveniordm.docs.cern.ch/reference/metadata for details
    resource_type = {
        "id": "dataset"
    }  # https://github.com/caltechlibrary/caltechdata/blob/main/app_data/vocabularies/resource_types.yaml
    creators = [
        {
            "person_or_org": {
                "type": "organizational",
                "name": "SXS Collaboration",
            },
        }
    ]
    title = f"Binary black-hole simulation {sxsid}"
    publication_date = datetime.datetime.utcnow().strftime("%Y-%M-%d")
    additional_titles = [
        {
            "title": f"{original_name}",
            "type": {
                "id": "alternative-title",
                "title": {"en": "Original simulation name"},
            },
        }
    ]
    description = 'Simulation of a black-hole binary system evolved by the <a href="https://www.black-holes.org/code/SpEC.html">SpEC code</a>.'
    rights = [{"id": "cc-by-4.0"}]
    subjects = [
        {"id": "http://www.oecd.org/science/inno/38235147.pdf?1.3"},
        {"subject": "astronomy"},
        {"subject": "astrophysics"},
        {"subject": "gravitational waves"},
        {"subject": "numerical relativity"},
        {"subject": "black holes"},
    ]
    version = "v10"
    default_preview_file = str(highest_lev / "metadata.json")

    metadata = {
        "access": {"files": "public", "record": "public"},
        "metadata": {
            "resource_type": resource_type,
            "creators": creators,
            "title": title,
            "publication_date": publication_date,
            "additional_titles": additional_titles,
            "description": description,
            "rights": rights,
            "subjects": subjects,
            "version": version,
        },
        "files": {
            "enabled": True,
            "default_preview": default_preview_file,
        },
    }

    metadata = {
        "access": {"record": "public", "files": "public"},
        "files": {"enabled": True},
        "metadata": {
            "creators": [
                {"person_or_org": {"family_name": "Brown", "given_name": "Troy", "type": "personal"}},
                {
                    "person_or_org": {
                        "family_name": "Collins",
                        "given_name": "Thomas",
                        "identifiers": [{"scheme": "orcid", "identifier": "0000-0002-1825-0097"}],
                        "name": "Collins, Thomas",
                        "type": "personal",
                    },
                    "affiliations": [{"id": "01ggx4157", "name": "Entity One"}],
                },
            ],
            "publication_date": "2020-06-01",
            "resource_type": {"id": "image-photo"},
            "title": "A Romans story",
        },
    }

    # Search for an existing record with this name
    existing = login.search(f'title:"{sxsid}"')

    if existing.get("hits", {}).get("total", 0) > 0:
        raise NotImplementedError()
        # find the latest version
        # update `metadata` with its content, and bump `publication_date` and `version`
        # check to see if any files have changed
        # if files have changed, edit the existing record and proceed
        # otherwise, exit early
    else:
        # create a new record and proceed
        # d = login.deposit(deposition_id=None)

        url = f"{login.base_url}api/records"
        r = login.session.post(url, data=json.dumps(metadata), headers={"Content-Type": "application/json"})
        if r.status_code != 201:
            print(f"Unable to create a new deposit on {url}.")
            try:
                print(r.json())
            except:
                pass
            r.raise_for_status()
            raise RuntimeError()  # Will only happen if the response was not strictly an error
        links = r.json()["links"]
        # url = "{0}api/records".format(login.base_url)
        # r = login.session.get(url, data=json.dumps(metadata))
        # assert r.status_code == 201, \
        #     f"Failed to create record (code: {r.status_code})"

    # d.update_metadata(metadata)
    # links = r.json()["links"]

    # Upload files
    for f in files:
        name = f.relative_to(directory)

        # Initiate the file
        r = login.session.post(links["files"], data=json.dumps([{"key": name}]))
        assert r.status_code == 201, f"Failed to create file {f} (code: {r.status_code})"
        file_links = r.json()["entries"][0]["links"]

        # Upload file content by streaming the data
        with open(f, "rb") as fp:
            r = login.session.put(file_links["content"], data=fp, headers={"Content-Type": "application/octet-stream"})
        assert r.status_code == 200, f"Failed to upload file contet {f} (code: {r.status_code})"

        # Commit the file.
        r = login.session.post(file_links["commit"])
        assert r.status_code == 200, f"Failed to commit file {f} (code: {r.status_code})"

    # Publish
    d.publish()

    return d
