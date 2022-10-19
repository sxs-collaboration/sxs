format_description = """New catalog format

{
    "sxs_format": "catalog_v2"
    "modified": "<YYYY-MM-DDThh:mm:ss.ssssss>",  # UTC time of last-modified record in this file

    "records": {  # Includes *all* versions
        "<sxs_id_versioned>": {  # SXS:(BBH|BHNS|NSNS):[0-9]{4,}v[0-9]{2,}
            "title": "<title>"  # Same as ["metadata"]["title"],
            "id": <id>,  # ~7-digit integer uniquely identifying this record
            "conceptrecid": "<conceptrecid>",  # ~7-digit integer (as string) collectively identifying all versions of this record
            "created": "<YYYY-MM-DDThh:mm:ss.ssssss>",  # UTC time of creation of this record on Zenodo
            "modified": "<YYYY-MM-DDThh:mm:ss.ssssss>",  # (UTC) Last modification of this record (possibly just Zenodo metadata modified)
            "bucket": "https://zenodo.org/api/files/<uuid>",  # Base URL for file uploads and downloads,
            "files": [
                {
                    "checksum": "<checksum>",  # MD5 checksum of file on Zenodo
                    "filename": "<filename>",  # Name of file; may contain slashes denoting directories
                    "filesize": <filesize>,  # Number of bytes in the file
                    "id": "<fileid>",  # A standard UUID (hexadecimal with characters in the pattern 8-4-4-4-12)
                },
                ...  # Other file descriptions in the order in which they were uploaded (not necessarily a meaningful order)
            ]
        },
        ...
    },

    "simulations": {  # Physical data (masses, spins, etc.) for all available SXS simulations, in the most recent version
        "<sxs_id>": {  # The SXS ID is a string like SXS:BHNS:0001 or SXS:BBH:1234
            "latest_version": "<sxs_id_versioned>",  # Entry in "records" containing the most recent version of this simulation's data
            "url": "<URL>",  # The URL of the Zenodo "conceptdoi" link, which *resolves to* the most-recent version
            #
            # NOTE: All of the following may be absent if this simulation is closed-access, or simply does not have metadata.
            #
            # Variable content describing (mostly) physical parameters of the system.  It's basically a
            # python-compatible version of the information contained in "metadata.txt" from the
            # highest-resolution run in the most-recent version of this simulation.  That file is meant to
            # be more-or-less as suggested in <https://arxiv.org/abs/0709.0093>.  The conversion to a
            # python-compatible format means that keys like "simulation-name" have had hyphens replaced by
            # underscores so that they can be used as variable names in python and any other sane language
            # (with apologies to Lisp).  As far as possible, values that are just strings in that file
            # have been converted into the relevant types -- like numbers, integers, and arrays.  Note
            # that some keys like eccentricity are sometimes numbers and sometimes the string "<number"
            # (meaning that the eccentricity is less than the number), which is necessarily a string.
            #
            # Below are just the first few keys that *may* be present.  Note that closed-access
            # simulations will have empty dictionaries here.
            #
            "simulation_name": "<directory_name>",  # This may be distinctly uninformative
            "alternative_names": "<some ugly thing>, ..., <sxs_id>",  # This may be a list of strings
            "initial_data_type": "<type>",  # Something like "BBH_CFMS"
            "object_types": "<type>",  # Currently "BHBH", "BHNS", or "NSNS"
            "number_of_orbits": <number>,  # This is a float, rather than an integer
            "reference_mass_ratio": <q>,  # Usually greater than 1 (exceptions are due to junk radiation)
            "reference_chi_eff": <chi_eff>,  # Dimensionless effective spin quantity
            "reference_chi1_perp": <chi1_perp>,  # Magnitude of component of chi1 orthogonal to "reference_orbital_frequency"
            "reference_chi2_perp": <chi2_perp>,  # Magnitude of component of chi2 orthogonal to "reference_orbital_frequency"
            "reference_mass1": <m2>,
            "reference_mass2": <m1>,
            "reference_dimensionless_spin1": [
                <chi1_x>,
                <chi1_y>,
                <chi1_z>
            ],
            "reference_dimensionless_spin2": [
                <chi2_x>,
                <chi2_y>,
                <chi2_z>
            ],
            "reference_eccentricity": <eccentricity>,  # A float or possibly a string containing "<" and a float
            "reference_orbital_frequency": [
                <omega_x>,
                <omega_y>,
                <omega_z>
            ],
            "reference_time": <time>,
            ...
        },
        ...
    }

}


"""
