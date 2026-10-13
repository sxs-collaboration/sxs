"""Documentation for metadata fields of individual simulations."""

import collections
import logging
import types
from math import inf
from textwrap import dedent
from enum import Flag, auto

# Possible groups for metadata fields
metadata_field_groups = {
    "id": "Identification",
    "refs": "References",
    "inputParam": "Input parameters for initial data",
    "IDmeasure": "Measurements of initial data",
    "referenceParam": "Reference quantities",
    "finalProperties": "Properties of merger/final quantities",
    "codeInfo": "Code information",
    "timeStamps": "Time stamps information",
}

# Defining class to distinguish between simulations.
class FieldPresence(Flag):
    BBH = auto()
    BHNS = auto()
    BNS = auto()
    SCATTERING = auto()
    MERGER = auto()
    ALL = BBH | BHNS | BNS | SCATTERING | MERGER  # All possible cases

FP = FieldPresence

# A field's description is just a dict; but this way we can provide
# some defaults, and later extend defaults, etc.
class MetadataField(collections.OrderedDict):
    """Information about a metadata field itself."""

    def __init__(self, *args, **kwargs):
        """Set initial values, and override with anything passed."""

        # Initialize as a dict
        super(MetadataField, self).__init__(*args, **kwargs)

        defaults = {
            "name":  "",           # Name of the field.
            "group": "",           # Which group this field is in.
            "intended_type": str,  # The intended type.
            "introduced": 2,       # The metadata version when this
                                   # field was introduced.
            "description": "",     # Description of the field, as markdown.
            "deprecated": inf,     # The metadata version when this
                                   # field was deprecated.  math.inf
                                   # means not deprecated.
            "computed": False,     # Whether field added by
                                   # sxs package
            "present_for": FP.ALL  # Whether the field is present
                                   # for vacuum mergers, matter mergers,
                                   # scattering or all.
        }

        # Update any missing values from defaults
        for k, v in defaults.items():
            self.setdefault(k, v)

# The metadata fields as a list
metadata_fields = [
    # group id
    MetadataField(
        group = "id",
        name  = "simulation_name",
        intended_type = str,
        introduced = 0,
        description = "A non-unique SXS-assigned identifier chosen "
        "before the simulation has been run. Useful only for SXS "
        "members building and debugging the catalog.",
    ),
    MetadataField(
        group = "id",
        name  = "alternative_names",
        intended_type = list[str],
        introduced = 0,
        description = "Comma-separated array of alternative names, "
        "longer, more descriptive, and/or indicating the specific "
        "series of simulations this configuration belongs to.  One "
        "of these alternative names is the `SXS:BBH:dddd` id-number, "
        "which is guaranteed to be unique."
    ),
    MetadataField(
        group = "id",
        name  = "keywords",
        intended_type = list[str],
        introduced = 0,
        description = "List of free-form keywords.  Presence of the "
        "keyword `deprecated` means that this simulation has been "
        "deprecated."
    ),
    MetadataField(
        group = "id",
        name  = "point_of_contact_email",
        intended_type = str,
        introduced = 0,
        description = "Contact information for questions.",
        deprecated = 2
    ),
    MetadataField(
        group = "id",
        name  = "job_archiver_email",
        intended_type = str,
        introduced = 2,
        description = "Email of person who archived this simulation into the "
        "catalog; usually the person who ran the simulation. Useful only for "
        "SXS members building and debugging the catalog."
    ),
    MetadataField(
        group = "id",
        name  = "authors_emails",
        intended_type = list[str],
        introduced = 0,
        deprecated = 2,
        description = "List of authors' emails."
    ),

    # group refs
    MetadataField(
        group = "refs",
        name  = "simulation_bibtex_keys",
        intended_type = list[str],
        introduced = 0,
        deprecated = 2,
        description = "References which should be cited if this "
        "simulation is used.",
    ),
    MetadataField(
        group = "refs",
        name  = "code_bibtex_keys",
        intended_type = list[str],
        introduced = 0,
        deprecated = 2,
        description = "List of bibtex keys which are references about "
        "the evolution code."
    ),
    MetadataField(
        group = "refs",
        name  = "initial_data_bibtex_keys",
        intended_type = list[str],
        introduced = 0,
        deprecated = 2,
        description = "List of bibtex keys which are references about "
        "the initial data code."
    ),
    MetadataField(
        group = "refs",
        name  = "quasicircular_bibtex_keys",
        intended_type = list[str],
        introduced = 0,
        deprecated = 2,
        description = ("List of bibtex keys which are references about "
                       "creating quasicircular initial data.")
    ),
    MetadataField(
        group = "refs",
        name  = "DOI_versions",
        intended_type = list[str],
        introduced = 2,
        computed = True,
        description = "Versions of simulation that are available via DOI.",
    ),
    MetadataField(
        group = "refs",
        name  = "citation_dois",
        intended_type = list[str],
        introduced = 2,
        description = "DOIs to cite when using this simulation."
    ),

    # group inputParam
    MetadataField(
        group = "inputParam",
        name  = "object_types",
        intended_type = str,
        introduced = 0,
        computed = True,
        description = "Keyword description to identify the types of both "
        "objects.  One of {'BHBH', 'BHNS', 'NSNS'}.",
    ),
    MetadataField(
        group = "inputParam",
        name  = "object1",
        intended_type = str,
        introduced = 0,
        description = "Keyword description to identify the object 1 "
        "type.  One of {'bh', 'ns'}."
    ),
    MetadataField(
        group = "inputParam",
        name  = "object2",
        intended_type = str,
        introduced = 0,
        description = "Keyword description to identify the object 2 "
        "type.  One of {'bh', 'ns'}."
    ),
    MetadataField(
        group = "inputParam",
        name  = "initial_data_type",
        intended_type = str,
        introduced = 0,
        description = dedent("""\
            Type of initial data.  One of

            - `BBH_CFMS` -- conformally flat, maximal slice;
            - `BBH_SKS` -- superposed Kerr-Schild;
            - `BBH_SHK` -- superposed harmonic Kerr-Schild [@Varma:2018sqd];
            - `BBH_SSphKS` -- superposed spherical Kerr-Schild [@Chen:2021rtb];
            - `BHNS`;
            - `NSNS`.""")
    ),
    MetadataField(
        group = "inputParam",
        name  = "initial_separation",
        intended_type = float,
        introduced = 0,
        description = r"Coordinate separation \(D_0\) between centers "
        "of compact objects, as passed to the initial data "
        "solver [@Cook:2004kt] [@Buonanno:2010yk] "
        "[@Ossokine:2015yla] (code units)."
    ),
    MetadataField(
        group = "inputParam",
        name  = "initial_orbital_frequency",
        intended_type = float,
        introduced = 0,
        description = r"Initial orbital frequency \(\Omega_0\) passed "
        "to the initial-data solver [@Buonanno:2010yk] "
        "[@Ossokine:2015yla] (code units)."
    ),
    MetadataField(
        group = "inputParam",
        name  = "initial_adot",
        intended_type = float,
        introduced = 0,
        description = r"Radial velocity parameter \(\dot{a}_0\) passed "
        "to the initial data solver [@Buonanno:2010yk] [@Pfeiffer:2007yz]."
    ),
    MetadataField(
        group = "inputParam",
        name  = "eos",
        intended_type = str,
        introduced = 1,
        present_for = FP.BNS | FP.BHNS | FP.MERGER,
        description = "Equation of state used for the NS during evolution."
    ),

    # group IDmeasure
    MetadataField(
        group = "IDmeasure",
        name  = "initial_ADM_energy",
        intended_type = float,
        introduced = 0,
        description = "ADM energy of the initial data (code units)."
    ),
    MetadataField(
        group = "IDmeasure",
        name  = "initial_ADM_linear_momentum",
        intended_type = list[float],
        introduced = 0,
        description = "ADM linear momentum of the initial data (code units)."
    ),
    MetadataField(
        group = "IDmeasure",
        name  = "initial_ADM_angular_momentum",
        intended_type = list[float],
        introduced = 0,
        description = "ADM angular momentum of the initial data (code units)."
    ),
    MetadataField(
        group = "IDmeasure",
        name  = "initial_mass1",
        intended_type = float,
        introduced = 0,
        description = "Christodoulou mass of apparent horizon 1 at "
        "initial data (code units)."
    ),
    MetadataField(
        group = "IDmeasure",
        name  = "initial_mass2",
        intended_type = float,
        introduced = 0,
        description = "Christodoulou mass of apparent horizon 2 at "
        "initial data (code units)."
    ),
    MetadataField(
        group = "IDmeasure",
        name  = "initial_mass_ratio",
        intended_type = float,
        introduced = 2,
        computed = True,
        description = "Mass ratio of the binary system at "
        "initial time.",
    ),
    MetadataField(
        group = "IDmeasure",
        name  = "initial_dimensionless_spin1",
        intended_type = list[float],
        introduced = 0,
        description = "Dimensionless spin of object 1 in the initial data."
    ),
    MetadataField(
        group = "IDmeasure",
        name  = "initial_dimensionless_spin2",
        intended_type = list[float],
        introduced = 0,
        description = "Dimensionless spin of object 2 in the initial data."
    ),
    MetadataField(
        group = "IDmeasure",
        name  = "initial_position1",
        intended_type = list[float],
        introduced = 0,
        description = "Initial coordinates of the center of body 1."
    ),
    MetadataField(
        group = "IDmeasure",
        name  = "initial_position2",
        intended_type = list[float],
        introduced = 0,
        description = "Initial coordinates of the center of body 2."
    ),
    MetadataField(
        group = "IDmeasure",
        name  = "initial_mass_withspin2",
        intended_type = float,
        introduced = 2,
        present_for = FP.BNS | FP.BHNS | FP.MERGER,
        description = "Mass calculated with rotational energy (i.e., using the "
        "Christodoulou formula, but with the measured NS spin).",
    ),

    # group referenceParam
    MetadataField(
        group = "referenceParam",
        name  = "relaxation_time",
        intended_type = float,
        introduced = 0,
        description = "Time at which we deem junk radiation to have "
        "sufficiently decayed (code units)."
    ),
    MetadataField(
        group = "referenceParam",
        name  = "reference_time",
        intended_type = float,
        introduced = 0,
        description = "Time at which reference quantities are extracted from "
        "the evolution (code units)."
    ),
    MetadataField(
        group = "referenceParam",
        name  = "reference_mass1",
        intended_type = float,
        introduced = 0,
        description = "Christodoulou mass of black hole 1 at "
        "reference time (code units)."
    ),
    MetadataField(
        group = "referenceParam",
        name  = "reference_mass2",
        intended_type = float,
        introduced = 0,
        description = "Christodoulou mass of black hole 2 at "
        "reference time (code units)."
    ),
    MetadataField(
        group = "referenceParam",
        name  = "reference_mass_ratio",
        intended_type = float,
        introduced = 2,
        description = "Mass ratio of the binary system at "
        "reference time.",
        computed = True
    ),
    MetadataField(
        group = "referenceParam",
        name  = "reference_dimensionless_spin1",
        intended_type = list[float],
        introduced = 0,
        description = "Dimensionless spin of object 1 at reference time."
    ),
    MetadataField(
        group = "referenceParam",
        name  = "reference_chi1_perp",
        intended_type = float,
        introduced = 2,
        description = "Magnitude of object 1's spin component in the "
        "orbital plane at reference time.",
        computed = True
    ),
    MetadataField(
        group = "referenceParam",
        name  = "reference_dimensionless_spin2",
        intended_type = list[float],
        introduced = 0,
        description = "Dimensionless spin of object 2 at reference time."
    ),
    MetadataField(
        group = "referenceParam",
        name  = "reference_chi2_perp",
        intended_type = float,
        introduced = 2,
        description = "Magnitude of object 2's spin component in the "
        "orbital plane at reference time.",
        computed = True
    ),
    MetadataField(
        group = "referenceParam",
        name  = "reference_chi_eff",
        intended_type = float,
        introduced = 2,
        description = "Effective spin of the binary system "
        "at reference time.",
        computed = True
    ),
    MetadataField(
        group = "referenceParam",
        name  = "reference_position1",
        intended_type = list[float],
        introduced = 0,
        description = "Coordinates of the center of body 1 at "
        "reference time."
    ),
    MetadataField(
        group = "referenceParam",
        name  = "reference_position2",
        intended_type = list[float],
        introduced = 0,
        description = "Coordinates of the center of body 2 at "
        "reference time."
    ),
    MetadataField(
        group = "referenceParam",
        name  = "reference_orbital_frequency",
        intended_type = list[float],
        introduced = 0,
        description = "Orbital angular frequency vector at reference time"
        " (code units)."
    ),
    MetadataField(
        group = "referenceParam",
        name  = "reference_mean_anomaly",
        intended_type = float,
        introduced = 0,
        description = "Mean anomaly at reference time."
    ),
    MetadataField(
        group = "referenceParam",
        name  = "reference_eccentricity",
        intended_type = float,
        introduced = 0,
        description = "Orbital eccentricity at reference time [@Mroue:2010re]."
    ),
    MetadataField(
        group = "finalProperties",
        name  = "number_of_orbits_from_reference_time",
        intended_type = float,
        introduced = 2,
        description = "Number of orbits from reference time until formation of "
        "a common apparent horizon."
    ),

    # group finalProperties
    MetadataField(
        group = "finalProperties",
        name  = "number_of_orbits",
        intended_type = float,
        introduced = 0,
        deprecated = 2,
        description = "Number of orbits until formation of a common apparent "
        "horizon.  Replaced by "
        "`number_of_orbits_from_reference_time` and "
        "`number_of_orbits_from_start`."
    ),
    MetadataField(
        group = "finalProperties",
        name  = "number_of_orbits_from_start",
        intended_type = float,
        introduced = 2,
        present_for = FP.BBH | FP.BHNS | FP.BNS | FP.MERGER,
        description = "Number of orbits from start of simulation until "
        "formation of a common apparent horizon."
    ),
    MetadataField(
        group = "finalProperties",
        name  = "common_horizon_time",
        intended_type = float,
        introduced = 0,
        present_for = FP.BBH | FP.MERGER,
        description = "Evolution time at which common horizon is first "
        "detected."
    ),
    MetadataField(
        group = "finalProperties",
        name  = "remnant_mass",
        intended_type = float,
        introduced = 0,
        present_for = FP.BBH | FP.MERGER,
        description = "Final Christodoulou mass of the remnant black hole after merger."
    ),
    MetadataField(
        group = "finalProperties",
        name  = "remnant_dimensionless_spin",
        intended_type = list[float],
        introduced = 0,
        present_for = FP.BBH | FP.MERGER,
        description = "Dimensionless spin of the remnant black hole after merger."
    ),
    MetadataField(
        group = "finalProperties",
        name  = "remnant_velocity",
        intended_type = list[float],
        introduced = 0,
        present_for = FP.BBH | FP.MERGER,
        description = "Linear velocity of the remnant black hole after merger."
    ),
    MetadataField(
        group = "finalProperties",
        name  = "disk_mass",
        intended_type = float,
        introduced = 1,
        present_for = FP.BNS | FP.BHNS | FP.MERGER,
        description = "Bound baryon mass on the disk at the end of the "
        "simulation (code units)."
    ),
    MetadataField(
        group = "finalProperties",
        name  = "ejecta_mass",
        intended_type = float,
        introduced = 1,
        present_for = FP.BNS | FP.BHNS | FP.MERGER,
        description = "Total amount of unbound mass generated by "
        "the end of the simulation (code units)."
    ),
    MetadataField(
        group = "finalProperties",
        name  = "end_of_trajectory_time",
        intended_type = float,
        introduced = 2,
        present_for = FP.BNS | FP.BHNS | FP.MERGER | FP.SCATTERING,
        description = "Time at end of trajectory before merger for bound "
        "systems or end of simulation for scattering systems."
    ),
    MetadataField(
        group = "finalProperties",
        name  = "final_time",
        intended_type = float,
        introduced = 1,
        description = "Time at end of simulation.",
        present_for = FP.BNS | FP.BHNS | FP.MERGER | FP.SCATTERING
    ),
    MetadataField(
        group = "finalProperties",
        name  = "merger_time",
        intended_type = float,
        introduced = 1,
        present_for = FP.BNS | FP.BHNS | FP.MERGER,
        description = "Time at which the density rises more than 3% above its "
        "original value."
    ),

    # group codeInfo
    MetadataField(
        group = "codeInfo",
        name  = "metadata_version",
        intended_type = int,
        introduced = 0,
        deprecated = 2,
        description = "This field has been replaced by the fields "
        "`metadata_format_revision` and "
        "`metadata_content_revision`.  The 2013 "
        "catalog [@Mroue:2013xna] implicitly carried metadata "
        "version 0. The 2019 catalog [@Boyle:2019kee] carried "
        "metadata version 1."
    ),
    MetadataField(
        group = "codeInfo",
        name  = "spec_revisions",
        intended_type = list[str],
        introduced = 0,
        description = "Array of git revisions of the evolution code."
    ),
    MetadataField(
        group = "codeInfo",
        name  = "spells_revision",
        intended_type = str,
        introduced = 0,
        description = "Git revision of initial data solver."
    ),
    MetadataField(
        group = "codeInfo",
        name  = "date_link_earliest",
        intended_type = str,
        introduced = 2,
        description = "Earliest link time of code used to perform this simulation."
    ),
    MetadataField(
        group = "codeInfo",
        name  = "internal_changelog",
        intended_type = dict,
        introduced = 2,
        description = "Text describing changes made in different "
        "`internal_minor_versions` of this local simulation. "
        "Always starts empty for new simulations."
    ),
    MetadataField(
        group = "codeInfo",
        name  = "internal_minor_version",
        intended_type = int,
        introduced = 2,
        description = "Incremented when anything changes in this local "
        "simulation that is not tracked by the fields "
        "`metadata_format_revision`, `metadata_content_revision`, or "
        "`postprocess_revision`. No relation to DOI revision numbers. "
        "Always starts at 0 for new simulations."
    ),
    MetadataField(
        group = "codeInfo",
        name  = "metadata_content_revision",
        intended_type = int,
        introduced = 2,
        description = "Incremented when values in the metadata "
        "change (which should seldom happen). "
        "Updated for all (non-deprecated) simulations at once. "
        "No relation to DOI revision numbers. "
        "At the time of this catalog release, all non-deprecated simulations "
        "carried `metadata_content_revision=1`."
    ),
    MetadataField(
        group = "codeInfo",
        name  = "metadata_format_revision",
        intended_type = int,
        introduced = 2,
        description = "Incremented when keys in the metadata "
        "change (which should seldom happen). "
        "Updated for all (non-deprecated) simulations at once. "
        "No relation to DOI revision "
        "numbers. At the time of this catalog release, all "
        "non-deprecated simulations carried "
        "`metadata_format_revision=2`."
    ),
    MetadataField(
        group = "codeInfo",
        name  = "postprocess_revision",
        intended_type = int,
        introduced = 2,
        description = "Incremented when anything changes that is not a "
        "raw SpEC output, such as re-computation of extrapolation, "
        "center-of-mass-correction, or memory-correction using newer "
        "algorithms or different parameters. Should not occur often. "
        "Updated for all (non-deprecated) simulations "
        "at once. No relation to DOI revision numbers. At the time of this "
        "catalog release, all non-deprecated simulations carried "
        "`postprocess_revision=1`."
    ),
    MetadataField(
        group = "codeInfo",
        name  = "t_relaxed_algorithm",
        intended_type = dict,
        introduced = 2,
        description = dedent("""\
            `t_relaxed_algorithm` is a dict. It contains fields:

            - `algorithm`: either `HHT` or `RMS`.
            - `reason`: present only for `RMS` method; text explaining why
              HHT method failed and we fell back to RMS.
            - `reference_time_method`: usually absent, but `set_by_hand` for
              certain simulations (almost head-on) where a reference time
              was explicitly set by hand and is not correlated with
              `relaxation_time`.""")
    ),
    MetadataField(
        group = "codeInfo",
        name  = "pbj_info",
        intended_type = dict,
        introduced = 2,
        description = dedent("""\
            This is a dict that contains { `base_lev` : `str`,
            `transition_time` : `float`, `base_lev_bitwise_identical` :
            `str` }

            - `base_lev` is the Lev that is shared between all PBandJ Levs.
              It is a string like `Lev3`.
            - `transition_time` is the time at which PBandJ happens. That
              is, before `transition_time`, all the Levs should be
              identical. If there is no PBandJ, then `transition_time` is
              0.0 and `base_lev` is the same as the Lev that the
              metadata.json file is in.
            - `base_lev_bitwise_identical` is either the string `true` or
              the string `false`. If `true`, this means that the current Lev
              (the one the metadata.json file is in) and its `base_lev`
              actually are bitwise identical up to (approx)
              `transition_time`. The `false` case occurs when someone ran an
              additional Lev at a later time, but `base_lev` had been
              deleted by the sysadmins, so the user reran `base_lev` and
              then did PBandJ to start the new lev. But the rerun of
              `base_lev` is not always bitwise identical to the original
              `base_lev` because something changed (timing-based stuff in
              SpEC, compiler version, libraries on the cluster, etc).""")
    ),

    # group timeStamps
    MetadataField(
        group = "timeStamps",
        name  = "date_postprocessing",
        intended_type = str,
        introduced = 2,
        description = "Timestamp of the most recent postprocessing of the raw "
        "simulation data to compute extrapolated, COM-corrected, "
        "memory-corrected waveforms."
    ),
    MetadataField(
        group = "timeStamps",
        name  = "date_run_earliest",
        intended_type = str,
        introduced = 2,
        description = "Timestamp of when this simulation was started."
    ),
    MetadataField(
        group = "timeStamps",
        name  = "date_run_latest",
        intended_type = str,
        introduced = 2,
        description = "Timestamp of when the last segment of this simulation "
        "started."
    ),

]

metadata_fields_dict = {f['name']: f for f in metadata_fields}
