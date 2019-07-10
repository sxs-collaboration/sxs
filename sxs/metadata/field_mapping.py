"""Dictionary mapping fields in the metadata to (short name, description) pairs

The keys of this dictionary should correspond to the keys present in the Metadata object (though not
all such keys must be here, and not all keys here must be present in a given Metadata object).  The
values should be two-tuples, the first element of which should be a short name that can fit in a
column header (and may include basic TeX commands); the second element of which should be a full
description of the item, to be displayed in tooltips, for example.

"""

metadata_field_mapping = {
    'name': ('Name', 'Primary short name by which this simulation is known'),
    'simulation_name': (r'Full name', 'Full simulation name'),
    'alternative_names': (r'Other names', 'Other names by which this run has been referenced'),
    'keywords': (r'Keywords', 'Keywords to qualitatively identify this simulation'),
    'point_of_contact_email': (r'Email', 'Point-of-contact email for this waveform.  Usually the person having placed the waveform into the repository'),
    'authors_emails': (r"Authors' emails", "Researchers who contributed to the generation of this waveform."),
    'simulation_bibtex_keys': (r'Sim bib keys', 'ADS BibTeX keys for this simulation'),
    'code_bibtex_keys': (r'Code bib keys', 'ADS BibTeX keys for the evolution code'),
    'initial_data_bibtex_keys': (r'ID bib keys', 'ADS BibTeX keys for the initial data'),
    'quasicircular_bibtex_keys': (r'QC bib keys', 'ADS BibTeX keys for eccentricity reduction'),
    'initial_data_type': (r'ID type', 'Initial-data type'),
    'initial_separation': (r'$\mathrm{sep}_{\mathrm{ini}}', 'Initial separation'),
    'initial_orbital_frequency': (r'$\Omega_{\mathrm{orb}\ \mathrm{ini}}$', 'Initial orbital frequency'),
    'initial_adot': (r'$\dot{a}_{\mathrm{ini}}$', 'initial rate of change of separation'),
    'object1': (r'Object 1', 'Type of object 1 (BH or NS)'),
    'object2': (r'Object 2', 'Type of object 2 (BH or NS)'),
    'initial_ADM_energy': (r'$E_{\mathrm{ADM}}$', 'ADM energy measured in the initial data'),
    'initial_ADM_linear_momentum': (r'$\vec{P}_{\mathrm{ADM}}$', 'ADM linear momentum measured in the initial data'),
    'initial_ADM_angular_momentum': (r'$\vec{L}_{\mathrm{ADM}}$', 'ADM angular momentum measured in the initial data'),
    'initial_mass1': (r'$M_{\mathrm{ini}, 1}$', 'Initial mass of object 1'),
    'initial_mass2': (r'$M_{\mathrm{ini}, 2}$', 'Initial mass of object 2'),
    'initial_dimensionless_spin1': (r'$\vec{chi}_{\mathrm{ini}, 1}$', 'Initial dimensionless spin of object 1'),
    'initial_dimensionless_spin2': (r'$\vec{chi}_{\mathrm{ini}, 2}$', 'Initial dimensionless spin of object 2'),
    'initial_position1': (r'$\vec{x}_{\mathrm{ini}, 1}$', 'Initial position of object 1'),
    'initial_position2': (r'$\vec{x}_{\mathrm{ini}, 2}$', 'Initial position of object 2'),
    'relaxation_time': (r'$t_{\mathrm{rel}}$', 'Time at which system has "relaxed" from initial transient stage'),
    'reference_time': (r'$t_{\mathrm{ref}}$', 'Time at which "reference" quantities are measured'),
    'reference_mass1': (r'$M_{\mathrm{ref}, 1}$', 'Mass of object 1 after reference time'),
    'reference_mass2': (r'$M_{\mathrm{ref}, 2}$', 'Mass of object 2 after reference time'),
    'reference_dimensionless_spin1': (r'$\vec{chi}_{\mathrm{ref}, 1}$', 'Dimensionless spin of object 1 after reference time'),
    'reference_dimensionless_spin2': (r'$\vec{chi}_{\mathrm{ref}, 2}$', 'Dimensionless spin of object 2 after reference time'),
    'reference_position1': (r'$\vec{x}_{\mathrm{ref}, 1}$', 'Position of object 1 after reference time'),
    'reference_position2': (r'$\vec{x}_{\mathrm{ref}, 2}$', 'Position of object 2 after reference time'),
    'reference_orbital_frequency': (r'$\Omega_{\mathrm{orb}\ \mathrm{ref}}$', 'Orbital frequency after reference time'),
    'reference_eccentricity': (r'$e_{\mathrm{ref}}$', 'Eccentricity after reference time'),
    'reference_mean_anomaly': (r'$M_{\mathrm{ref}}$', 'Mean anomaly after reference time'),
    'merger_time': (r'$t_{\mathrm{merg}}$', 'Merger time'),
    'number_of_orbits': (r'$N_{\mathrm{orb}}$', 'Number of orbits from beginning of simulation to merger'),
    'final_time': (r'$t_{\mathrm{fin}}$', 'Final time'),
    'remnant_mass': (r'$M_{\mathrm{rem}}$', 'Mass of remnant object'),
    'remnant_dimensionless_spin': (r'$\vec{\chi}_{\mathrm{rem}}$', 'Dimensionless spin of remnant object'),
    'remnant_velocity': (r'$v{\mathrm{rem}}$', 'Velocity of remnant object'),
    'disk_mass': (r'$M_{\mathrm{disk}}$', 'Mass of any remaining disk'),
    'ejecta_mass': (r'$M_{\mathrm{ejecta}}', 'ejecta mass'),
    'spec_revisions': (r'SpEC revisions', 'SpEC revisions'),
    'spells_revision': (r'SPELLS revision', 'SPELLS revision'),
}
