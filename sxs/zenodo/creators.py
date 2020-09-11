known_creators = {
    'Michael Boyle': {
        'name': 'Boyle, Michael',
        'affiliation': 'Cornell University',
        'orcid': '0000-0002-5075-5116',
    },
    'Mike Boyle': {
        'name': 'Boyle, Michael',
        'affiliation': 'Cornell University',
        'orcid': '0000-0002-5075-5116',
    },
    'Luisa Buchman': {
        'name': 'Buchman, Luisa',
        'affiliation': 'Caltech',
    },
    'Matthew Duez': {
        'name': 'Duez, Matthew',
        'affiliation': 'Washington State',
    },
    'Francois Foucart': {
        'name': 'Foucart, Francois',
        'affiliation': 'University of New Hampshire',
    },
    'Michael Grudich': {
        'name': 'Grudich, Michael',
        'affiliation': 'Caltech',
    },
    'Dan Hemberger': {
        'name': 'Hemberger, Dan',
    },
    'Larry Kidder': {
        'name': 'Kidder, Larry',
        'affiliation': 'Cornell University',
    },
    'Geoffrey Lovelace': {
        'name': 'Lovelace, Geoffrey',
        'affiliation': 'Cal State Fullerton',
        'orcid': '0000-0002-7084-1070',
    },
    'Ilana MacDonald': {
        'name': 'MacDonald, Ilana',
    },
    'Abdul Mroue': {
        'name': 'Mroue, Abdul',
    },
    'Harald Pfeiffer': {
        'name': 'Pfeiffer, Harald',
        'affiliation': 'AEI Potsdam',
        'orcid': '0000-0001-9288-519X',
    },
    'Mark Scheel': {
        'name': 'Scheel, Mark',
        'affiliation': 'Caltech',
    },
    'Bela Szilagyi': {
        'name': 'Szilagyi, Bela',
    },
    'Nicholas Taylor': {
        'name': 'Taylor, Nicholas',
    },
    'Anil Zenginoglu': {
        'name': 'Zenginoglu, Anil',
    }
}

creators_emails = {
    'mgiesler@tapir.caltech.edu': {'name': 'Matt Giesler'},
    # 'glovelace@fullerton.edu': { 'name': 'Geoffrey Lovelace' },  # in the default list, so don't worry about it
    'kevinbarkett@yahoo.com': {'name': 'Kevin Barkett'},
    'jonathan.blackman.0@gmail.com': {'name': 'Jonathan Blackman'},
    'ian.hinder@aei.mpg.de': {'name': 'Ian Hinder'},
    # 'pfeiffer@cita.utoronto.ca': { 'name': 'Harald Pfeiffer' },  # in the default list, so don't worry about it
    'tonyc@cita.utoronto.ca': {'name': 'Tony Chu'},
    'kbarkett@caltech.edu': {'name': 'Kevin Barkett'},
    'serguei.ossokine@aei.mpg.de': {'name': 'Serguei Ossokine'},
    'ossokine@aei.mpg.de': {'name': 'Serguei Ossokine'},
}

default_creators = [
    known_creators[person]
    for person in [
        'Larry Kidder', 'Harald Pfeiffer', 'Mark Scheel',
        'Michael Boyle', 'Dan Hemberger', 'Geoffrey Lovelace', 'Bela Szilagyi'
    ]
]
