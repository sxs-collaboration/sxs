from .api import Login, Deposition
from .creators import creators


def publish_sxs_bbh_simulation(sxs_bbh_directory_name, sandbox=False, deposition_id=None,
                               access_token=None, access_token_path=None, session=None):
    import re
    import os

    if not os.path.isdir(sxs_bbh_directory_name):
        print('The input directory name "{0}" does not appear to be a directory.'.format(sxs_bbh_directory_name))
        raise ValueError(sxs_bbh_directory_name)
    sxs_bbh = os.path.basename(sxs_bbh_directory_name)
    if not sxs_bbh.startswith('SXS:BBH:') or len(sxs_bbh) <= 8:
        print('Input directory "{0}" does not start with "SXS:BBH:" followed by some identifier.'.format(sxs_bbh))
        print('This function is only designed to handle SXS:BBH directories.')
        print('Since the input directory does not conform to this specification,')
        print('this function cannot construct the data with any robustness.')
        print('Feel free to copy this function\'s source to publish to Zenodo')
        print('as you wish, and submit a Pull Request if you feel that it would')
        print('make a useful addition to this module.')
        raise ValueError(sxs_bbh_directory_name)

    # Log in to zenodo
    l = Login(sandbox, access_token, access_token_path, session)

    # Get this deposition
    if deposition_id is not None:
        d = l.deposition(deposition_id)
        title = d.title
    else:
        # Check to see if this simulation exists in the list of the user's depositions or in the sxs community
        title = 'Binary black-hole simulation {0}'.format(sxs_bbh)
        matching_depositions = l.list_depositions(q='title: "{0}"'.format(title))
        if len(matching_depositions) == 1:
            deposition_id = matching_depositions[0]['id']
            print('A deposition with title "{0}"'.format(title))
            print('has already been started with this login.  Opening it for editing.')
            d = l.deposition(deposition_id)
        elif len(matching_depositions) > 1:
            print('Multiple depositions titled "{0}" have been found.'.format(title))
            raise ValueError(title)
        elif len(matching_depositions) == 0:
            # Check to see if this simulation is already in sxs
            r = l.session.get("https://zenodo.org/api/records/", params={'q': 'title: "{0}"'.format(title)})
            if r.status_code != 200:
                print('An unknown error occurred when trying to access https://zenodo.org/api/records/.')
                try:
                    print(r.json())
                except:
                    pass
                r.raise_for_status()
                raise RuntimeError()  # Will only happen if the response was not strictly an error
            communities = [community['identifier'] for representation in r.json()
                           for community in representation['metadata']['communities']]
            if 'sxs' in communities:
                print('')
                raise ValueError(title)
            d = l.new_deposition
        
    
    
    # Get the list of files we'll be uploading and compare to files already in the deposition


    # Convert each metadata.txt file to a metadata.json file


    # Get list of creators and keywords


    # Construct the Zenodo metadata
    new_metadata = {
        'title': title,
        'description': description,
        'keywords': keywords,
        'creator': creators,
        'upload_type': 'dataset',
        'access_right': 'open',
        'license': 'cc-by',
        'communities': [{'identifier': 'sxs'}],
    }
    metadata = d.metadata
    metadata.update(new_metadata)

    
        



    
    OrderedDict(sorted(foo.iteritems(), key=lambda x: x[1]['depth']))


