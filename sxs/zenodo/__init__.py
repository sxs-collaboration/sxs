from .api import Login, Deposit, Records


def publish_sxs_bbh_simulation(sxs_bbh_directory_name, sandbox=False, deposition_id=None,
                               access_token=None, access_token_path=None, session=None,
                               creators=None, description=None, keywords=None,
                               access_right='open', license='cc-by'):
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

    # Get this deposit
    if deposition_id is not None:
        d = l.deposit(deposition_id)
        title = d.title
    else:
        # Check to see if this simulation exists in the list of the user's deposits or in the sxs community
        title = 'Binary black-hole simulation {0}'.format(sxs_bbh)
        matching_deposits = l.list_deposits(q='title: "{0}"'.format(title))
        if len(matching_deposits) == 1:
            deposition_id = matching_deposits[0]['id']
            print('A deposit with title "{0}"'.format(title))
            print('has already been started with this login.  Opening it for editing.')
            d = l.deposit(deposition_id)
        elif len(matching_deposits) > 1:
            print('Multiple deposits titled "{0}" have been found.'.format(title))
            raise ValueError(title)
        elif len(matching_deposits) == 0:
            # Check to see if this simulation is already in sxs but not owned by this login
            records = Records.search('title: "{0}"'.format(title))
            records = [r for r in records if r['title'] == title]  # Ensure *exact* match
            communities = [community['identifier']
                           for representation in records
                           for community in representation['metadata']['communities']]
            if 'sxs' in communities:
                print('A record exists on Zenodo with exactly the name "{0}",'.format(title))
                print('but is apparently owned by someone else in the "sxs" community.')
                print('Please contact the owner of that record if you wish to update it.')
                print('Web link: {0}'.format(records[0]['links']['html']))
                raise ValueError(title)
            elif len(communities) > 0:
                print('A record exists on Zenodo with exactly the name "{0}",'.format(title))
                print('but is apparently not owned by you, nor is it in the "sxs" community.')
                print('Web link: {0}'.format(records[0]['links']['html']))
                raise ValueError(title)
            d = l.new_deposit


    # Get the list of files we'll be uploading and compare to files already in the deposit to see if
    # any have changed.  If so, we need to create a new version.  Otherwise, we can just edit this
    # version.
    
    if len(new_files) == 0:
        print('Zenodo requires that there be at least one file')
        raise ValueError('No files found')


    # Convert each metadata.txt file to a metadata.json file sorted with interesting stuff at the
    # top of the file, so it appears on Zenodo's preview.

    ### OrderedDict(sorted(foo.iteritems(), key=lambda x: x[1]['depth']))


    # Get list of creators, keywords, and description

    if description is None:
        description = dedent("""\
        Simulation of a black-hole binary system evolved by the <a href="">SpEC code</a>.
        """)


    # Construct the Zenodo metadata
    new_metadata = {
        'title': title,
        'description': description,
        'keywords': keywords,
        'creator': creators,
        'upload_type': 'dataset',
        'access_right': access_right,
        'license': 'cc-by',
        'communities': [{'identifier': 'sxs'}],
    }
    metadata = d.metadata
    metadata.update(new_metadata)  # Ensure that fields we haven't changed are still present
    d.update_metadata(metadata)

    # Publish this version
    d.publish()



    

