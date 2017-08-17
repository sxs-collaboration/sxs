

def parse_metadata(filename):
    """Read a file of metadata and return a python dictionary object
    
    """
    if filename.endswith('.json'):
        import json
        with open(filename) as data_file:
            return json.load(data_file)
    else if filename.endswith('.txt'):
        return parse_metadata_txt(filename)
    else:
        raise ValueError("Don't understand ending of filename '{0}'; should be '.json' or '.txt'.".format(filename))

