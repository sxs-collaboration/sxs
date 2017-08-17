import os.path
import collections
import re
import json


_valid_identifier_pattern = re.compile('\W|^(?=\d)')
def _valid_identifier(key):
    return _valid_identifier_pattern.sub('_', key)

_metadata_key_map = {
    # This should be a dictionary of `valid_identifier: metadata_key` pairs for any pair that isn't
    # adequately covered by replacing underscores with dashes.  For example, this would be a working
    # but unnecessary pair:
    #   'simulation_name': 'simulation-name',
}
def _valid_identifier_to_metadata_key(key):
    return _metadata_key_map.get(key, key.replace('_', '-'))



class Metadata(collections.OrderedDict):
    """Object to collect metadata

    This object is essentially a `collections.OrderedDict`, with a few extra features:
      1) Keys are always forced to be valid python identifiers, whether setting or getting.
      2) There are a few extra methods for constructing these objects from json data or files or txt
         files, and outputting them to json or txt files.
      3) If an attribute does not exist on the object, the keys are searched.  This allows for
         tab-completion on the key names.
      4) Properties (methods without parentheses) that attempt to automatically determine the
         resolution and lev number from the 'simulation-name'.

    """

    @classmethod
    def from_file(cls, filename):
        """Read a file into a Metadata object

        The input file may either end in `.txt`, for an old-style metadata.txt file, or in `.json`,
        for a JSON file.  If neither ending is present, the function searches for files with the
        given `filename` plus both `.txt` and `.json`, and reads the newer one or at least the
        existing one.

        """
        if filename.endswith('.json'):
            return cls.from_json_file(filename)
        elif filename.endswith('.txt'):
            return cls.from_txt_file(filename)
        else:
            json_present = os.path.isfile(filename + '.json')
            txt_present = os.path.isfile(filename + '.txt')
            if json_present and not txt_present:
                return cls.from_json_file(filename + '.json')
            elif txt_present and not json_present:
                return cls.from_txt_file(filename + '.txt')
            elif json_present and txt_present:
                json_time = os.path.getmtime(filename + '.json')
                txt_time = os.path.getmtime(filename + '.txt')
                if txt_time > json_time:
                    return cls.from_txt_file(filename + '.txt')
                else:
                    return cls.from_json_file(filename + '.json')
            else:
                raise ValueError('Could not find file named "{0}", "{0}.txt", or "{0}.json"'.format(filename))

    @classmethod
    def from_json_data(cls, json_data):
        return json.load(json_data, object_pairs_hook=cls)

    @classmethod
    def from_json_file(cls, json_file):
        with open(json_file, 'r') as metadata_file:
            metadata = json.load(metadata_file, object_pairs_hook=cls)
        return metadata

    @classmethod
    def from_txt_file(cls, txt_file, cache_json=True):
        """Read metadata.txt file into Metadata object with valid python identifiers for keys

        A standard metadata.txt file is close to being an executable python script that just defines a bunch of
        constants.  The three problems with the metadata.txt format are:

          1) variable names contain dashes, which is the subtraction operator in python,
          2) strings are not enclosed in quotes, and
          3) lists are not enclosed in brackets

        It is easy to correct these problems.  In particular, (1) is resolved by changing dashes to underscores in the
        identifiers.  A bug in SpEC's metadata.txt files -- whereby some comment lines are missing the initial `#` -- is
        also fixed.

        Note that this function is not very flexible when it comes to generalizing the syntax of the metadata.txt files.
        In particular, it assumes that the right-hand sides are either numbers or strings (or lists of either numbers or
        strings).  For example, I think I've seen cases where the eccentricity is given as something like "<1e-5".  Since
        python has no "less-than" type, this is converted to a string.  But generally, this does seem to work on
        metadata.txt files in the SXS waveform repository.

        """
        from ast import literal_eval
        assignment_pattern = re.compile(r"""([-A-Za-z0-9]+)\s*=\s*(.*)""")
        string_pattern = re.compile(r"""[A-DF-Za-df-z<>@]""")  # Ignore 'e' and 'E' because they may appear in numbers
        metadata = cls()

        with open(txt_file, 'r') as metadata_file:
            for line in metadata_file:
                # Work around bug where some lines of dashes are missing the comment character
                if line.startswith('-'):
                    continue

                # Process each assignment line, skipping comments and unrecognized lines
                match = assignment_pattern.match(line)
                if match:
                    variable, quantity = match.groups()

                    # It was a stupid choice to make variables contain dashes
                    variable = _valid_identifier(variable)

                    # If `quantity` is an empty string, we should just replace it with an empty list
                    if not quantity or quantity == '\n':
                        quantity = []
                    else:
                        if string_pattern.search(quantity):
                            # If this is a string, strip whitespace from it, split lists and place
                            # brackets around them, and place quotation marks around each element
                            quantities = [q.strip() for q in quantity.split(",")]
                            if "," in quantity:
                                quantity = "['" + "', '".join(quantities) + "']"
                            else:
                                quantity = "'" + quantities[0] + "'"
                        else:
                            # Otherwise, just place brackets around lists of strings
                            if "," in quantity:
                                quantity = "[" + quantity + "]"

                    # Add this line to the metadata, whether or not it's been modified
                    metadata[variable] = literal_eval(quantity)

        if cache_json:
            # Skip the text processing next time, and just go straight to json
            txt_index = txt_file.rfind('.txt')
            if txt_index == -1:
                json_file = txt_file + '.json'
            else:
                json_file = txt_file[:txt_index] + '.json'
            metadata.to_json_file(json_file)

        return metadata

    def to_json(self, indent=4, separators=(',', ': ')):
        return json.dumps(self, indent=indent, separators=separators)

    def to_json_file(self, json_file, indent=4, separators=(',', ': ')):
        with open(json_file, 'w') as f:
            f.write(self.to_json(indent=indent, separators=separators))

    def to_txt(self):
        raise NotImplementedError()

    def to_txt_file(self, txt_file):
        with open(txt_file, 'w') as f:
            f.write(self.to_txt())
        
    def __init__(self, *args, **kwargs):
        """Initialize the OrderedDict, ensuring that all keys have been converted to valid identifiers

        This function intercepts the allowed args and kwargs and converts any keys before simply
        calling the base class's initialization function.

        """
        if len(args) > 0:
            args = list(args)
            if isinstance(args[0], collections.abc.Mapping):
                mapping = args[0]
                args[0] = OrderedDict([(_valid_identifier(key), mapping[key]) for key in mapping])
            else:
                iterable = args[0]
                args[0] = [(_valid_identifier(k), v) for k, v in iterable]
        if len(kwargs) > 0:
            for key in list(kwargs):
                kwargs[_valid_identifier(key)] = kwargs.pop(key)
        super(Metadata, self).__init__(*args, **kwargs)

    @classmethod
    def fromkeys(cls, iterable):
        iterable = [(_valid_identifier(k), v) for k, v in iterable]
        return super(Metadata, self).fromkeys(iterable)

    @property
    def resolution(self):
        """Try to determine the resolution from the 'simulation-name' field"""
        simulation_name = self['simulation_name']
        last_slash_index = simulation_name.rindex('/')
        return simulation_name[last_slash_index+1:]

    @property
    def lev(self):
        """Try to determine a numeric "Lev" number from the 'simulation-name' field"""
        resolution = self.resolution
        return int(resolution.replace('Lev', ''))

    def __contains__(self, key):
        return super(Metadata, self).__contains__(_valid_identifier(key))

    def __delattr__(self, name):
        super(Metadata, self).__delattr__(_valid_identifier(name))

    def __delitem__(self, key):
        super(Metadata, self).__delitem__(_valid_identifier(key))

    def __getattribute__(self, name):
        """Include keys as attributes

        This allows retrieval of a key like `md['simulation_name']` as `md.simulation_name`.

        """
        try:
            return super(Metadata, self).__getattribute__(name)
        except AttributeError as e:
            try:
                return self[_valid_identifier(name)]
            except KeyError:
                raise e

    def __dir__(self):
        """Ensure that the keys are included in tab completion"""
        return sorted(set(super(Metadata, self).__dir__()) | set(self.keys()))

    def __getitem__(self, key):
        return super(Metadata, self).__getitem__(_valid_identifier(key))

    def __setattr__(self, name, value):
        name = _valid_identifier(name)
        super(Metadata, self).__setattr__(name, value)

    def __setitem__(self, key, value):
        super(Metadata, self).__setitem__(_valid_identifier(key), value)

    def get(self, key, default=None):
        return super(Metadata, self).get(_valid_identifier(key), default)

    def pop(self, key, default=None):
        return super(Metadata, self).pop(_valid_identifier(key), default)

    def setdefault(self, key, default=None):
        return super(Metadata, self).setdefault(_valid_identifier(key), default)

    def update(self, mapping_or_iterable=None):
        if isinstance(mapping_or_iterable, collections.abc.Mapping):
            mapping_or_iterable = OrderedDict([(_valid_identifier(key), mapping_or_iterable[key]) for key in mapping_or_iterable])
        elif isinstance(mapping_or_iterable, collections.abc.Iterable):
            mapping_or_iterable = [(_valid_identifier(k), v) for k, v in mapping_or_iterable]
        return super(Metadata, self).update(mapping_or_iterable)


# class Metadata(object):

#     def __init__(self, metadata):
#         from collections import OrderedDict
#         if isinstance(metadata, dict):
#             # Just make sure that this is specifically an *ordered* dictionary
#             metadata = OrderedDict(metadata)
#         else if isinstance(metadata, basestr):
#             if metadata.endswith('.json'):
#                 with open(metadata) as file:
#                     metadata = json.load(metadata, object_pairs_hook=OrderedDict)
#             else if metadata.endswith('.txt'):
#                 pass
#             else:
#                 raise ValueError('Cannot understand metadata file "{0}" format from extension'.format(metadata))
#         else:
#             raise ValueError('Unknown input of type "{0}"'.format(type(metadata)))
        
#         self.__dict__.update(metadata)
#         self.all = metadata

#     def __len__(self):
#         return len(self.metadata)

#     def __getitem__(self, key):
#         """Ensure that the key has all dashes replaced by underscores"""
#         return self.metadata[key.replace('-', '_')]

#     def __setitem__(self, key, value):
#         """Ensure that the key has all dashes replaced by underscores"""
#         key = key.replace('-', '_')
#         self.metadata[key] = value
#         self.__dict__.update({key: value})

#     def __delitem__(self, key):
#         del self.metadata[key.replace('-', '_')]

#     def __iter__(self):
#         for key in self.metadata:
#             yield key

#     def __reversed__(self):
#         for key in reversed(self.metadata):
#             yield key

#     def items(self):
#         """Return iterator over (key, value) pairs"""
#         for key, value in self.metadata.items():
#             yield key, value

#     def __contains__(self, key):
#         return (key.replace('-', '_') in self.metadata)

#     def get(self, key, default=''):
#         return self.get(key.replace('-', '_'), default)

#     def __repr__(self):
#         """Return an unambiguous string representation"""
#         return 'Metadata({0})'.format(repr(self.metadata))
#         # return '{' + ',\n '.join(["'{0}': {1}".format(key, val) for key, val in self.__dict__.items()]) + ',}'

#     def __str__(self):
#         """Return a compact (but hopefully descriptive and unique) string representation"""
#         string = ("Simulation name: {0}\n"
#                   + "Alternative names: {1}\n"
#                   + "Masses: {2}, {3}\n"
#                   + "Spins: {4},\n       {5}\n")
#         return string.format(self.get('simulation_name'), self.get('alternative_names'),
#                              self.get('relaxed_mass1'), self.get('relaxed_mass2'),
#                              self.get('relaxed_spin1'), self.get('relaxed_spin2'))

#     def _repr_html_(self):
#         """Return a string representation of this object to appear in jupyter notebooks"""
#         string = ("Simulation name: {0}<br/>\n"
#                   + "Alternative names: {1}<br/>\n"
#                   + "Masses: {2}, {3}<br/>\n"
#                   + "Spins:<br/>\n&nbsp;&nbsp;&nbsp;&nbsp;{4},<br/>\n&nbsp;&nbsp;&nbsp;&nbsp;{5}\n")
#         return string.format(self.get('simulation_name'), self.get('alternative_names'),
#                              self.get('relaxed_mass1'), self.get('relaxed_mass2'),
#                              self.get('relaxed_spin1'), self.get('relaxed_spin2'))
