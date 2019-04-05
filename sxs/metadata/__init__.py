from __future__ import absolute_import, division, print_function
import os.path
import collections
import warnings
import re
import json

from .catalog import (read_catalog, drop_all_but_highest_levs, drop_all_but_selected_resolutions,
                      key_by_alternative_name, symlink_runs)
from .field_mapping import metadata_field_mapping
from .fields import metadata_fields
from .web import create_web_files
from .initial_data import uses_new_initial_data, get_initial_data_details


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


def _mkdir_recursively(path):
    import errno    
    import os
    import os.path
    try:
        os.makedirs(os.path.abspath(path))
    except OSError as e:
        if e.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


class Metadata(collections.OrderedDict):
    """Object to collect metadata

    Note that the constructor is not generally useful from outside this class.  See the
    Metadata.from_file, etc., functions for more useful factory functions.

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
    def from_file(cls, filename, ignore_invalid_lines=False):
        """Read a file into a Metadata object

        The input file may either end in `.txt`, for an old-style metadata.txt file, or in `.json`,
        for a JSON file.  If neither ending is present, the function searches for files with the
        given `filename` plus both `.txt` and `.json`, and reads the newer one or at least the
        existing one.

        """
        if filename.endswith('.json'):
            return cls.from_json_file(filename)
        elif filename.endswith('.txt'):
            return cls.from_txt_file(filename, ignore_invalid_lines=ignore_invalid_lines)
        else:
            json_present = os.path.isfile(filename + '.json')
            txt_present = os.path.isfile(filename + '.txt')
            if json_present and not txt_present:
                return cls.from_json_file(filename + '.json')
            elif txt_present and not json_present:
                return cls.from_txt_file(filename + '.txt', ignore_invalid_lines=ignore_invalid_lines)
            elif json_present and txt_present:
                json_time = os.path.getmtime(filename + '.json')
                txt_time = os.path.getmtime(filename + '.txt')
                if txt_time > json_time:
                    return cls.from_txt_file(filename + '.txt', ignore_invalid_lines=ignore_invalid_lines)
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
        metadata['metadata_path'] = json_file
        return metadata

    @classmethod
    def from_txt_file(cls, txt_file, ignore_invalid_lines=False, cache_json=True):
        """Read metadata.txt file into Metadata object with valid python identifiers for keys

        A standard metadata.txt file is close to being an executable python script that just defines
        a bunch of constants.  The three main problems with the metadata.txt format are:

          1) variable names contain dashes, which is the subtraction operator in python,
          2) strings are not enclosed in quotes, and
          3) lists are not enclosed in brackets

        It is easy to correct these problems.  In particular, (1) is resolved by changing dashes to
        underscores in the identifiers.  A bug in SpEC's metadata.txt files -- whereby some comment
        lines are missing the initial `#` -- is also fixed.  There are also occasional other
        problems, like commas missing from lists.  All syntax errors as of this writing are fixed in
        this function.

        Note that this function is not very flexible when it comes to generalizing the syntax of the
        metadata.txt files.  In particular, it assumes that the right-hand sides are either numbers
        or strings (or lists of either numbers or strings).  For example, I think I've seen cases
        where the eccentricity is given as something like "<1e-5".  Since python has no "less-than"
        type, this is converted to a string.  But generally, this does seem to work on metadata.txt
        files in the SXS waveform repository.

        """
        # This is considered the safe way to evaluate strings "containing Python values from
        # untrusted sources without the need to parse the values oneself".  We will use it to parse
        # the right-hand side expressions from the metadata.txt file, once we've appropriately
        # modified them, so that python will get to decide what should be a string, float, int,
        # list, etc.
        from ast import literal_eval

        assignment_pattern = re.compile(r"""([-A-Za-z0-9]+)\s*=\s*(.*)""")
        string_pattern = re.compile(r"""[A-DF-Za-df-z<>@]""")  # Ignore 'e' and 'E' because they may appear in numbers
        multispace_pattern = re.compile(r"""\s+""")
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
                        quantity = '[]'
                    else:
                        q = quantity.strip()
                        if '[unknown]' in q.lower():
                            metadata[variable] = 'NaN'
                            continue
                        elif ((q.startswith('"') and q.endswith('"'))
                            or (q.startswith("'") and q.endswith("'"))):
                            # If the whole thing is quoted, just leave it as is
                            quantity = q
                        elif string_pattern.search(quantity):
                            # If this is a string, strip whitespace from it, split lists and place
                            # brackets around them, and place quotation marks around each element
                            quantities = [q.strip() for q in quantity.split(",")]
                            if "," in quantity:
                                quantity = "['" + "', '".join(quantities) + "']"
                            else:
                                quantity = "'" + quantities[0] + "'"
                        else:
                            # Otherwise, just place brackets around lists of non-strings
                            quantity = quantity.strip()
                            if "," in quantity:
                                quantity = "[" + quantity + "]"
                            elif " " in quantity:
                                quantity = "[" + multispace_pattern.sub(',', quantity) + "]"

                    # Add this line to the metadata, whether or not it's been modified
                    try:
                        metadata[variable] = literal_eval(quantity)
                    except SyntaxError as e:
                        message = ("\nWhile parsing {0}, transformed input text:\n".format(txt_file)
                                   + "    " + line.rstrip()
                                   + "\ninto quantity\n"
                                   + "    " + variable + " = " + quantity
                                   + "\nParsing this using `ast.literal_eval` resulted in a SyntaxError.\n")
                        if cache_json and ignore_invalid_lines:
                            cache_json = False
                            message += "JSON caching will be turned off for this file until the error is fixed.\n"
                        warnings.warn(message)
                        if not ignore_invalid_lines:
                            raise e

        if cache_json:
            # Skip the text processing next time, and just go straight to json
            txt_index = txt_file.rfind('.txt')
            if txt_index < 1:
                json_file = txt_file + '.json'
            else:
                json_file = txt_file[:txt_index] + '.json'
            metadata.to_json_file(json_file)

        metadata['metadata_path'] = txt_file
        return metadata

    def to_json(self, indent=4, separators=(',', ': ')):
        return json.dumps(self, indent=indent, separators=separators)

    def to_json_file(self, json_file, indent=4, separators=(',', ': ')):
        with open(json_file, 'w') as f:
            f.write(self.to_json(indent=indent, separators=separators))

    def to_txt(self):
        def deformat(value):
            """Basically undo the nice formatting of `from_txt_file`"""
            if isinstance(value, list):
                return ','.join(['{0}'.format(item) for item in value])
            else:
                return '{0}'.format(value)
        return '\n'.join(['{0} = {1}'.format(_valid_identifier_to_metadata_key(key), deformat(self[key])) for key in self])

    def to_txt_file(self, txt_file):
        with open(txt_file, 'w') as f:
            f.write(self.to_txt() + '\n')

    def __init__(self, *args, **kwargs):
        """Initialize the OrderedDict, ensuring that all keys have been converted to valid identifiers

        Note that this constructor is not frequently used directly from outside this code.  See
        Metadata.from_file, etc., for more useful factory functions.

        This function intercepts the allowed args and kwargs and converts any keys before simply
        calling the base class's initialization function.

        """
        if len(args) > 0:
            args = list(args)
            if isinstance(args[0], collections.Mapping):
                mapping = args[0]
                args[0] = OrderedDict([(_valid_identifier(key), mapping[key]) for key in mapping])
            else:
                iterable = args[0]
                args[0] = [(_valid_identifier(k), v) for k, v in iterable]
        if len(kwargs) > 0:
            for key in list(kwargs):
                kwargs[_valid_identifier(key)] = kwargs.pop(key)
        super(Metadata, self).__init__(*args, **kwargs)

    def add_com_parameters(self, raise_on_errors=False):
        """Add any translation and boost parameters found in all CoM files in this directory"""
        import os.path
        import glob
        import h5py

        path = os.path.dirname(self.get('metadata_path', '.'))
        com_parameters = self.get('com_parameters', {})
        for file_name in glob.glob(os.path.join(path, '*_CoM.h5')):
            try:
                with h5py.File(file_name, 'r') as f:
                    for g in f:
                        g_parameters = {}
                        if hasattr(f[g], 'attrs'):
                            if 'space_translation' in f[g].attrs:
                                g_parameters['space_translation'] = list(f[g].attrs['space_translation'])
                            if 'boost_velocity' in f[g].attrs:
                                g_parameters['boost_velocity'] = list(f[g].attrs['boost_velocity'])
                        if g_parameters:
                            com_parameters['{0}/{1}'.format(os.path.basename(file_name), g)] = g_parameters
            except:
                if raise_on_errors:
                    raise
        if com_parameters:
            self['com_parameters'] = com_parameters
        return self

    def add_standard_parameters(self, raise_on_errors=False):
        """Add standard parameters that aren't included in the default metadata fields

        New parameters include 'object_types', 'initial_mass_ratio', and 'relaxed_mass_ratio'.  If
        'relaxed_dimensionless_spin*' are not present, but the parameters necessary to compute them
        are, they are also added.  Finally, we also add 'relaxed_chi_eff', 'relaxed_chi1_perp', and
        'relaxed_chi2_perp'.

        """
        import numpy as np
        if 'object1' in self and 'object2' in self:
            self['object_types'] = ''.join(sorted([self['object1'].upper(), self['object2'].upper()]))
        if 'initial_mass1' in self and 'initial_mass2' in self:
            try:
                mass_ratio = float(self['initial_mass1']) / float(self['initial_mass2'])
                self['initial_mass_ratio'] = mass_ratio
            except:
                if raise_on_errors:
                    raise
        if 'relaxed_mass1' in self and 'relaxed_mass2' in self:
            try:
                mass_ratio = float(self['relaxed_mass1']) / float(self['relaxed_mass2'])
                self['relaxed_mass_ratio'] = mass_ratio
            except:
                if raise_on_errors:
                    raise
        if 'relaxed_dimensionless_spin1' not in self:
            if 'relaxed_spin1' in self and 'relaxed_mass1' in self:
                try:
                    self['relaxed_dimensionless_spin1'] = list(np.array(self['relaxed_spin1']) / self['relaxed_mass1']**2)
                except:
                    if raise_on_errors:
                        raise
        if 'relaxed_dimensionless_spin2' not in self:
            if 'relaxed_spin2' in self and 'relaxed_mass2' in self:
                try:
                    self['relaxed_dimensionless_spin2'] = list(np.array(self['relaxed_spin2']) / self['relaxed_mass2']**2)
                except:
                    if raise_on_errors:
                        raise
        if 'initial_dimensionless_spin1' not in self:
            if 'initial_spin1' in self and 'initial_mass1' in self:
                try:
                    self['initial_dimensionless_spin1'] = list(np.array(self['initial_spin1']) / self['initial_mass1']**2)
                except:
                    if raise_on_errors:
                        raise
        if 'initial_dimensionless_spin2' not in self:
            if 'initial_spin2' in self and 'initial_mass2' in self:
                try:
                    self['initial_dimensionless_spin2'] = list(np.array(self['initial_spin2']) / self['initial_mass2']**2)
                except:
                    if raise_on_errors:
                        raise
        if ('relaxed_mass1' in self and 'relaxed_mass2' in self and 'relaxed_orbital_frequency' in self
            and 'relaxed_dimensionless_spin1' in self and 'relaxed_dimensionless_spin2' in self):
                try:
                    m1 = float(self['relaxed_mass1'])
                    m2 = float(self['relaxed_mass2'])
                    chi1 = np.array(self['relaxed_dimensionless_spin1'], dtype=float)
                    chi2 = np.array(self['relaxed_dimensionless_spin2'], dtype=float)
                    L = np.array(self['relaxed_orbital_frequency'], dtype=float)
                    L /= np.linalg.norm(L)
                    chi1L = np.dot(chi1, L)
                    chi2L = np.dot(chi2, L)
                    chi1perp = np.cross(chi1, L)
                    chi2perp = np.cross(chi2, L)
                    self['relaxed_chi_eff'] = (m1*chi1L+m2*chi2L)/(m1+m2)
                    self['relaxed_chi1_perp'] = np.linalg.norm(chi1perp)
                    self['relaxed_chi2_perp'] = np.linalg.norm(chi2perp)
                except:
                    if raise_on_errors:
                        raise
        return self

    def add_extras(self, raise_on_errors=False):
        """Add information to the metadata from other files in its directory"""
        self.add_com_parameters(raise_on_errors=raise_on_errors)
        self.add_standard_parameters(raise_on_errors=raise_on_errors)
        return self

    def reorder_keys(self, order=None):
        """Return a copy of this object with keys reordered

        It is sometimes nice to reorder the keys of the metadata to display the most interesting
        quantities first.  The usual order output by SpEC, for example, hides crucial quantities
        like the masses and spins after lots of uninteresting keys like the author list and various
        bibtex groups.  This function allows the keys to be reordered using exact matches and
        regular expressions.

        """
        import re
        if order is None:
            order = [
                'simulation_name',
                'alternative_names',
                'initial_data_type',
                'eos',
                'object_types',
                'number_of_orbits',
                'relaxed_mass_ratio',
                'relaxed_chi_eff',
                'relaxed_chi1_perp',
                'relaxed_chi2_perp',
                'relaxed_eccentricity',
                'relaxed_dimensionless_spin1',
                'relaxed_dimensionless_spin2',
                'relaxed_orbital_frequency',
                'relaxed_mass1',
                'relaxed_mass2',
                'relaxed.*',
            ]
        original_keys = list(self)
        new = type(self)()
        for ordered_key in order:
            if ordered_key in original_keys:
                new[ordered_key] = self[ordered_key]
                original_keys.remove(ordered_key)
            else:
                key_pattern = re.compile(ordered_key)
                for key in list(original_keys):  # Iterate over a *copy* of the original_keys list
                    if key_pattern.match(key):
                        new[key] = self[key]
                        original_keys.remove(key)
        for key in original_keys:
            new[key] = self[key]
        return new

    @classmethod
    def fromkeys(cls, iterable):
        iterable = [(_valid_identifier(k), v) for k, v in iterable]
        return super(Metadata, self).fromkeys(iterable)

    @property
    def resolution(self):
        """Try to determine the resolution from the 'simulation-name' field"""
        simulation_name = self['simulation_name']
        last_slash_index = simulation_name.rindex(os.sep)
        return simulation_name[last_slash_index+1:]

    @property
    def lev(self):
        """Try to determine an integer "Lev" number from the 'simulation-name' field"""
        resolution = self.resolution
        return int(resolution.replace('Lev', ''))

    @property
    def simulation_group(self):
        """Remove any trailing '/LevN' part of the simulation-name"""
        simulation_name = self['simulation_name']
        last_slash_index = simulation_name.rindex(os.sep + 'Lev')
        if last_slash_index < 1:
            last_slash_index = len(simulation_name)
        return simulation_name[:last_slash_index]

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
        return sorted(set(super(Metadata, self).__dir__() + list(self.keys())))

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
        if isinstance(mapping_or_iterable, collections.Mapping):
            mapping_or_iterable = collections.OrderedDict([(_valid_identifier(key), mapping_or_iterable[key]) for key in mapping_or_iterable])
        elif isinstance(mapping_or_iterable, collections.Iterable):
            mapping_or_iterable = [(_valid_identifier(k), v) for k, v in mapping_or_iterable]
        return super(Metadata, self).update(mapping_or_iterable)


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
