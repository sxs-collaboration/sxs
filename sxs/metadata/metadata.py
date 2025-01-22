"""Container for metadata of individual simulations"""

import re
import collections


_valid_identifier_pattern = re.compile(r'\W|^(?=\d)')


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
    """Interface to metadata

    Note that the constructor is not generally useful from outside this class.  See
    Metadata.from_file, etc., for more useful initialization functions.

    This object is essentially a `collections.OrderedDict`, with a few extra features:
      1) Keys are always forced to be valid python identifiers, whether setting or
         getting.
      2) There are a few extra methods for constructing these objects from json data or
         files or txt files, and outputting them to json or txt files.
      3) If an attribute does not exist on the object, the keys are searched.  This
         allows for tab-completion on the key names.
      4) Properties (methods without parentheses) that attempt to automatically
         determine the resolution and lev number from the 'simulation-name'.

    """

    @classmethod
    def from_file(cls, file_name, ignore_invalid_lines=False, cache_json=True):
        """Read a file into a Metadata object

        Parameters
        ----------
        file_name : str
            The input `file_name` may either end in `.txt`, for an old-style
            metadata.txt file, or in `.json`, for a JSON file.  If neither ending is
            present, the function searches for files with the given `file_name` plus
            both `.txt` and `.json`, and reads the newer one or at least the existing
            one.
        ignore_invalid_lines : bool, optional
            If True, invalid lines will be ignored; otherwise, an error will be raised.
        cache_json : bool, optional
            If True, and the `.json` file does not exist or is older than the `.txt`
            file, a new `.json` file will be created with contents equivalent to the
            `.txt` file.

        Raises
        ------
        SyntaxError : on  `.txt` parse errors when `ignore_invalid_lines` is  False.

        """
        from pathlib import Path

        path = Path(file_name).expanduser().resolve()
        json_path = path.with_suffix(".json")
        txt_path = path.with_suffix(".txt")

        if path.suffix == ".json":
            return cls.from_json_file(json_path)
        elif path.suffix == ".txt":
            return cls.from_txt_file(txt_path, ignore_invalid_lines=ignore_invalid_lines, cache_json=cache_json)
        else:
            json_exists = json_path.exists()
            txt_exists = txt_path.exists()
            if json_exists and not txt_exists:
                return cls.from_json_file(json_path)
            elif txt_exists and not json_exists:
                return cls.from_txt_file(txt_path, ignore_invalid_lines=ignore_invalid_lines, cache_json=cache_json)
            elif json_exists and txt_exists:
                json_time = json_path.stat().st_mtime
                txt_time = txt_path.stat().st_mtime
                if txt_time > json_time:
                    return cls.from_txt_file(txt_path, ignore_invalid_lines=ignore_invalid_lines, cache_json=cache_json)
                else:
                    return cls.from_json_file(json_path)
            else:
                raise ValueError(f"Could not find file named '{json_path}' or '{txt_path}'")

    load = from_file

    @classmethod
    def from_json_data(cls, json_data):
        """Read metadata from a file-like object (not a file)

        Parameters
        ----------
        json_data : file-like object
            Note that this cannot be just a file name; it must be a file-like object
            (such as an open file handle).  See the `from_json_file` function if you
            just want to pass the path to a file.

        See Also
        --------
        sxs.Metadata.from_json_file : `.json` files
        sxs.Metadata.from_file : reads `.txt` or `.json` files

        """
        import json
        # noinspection PyTypeChecker
        return json.load(json_data, object_pairs_hook=cls)

    @classmethod
    def from_json_file(cls, json_file):
        """Read metadata.json file

        Parameters
        ----------
        json_file : str
            The path to an `metadata.json` file.

        See Also
        --------
        sxs.Metadata.from_file : reads `.txt` or `.json` files

        """
        import json
        from pathlib import Path
        path = Path(json_file).expanduser().resolve().with_suffix(".json")
        with path.open(mode="r") as metadata_file:
            # noinspection PyTypeChecker
            metadata = json.load(metadata_file, object_pairs_hook=cls)
        metadata["metadata_path"] = str(json_file)
        return metadata

    @classmethod
    def from_txt_file(cls, txt_file, ignore_invalid_lines=False, cache_json=True):
        """Read metadata.txt file

        Parameters
        ----------
        txt_file : str
            The path to an old-style `metadata.txt` file.
        ignore_invalid_lines : bool, optional
            If True, invalid lines will be ignored; otherwise, an error will be raised.
        cache_json : bool, optional
            If True, a new `.json` file will be created with contents equivalent to the
            `.txt` file.

        See Also
        --------
        sxs.Metadata.from_file : reads `.txt` or `.json` files

        Notes
        -----
        A standard metadata.txt file is close to being an executable python script that
        just defines a bunch of constants.  The three main problems with the
        metadata.txt format are

          1) variable names contain dashes, which is the subtraction operator in python,
          2) strings are not enclosed in quotes, and
          3) lists are not enclosed in brackets

        It is easy to correct these problems.  In particular, (1) is resolved by
        changing dashes to underscores in the identifiers.  A bug in SpEC
        metadata.txt files -- whereby some comment lines are missing the initial `#` --
        is also fixed.  There are also occasional other problems, like commas missing
        from lists.  All syntax errors as of this writing are fixed in this function.

        Note that this function is not very flexible when it comes to generalizing the
        syntax of the metadata.txt files.  In particular, it assumes that the
        right-hand sides are either numbers or strings (or lists of either numbers or
        strings).  For example, I think I've seen cases where the eccentricity is given
        as something like "<1e-5".  Since python has no "less-than" type, this is
        converted to a string.  But generally, this does seem to work on metadata.txt
        files in the SXS waveform repository.

        """
        # This is considered the safe way to evaluate strings "containing Python values
        # from untrusted sources without the need to parse the values oneself".  We
        # will use it to parse the right-hand side expressions from the metadata.txt
        # file, once we've appropriately modified them, so that python will get to
        # decide what should be a string, float, int, list, etc.
        import warnings
        from ast import literal_eval

        from pathlib import Path
        path = Path(txt_file).expanduser().resolve().with_suffix(".txt")

        assignment_pattern = re.compile(r"""([-A-Za-z0-9]+)\s*=\s*(.*)""")
        string_pattern = re.compile(r"""[A-DF-Za-df-z<>@]""")  # Ignore "e" and "E" because they may appear in numbers
        multi_space_pattern = re.compile(r"""\s+""")
        metadata = cls()

        with path.open("r") as metadata_file:
            for line in metadata_file:
                # Work around bug where some lines of dashes are missing the comment character
                if line.startswith("-"):
                    continue

                # Process each assignment line, skipping comments and unrecognized lines
                match = assignment_pattern.match(line)
                if match:
                    variable, quantity = match.groups()

                    # It was a stupid choice to make variables contain dashes
                    variable = _valid_identifier(variable)

                    # If `quantity` is an empty string, we should just replace it with an empty list
                    if not quantity or quantity == "\n":
                        quantity = "[]"
                    else:
                        q = quantity.strip()
                        if "[unknown]" in q.lower():
                            metadata[variable] = "NaN"
                            continue
                        elif (q.startswith('"') and q.endswith('"')) or (q.startswith("'") and q.endswith("'")):
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
                                quantity = "[" + multi_space_pattern.sub(',', quantity) + "]"

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
            metadata.to_json_file(path.with_suffix(".json"))

        metadata["metadata_path"] = str(txt_file)
        return metadata

    def to_json(self, indent=4, separators=(",", ": ")):
        """Export to JSON string"""
        import json
        return json.dumps(self, indent=indent, separators=separators)

    def to_json_file(self, json_file, indent=4, separators=(",", ": ")):
        """Write to JSON file"""
        from pathlib import Path
        path = Path(json_file).expanduser().resolve().with_suffix(".json")
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as f:
            f.write(self.to_json(indent=indent, separators=separators))

    save = to_json_file

    def to_txt(self):
        """Export to string like metadata.txt contents"""
        def deformat(value):
            """Basically undo the nice formatting of `from_txt_file`"""
            if isinstance(value, list):
                return ",".join(["{0}".format(item) for item in value])
            else:
                return f"{value}"
        return "\n".join([f"{_valid_identifier_to_metadata_key(key)} = {deformat(self[key])}" for key in self])

    def to_txt_file(self, txt_file):
        """Write to file in metadata.txt format"""
        from pathlib import Path
        path = Path(txt_file).expanduser().resolve().with_suffix(".txt")
        path.mkdir(parents=True, exist_ok=True)
        with path.open("w") as f:
            f.write(self.to_txt() + "\n")

    def __init__(self, *args, **kwargs):
        """Initialize the OrderedDict, converting all keys to valid identifiers

        Note that the constructor is not generally useful from outside this class.  See
        Metadata.from_file, etc., for more useful initialization functions.

        This function intercepts the allowed args and kwargs and converts any keys
        before simply calling the base class's initialization function.

        """
        import collections
        if len(args) > 0:
            args = list(args)
            if isinstance(args[0], collections.abc.Mapping):
                mapping = args[0]
                args[0] = collections.OrderedDict([(_valid_identifier(key), mapping[key]) for key in mapping])
            else:
                iterable = args[0]
                args[0] = [(_valid_identifier(k), v) for k, v in iterable]
        if len(kwargs) > 0:
            for key in list(kwargs):
                kwargs[_valid_identifier(key)] = kwargs.pop(key)
        super(Metadata, self).__init__(*args, **kwargs)

    def add_com_parameters(self, first_only=True, raise_on_errors=False):
        """Add any translation and boost parameters found in all CoM files in this directory

        Adds a new key `com_parameters` to the top level of the metadata dictionary
        containing the `space_translation` and `boost_velocity` parameters for the COM
        correction.

        Parameters
        ----------
        first_only : bool, optional
            If True, add the first set of parameters directly under the top-level key
            `com_parameters`; otherwise, add separate entries under that key for each
            file ending in `_CoM.h5`.
        raise_on_errors : bool, optional
            If False, suppress any exceptions that happen in the core loop of this
            function; otherwise, raise.

        """
        import os.path
        import glob
        import h5py
        import re
        import ast

        path = os.path.dirname(self.get("metadata_path", "."))
        com_parameters = self.get("com_parameters", {})
        for file_name in reversed(sorted(glob.glob(os.path.join(path, "*_CoM.h5")))):
            try:
                with h5py.File(file_name, "r") as f:
                    for g in f:
                        g_parameters = {}
                        if hasattr(f[g], "attrs") and "space_translation" in f[g].attrs:
                            g_parameters["space_translation"] = list(f[g].attrs["space_translation"])
                        if hasattr(f[g], "attrs") and "boost_velocity" in f[g].attrs:
                            g_parameters["boost_velocity"] = list(f[g].attrs["boost_velocity"])
                        if "History.txt" in f[g]:
                            history = f[g]["History.txt"][()]
                            if hasattr(history, "decode"):
                                history = history.decode()
                            for parameter_name in ["boost_velocity", "space_translation"]:
                                if parameter_name not in g_parameters:
                                    pattern = rf'"{parameter_name}": array\((.*?)\)'
                                    matches = re.search(pattern, history)
                                    if matches:
                                        g_parameters[parameter_name] = ast.literal_eval(matches.group(1))
                        if first_only and "space_translation" in g_parameters and "boost_velocity" in g_parameters:
                            self["com_parameters"] = g_parameters
                            return self
                        if g_parameters:
                            com_parameters["{0}/{1}".format(os.path.basename(file_name), g)] = g_parameters
            except Exception:
                if raise_on_errors:
                    raise
        if com_parameters:
            self["com_parameters"] = com_parameters
        return self

    def add_standard_parameters(self, raise_on_errors=False):
        """Add standard parameters that aren't included in the default metadata fields

        New parameters include 'object_types', 'initial_mass_ratio', and
        'reference_mass_ratio'.  If 'reference_dimensionless_spin*' are not present,
        but the parameters necessary to compute them are, they are also added.
        Finally, we also add 'reference_chi_eff', 'reference_chi1_perp', and
        'reference_chi2_perp'.

        """
        import math
        import numpy as np

        def stringify_nan(number):
            if math.isnan(number):
                return "NaN"
            else:
                return number

        def stringify_nans(array):
            return [stringify_nan(number) for number in array]

        if "object1" in self and "object2" in self:
            self["object_types"] = "".join(sorted([self["object1"].upper(), self["object2"].upper()]))
        if "initial_mass1" in self and "initial_mass2" in self:
            try:
                mass_ratio = float(self["initial_mass1"]) / float(self["initial_mass2"])
                self["initial_mass_ratio"] = stringify_nan(mass_ratio)
            except Exception:
                if raise_on_errors:
                    raise
        if "reference_mass1" in self and "reference_mass2" in self:
            try:
                mass_ratio = float(self["reference_mass1"]) / float(self["reference_mass2"])
                self["reference_mass_ratio"] = stringify_nan(mass_ratio)
            except Exception:
                if raise_on_errors:
                    raise
        if "reference_dimensionless_spin1" not in self:
            if "reference_spin1" in self and "reference_mass1" in self:
                try:
                    self["reference_dimensionless_spin1"] = stringify_nans(
                        np.array(self["reference_spin1"]) / self["reference_mass1"]**2
                    )
                except Exception:
                    if raise_on_errors:
                        raise
        if "reference_dimensionless_spin2" not in self:
            if "reference_spin2" in self and "reference_mass2" in self:
                try:
                    self["reference_dimensionless_spin2"] = stringify_nans(
                        np.array(self["reference_spin2"]) / self["reference_mass2"]**2
                    )
                except Exception:
                    if raise_on_errors:
                        raise
        if "initial_dimensionless_spin1" not in self:
            if "initial_spin1" in self and "initial_mass1" in self:
                try:
                    self["initial_dimensionless_spin1"] = stringify_nans(
                        np.array(self["initial_spin1"]) / self["initial_mass1"]**2
                    )
                except Exception:
                    if raise_on_errors:
                        raise
        if "initial_dimensionless_spin2" not in self:
            if "initial_spin2" in self and "initial_mass2" in self:
                try:
                    self["initial_dimensionless_spin2"] = stringify_nans(
                        np.array(self["initial_spin2"]) / self["initial_mass2"]**2
                    )
                except Exception:
                    if raise_on_errors:
                        raise
        if ("reference_mass1" in self and "reference_mass2" in self and "reference_orbital_frequency" in self
                and "reference_dimensionless_spin1" in self and "reference_dimensionless_spin2" in self):
            try:
                m1 = float(self["reference_mass1"])
                m2 = float(self["reference_mass2"])
                chi1 = np.array(self["reference_dimensionless_spin1"], dtype=float)
                chi2 = np.array(self["reference_dimensionless_spin2"], dtype=float)
                L = np.array(self["reference_orbital_frequency"], dtype=float)
                L /= np.linalg.norm(L)
                chi1L = np.dot(chi1, L)
                chi2L = np.dot(chi2, L)
                chi1perp = np.cross(chi1, L)
                chi2perp = np.cross(chi2, L)
                self["reference_chi_eff"] = stringify_nan((m1*chi1L+m2*chi2L)/(m1+m2))
                self["reference_chi1_perp"] = stringify_nan(np.linalg.norm(chi1perp))
                self["reference_chi2_perp"] = stringify_nan(np.linalg.norm(chi2perp))
            except Exception:
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

        It is sometimes nice to reorder the keys of the metadata to display the most
        interesting quantities first.  The usual order output by SpEC, for example,
        hides crucial quantities like the masses and spins after lots of uninteresting
        keys like the author list and various bibtex groups.  This function allows the
        keys to be reordered using exact matches and regular expressions.

        """
        import re
        if order is None:
            order = [
                "url",
                "simulation_name",
                "alternative_names",
                "initial_data_type",
                "eos",
                "object_types",
                "number_of_orbits",
                "reference_mass_ratio",
                "reference_chi_eff",
                "reference_chi1_perp",
                "reference_chi2_perp",
                "reference_eccentricity",
                "reference_dimensionless_spin1",
                "reference_dimensionless_spin2",
                "reference_orbital_frequency",
                "reference_mass1",
                "reference_mass2",
                "reference.*",
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

    # @classmethod
    # def fromkeys(cls, iterable):
    #     iterable = [(_valid_identifier(k), v) for k, v in iterable]
    #     return super(Metadata, cls).fromkeys(iterable)

    @property
    def resolution(self):
        """Try to determine the resolution from the "simulation-name" field"""
        import os
        simulation_name = self["simulation_name"]
        last_slash_index = simulation_name.rindex(os.sep)
        return simulation_name[last_slash_index+1:]

    @property
    def lev(self):
        """Try to determine an integer "Lev" number from the "simulation-name" field"""
        resolution = self.resolution
        return int(resolution.replace("Lev", ""))

    @property
    def simulation_group(self):
        """Remove any trailing "/LevN" part of the simulation-name"""
        import os
        simulation_name = self["simulation_name"]
        last_slash_index = simulation_name.rindex(os.sep + "Lev")
        if last_slash_index < 1:
            last_slash_index = len(simulation_name)
        return simulation_name[:last_slash_index]

    ###############################################################################
    ###############################################################################
    #                                                                             #
    # The following methods mirror equivalent methods in OrderedDict, but also    #
    # ensure that any keys are converted to valid identifiers first.  This        #
    # enables this object to be used more like the original metadata.txt format.  #
    #                                                                             #
    ###############################################################################
    ###############################################################################

    def __contains__(self, key):
        return super(Metadata, self).__contains__(_valid_identifier(key))

    def __delattr__(self, name):
        super(Metadata, self).__delattr__(_valid_identifier(name))

    def __delitem__(self, key):
        super(Metadata, self).__delitem__(_valid_identifier(key))

    def __getattribute__(self, name):
        """Include keys as attributes

        This allows retrieval of a key like `md["simulation_name"]` as `md.simulation_name`.

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
        return list(sorted(set(super(Metadata, self).__dir__()))) + list(self.keys())

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

    def update(self, mapping_or_iterable=None, **kwargs):
        if isinstance(mapping_or_iterable, collections.abc.Mapping):
            mapping_or_iterable = collections.OrderedDict(
                [(_valid_identifier(key), mapping_or_iterable[key]) for key in mapping_or_iterable]
            )
        elif isinstance(mapping_or_iterable, collections.abc.Iterable):
            mapping_or_iterable = [(_valid_identifier(k), v) for k, v in mapping_or_iterable]
        return super(Metadata, self).update(mapping_or_iterable)
