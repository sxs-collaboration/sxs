#!/usr/bin/env python

# Copyright (c) 2018, Michael Boyle
# See LICENSE file for details: <https://github.com/moble/quaternion/blob/master/LICENSE>

# Construct the version number from the date and time this python version was created.
from os import environ
if "package_version" in environ:
    version = environ["package_version"]
    print("Setup.py using environment version='{0}'".format(version))
else:
    print("The variable 'package_version' was not present in the environment")
    try:
        # For cases where this is being installed from git.  This gives the true version number.
        from sys import platform
        from subprocess import check_output
        on_windows = ('win' in platform.lower() and not 'darwin' in platform.lower())
        use_shell = not on_windows
        version = check_output("""git log -1 --format=%cd --date=format:'%Y.%-m.%-d.%-H.%-M.%-S'""", shell=use_shell).decode('ascii').rstrip()
        print("Setup.py using git log version='{0}'".format(version))
    except:
        # For cases where this isn't being installed from git.  This gives the wrong version number,
        # but at least it provides some information.
        try:
            from time import strftime, gmtime
            try:
                version = strftime("%Y.%-m.%-d.%-H.%-M.%-S", gmtime())
            except ValueError:  # because Windows
                version = strftime("%Y.%m.%d.%H.%M.%S", gmtime())
            print("Setup.py using strftime version='{0}'".format(version))
        except:
            version = '0.0.0'
            print("Setup.py failed to determine the version; using '{0}'".format(version))
with open('_version.py', 'w') as f:
    f.write('__version__ = "{0}"'.format(version))


# long_description = """"""

if __name__ == "__main__":
    from os import getenv
    from setuptools import setup
    setup(name='sxs',
          packages = ['sxs', 'sxs.metadata', 'sxs.doxygen', 'sxs.references', 'sxs.zenodo', 'sxs.zenodo.api'],
          scripts = ['scripts/sxs'],
          version=version,
          description = 'A collection of python code used by the SXS collaboration',
          url='https://github.com/moble/sxs',
          author='Michael Boyle',
          author_email='mob22@cornell.edu',
          # long_description=long_description,
          # download_url = 'https://github.com/moble/sxs/archive/master.tar.gz',
    )
