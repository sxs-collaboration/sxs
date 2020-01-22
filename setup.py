#!/usr/bin/env python

# Copyright (c) 2019, Simulating eXtreme Spacetimes Collaboration
# See LICENSE file for details: <https://github.com/sxs-collaboration/sxs/blob/master/LICENSE>

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
with open('sxs/_version.py', 'w') as f:
    f.write('__version__ = "{0}"'.format(version))


long_description = """\
This package provides a number of utilities for use by the SXS collaboration, and others who use our data.
For example, the `metadata` subpackage provides functions for reading and analyzing the metadata files
provided with SXS simulations for describing the physics they represent.  Another important subpackage is
`zenodo`, which provides handy functions for interacting with zenodo.org in general, and in particular our
collection of simulation data on that site.

A handy command-line interface is also installed when you install this package.  It enables direct calling to
any function in this package.  For example, if you want to upload a directory to Zenodo, you could use the
`sxs.zenodo.upload` function through python, or you could just run the command `sxs zenodo upload`.  Any
arguments to these functions should be passed immediately following the function name, followed by any keyword
arguments in the usual `--keyword=value` format.  Run `sxs --help` for more information.

"""

if __name__ == "__main__":
    from os import getenv
    from setuptools import setup
    setup(name='sxs',
          packages = ['sxs', 'sxs.metadata', 'sxs.doxygen',
                      'sxs.format', 'sxs.format.lvc', 'sxs.references',
                      'sxs.utilities', 'sxs.utilities.decimation',
                      'sxs.validate', 'sxs.zenodo', 'sxs.zenodo.api'],
          scripts = ['scripts/sxs'],
          include_package_data=True,
          version=version,
          description = 'A collection of python code used by the SXS collaboration',
          url='https://github.com/sxs-collaboration/sxs',
          author='Michael Boyle',
          author_email='mob22@cornell.edu',
          long_description=long_description,
          install_requires=[
              'numpy',
              'scipy',
              'h5py',
              'requests',  # For interacting over HTTP
              'requests_toolbelt',  # For dumping information about requests/responses
              'ads',  # For searching ADS via the API
              'pylatexenc',  # For converting unicode to latex
              'lxml',  # For parsing doxygen
              'feedparser',  # For parsing arxiv responses
              'tqdm',  # For nice progress bars
              'pytz',  # For timezone information
          ],
          # download_url = 'https://github.com/moble/sxs/archive/master.tar.gz',
    )
