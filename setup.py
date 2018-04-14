from setuptools import setup
from distutils.command.build_py import build_py


def calculate_version():
    """Construct a version string from date, time, and short git hash

    The resulting version looks something like "2018.4.11.12.56.27.dev71379535", where the first six
    numbers are the date and time, while the part after "dev" is the short git hash converted to
    decimal.  This version string has a few nice features:
      * This is a valid PEP 440 version number
      * It's automatic
      * It's serial
      * It contains all the information I can think of that might be needed

    Unfortunately, in order to comply with PEP 440 (and more specifically, not break various python
    tools that analyze version numbers), the short git hash must be converted from hexadecimal to
    decimal.  To convert it back, you can run a python command like

        >>> hex(71379535)[2:]
        '4412a4f'

    """
    import subprocess
    from warnings import warn
    from sys import platform
    on_windows = ('win' in platform.lower() and not 'darwin' in platform.lower())
    use_shell = not on_windows
    try:
        git_revision = subprocess.check_output("""git show -s --format="%ci %h" HEAD""", shell=use_shell).decode('ascii').rstrip()
        date, time, utc_offset, short_hash = git_revision.split(' ')
        date = date.replace('-', '.').strip()  # make date an acceptable version string
        time = time.replace(':', '.').strip()  # make time an acceptable version string
        short_hash = short_hash.strip()  # remove newline and any other whitespace
        short_hash = int(short_hash, 16)  # So that it's a valid PEP 440 version identifier
        dirty = bool(subprocess.call("git diff-files --quiet --", shell=use_shell))
        dirty = dirty or bool(subprocess.call("git diff-index --cached --quiet HEAD --", shell=use_shell))
        version = '{0}.{1}.dev{2}'.format(date, time, short_hash)
        if dirty:
            version += '+dirty'
        exec('putative__version__ = "{0}"'.format(version))  # see if this will raise an error for some reason
    except Exception as e:
        # If any of the above failed for any reason whatsoever, fall back on this dumb version
        from os import getenv
        from datetime import datetime
        warning = ('\nThe `calculate_version` function failed to get the git version.'
                   'Maybe your version of python (<2.7?) is too old.  Here\'s the exception:\n' + str(e) + '\n'
                   'This should not be a problem, unless you need an accurate version number.'
                   'Continuing on, in spite of it all...\n')
        warn(warning)
        date, time = datetime.now().isoformat().split('T')
        date = date.replace('-', '.').strip()
        time = time[:8].replace(':', '.').strip()
        version = '0.0.0.dev' + date + '.' + time
    return version


class build_py_and_copy_version(build_py):
    """Add version-copying step to standard build_py"""
    def run(self):
        build_py.run(self)
        version = calculate_version()
        print('build_py_copy_version using __version__ = "{0}"'.format(version))
        if not self.dry_run:
            import os.path
            for package in self.packages:
                with open(os.path.join(self.build_lib, os.path.join(*package.split('.')), '_version.py'), 'w') as fobj:
                    fobj.write('__version__ = "{0}"'.format(version))


setup(
    name = 'sxs',
    packages = ['sxs', 'sxs.metadata', 'sxs.doxygen', 'sxs.references', 'sxs.zenodo', 'sxs.zenodo.api'],
    scripts = ['sxs/doxygen/validate_xml.py'],
    version=calculate_version(),
    cmdclass={'build_py': build_py_and_copy_version},
    description = 'A collection of python code used by the SXS collaboration',
    author = 'SXS Collaboration',
    author_email = 'web-admin@black-holes.org',
    url = 'https://www.black-holes.org/',
    # download_url = 'https://github.com/moble/sxs/archive/master.tar.gz',
)
