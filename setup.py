from distutils.core import setup
from setup_helpers import calculate_version, build_py_and_copy_version

setup(
    name = 'sxs',
    packages = ['sxs', 'sxs.metadata'],
    version=calculate_version(),
    cmdclass={'build_py': build_py_and_copy_version},
    description = 'A collection of python code used by the SXS collaboration',
    author = 'SXS Collaboration',
    author_email = 'web-admin@black-holes.org',
    url = 'https://www.black-holes.org/',
    download_url = 'https://github.com/moble/sxs/archive/master.tar.gz',
)
