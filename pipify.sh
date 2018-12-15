/bin/rm -rf build dist sxs.egg-info

# # This is deprecated
# python setup.py bdist_wheel --universal upload

# This is the new way
pip install --quiet --upgrade twine
python setup.py sdist bdist_wheel --universal
twine upload dist/*
