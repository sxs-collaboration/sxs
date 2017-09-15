#! /usr/bin/env python

from __future__ import division, print_function
from warnings import warn
from lxml.etree import parse, XMLSyntaxError, XMLParser

def validate_xml(source, parser=None, skip_entity_errors=True):
    if parser is None:
        parser = XMLParser()
        # parser._parser.UseForeignDTD(True)
    try:
        parse(source, parser)
    except XMLSyntaxError as e:
        e = str(e)
        if skip_entity_errors and e.startswith('Entity'):
            pass
        else:
            if isinstance(source, str):
                warning = '\n\tParsing failed for {0}:\n\t\t'.format(source) + str(e)
            else:
                warning = '\n\tParsing failed:\n\t\t' + str(e)
            warn(warning)
            return False
    return True

def validate_directory(directory='.', skip_entity_errors=True):
    import glob, os
    os.chdir(directory)
    total_file_count = 0
    bad_file_count = 0
    for file_name in glob.glob('*.xhtml'):
        total_file_count += 1
        if not validate_xml(file_name, skip_entity_errors=skip_entity_errors):
            bad_file_count += 1
    return total_file_count, bad_file_count

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", help="Directory containing all the files", default=".")
    parser.add_argument("--show_entity_errors", help="Show errors caused by HTML entities", default=False)
    args = parser.parse_args()
    skip = not args.show_entity_errors
    total_file_count, bad_file_count = validate_directory(directory=args.directory, skip_entity_errors=skip)
    print('Bad files / total files = {0} / {1} = {2}'.format(bad_file_count, total_file_count, bad_file_count/total_file_count))
