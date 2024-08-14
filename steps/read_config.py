#! /usr/bin/env python
#
# read_config.py
#
# Nick Barnes, Ravenbrook Limited, 2010-01-16
# Avi Persin, Revision 2016-01-06

from settings import *

"""
Python code to read the various config and station files used by
GISTEMP:
"""


def get_changes_dict():
    """Reads the file input/Ts.strange.v4.list.IN_full and returns a
    dict result.  Each line in that file begins with a 12-digit
    station ID - actually the tuple (country-code, WMO station,
    modifier, duplicate) - and ends with either yyyy/mm, specifying a
    month datum to omit or with xxxx-yyyy, specifying years to omit.
    xxxx can be 0, meaning from the beginning. yyyy can be 9999,
    meaning to the end.  The dict is a map from ID to
    ('month',yyyy,mm) or ('years',xxxx,yyyy).
    """

    dict = {}
    for line in open(INPUT_DIR + 'Ts.strange.v4SCAR.list.IN_full', 'r'):
        split_line = line.split()
        id = split_line[0]
        try:
            year1, year2 = map(int, split_line[-1].split("-"))
            val = ("years", year1, year2)
        except ValueError:
            year, month = map(int, split_line[-1].split("/"))
            val = ("month", year, month)
        dict[id] = dict.get(id, [])
        dict[id].append(val)
    return dict
