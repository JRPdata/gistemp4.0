#! /usr/bin/env python
# step0.py
#
# Nick Barnes and David Jones.
# Copyright (C) Ravenbrook Limited, 2008-2010.
# Copyright (C) Climate Code Fonudation, 2010-2012.
# Avi Persin, Revision 2016-01-06
#
# BSD license, see license.txt
from settings import *

"""
Python code for the STEP0 part of the GISTEMP algorithm: combining
diverse inputs into a single dataset.
"""


def earthly(records):
    """
    `records` is a dict of records. records for stations that
    have no lat/lon metadata ("not on Earth") are removed. The
    input dict is not modified, a fresh sequence of (id, record)
    pairs dict is returned.
    """
    for id, record in records.items():
        station = record.station
        if (-90.0 <= station.lat <= 90.0) and (-180.0 <= station.lon <= 180.0):
            yield id, record
        else:
            print("%s has invalid latitude/longitude" % station.uid)


def append_scar():
    scar_data = open(INPUT_DIR + "antarc1.list").readlines() + open(INPUT_DIR + "antarc2.list").readlines() + open(
        INPUT_DIR + "antarc3.list").readlines()
    with open(INPUT_DIR + "v4.inv", "a") as file:
        for station in scar_data:
            id = station[:11].strip()
            name = station[12:41]
            lat = station[43:49].strip().ljust(8, '0')
            lon = station[49:57].strip().ljust(8, '0')
            file.write(id + " " + lat + "   " + lon + " -999.0 " + name + "\n")
    file.close()


def step0(input):
    """
    An iterator for Step 0.  Produces a stream of `giss_data.Series`
    instances.  *input* should be an instance that has an open()
    method.  input.open(x) is called for each data source x.
    (typically, this input object is made by the tool.io.step0_input()
    function).
    """
    if len(open(INPUT_DIR + "v4.inv", 'r').readlines()[0].split()) > 5:
        pass
    else:
        # append scar data to inv file
        # append_scar()

        # generate BI for inv file
        from tool import generate_brightness
        generate_brightness.run()

    # Read each data input into dictionary form.
    data = {}
    for source in input.sources:
        print("Load %s records" % source.upper())

        records = input.open(source)
        data[source] = dict((record.uid, record) for record in records)

    # Join all data sources together; and remove stations with
    # no lat/lon metadata.
    records = {}
    for source in input.sources:
        records.update(earthly(data[source]))

    # We sort here - as does GISTEMP - so that all the records for a
    # given 11-digit station ID are grouped together, ready for
    # combining in the next step.
    for _, record in sorted(records.items()):
        if record:
            yield record
