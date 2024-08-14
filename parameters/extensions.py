#! /usr/bin/env python
#
# parameters/extensions.py
#
# Nick Barnes, Ravenbrook Limited, 2010-02-15
# Avi Persin, Revision 2016-01-06

"""Parameters controlling cccgistemp extensions to the standard
GISTEMP algorithm.

Parameters controlling the standard GISTEMP algorithm, or obsolete
features of GISTEMP, are in other parameter files.
"""

data_sources = "ghcn"
"""
Data sources that are used for the analysis (space separated string).
'ghcn' is the Global Historical Climate Network (NOAA/NCDC), version 4,
which includes antarctic data.

These sources are no longer used by GISTEMP but should still work with
cccgistemp:

'ghcn.v2' is GHCN version 2;
'ushcn' is United Stated Historical Climate Network (NOAA), version 2;
"""

ocean_source = "ersstv5"
"""
Which (single) source is used for ocean data. The file
input/SBBX.xxxx is used (case insensitive search).
"""

element = ''
"""
Which meteorological element is analysed when a (GHCN-M v3
format) data file contains more than one.

When using ISTI files, set this to "TAVG", "TMIN", "TMAX" as
appropriate.
"""

augment_metadata = ''
"""
(In the usual analysis this parameter is empty) This parameter enables
additional metadata fields to be read from a file.  The format is
"name=colA,colB,colC" (with an arbitrary number of comma separated
columuns).  The file called *name* is opened; each row is a comma
separated sequence of field values, with the fields being 'colA',
'colB', and so on.  There must be exactly one column called 'uid'.
If a station with the same uid appears in the ordinary metadata
(usually sourced from the v2.inv file) then the extra fields are
associated with the station.
"""

work_file_format = "v4"
"""
The format of the intermediate files written to the 'work' directory:
'v3' for GHCN v3.
"""
