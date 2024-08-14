#! /usr/bin/env python
#
# parameters/obsolete.py
#
# Nick Barnes, Ravenbrook Limited, 2010-02-15
# Avi Persin, Revision 2016-01-06

"""Parameters controlling obsolete features of the GISTEMP algorithm.

Parameters controlling the standard GISTEMP algorithm, or cccgistemp
extensions to the GISTEMP algorithm, are in other parameter files.
"""

# STEP 1
#
# The following parameters control behaviour which used to be a part
# of GISTEMP step 1.

combine_records = False
"""Whether to attempt to combine station records.  The first part of
STEP 1 used to attempt to combine any non-overlapping records from the
same station into a combined station record, and then to combine
overlapping records.  This was dropped by GISTEMP v3 in 2011-12, when
GHCN-M 3 was adopted as the source data set: this sort of combination
is done in the preparation of GHCN-M 3.
"""

station_combine_min_overlap = 4
"""The minimum number of years of overlap, between a combined record
and a candidate record, to allow the candidate to be added into the
combined record."""

station_combine_bucket_radius = 10
"""Used when deciding whether to add a non-overlapping station record
to a combined record.  This number of years is considered, either side
of the centre of the potential new combined record.  If there are
enough valid years in each record (see *station_combine_min_mid_years*),
and the difference between the average anomaly of the combined record
and the average anomaly of the candidate record, over those years, is
less than the standard deviation of the combined record, then the
record is combined."""

station_combine_min_mid_years = 5
"""The minimum number of acceptable year records in the central
"bucket", for both the combined record and the candidate record, when
combining non-overlapping station records."""
