#!/usr/local/bin/python3.4
#
# gio.py
#
# Paul Ollis and David Jones, 2010-03-10
# Avi Persin, Revision 2016-01-06
#
# (was previously called giss_io.py, then io.py)

"""GISTEMP Input/Output.  Readers and writers for datafiles used by NASA
GISS GISTEMP.

Some of these file formats are peculiar to GISS, others are defined and
used by other bodies (such as NOAA's v2.mean format).
"""
# Clear Climate Code
import copy
import itertools
import math
import re
import struct
import csv

import warnings

import fort
import numpy as np

import parameters
from settings import *
from steps import giss_data

#: Integer code used to indicate missing data.
#:
#: This is units of 0.1 celsius. This code is only used when
#: reading or writing input and working files.
MISSING = 9999


# For all plausible integers, converting them from external to internal
# to external again will preserve their value.
def tenths_to_float(t):
    if t == MISSING:
        return giss_data.MISSING
    return t * 0.1


def internal_to_external(series, scale=0.1):
    """Convert a series of values to external representation by
    converting to integer tenths (or other scale).  Normally
    this is used to convert a series from degrees Celcius to tenths
    of a degree.

    :Param series:
        A list or iterable of floating point value; usually each value
        represents a temperature in Celsius.

    :Return:
        A new list of values (ints).

    """

    # Note: 1/0.1 == 10.0; 1/0.01 == 100.0 (in other words even though
    # 0.1 and 0.01 are not stored exactly, their reciprocal is exactly
    # an integer)
    scale = 1.0 / scale

    def toint(f):
        if abs(f - giss_data.MISSING) < 0.01:
            return MISSING
        return int(round(f * scale))

    return [toint(v) for v in series]


def convert_tenths_to_float(tenths_series):
    """The inverse of `internal_to_external`."""
    return [tenths_to_float(v) for v in tenths_series]


def open_or_uncompress(filename):
    """Opens the text file `filename` for reading.  If this fails then
    it attempts to find a compressed version of the file by appending
    '.gz' to the name and opening that (uncompressing it on
    the fly).

    """
    try:
        return open(filename)
    except IOError:
        # When none of filename, nor filename.gz exists we
        # want to pretend that the exception comes from the original
        # call to open, above.  Otherwise the user can be confused by
        # claims that "foo.gz" does not exist when they tried to open
        # "foo".  In order to maintain this pretence, we have to get
        # the exception info and save it. See
        # http://blog.ianbicking.org/2007/09/12/re-raising-exceptions/
        import sys

        exception = sys.exc_info()
        try:
            import gzip

            return gzip.open(filename + '.gz')
        except IOError:
            pass
        raise (exception[0], exception[1], exception[2])


class SubboxWriter(object):
    """Produces a GISTEMP SBBX (subbox); typically the output of
    step3 (and 4), and the input to step 5.
     
     This is a .npz file containing numpy arrays.
    """

    def __init__(self, file):
        self.file = open(file + '.npz', 'wb')
        self.meta = None
        self.buf_record = None
        self.result = []

    def _flush(self, record):
        # Bounding box; to be converted to integer hundredths.
        # Conventionally the 4 elements of the box are southern
        # latitude, northern latitude, western longitude, eastern
        # longitude (but the code doesn't care).
        box = record.box
        box = [int(round(x * 100)) for x in box]
        if record.station_months == 0:
            # Write as trimmed record.
            rec = [1, box[0], box[1], box[2], box[3], record.stations, record.station_months, record.d]
        else:
            rec = [len(record), box[0], box[1], box[2], box[3], record.stations, record.station_months, record.d]
        series = np.asarray(record.series)
        return [rec, series]

    def write(self, record):
        if self.meta is None:
            assert hasattr(record, "precipitation_flag"), "First record must be SubboxMetaData"
            rec = np.array([record.mo1, record.kq, record.mavg, record.monm, record.monm4, record.yrbeg,
                            record.missing_flag, record.precipitation_flag, record.title], dtype=object)
            self.meta = rec
        else:
            self.result.append(self._flush(record))

    def close(self):
        self.file.close()


class SubboxReader(object):
    """Reads GISS subbox files (SBBX).  These files are output by Step
    3, and consumed by Step 5.  Step 4 both reads and writes a subbox
    file.
    """

    def __init__(self, rawfile, bos='>', celltype=None):
        self.bos = bos
        self.f = fort.File(rawfile, bos=self.bos)
        rec = self.f.readline()
        (self.mo1, kq, mavg, monm, monm4, yrbeg, missing_flag,
         precipitation_flag, title) = struct.unpack(self.bos + '8i80s', rec)

        self.meta = giss_data.SubboxMetaData(self.mo1, kq, mavg, monm,
                                             monm4, yrbeg, missing_flag, precipitation_flag, title)

        assert self.meta.mavg == 6, "Only monthly averages supported"

        if celltype is None:
            if "sea" in title.lower().split():
                celltype = 'P'
            else:
                celltype = 'C'
        self.celltype = celltype

        # Synthesize a gridding radius by parsing it out of the title.
        radiusre = r'CR (\d+) *KM'
        m = re.search(radiusre, title.decode("utf-8"))
        if m:
            radius = int(m.group(1))
            self.meta.gridding_radius = radius

    def info(self):
        """Return a length 8 sequence corresponding to the INFO array
        record in the binary file.
        """
        m = self.meta
        return [self.mo1, m.kq, m.mavg, m.monm,
                m.monm4, m.yrbeg, m.missing_flag, m.precipitation_flag]

    def __iter__(self):
        yield self.meta

        rec = self.f.readline()
        while rec:
            mo1 = self.mo1
            fmt = "iiiiiiif%df" % mo1
            fields = list(struct.unpack(self.bos + fmt, rec))
            series = fields[8:]
            self.mo1 = fields[0]
            # Make an attributes dictionary.
            # The box boundaries are fields[1:5], but we need to scale
            # them to fractional degrees first:
            for i in range(1, 5):
                fields[i] /= 100.0
            attr = dict(zip(
                ['lat_S', 'lat_N', 'lon_W', 'lon_E',
                 'stations', 'station_months', 'd'],
                fields[1:8]))
            attr['box'] = fields[1:5]
            subbox = giss_data.Series(series=series,
                                      celltype=self.celltype, **attr)
            rec = self.f.readline()
            yield subbox

    def __getattr__(self, name):
        return getattr(self.meta, name)


class SubboxReaderNpz(object):
    """Reads GISS subbox files (SBBX).  These files are output by Step
    3, and consumed by Step 5.  Step 4 both reads and writes a subbox
    file.
    
    Reads a .npz file consisting of numpy arrays.
    """

    def __init__(self, file, celltype=None):
        self.f = np.load(file + '.npz')
        meta = self.f['meta']
        title = meta[-1]
        self.meta = giss_data.SubboxMetaData(*meta)
        self.f.files.remove('meta')
        self.f.files = sorted(self.f.files, key=lambda t: int(t.split('_')[1]))
        assert self.meta.mavg == 6, "Only monthly averages supported"
        if type(title) is bytes:
            title = title.decode()

        if celltype is None:
            if "sea" in title.lower().split():
                celltype = 'P'
            else:
                celltype = 'C'
        self.celltype = celltype
        # Synthesize a gridding radius by parsing it out of the title.

        import parameters
        radius = getattr(parameters, "gridding_radius")
        self.meta.gridding_radius = radius

    def info(self):
        """Return a length 8 sequence corresponding to the INFO array
        record in the binary file.
        """
        m = self.meta
        return [self.mo1, m.kq, m.mavg, m.monm,
                m.monm4, m.yrbeg, m.missing_flag, m.precipitation_flag]

    def __iter__(self):
        yield self.meta
        for rec in self.f:
            rec = self.f[rec]
            fields = rec[0]
            series = rec[1]
            # Make an attributes dictionary.
            # The box boundaries are fields[1:5], but we need to scale
            # them to fractional degrees first:
            for i in range(1, 5):
                fields[i] /= 100.0
            attr = dict(zip(['lat_S', 'lat_N', 'lon_W', 'lon_E', 'stations', 'station_months', 'd'], fields[1:8]))
            attr['box'] = fields[1:5]
            subbox = giss_data.Series(series=series, celltype=self.celltype, **attr)
            yield subbox

    def __getattr__(self, name):
        return getattr(self.meta, name)


def GHCNV4Reader(path=None, file=None, meta=None,
                 year_min=None, scale=None, element=None):
    """Reads a file in GHCN V4 .dat format and yields each station
    record (as a giss_data.Series instance).  For now, this treats
    all the data for a station as a single record (contrast with GHCN V2
    which could have several "duplicates" for a single station).

    If a *meta* dict is supplied then the Series instance will have its
    "station" attribute set to value corresponding to the 11-digit ID in
    the *meta* dict.

    If `year_min` is specified, then only years >= year_min are kept
    (the default, None, keeps all years).

    If *scale* is specified then the (integer) values in the file are
    multiplied by *scale* before being returned.  When it is not
    specified (the normal case), the scale is derived from the element
    specified in the file (normally for monthly means this is "TAVG" and
    the scale implied by that is 0.01 (degrees C)).

    The default behaviour for *element* (when None) is to assume
    that the input file contains exactly one sort of element and the
    values for this element are returned. In the usual analysis this
    is a GHCN-M v4 file containing TAVG. If the input file
    contains more than one sort of element (such as an ISTI
    file) then an exception will be raised. In order to read
    from a file containing more than one sort of element,
    specify the element with this argument. For example,
    element='TMIN'.

    See ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v4/readme.txt for format
    of this file.
    """

    if path:
        inp = open(path)
    else:
        inp = file

    def id11(l):
        """Extract the 11-digit station identifier."""
        return l[:11]

    noted_element = {}

    def note_element(element):
        """
        Print the meteorological element we are reading (the
        first time we see it). Also, abort if we see more than
        one sort of element (for example, we try reading an ISTI
        file that has TAVG, TMIN, TMAX, without specifying which
        one we want).
        """

        if element in noted_element:
            return
        noted_element[element] = True

        if len(noted_element) > 1:
            raise Exception("File contains more than one sort of element: %r" % noted_element.keys())

        friendly = dict(TAVG='average temperature',
                        TMIN='mean minimum temperature',
                        TMAX='mean maximum temperature')
        print("(Reading %s)" % friendly[element])

    element_scale = dict(TAVG=0.01, TMIN=0.01, TMAX=0.01)
    reject = 'DKOSTW'

    def convert(s):
        """Convert single value. *s* is the 8 character string: 5
        characters for value, 3 for flags."""

        # This function captures *multiplier* which can, in principle,
        # change for each line.

        v = int(s[:5])
        # Flags for Measurement (missing days), Quality, and
        # Source.
        m, q, s = s[5:8]
        if q in reject or v == -9999:
            v = MISSING
        else:
            v *= multiplier
        return v

    all_missing = [MISSING] * 12

    for id, lines in itertools.groupby(inp, id11):
        key = dict(uid=id, first_year=year_min)
        if meta and meta.get(id):
            key['station'] = meta[id]
        record = giss_data.Series(**key)
        for line in lines:
            year = int(line[11:15])
            found_element = line[15:19]
            if element and element != found_element:
                continue
            note_element(found_element)
            if scale:
                multiplier = scale
            else:
                multiplier = element_scale[found_element]
            values = [convert(line[i:i + 8]) for i in range(19, 115, 8)]
            if values != all_missing:
                record.add_year(year, values)
        if len(record) != 0:
            yield record


class GHCNV3Writer(object):
    """Write a file in GHCN v3 format. See also GHCNV4Reader.  The
    format is documented in
    ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v3/README .  If the records
    have an 'element' property, then that is used for the 'element'
    field in the GHCN V3 file, otherwise 'TAVG' is used.
    """

    def __init__(self, path=None, file=None, scale=0.01, **k):
        if path is not None:
            self.f = open(path, "w")
        else:
            self.f = file
        self.scale = scale

    def to_text(self, t):
        if t == MISSING:
            return "-9999"
        else:
            return "%5d" % t

    def write(self, record):
        """Write an entire record out."""
        for year in range(record.first_year, record.last_year + 1):
            if not record.has_data_for_year(year):
                continue
            element = getattr(record, 'element', 'TAVG')
            self.writeyear(record.uid, element, year, record.get_a_year(year))

    def writeyear(self, uid, element, year, temps):

        """Write a single year's worth of data out.  *temps* should
        contain 12 monthly values."""

        if len(uid) > 11:
            # Convert GHCN v2 style identifier into 11-digit v3
            # identifier; use 12th digit for the source flag.
            uid = uid[:12]
            assert len(uid) == 12
            sflag = uid[11]
        elif len(uid) == 6:
            # Assume it's a 6 digit identifier from USHCN.
            uid = '42500' + uid
            sflag = 'U'
        else:
            sflag = ' '
        id11 = "%-11.11s" % uid
        assert len(element) == 4

        tstrings = [self.to_text(t)
                    for t in internal_to_external(temps, scale=self.scale)]
        flags = ['  ' + sflag] * 12

        self.f.write('%s%04d%s%s\n' % (uid, year, element,
                                       ''.join(t + flag for t, flag in zip(tstrings, flags))))

    def close(self):
        self.f.close()


antarc_discard_re = re.compile(r'^$|^Get |^[12A-Z]$')
antarc_temperature_re = re.compile(r'^(.*) .* *temperature')


def read_antarctic(path, station_path, discriminator,
                   meta=None, year_min=None):
    stations = read_antarc_station_ids(station_path, discriminator)
    record = None
    for line in open(path):
        if antarc_discard_re.search(line):
            continue
        station_line = antarc_temperature_re.match(line)
        if station_line:
            station_name = station_line.group(1)
            station_name = station_name.replace('\\', '')
            id12 = stations[station_name]
            if record is not None:
                yield record
            key = dict(uid=id12, first_year=year_min)
            id11 = id12[:11]
            if meta and meta.get(id11):
                key['station'] = meta[id11]
            record = giss_data.Series(**key)
            continue
        line = line.strip()
        if line.find('.') >= 0 and line[0] in '12':
            year, data = read_antarc_line(line)
            if year >= giss_data.BASE_YEAR:
                record.add_year(year, data)
    if record is not None:
        yield record


austral_discard_re = re.compile(r'^$|:')
austral_header_re = re.compile(r'^\s*(.+?)  .*(E$|E )')


def read_australia(path, station_path, discriminator,
                   meta=None, year_min=None):
    stations = read_antarc_station_ids(station_path, discriminator)
    record = None
    for line in open(path):
        if austral_discard_re.search(line):
            continue
        station_line = austral_header_re.match(line)
        if station_line:
            station_name = station_line.group(1).strip()
            id12 = stations[station_name]
            if record is not None:
                yield record
            key = dict(uid=id12, first_year=year_min)
            id11 = id12[:11]
            if meta and meta.get(id11):
                key['station'] = meta[id11]
            record = giss_data.Series(**key)
            continue
        line = line.strip()
        if line.find('.') >= 0 and line[0] in '12':
            year, data = read_antarc_line(line)
            if year >= giss_data.BASE_YEAR:
                record.add_year(year, data)

    if record is not None:
        yield record


def read_antarc_line(line):
    """Convert a single line from the Antarctic/Australasian dataset
    files into a year and a 12-tuple of floats (the temperatures in
    Centigrade).
    """

    year = int(line[:4])
    line = line[4:]
    tuple = []
    if line[6] == '.' or line[7] == '-':
        # Some of the datasets are 12f8.1 with missing values as '       -'.
        for i in range(0, 12):
            tuple.append(read_float(line[i * 8:i * 8 + 8]))
    else:
        # Others are xx12f7.1 or xxx12f7.1 with missing values as '       '.
        np = line.find('.')
        if np < 0:
            raise (ValueError, "Non-data line encountered: '%s'" % line)
        position = (np % 7) + 2
        for i in range(0, 12):
            tuple.append(read_float(line[i * 7 + position:i * 7 + 7 + position]))
    return year, tuple


def read_antarc_station_ids(path, discriminator):
    """Reads a SCARs station ID files and returns a dictionary
    mapping station name to the 12-digit station ID.
    """
    dict = {}
    for line in open(path):
        id11 = line[:11]
        station = line[12:42].strip()
        dict[station] = id11 + discriminator
    return dict


def station_metadata(path=None, file=None, format='giss_v4'):
    """Read station metadata from file, return it as a dictionary.
    *format* specifies the format of the metadata can be:
    'giss_v4' for GHCN v4 ;

    Here are two typical lines, with a record diagram

    id---------xname--------------------------xlat---xlon----x1---2----34----5-6-7-8-910grveg-----------GU--11
    0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
    40371148001 ALMASIPPI,MA                    49.55  -98.20  274  287R   -9FLxxno-9x-9COOL FIELD/WOODSA1   0
    42572530000 CHICAGO/O'HARE, ILLINOIS        42.00  -87.90  205  197U 6216FLxxno-9A 1COOL CROPS      C3 125

       uid                 40371148001          42572530000
          The unique ID of the station. This is held as an 11 digit string.
       name                ALMASIPPI,MA         CHICAGO/O'HARE, ILLINOIS
        The station's name.
       lat                 49.55                42.00
        The latitude, in degrees (two decimal places).
       lon                 -98.20               -87.90
        The longitude, in degrees (two decimal places).
    1  stelev              274                  205
        The station elevation in metres.
    2  grelev              287                  197
        The grid elevation in metres (value taken from gridded dataset).
    3  popcls              R                    U
        'R' for rural,  'S' for semi-urban, 'U' for urban.
    4  popsiz              -9                   6216
        Population of town in thousands.
    5  topo                FL                   FL
        The topography.
    6  stveg               xx                   xx
    7  stloc               no                   no
        Whether the station is near a lake (LA) or ocean (OC).
    8  ocndis              -9                   -9
    9  airstn              x                    A
    10 towndis             -9                   1
       grveg               COOL FIELD/WOODS     COOL CROPS
        An indication of vegetation, from a gridded dataset. For example,
        'TROPICAL DRY FOR'.
    G  popcss              A                    C
        Population class based on satellite lights (GHCN value).
    U  us_light            1                    3
        Urban/Rural flag based on satellite lights for US stations
        (' ' for non-US stations).  '1' is dark, '3' is bright.
    11 global_light        0                    125
    Global satellite nighttime light value.  Range 0-186 (at
    least).
    """
    # Do not supply both arguments!
    assert not (file and path)
    assert format in ('giss_v3', 'v3', 'giss_v4')
    if path:
        try:
            file = open(path)
        except IOError:
            warnings.warn("Could not load %s metadata file: %s" % (format, path))
            return {}
    assert file

    def blank_int(s):
        """
        Convert a field to int, or if it is blank, convert to None.
        """

        if s == '' or s.isspace():
            return None
        return int(s)

    # Fields are named after the designators used in the GHCN v4
    # documentation
    # except for:
    # uid (GHCN: ID), lat (GHCN: latitude), lon (GHCN: longitude),
    # us_light (GISTEMP specific field for nighttime satellite
    # brightness over the US, see Hansen et al 2001), global_light
    # (GISTEMP specific field for global nighttime satellite
    # brightness).


    v3fields = dict(
        uid=(0, 11, str),
        name=(12, 42, str),
        lat=(43, 49, float),
        lon=(50, 57, float),
        stelev=(58, 62, int),
        grelev=(62, 67, blank_int),
        popcls=(67, 68, str),
        popsiz=(68, 73, blank_int),
        topo=(73, 75, str),
        stveg=(75, 77, str),
        stloc=(77, 79, str),
        ocndis=(79, 81, blank_int),
        airstn=(81, 82, str),
        towndis=(82, 84, blank_int),
        grveg=(84, 100, str),
        popcss=(100, 101, str),
        global_light=(101, 106, blank_int),  # GISTEMP only
        berkeley=(106, 109, str),  # GISTEMP only; comment suggests derived from Berkeley Earth.
    )

    # See ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v4/readme.txt for format
    # of GHCN's original metadata file.
    v3_ghcn_fields = dict(
        uid=(0, 11, str),
        lat=(12, 20, float),
        lon=(21, 30, float),
        stelev=(31, 37, float),
        name=(38, 68, str),
        grelev=(69, 73, blank_int),
        popcls=(73, 74, str),
        popsiz=(75, 79, blank_int),
        topo=(79, 81, str),
        stveg=(81, 83, str),
        stloc=(83, 85, str),
        ocndis=(85, 87, blank_int),
        airstn=(87, 88, str),
        towndis=(88, 90, blank_int),
        grveg=(90, 106, str),
        popcss=(106, 107, str),
    )
    v4fields = dict(
        uid=(0, 11, str),
        lat=(12, 20, float),
        lon=(21, 30, float),
        stelev=(31, 37, float),
        name=(37, 69, str),
        global_light=(69, 74, blank_int), # was (69, 71, blank_int)
    )

    if 'giss_v3' == format:
        fields = v3fields
    elif 'giss_v4' == format:
        fields = v4fields
    elif 'v3' == format:
        fields = v3_ghcn_fields

    result = {}
    for line in file:
        d = dict((field, convert(line[a:b])) for field, (a, b, convert) in fields.items())
        result[d['uid']] = giss_data.Station(**d)
    return result


def augmented_station_metadata(path=None, file=None, format='v3'):
    """Reads station metadata just like station_metadata() but
    additionally augments records with metadata obtained from another
    file, specified by parameters.augment_metadata.
    """
    meta = station_metadata(path=path, file=file, format=format)
    augments = parameters.augment_metadata

    if augments:
        path, columns = augments.split('=')
        columns = columns.split(',')
        assert 'uid' in columns
        for row in open(path):
            row = row.strip().split(',')
            d = dict(zip(columns, row))
            # Convert things that look like numbers, to numbers.
            # (except for uid, which is always a string)
            for k, v in d.items():
                if k == 'uid':
                    continue
                try:
                    v = float(v)
                except ValueError:
                    pass
                d[k] = v
            uid = d['uid']
            if uid in meta:
                meta[uid].__dict__.update(d)
    return meta


def read_generic_v3(name):
    """
    Reads a "generic" source in GHCN-M v3 format.

    *name* should be the name of the .dat file, in the input/
    directory.  The .inv file found alongside will be used for
    metadata.
    """

    filename = os.path.join(BASE_PATH + "input", name)
    if not name.endswith('.dat'):
        raise Exception("Don't know where to look for .inv file for GHCN-M v2 format file: %r" % name)

    if parameters.element:
        element = parameters.element
    else:
        element = None

    invfile = name[:-4] + '.inv'
    invfile = os.path.join(BASE_PATH + "input", invfile)
    return GHCNV4Reader(file=open(filename),
                        meta=augmented_station_metadata(invfile, format='v3'),
                        year_min=giss_data.BASE_YEAR,
                        element=element)


def read_float(s):
    """Returns the float converted from the argument string.
    If float conversion fails, returns MISSING.
    """

    try:
        return float(s)
    except:
        return giss_data.MISSING


_v3meta = None


def v3meta():
    """Return the GHCN v3 metadata.  Loading it (from the modified
    version of the file supplied by GISS) if necessary.
    """

    # It's important that this file be opened lazily, and not at module
    # load time (if "input/" hasn't been populated yet, it won't be
    # found).

    global _v3meta

    v3inv = os.path.join(INPUT_DIR, 'v4.inv')
    if not _v3meta:
        _v3meta = augmented_station_metadata(v3inv, format='giss_v4')
    return _v3meta


def maskboxes(inp, grid):
    """Read a step5mask file box by box.  Yield (value, box) pair.
    """
    for row, box in zip(inp, grid):
        lat = float(row[:5])
        lon = float(row[5:11])
        s, n, w, e = box
        # If either of these fail, the input mask is in wrong sequence.
        assert s < lat < n
        assert w < lon < e
        v = float(row[16:21])
        yield v, box


class Input:
    """Generally one instance is created: the result of
    step0_input()."""

    def __init__(self):
        self.sources = parameters.data_sources.split()

    def open(self, source):
        """Open the source (specified as a string), and return an
        iterator."""
        if source == 'ghcn' or re.match('ghcnm.(tavg|tmax|tmin)', source):
            if source == 'ghcn':
                ghcn4file = INPUT_DIR + 'ghcnm.tavg.qcf.dat'
            else:
                if source.endswith('.dat'):
                    pass
                else:
                    source += '.qca.dat'
                ghcn4file = os.path.join(INPUT_DIR, source)
            invfile = INPUT_DIR + 'v4.inv'

            return GHCNV4Reader(file=open(ghcn4file),
                                meta=augmented_station_metadata(invfile, format='giss_v4'),
                                year_min=giss_data.BASE_YEAR)
        if source == 'scar':
            return itertools.chain(
                read_antarctic(INPUT_DIR + "antarc1.txt", INPUT_DIR + "antarc1.list", '8',
                               meta=v3meta(), year_min=giss_data.BASE_YEAR),
                read_antarctic(INPUT_DIR + "antarc3.txt", INPUT_DIR + "antarc3.list", '9',
                               meta=v3meta(), year_min=giss_data.BASE_YEAR),
                read_australia(INPUT_DIR + "antarc2.txt", INPUT_DIR + "antarc2.list", '7',
                               meta=v3meta(), year_min=giss_data.BASE_YEAR))
        if source.endswith('.dat'):
            return read_generic_v3(source)
        raise Exception("Cannot open source %r" % source)


# Each of the stepN_input functions below produces an iterator that
# yields that data for that step feeding from data files.
# Each of the stepN_output functions below is effectively a "tee" that
# writes the data to a file; they each take a data object (an
# iterator), write each item to a file, and yield each item.
def step0_input():
    input = Input()
    return input


def choose_writer():
    """Choose a record writer function, according to
    parameters.work_file_format, and return (function,filext) pair."""

    format = parameters.work_file_format

    if format == 'v3' or format == 'v4':
        writer = GHCNV3Writer
    return writer, format


def generic_output_step(n):
    """Return a generic output routine for step *n*."""

    def output(data):
        writer, ext = choose_writer()
        path = os.path.join(WORK_DIR, 'step%d.%s' % (n, ext))
        out = writer(path=path)
        for thing in data:
            out.write(thing)
            yield thing
        print("Step %d: closing output file." % n)
        out.close()
        progress = open(PROGRESS_DIR + 'progress.txt', 'a')
        progress.write("\nStep %d: closing output file.\n" % n)

    return output


step0_output = generic_output_step(0)


def step1_input():
    return GHCNV4Reader(WORK_DIR + "step0.v4",
                        meta=v3meta(),
                        year_min=giss_data.BASE_YEAR)


step1_output = generic_output_step(1)


def step2_input():
    return GHCNV4Reader(WORK_DIR + "step1.v4", meta=v3meta())


step2_output = generic_output_step(2)


def step3_input():
    return GHCNV4Reader(WORK_DIR + "step2.v4", meta=v3meta())


STEP3_OUT = os.path.join(RESULT_DIR, 'SBBX1880.Ts.GHCN.CL.PA.1200')


def step3_output(data):
    out = SubboxWriter(STEP3_OUT)
    writer, ext = choose_writer()
    textout = writer(path=(WORK_DIR + 'step3.%s' % ext), scale=0.01)
    gotmeta = False
    for thing in data:
        out.write(thing)
        if gotmeta:
            textout.write(thing)
        gotmeta = True
        yield thing
    np.savez_compressed(out.file, *out.result, meta=out.meta)
    print("Step 3: closing output file")
    out.close()
    textout.close()
    progress = open(PROGRESS_DIR + 'progress.txt', 'a')
    progress.write("\nStep3: closing output file\n")


def step3c_input():
    """Use the output from the ordinary Step 3."""

    land = SubboxReaderNpz(STEP3_OUT)
    return iter(land)


def make_3d_array(a, b, c):
    """Create an array with three dimensions.

    Actually a list-of-lists-of-lists, but the result can be treated as an
    array with dimensions ``[a][b][c]``.

    """
    x = [0.0] * c
    arr = []
    for i in range(a):
        arr.append([list(x) for j in range(b)])

    return arr


def step4_find_monthlies(latest_year, latest_month):
    dates = {}
    filename_re = re.compile('^oiv2mon\.([0-9][0-9][0-9][0-9])([0-9][0-9])(\.gz)?$')
    for f in os.listdir(INPUT_DIR):
        m = filename_re.match(f)
        if m:
            year = int(m.group(1))
            month = int(m.group(2))
            if (year, month) > (latest_year, latest_month):
                if m.group(3):
                    f = f[:-3]
                dates[(year, month)] = 'input/' + f
    l = list(dates.items())
    l.sort()
    return l


def step4_load_sst_monthlies(latest_year, latest_month):
    files = step4_find_monthlies(latest_year, latest_month)
    if not files:
        print("No more recent sea-surface data files.\n")
        return None

    first_year = files[0][0][0]
    last_year = files[-1][0][0]
    n_years = last_year - first_year + 1

    # Read in the SST data for recent years
    sst = make_3d_array(360, 180, 12 * n_years)

    dates = []
    for (date, file) in files:
        dates.append(date)
        (year, month) = date
        f = open_or_uncompress(file)
        print("reading", file)
        f = fort.File(f, bos=">")
        f.readline()  # discard first record
        data = f.readline()
        f.close()
        month = 12 * (year - first_year) + month - 1
        p = 0
        for lat in range(180):
            for long in range(360):
                v, = struct.unpack(">f", data[p:p + 4])
                p += 4
                sst[long][lat][month] = v

    return sst, dates


# This is used to extract the end month/year from the title of the
# SBBX file. This file can either be the usual ERSST file, or a
# variety of alternatives, including the traditional HadR2 file.
alternatives = ("(?:" +
                "|".join([
                    "Had: 1880-11/1981, oi2: 12/1981-",  # traditional
                    "ERSST 1880-1981  OISSTadj 1982-",
                    "HadISST 01/1880 -",
                    "ERSST 01/1880 -",  # since 2013, the usual analysis
                    "ERSSTv4 01/1880 -",  # since 08/2015, the usual analysis
                    "ERSSTv5 01/1880 -",  # for ERSST v5
                ]) + ")")
rTitle = re.compile(r"Monthly Sea Surface Temperature anom \(C\) " +
                    alternatives +
                    " *(\d+)/(\d+)")


def find_ocean_file():
    """
    From parameters.ocean_source, find the filename for ocean
    data and return it.
    """

    source = parameters.ocean_source.upper()
    dir = INPUT_DIR
    for name in os.listdir(dir):
        if name.upper() == "SBBX." + source:
            return os.path.join(dir, name)

    raise Exception("Can't find ocean source %r." %
                    parameters.ocean_source)


def step4_input(land):
    # The "land is None" check allows Step 4 to be run on its
    # own, loading the land data from work files in that case.
    if land is None:
        land = SubboxReaderNpz(STEP3_OUT)
    ocean_file = find_ocean_file()
    ocean = SubboxReader(open(ocean_file, 'rb'))
    ocean.meta.ocean_source = parameters.ocean_source

    m = rTitle.match(ocean.meta.title.decode("utf-8"))
    if m is None:
        print("The title in %s does not look right\n" % ocean_file)
        print("Unable to determine end month/year from:\n")
        print("  %r\n" % ocean.meta.title)
        sys.exit(1)
    end_month = int(m.group(1))
    end_year = int(m.group(2))
    monthlies = step4_load_sst_monthlies(end_year, end_month)
    return land, ocean, monthlies


def step4_output(data):
    # The Step 4 output is slightly unusual, it is an iterable of pairs.
    # We only want to write the records from the right-hand item (the
    # ocean data).  The left-hand items are land data, already written
    # by Step 3.
    out = SubboxWriter(RESULT_DIR + "SBBX.SST")
    for land, ocean in data:
        out.write(ocean)
        yield land, ocean
    np.savez_compressed(out.file, *out.result, meta=out.meta)
    print("Step4: closing output file")
    out.close()
    progress = open(PROGRESS_DIR + 'progress.txt', 'a')
    progress.write("\nStep4: closing output file\n")


def step5_input(data):
    if not data:
        land = SubboxReaderNpz(STEP3_OUT)
        try:
            ocean = SubboxReaderNpz(RESULT_DIR + 'SBBX.SST')
            ocean.meta.ocean_source = parameters.ocean_source
        except IOError:
            data = ensure_landocean(iter(land))
        else:
            data = zip(land, ocean)
    else:
        data = ensure_landocean(data)

    # Add optional mask.
    try:
        p = os.path.join(BASE_PATH + 'input', 'step5mask')
        mask = open(p)
        print("Using mask from", p)
    except IOError:
        mask = None
    meta = next(data)
    if mask is None:
        yield (None,) + tuple(meta)
        for land, ocean in data:
            yield None, land, ocean
    else:
        yield ('mask from %s' % mask.name,) + tuple(meta)
        for maskrow, (land, ocean) in zip(mask, data):
            maskv = float(maskrow[16:21])
            yield maskv, land, ocean


def ensure_landocean(data):
    """Ensure that the data stream has a land and an ocean series.  If
    we are piping Step 3 straight into Step 5 then we only have a land
    series.  In that case we synthesize missing ocean data.
    """

    # First item from iterator is normally a pair of metadata objects,
    # one for land, one for ocean.  If we are piping step3 straight into
    # step5 then it is not a pair.

    meta = next(data)
    try:
        land_meta, ocean_meta = meta
    except (TypeError, ValueError):
        # Use the land meta object for both land and ocean data
        land_meta = meta
        ocean_meta = copy.copy(meta)
        print("No ocean data; using land data only")
        data = add_blank(data, 'ocean')

    if land_meta is None:
        # Synthesize land data
        land_meta = copy.copy(ocean_meta)
        print("No land data; using ocean data only")
        data = add_blank(extract_ocean(data), 'land')

    yield land_meta, ocean_meta
    for series in data:
        yield series


def extract_ocean(data):
    for land, ocean in data:
        yield ocean


def add_blank(data, required):
    """Augment a single data series with blank data to make a data
    series pair.  *required* should be 'land' to synthesize the first of
    the pair; or 'ocean' to synthesize the second of the pair.
    """

    assert required in ('land', 'ocean')

    for this_box in data:
        other_series = [MISSING] * len(this_box.series)
        other_box = giss_data.Series(series=other_series,
                                     box=this_box.box,
                                     stations=0, station_months=0,
                                     d=MISSING)
        if required == 'land':
            yield other_box, this_box
        else:
            yield this_box, other_box


def step5_bx_output(meta, data):
    """Write box (BX) output file."""
    title = meta.title
    # Usually one of 'land', 'ocean', 'mixed'.
    mode = meta.mode
    result = []
    boxf = open(os.path.join(RESULT_DIR, make_filename(meta, 'BX') + '.npz'), 'wb')
    info = info_from_meta(meta)
    info.append(title)
    info = np.array(info, dtype=object)

    for record in data:
        avgr, wtr, ngood, box = record
        rec = [np.asarray(avgr), np.asarray(wtr), [ngood, box]]
        result.append(rec)
        yield record

    np.savez_compressed(boxf, *result, meta=info)
    print("Step 5: Closing box file:", boxf.name)
    boxf.close()
    progress = open(PROGRESS_DIR + 'progress.txt', 'a')
    progress.write("\nStep 5: Closing box file:" + boxf.name + '\n')


def make_filename(meta, kind):
    """Using the metadata in *meta* make a filename for an output file
    of type *kind*.  *kind* is usually one of 'ZON' or 'BX'.

    A typical filename is
      result/mixedBX.Ts.ERSST.GHCN.CL.PA.1200
    """

    if hasattr(meta, 'gridding_radius'):
        radius = ".%.0f" % meta.gridding_radius
    else:
        # No gridding radius specified in *meta*; typically, an ocean
        # file.
        radius = ''

    mode = meta.mode

    return (
        meta.mode + kind + '.Ts' + name_sources(meta, mode) +
        '.CL.PA' + radius)


def name_sources(meta, mode):
    """
    From the meta data and the mode, create a filename fragment
    that names the sources used.
    """

    land_source = ''
    ocean_source = ''

    if mode in ('land', 'mixed'):
        land_source = '.GHCN'
    if mode in ('ocean', 'mixed'):
        ocean_source = '.' + meta.ocean_source.upper()

    return ocean_source + land_source


def make_text_filename(meta, mode, part):
    return mode + part + '.Ts' + name_sources(meta, mode) + '.CL.PA.txt'


def info_from_meta(meta):
    """Take a metadata object (any object with certain fields) and
    return an 8 element "info" list; the 8 elements form the header of
    various GISS specific Fortran binary files.
    """

    return [meta.mo1, meta.kq, meta.mavg, meta.monm,
            2 * meta.monm + 5, meta.yrbeg, meta.missing_flag,
            meta.precipitation_flag]


def step5_mask_output(data):
    """Output the landmask used by Step 5."""
    # metadata
    yield next(data)

    out = open(os.path.join(WORK_DIR, 'step5mask'), 'w')

    for datum in data:
        mask, land, ocean = datum
        assert boxid_eq(land.uid, ocean.uid)
        out.write("%sMASK%.3f\n" % (land.uid, mask))
        yield datum
    out.close()


def boxid_eq(uid1, uid2):
    """Compare two IDs both of which are from subboxes.  They should be
    of the form -SS.S+EEE.EC.  They should be the same, but due to
    rounding differences, they can differ by 1 in the last digit."""

    lat1 = float(uid1[:5])
    lon1 = float(uid1[5:11])
    lat2 = float(uid2[:5])
    lon2 = float(uid2[5:11])
    return abs(lat1 - lat2) < 0.15 and abs(lon1 - lon2) < 0.15


def step5_zone_titles():
    """Return the titles used for the 16 zones."""

    # Boundaries (degrees latitude, +ve North) of the 8 basic belts.
    band = ['90.0 N',
            '64.2 N',
            '44.4 N',
            '23.6 N',
            'EQUATOR',
            '23.6 S',
            '44.4 S',
            '64.2 S',
            '90.0 S']
    # Accumulate the titles here.
    titles = ['  LATITUDE BELT FROM %7s TO %7s' % (band[j + 1], band[j])
              for j in range(8)]

    titles += [
        '  NORTHERN LATITUDES: 23.6 N TO  90.0 N',
        '       LOW LATITUDES: 23.6 S TO  23.6 N',
        '  SOUTHERN LATITUDES: 90.0 S TO  23.6 S',
        '  NHEM.MID LATITUDES: 23.6 N TO  64.2 N',
        '  SHEM.MID LATITUDES: 64.2 S TO  23.6 S',
        'NORTHERN HEMISPHERE',
        'SOUTHERN HEMISPHERE',
        'GLOBAL MEANS']

    # Ensure all titles are 80 characters long.
    titles = list(map(lambda s: ('%-80s' % s)[:80], titles))
    return titles


def step5_output(results):
    """Generate final Step 5 output files.  *results* is a sequence of
    tuples, each tuples corresponding to the zonal results for
    an entire analysis (typically 3 analyses: land, ocean, mixed).  The
    contents of the tuple itself are a bit baroque, see `annzon` for
    details.
    """
    for item in results:
        step5_output_one(item)
    to_csv(['mixedGLB.Ts.ERSSTV5.GHCN.CL.PA.txt', 'mixedNH.Ts.ERSSTV5.GHCN.CL.PA.txt',
            'mixedSH.Ts.ERSSTV5.GHCN.CL.PA.txt', 'mixedZonAnn.Ts.ERSSTV5.GHCN.CL.PA.txt',
            'landGLB.Ts.GHCN.CL.PA.txt', 'landNH.Ts.GHCN.CL.PA.txt', 'landSH.Ts.GHCN.CL.PA.txt',
            'landZonAnn.Ts.GHCN.CL.PA.txt'])
    return "Step 5 Completed"


def to_csv(filenames=''):
    for filename in filenames:
        create_csv(RESULT_DIR + filename)


def create_csv(dir_name, txt_title=""):
    if "/graph.txt" in dir_name:
        file = open(dir_name.replace('/graph.txt', '/graph.csv'), 'w')
    else:
        file = open(dir_name.replace('.txt', '.csv'), 'w')
    wr = csv.writer(file, lineterminator='\n')
    if txt_title == "":
        txt_title = set_display_name(dir_name.split('/')[-1][:-4])
    wr.writerow([txt_title])
    for line in open(dir_name).readlines():
        if re.match(r'\d{4}', line) or line[:4] == "Year" or line[:4] == "Year":
            csv_line = line.split()
            if "|" in csv_line:
                csv_line.remove("|")
            if csv_line == ['Year', 'Glob', 'NHem', 'SHem', '-90N', '-24N', '-24S', '-90N', '-64N', '-44N', '-24N',
                            '-EQU', '-24S', '-44S', '-64S']:
                file.seek(0)
                csv_line = ['Year', 'Glob', 'NHem', 'SHem', '24N-90N', '24S-24N', '90S-24S', '64N-90N', '44N-64N',
                            '24N-44N', 'EQU-24N', '24S-EQU', '44S-24S', '64S-44S', '90S-64S']
            if len(csv_line) == 16:
                csv_line = csv_line[:-1]
            wr.writerow(csv_line)


def set_display_name(filename):
    display_names = {"landGLB.Ts.GHCN.CL.PA": "Station: Global Means",
                     "landNH.Ts.GHCN.CL.PA": "Station: Northern Hemispheric Means",
                     "landSH.Ts.GHCN.CL.PA": "Station: Southern Hemispheric Means",
                     "mixedGLB.Ts.ERSSTV5.GHCN.CL.PA": "Land-Ocean: Global Means",
                     "mixedNH.Ts.ERSSTV5.GHCN.CL.PA": "Land-Ocean: Northern Hemispheric Means",
                     "mixedSH.Ts.ERSSTV5.GHCN.CL.PA": "Land-Ocean: Southern Hemispheric Means",
                     "landZonAnn.Ts.GHCN.CL.PA": "Station: Annual Zonal Means",
                     "mixedZonAnn.Ts.ERSSTV5.GHCN.CL.PA": "Land-Ocean: Annual Zonal Means",
                     }
    if filename in display_names:
        return display_names[filename]
    else:
        return ""


def step5_output_one(item):
    (meta, data, wt, ann, monmin) = item
    title = meta.title

    try:
        title = title.decode()
    except:
        pass

    XBAD = 9999
    iy1tab = 1880
    zone_titles = step5_zone_titles()
    months_data = meta.months_data
    cur_year, m = divmod(months_data - 1, 12)
    iyrs = cur_year - 1879
    m += 1
    cur_month = str(m).zfill(2)
    # If we have the complete year we can add the current year to the zonal means files
    if m == 12:
        iyrsp = iyrs
    else:
        iyrsp = iyrs - 1
    cur_year = str(cur_year)
    titl2 = ' zones:  90->64.2->44.4->23.6->0->-23.6->-44.4->-64.2->-90                      '
    iyrbeg = meta.yrbeg
    jzm = len(ann)
    monm = iyrs * 12

    mode = meta.mode
    out = open_step5_outputs(meta, mode)
    if mode == 'mixed':
        # Check that land and ocean have same range, otherwise, divert
        # output.
        if meta.land_month_range != meta.ocean_month_range:
            # Send output to a set of files starting with "tainted".
            # Note that the original, "mixed", files will have been
            # truncated: This stops anyone using their contents.
            out = open_step5_outputs(meta, 'tainted')

    # Create and write out the header record of the output files.
    # Remove everything up to the first ')' of the title.
    if mode == 'mixed':
        data_category = 'Land-Ocean'
        sources = 'sources:  GHCN-v4 1880-' + cur_month + '/' + cur_year + ' + SST: ERSST v5 1880-' + cur_month + '/' + cur_year + '\n'
    elif mode == 'ocean':
        data_category = 'Ocean'
        sources = 'sources:  SST: ERSST v5 1880-' + cur_month + '/' + cur_year + '\n'
    else:
        data_category = 'Station'
        sources = 'sources:  GHCN-v4 1880-' + cur_month + '/' + cur_year + ' (meteorological stations only)\n'

    header = ' ' * 18 + 'Annual mean ' + data_category + ' Temperature Index in degrees Celsius\n' + \
             ' ' * 38 + 'selected zonal means\n' + ' ' * 38 + '--------------------\n' + ' ' * 20 + \
             sources + ' ' * 20 + \
             'using elimination of outliers and homogeneity adjustment\n' + ' ' * 25 + \
             'Note: ***** = missing - base period: 1951-1980\n'
    print(header, file=out[0])
    # iord literal borrowed exactly from Fortran...
    iord = [16, 14, 15, 9, 10, 11, 1, 2, 3, 4, 5, 6, 7, 8]
    # ... and then adjusted for Python index convention.
    iord = list(map(lambda x: x - 1, iord))

    # Display the annual means.

    def annasstr(z):
        """Helper function that returns the annual anomaly for zone *z*
        as a string representation of an integer (the integer is the
        anomaly scaled by 100 to convert to centikelvin).

        The returned value is a string that is 5 characters long.  If
        the integer will not fit into a 5 character string, '*****' is
        returned (this emulates the Fortran convention of formatting
        999900 (which is the XBAD value in centikelvin) as a '*****'.

        The year, *iy*, is lexically captured which is a bit horrible.
        """
        x = int(math.floor(100 * ann[z][iy] + 0.5))
        x = '%5d' % x
        if len(x) > 5:
            return '*****'
        return x

    banner = """
                           24N   24S   90S     64N   44N   24N   EQU   24S   44S   64S   90S
Year  Glob  NHem  SHem    -90N  -24N  -24S    -90N  -64N  -44N  -24N  -EQU  -24S  -44S  -64S
""".strip('\n')
    print(banner, file=out[0])
    for iy in range(iy1tab - iyrbeg, iyrsp):
        iyr = iyrbeg + iy
        row_data = [annasstr(zone) for zone in iord]
        row_data = [str(int(x) / 100) if "*" not in x else x for x in row_data]
        row_data = ["-" + x[2:] if x[:2] == "-0" else x for x in row_data]
        row_data = [x[1:] if x[0] == "0" else x for x in row_data]
        for i in range(len(row_data)):
            if str(row_data[i])[-2] == '.':
                row_data[i] = (str(row_data[i]) + "0").rjust(5)
            else:
                row_data[i] = str(row_data[i]).rjust(5)

        formatting = '%4d' + ' %s' * 3 + '  ' + ' %s' * 3 + '  ' + ' %s' * 8
        values = tuple([iyr] + row_data)
        print(formatting % values, file=out[0])
    # The trailing banner is just like the repeated banner, except that
    # "Year  Glob  NHem  SHem" appears on on the first line, not the
    # second line (and the same for the "Year" that appears at the end
    # of the line).  *sigh*.
    banner = banner.split('\n')

    banner[0] = banner[1][:24] + banner[0][24:] + ' Year'
    banner[1] = ' ' * 24 + banner[1][24:-5]
    banner = '\n'.join(banner)

    tit = ['GLOBAL', 'N.HEMISPH.', 'S.HEMISPH.']
    # Shift the remaining 3 output files so that the indexing works out.
    out = out[1:]
    banner = 'Year    Jan   Feb   Mar   Apr   May   Jun   Jul   Aug   Sep   Oct   Nov   Dec     J-D   D-N     DJF   MAM   JJA   SON'
    # All the "WRITE(96+J" stuff in the Fortran is replaced with this
    # enumeration into the *out* array (an array of file descriptors).
    for j, outf in enumerate(out):
        header = ' ' * 18 + tit[j] + ' ' + data_category + \
                 ' Temperature Index in degrees Celsius   base period: 1951-1980\n\n' + \
                 ' ' * 30 + sources + \
                 ' ' * 30 + 'using elimination of outliers and homogeneity adjustment\n' + \
                 ' ' * 30 + 'Notes: 1950 DJF = Dec 1949 - Feb 1950 ;  ***** = missing\n\n' + \
                 ' ' * 83 + 'AnnMean'
        print(header, file=outf)
        print(banner, file=outf)
        for iy in range(iy1tab - iyrbeg, iyrs):
            # Each year is formatted as a row of 18 numbers (12 months,
            # 2 different annual anomalies, and 4 seasonal).
            row = [100 * XBAD] * 18

            # *data* for this zone, avoids some duplication of code.

            zdata = data[iord[j]]

            # 4 seasons.
            season = [9999] * 4
            if iy > 0:
                season[0] = zdata[iy - 1][11] + zdata[iy][0] + zdata[iy][1]
            for s in range(1, 4):
                season[s] = sum(zdata[iy][s * 3 - 1:s * 3 + 2])
            # Each season slots into slots 14 to 17 of *row*.
            for s, x in enumerate(season):
                if x < 8000:
                    row[14 + s] = int(round(100.0 * x / 3))
            # Meteorological Year is average of 4 seasons.
            metann = sum(season)
            if metann < 8000:
                row[13] = int(round(100.0 * metann / 12))
            # Calendar year as previously computed.
            calann = ann[iord[j]][iy]

            # For final year of data, suppress annual anomaly value unless
            # December is present (assuming a full year of data?).
            if iy == iyrs - 1 and zdata[iy][-1] > 8000:
                calann = 9999

            if calann < 8000:
                row[12] = int(round(100.0 * ann[iord[j]][iy]))
            # Fill in the months.
            for m in range(12):
                row[m] = int(round(100.0 * zdata[iy][m]))

            # Convert each of *row* to a string, storing the results in
            # *sout*.
            formatted_row = [None] * len(row)
            for i, x in enumerate(row):
                x = '%5d' % x
                if len(x) > 5:
                    x = '  ***'
                formatted_row[i] = x
            year = iyrbeg + iy

            formatted_row = [str(int(x) / 100) if "*" not in x else x for x in formatted_row]
            formatted_row = ["-" + x[2:] if x[:2] == "-0" else x for x in formatted_row]
            formatted_row = [x[1:] if x[0] == "0" else x for x in formatted_row]
            for i in range(len(formatted_row)):
                if str(formatted_row[i])[-2] == '.':
                    formatted_row[i] = (str(formatted_row[i]) + "0").rjust(6)
                else:
                    formatted_row[i] = str(formatted_row[i]).rjust(6)
            print(('%4d ' + '%s' * 12 + '  %s%s  ' + '%s' * 4) % tuple([year] + formatted_row),
                  file=outf)

    # Save monthly means on disk.
    zono = open(os.path.join(RESULT_DIR, make_filename(meta, 'ZON') + '.npz'), 'wb')
    result = []

    titl2 = titl2.encode()
    if type(title) is str:
        title = title.encode()
    meta_data = info_from_meta(meta) + [title, titl2]
    meta_data = np.array(meta_data, dtype=object)
    for jz in range(jzm):
        result.append([[zone_titles[jz].encode()], np.array(data[jz]).ravel()])
    np.savez_compressed(zono, *result, meta=meta_data)
    zono.close()


def open_step5_outputs(meta, mode):
    """
    Open a set of Step 5 output files (there are 4) and return a list of
    the open file objects.  *meta* is the metadata for this
    analysis.  *mode* is a prefix used for the names, it is
    usually one of 'land', 'ocean', or 'mixed'.
    """

    parts = ['ZonAnn', 'GLB', 'NH', 'SH']
    files = [
        open(os.path.join(RESULT_DIR, make_text_filename(meta, mode, part)), 'w')
        for part in parts]
    return files
