"""
Minimal reader/writer for gfortran-style sequential unformatted files
compiled with -fconvert=big-endian -frecord-marker=4 (i.e. each record is
wrapped with a 4-byte big-endian record-length marker before and after).
"""
import struct


def read_fortran_record(f):
    """Read one record, return its raw payload bytes, or None at EOF."""
    marker = f.read(4)
    if len(marker) == 0:
        return None
    if len(marker) != 4:
        raise IOError("truncated record marker")
    n = struct.unpack('>i', marker)[0]
    data = f.read(n)
    if len(data) != n:
        raise IOError("truncated record payload")
    end = f.read(4)
    n2 = struct.unpack('>i', end)[0]
    if n != n2:
        raise IOError(f"record marker mismatch: {n} != {n2}")
    return data


def write_fortran_record(f, payload: bytes):
    n = len(payload)
    marker = struct.pack('>i', n)
    f.write(marker)
    f.write(payload)
    f.write(marker)
