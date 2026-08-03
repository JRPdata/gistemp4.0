# Compares two subbox files (paths as args) to verify they are near equivalent
# (within float32 tolerance=1e-5)
# run from ./tool
import sys
import struct
import math
import fort

TOLERANCE = 1e-5
MISSING = 9999


class SubboxReader:

    def __init__(self, rawfile, bos='>'):
        self.bos = bos
        self.file_obj = open(rawfile, 'rb')
        self.f = fort.File(self.file_obj, bos=bos)

        rec = self.f.readline()

        (
            self.mo1,
            self.kq,
            self.mavg,
            self.monm,
            self.monm4,
            self.yrbeg,
            self.missing_flag,
            self.precipitation_flag,
            self.title
        ) = struct.unpack(self.bos + '8i80s', rec)

        self.title = self.title.decode('latin-1').strip()

        assert self.mavg == 6

    def records(self):

        fmt = "iiiiiiif%df" % self.mo1
        expected_size = struct.calcsize(self.bos + fmt)

        while True:

            rec = self.f.readline()

            if not rec:
                break

            if len(rec) != expected_size:
                continue

            fields = struct.unpack(self.bos + fmt, rec)

            lat_S, lat_N, lon_W, lon_E = [
                f / 100.0 for f in fields[1:5]
            ]

            series = fields[8:]

            yield (
                lat_S,
                lat_N,
                lon_W,
                lon_E,
                series
            )

    def close(self):
        self.file_obj.close()


def compare_sbbx(file1, file2, tolerance=0.0):

    a = SubboxReader(file1)
    b = SubboxReader(file2)

    print("Comparing:")
    print(file1)
    print(file2)

    print("Years:", a.yrbeg, b.yrbeg)
    print("Months:", a.monm, b.monm)

    if a.yrbeg != b.yrbeg or a.monm != b.monm:
        raise RuntimeError("SBBX files have different time ranges")

    n_diff = 0
    n_values = 0
    max_diff = 0
    max_info = None

    for irec, (ra, rb) in enumerate(zip(a.records(), b.records())):

        lat_Sa, lat_Na, lon_Wa, lon_Ea, sa = ra
        lat_Sb, lat_Nb, lon_Wb, lon_Eb, sb = rb

        if (
            lat_Sa != lat_Sb or
            lat_Na != lat_Nb or
            lon_Wa != lon_Wb or
            lon_Ea != lon_Eb
        ):
            raise RuntimeError(
                f"Record {irec}: spatial mismatch"
            )

        for i, (va, vb) in enumerate(zip(sa, sb)):

            n_values += 1

            if va == MISSING and vb == MISSING:
                continue

            diff = va - vb

            if abs(diff) > tolerance:

                year = a.yrbeg + i // 12
                month = i % 12 + 1
                """
                print(
                    f"DIFF record={irec} "
                    f"lat={lat_Sa:.2f}:{lat_Na:.2f} "
                    f"lon={lon_Wa:.2f}:{lon_Ea:.2f} "
                    f"{year}-{month:02d} "
                    f"A={va:.10f} "
                    f"B={vb:.10f} "
                    f"diff={diff:.10f}"
                )
                """

                n_diff += 1

            if abs(diff) > max_diff:
                max_diff = abs(diff)
                max_info = (
                    irec,
                    lat_Sa,
                    lat_Na,
                    lon_Wa,
                    lon_Ea,
                    i,
                    va,
                    vb,
                    diff
                )

    print()
    print("Compared values:", n_values)
    print("Differences:", n_diff)

    print("MAX DIFFERENCE:")
    print(max_info)
    print("max diff =", max_diff)

    a.close()
    b.close()


if __name__ == "__main__":

    if len(sys.argv) != 3:
        print(
            "Usage: python compare_sbbx.py file1 file2"
        )
        sys.exit(1)

    compare_sbbx(
        sys.argv[1],
        sys.argv[2],
        tolerance=TOLERANCE
    )
