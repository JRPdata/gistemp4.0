import sys
import struct
# fort is the fort.py file in the gistemp directory tool/
import fort
import math

#: Magic number used to indicate missing data.
MISSING = 9999


class SubboxReader:
    """Reads GISS subbox files (SBBX) as per the original Fortran program."""

    def __init__(self, rawfile, bos='>'):
        self.bos = bos
        # Open the file using the built-in open() in binary mode and pass the file object to fort.File
        self.file_obj = open(rawfile, 'rb')
        self.f = fort.File(self.file_obj, bos='>')  # Adjust 'bos' as necessary
        # Proceed with the rest of the initialization code
        rec = self.f.readline()
        (self.mo1, self.kq, self.mavg, self.monm, self.monm4, self.yrbeg,
         self.missing_flag, self.precipitation_flag, self.title) = struct.unpack(self.bos + '8i80s', rec)

        self.title = self.title.decode('latin-1').strip()

        assert self.mavg == 6, "Only monthly averages supported"

    def inspect_record(self, rec):
        if len(rec) == 36:
            print("Record size is 36 bytes.")
            print("Content of the record:", rec)
        elif len(rec) != 6992:
            print(f"Unexpected record size: {len(rec)} bytes.")
        else:
            print("Record size is as expected.")

        # Further checks can be implemented based on expected content

    def read_data(self):
        """Reads the data part of the SBBX file."""
        tin = [self.missing_flag] * self.monm
        output_data = []

        for N in range(8000):
            try:
                #print(N)
                rec = self.f.readline()
            except:
                continue
            if not rec:
                break
            fmt = "iiiiiiif%df" % self.mo1

            # Define the format string based on expected record structure
            fmt = "iiiiiiif%df" % self.mo1
            expected_size = struct.calcsize(self.bos + fmt)

            if len(rec) != expected_size:
                #print(f"Error: Expected {expected_size} bytes but got {len(rec)} bytes. Skipping record")
                continue
            else:
                #print("Expected {expected_size:", expected_size, len(rec))
                pass

            fields = list(struct.unpack(self.bos + fmt, rec))
            series = fields[8:]
            lat_S, lat_N, lon_W, lon_E = [f / 100.0 for f in fields[1:5]]
            end_year_plus_one = int(math.ceil(self.yrbeg + self.monm / 12))
            for iyr in range(self.yrbeg, end_year_plus_one):
                iok = sum(1 for x in series[:12] if x != self.missing_flag)
                if iok > 0:
                    year_data = f"{lat_S:7.2f}{lat_N:7.2f}{lon_W:8.2f}{lon_E:8.2f}{iyr:5d}" + \
                                ''.join(f"{v:8.2f}" for v in series[:12])
                    yield year_data
                series = series[12:]

def read_sbbx_to_txt(filename):
    reader = SubboxReader(filename)
    end_year = int(math.ceil(reader.yrbeg + reader.monm / 12 - 1))
    title_str = f"{reader.title: <80} {reader.yrbeg}-{end_year} varying base period"
    print(title_str)
    print(f"Missing data flag: {reader.missing_flag}")

    output_filename = filename.strip() + '.txt'
    with open(output_filename, 'w') as outfile:
        outfile.write(f"{title_str}\n")
        outfile.write(" latitude-rnge longitude-range year     Jan     Feb     Mar     Apr     May     Jun     Jul     Aug     Sep     Oct     Nov     Dec\n")
        for line in reader.read_data():
            outfile.write(line + '\n')

    reader.file_obj.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py sbbx-file_name")
        sys.exit(1)

    filename = sys.argv[1]
    read_sbbx_to_txt(filename)
