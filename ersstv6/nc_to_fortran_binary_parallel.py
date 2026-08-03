# nc_to_fortran_binary_parallel.py
# generate the fortran binaries into input_files/

# optional -jN (optional parallelization, not really needed for computers with SSDs)
import netCDF4 as nc
import numpy as np
import struct
from datetime import datetime
import sys
import os
import glob
import re
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

GISTEMP_PATH = toolconfig.GISTEMP_PATH
SST_NCEI_PATH = toolconfig.SST_NCEI_PATH
TOOL_PATH = GISTEMP_PATH / "tool"
sys.path.insert(0, TOOL_PATH)
# this folder
ERSSTV6_PROGRAM_PATH = GISTEMP_PATH / "ersstv6"

input_pattern = "ersst.v*.nc"

from fort import File

missing_flag = 9999.0

# skip sst temps earlier than this year
min_year = 1880

def parse_date_from_filename(filename):
    """
    Extract YYYYMMDD from:
    ersst.v{digit}.{YYYYMMDD}.nc
    """
    m = re.search(r'ersst\.v\d+\.(\d{6})\.nc$', filename)
    if not m:
        raise ValueError(f"Could not parse date from {filename}")

    date = datetime.strptime(f"{m.group(1)}01", "%Y%m%d")
    return date.year, date.month


def verify_contiguous_months(files):
    """
    Ensure files contain every month from first to last without gaps.
    """

    dates = [parse_date_from_filename(f) for f in files]

    for (y1, m1), (y2, m2) in zip(dates[:-1], dates[1:]):

        # advance one month
        if m1 == 12:
            expected = (y1 + 1, 1)
        else:
            expected = (y1, m1 + 1)

        if (y2, m2) != expected:
            raise RuntimeError(
                f"Missing month between {y1}-{m1:02d} and {y2}-{m2:02d}"
            )


def process_one_file(netcdf_file):

    import os
    import threading
    import os
    import time
    import multiprocessing as mp

    year, month = parse_date_from_filename(netcdf_file)

    if year < min_year:
        return None

    with nc.Dataset(netcdf_file, 'r') as ds:

        if 'sst' not in ds.variables:
            raise ValueError(f"No sst variable in {netcdf_file}")

        sst = ds.variables['sst']

        # Expected:
        # float32 sst(time, lev, lat, lon)
        if sst.ndim != 4:
            raise ValueError(
                f"Unexpected sst dimensions in {netcdf_file}: {sst.shape}"
            )

        # take first time and first level
        data_var = sst[0, 0, :, :]

        # dynamically determine latitude ordering
        lat_var_name = 'lat' if 'lat' in ds.variables else 'latitude'
        lats = ds.variables[lat_var_name][:]

        is_descending = lats[0] > lats[-1]

        # Need mask before filling
        fraction_unmasked = np.float32(
            1 - np.mean(np.ma.getmaskarray(data_var))
        )

        data = np.ma.filled(data_var, missing_flag)

        if is_descending:
            data = data[::-1, :]

    if data.shape != (89, 180):
        raise ValueError(
            f"{netcdf_file}: expected (89,180), got {data.shape}"
        )

    title = (
        f'Sea Surface Temperatures\t(C)\t'
        f'{year}\t{month:02d} - {year}\t{month:02d} {fraction_unmasked}'
    )

    title = title.ljust(80)

    header = title.encode('ascii')

    output_file = str(ERSSTV6_PROGRAM_PATH / 'input_files' / f"ersst.v6.{year}{month:02d}.bin")

    data = data.astype(np.float32)

    with open(output_file, 'wb') as f:
        s = File(f, bos='>')

        arr = data.flatten().tolist()

        header_and_data_bytes = struct.pack(
            f'>80s{len(arr)}f',
            header,
            *arr
        )

        s.writeline(header_and_data_bytes)

    return output_file


def convert_netcdf_to_fortran_binaries(jobs):

    files = sorted(glob.glob(input_pattern))

    if not files:
        raise RuntimeError("No NetCDF files found")

    print(f"Found {len(files)} files")

    verify_contiguous_months(files)

    print(
        "Verified continuous sequence:",
        parse_date_from_filename(files[0]),
        "to",
        parse_date_from_filename(files[-1])
    )
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from tqdm import tqdm

    with ProcessPoolExecutor(max_workers=jobs) as executor:

        futures = {
            executor.submit(process_one_file, f): f
            for f in files
        }

        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Converting ERSSTv6",
        ):
            src = futures[future]

            try:
                result = future.result()

            except Exception:
                print(f"FAILED: {src}")
                raise


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Convert ERSSTv6 monthly files to Fortran binaries"
    )
    
    parser.add_argument(
        "-i",
        "--inputpath",
        type=str,
        default=str(SST_NCEI_PATH),
        help="Path containing monthly (NCEI) ERSSTv6 files, i.e. ersst.v*.nc files"
    )

    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="Number of parallel workers (default: 1)"
    )

    args = parser.parse_args()
    
    convert_netcdf_to_fortran_binaries(args.jobs)
