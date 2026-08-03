# ERSSTv6_to_SBBX.py
# (parallelized) coordinator:
#  1. converts netCDF monthly files to fortran binaries
#  2. runs the fortran programs in order to create the ERSSTv6 SBBX file
#     it is an alternative to running the standalone nc_to_fortran_binary_parallel.py followed by updateERSST.sh

import netCDF4 as nc
import numpy as np
import struct
from datetime import datetime
import sys
import os
import re
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import subprocess
import toolconfig
from pathlib import Path
import glob

GISTEMP_PATH = toolconfig.GISTEMP_PATH
SST_NCEI_PATH = toolconfig.SST_NCEI_PATH
TOOL_PATH = GISTEMP_PATH / "tool"
sys.path.insert(0, str(TOOL_PATH))
# this folder
ERSSTV6_PROGRAM_PATH = GISTEMP_PATH / "ersstv6"

INPUT_PATTERN = "ersst.v*.nc"

# Paths to compiled fortran binaries from NASA/GISS GISTEMP (in this folder)
MASK_REGRID = str(ERSSTV6_PROGRAM_PATH / "MaskRegrid")
REARRANGE = str(ERSSTV6_PROGRAM_PATH / "rearrange.ERSST")
TRIM = str(ERSSTV6_PROGRAM_PATH / "trimSBBX")

from fort import File
# GISTEMP specifics
missing_flag = 9999.0
# skip sst temps earlier than this year
min_year = 1880

require_min_year = True # Make sure we have all the data going back to at least 1880 (even olderer is fine if present as we will skip it)

def run_pipeline(yr1, month1, yr2, month2):

    print("Removing old files")

    for f in [
        "ERdSST_monthly",
        "SBBX.ERSST.upd",
        "SBBX.ERSSTv6"
    ]:
        if os.path.exists(f):
            os.remove(f)


    print("Running MaskRegrid")
    result = subprocess.run(
        [
            MASK_REGRID,
            str(yr1),
            str(month1),
            str(yr2),
            str(month2),
            "0.5"
        ],
        cwd=str(ERSSTV6_PROGRAM_PATH),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )

    print(result.stdout)

    if result.returncode != 0:
        raise RuntimeError(
            f"MaskRegrid failed (exit code {result.returncode})"
        )

    if "could not find file" in result.stdout.lower():
        raise RuntimeError(
            "MaskRegrid could not find an input file"
        )
            
    print("Created:")
    
    subprocess.run(["ls", "-l", "ERdSST_monthly"], cwd=str(ERSSTV6_PROGRAM_PATH))


    print("Running rearrange.ERSST")

    subprocess.run(
        [
            REARRANGE,
            str(yr1),
            str(month1),
            str(yr2),
            str(month2)
        ],
        cwd=str(ERSSTV6_PROGRAM_PATH),
        check=True
    )


    print("Created:")
    subprocess.run(["ls", "-l", "SBBX.ERSST.upd"], cwd=str(ERSSTV6_PROGRAM_PATH),)


    print("Running trimSBBX")

    subprocess.run(
        [
            TRIM,
            "SBBX.ERSST.upd",
            "SBBX.ERSSTv6"
        ],
        cwd=str(ERSSTV6_PROGRAM_PATH),
        check=True
    )

    print("Created:")
    subprocess.run(["ls", "-l", "SBBX.ERSSTv6"], cwd=str(ERSSTV6_PROGRAM_PATH))

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

    output_file = f'input_files/ersst.v6.{year}{month:02d}.bin'

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
    pattern = os.path.join(str(SST_NCEI_PATH), INPUT_PATTERN)
    files = sorted(glob.glob(pattern))

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

    first_year, first_month = parse_date_from_filename(files[0])
    last_year, last_month = parse_date_from_filename(files[-1])
    
    if require_min_year:
        if first_year != min_year and first_month != 1:
            raise RuntimeError(f"As required min_year == {required_min_year}, monthly .nc data is required going back to {f}-01. Only have from {first_year}-{first_month}.")

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

    print("All ERSST fortran monthly binaries successfully created.")

    return first_year, first_month, last_year, last_month

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Convert ERSSTv6 monthly files to Fortran binaries"
    )

    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="Number of parallel workers (default: 1)"
    )

    args = parser.parse_args()

    yr1, month1, yr2, month2 = (
        convert_netcdf_to_fortran_binaries(args.jobs)
    )

    print(f"Running fortran programs to combine monthly fortran binaries from {min_year}-{month1:02d} to {yr2}-{month2:02d}")
    run_pipeline(
        min_year,
        month1,
        yr2,
        month2
    )
