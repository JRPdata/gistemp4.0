#!/usr/bin/env python3
"""
ersstv6_to_sbbx.py

Consolidated, pure-Python reproduction of the GISTEMP ERSSTv6 -> SBBX
pipeline that is normally done by three separate Fortran programs
(MaskRegrid.f, rearrange_ERSST.f, trimSBBX.f) glued together by
ERSSTv6_to_SBBX.py / nc_to_fortran_binary_parallel.py.

Differences from the original pipeline (by design, per requirements):
  * No intermediate Fortran monthly binaries (ersst.v6.yyyymm.bin) and no
    intermediate ERdSST_monthly / SBBX.ERSST.upd files are written -- .nc
    files are read directly and regridded/rearranged/trimmed in memory.
  * The open-ocean mask is read from a .npz file (produced once via
    convert_ersst_mask_to_npz.py) instead of the compiled MaskRegrid.exe
    reading ERSST_open_ocean_mask.bin.
  * This script only supports a full "from scratch" rebuild starting at
    1880-01 (like the original when SBBX.ERSSTv6 does not yet exist / is
    not being "updated") -- it does not implement rearrange_ERSST.f's
    incremental-update-of-an-existing-SBBX.ERSSTv6 code path.

Everything else -- grid geometry, box ordering, trimming rules, and the
Fortran unformatted binary layout of the output SBBX file -- is translated
as literally as possible from the .f sources so the output can be
byte-compared against a real Fortran-produced SBBX.ERSSTv6.

Optionally also writes a .npz form of the (untrimmed, fixed-length,
9999-padded) subbox data, in the same "arr"/"meta" layout used elsewhere
in this project for SBBX-derived npz files.
"""
import argparse
import glob
import os
import re
import sys
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import netCDF4 as nc

from ea_grid import EAGrid
from fortran_io import write_fortran_record

import toolconfig

try:
    from tqdm import tqdm
except ImportError:  # tqdm is a nicety, not a hard requirement
    def tqdm(iterable, **kwargs):
        return iterable

MIN_YEAR = 1880       # GISTEMP iyrbeg -- fixed start of the SBBX time series
MISSING_FLAG = np.float32(9999.0)
NUM_SUBBOXES = 8000
FILENAME_RE = re.compile(r'ersst\.v(\d+)\.(\d{6})\.nc$')

from pathlib import Path
GISTEMP_PATH = toolconfig.GISTEMP_PATH
SST_NCEI_PATH = toolconfig.SST_NCEI_PATH
TOOL_PATH = GISTEMP_PATH / "tool"
sys.path.insert(0, TOOL_PATH)
# this folder
ERSSTV6_PROGRAM_PATH = GISTEMP_PATH / "ersstv6py"

# ----------------------------------------------------------------------
# NetCDF reading (mirrors process_one_file() in ERSSTv6_to_SBBX.py, but
# returns the array instead of writing a Fortran binary file to disk)
# ----------------------------------------------------------------------
def parse_date_from_filename(filename):
    m = FILENAME_RE.search(filename)
    if not m:
        raise ValueError(f"Could not parse date from {filename}")
    date = datetime.strptime(f"{m.group(2)}01", "%Y%m%d")
    return date.year, date.month


def verify_contiguous_months(files):
    dates = [parse_date_from_filename(f) for f in files]
    for (y1, m1), (y2, m2) in zip(dates[:-1], dates[1:]):
        expected = (y1 + 1, 1) if m1 == 12 else (y1, m1 + 1)
        if (y2, m2) != expected:
            raise RuntimeError(f"Missing month between {y1}-{m1:02d} and {y2}-{m2:02d}")


def read_one_nc(netcdf_file):
    """Returns (year, month, data) with data shaped (89,180) float32, ascending
    latitude, missing values filled with MISSING_FLAG -- or None if the file's
    year predates MIN_YEAR."""
    year, month = parse_date_from_filename(netcdf_file)
    if year < MIN_YEAR:
        return None

    with nc.Dataset(netcdf_file, 'r') as ds:
        if 'sst' not in ds.variables:
            raise ValueError(f"No sst variable in {netcdf_file}")

        sst = ds.variables['sst']
        if sst.ndim != 4:
            raise ValueError(f"Unexpected sst dimensions in {netcdf_file}: {sst.shape}")

        data_var = sst[0, 0, :, :]

        lat_var_name = 'lat' if 'lat' in ds.variables else 'latitude'
        lats = ds.variables[lat_var_name][:]
        is_descending = lats[0] > lats[-1]

        data = np.ma.filled(data_var, MISSING_FLAG)
        if is_descending:
            data = data[::-1, :]

    if data.shape != (89, 180):
        raise ValueError(f"{netcdf_file}: expected (89,180), got {data.shape}")

    return year, month, data.astype(np.float32)


def load_all_months(input_folder, jobs):
    pattern = os.path.join(input_folder, "ersst.v*.nc")
    in_files = sorted(glob.glob(pattern))
    if not in_files:
        raise RuntimeError(f"No NetCDF files found matching {pattern}")

    # Some source directories contain files predating MIN_YEAR (1880); skip
    # those entirely rather than feeding them into the contiguity check.
    files = []
    for file in in_files:
        year, month = parse_date_from_filename(file)
        if year >= MIN_YEAR:
            files.append(file)
    if not files:
        raise RuntimeError(f"No NetCDF files with year >= {MIN_YEAR} found matching {pattern}")

    print(f"Found {len(in_files)} files total ({len(files)} from {MIN_YEAR} onwards)")
    verify_contiguous_months(files)
    first_year, first_month = parse_date_from_filename(files[0])
    last_year, last_month = parse_date_from_filename(files[-1])
    print(f"Verified continuous sequence: ({first_year}, {first_month}) to ({last_year}, {last_month})")

    if first_year != MIN_YEAR or first_month != 1:
        raise RuntimeError(
            f"This tool requires monthly .nc data starting at {MIN_YEAR}-01; "
            f"the earliest available file is {first_year}-{first_month:02d}."
        )

    months = {}
    if jobs and jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            futures = {executor.submit(read_one_nc, f): f for f in files}
            for future in tqdm(as_completed(futures), total=len(futures), desc="Reading ERSSTv6 .nc files"):
                src = futures[future]
                try:
                    result = future.result()
                except Exception:
                    print(f"FAILED: {src}")
                    raise
                if result is not None:
                    yr, mo, data = result
                    months[(yr, mo)] = data
    else:
        for f in tqdm(files, desc="Reading ERSSTv6 .nc files"):
            result = read_one_nc(f)
            if result is not None:
                yr, mo, data = result
                months[(yr, mo)] = data

    return months, last_year, last_month


# ----------------------------------------------------------------------
# Title construction (mirrors the fixed-format CHARACTER*80 title built
# in rearrange_ERSST.f)
# ----------------------------------------------------------------------
def build_title(last_year, last_month):
    part1 = 'Monthly Sea Surface Temperature anom (C)'   # title(1:40)
    assert len(part1) == 40
    part2 = ' ERSSTv6 01/1880 - '                          # title(41:59)
    assert len(part2) == 19
    part3 = f'{last_month:02d}/{last_year}'                # title(60:66)
    assert len(part3) == 7
    title = (part1 + part2 + part3).ljust(80)
    assert len(title) == 80
    return title


# ----------------------------------------------------------------------
# Main pipeline
# ----------------------------------------------------------------------
def run(args):
    grid = EAGrid()

    mask_npz = np.load(args.open_ocean_mask_npz_path)
    mask = mask_npz['mask']  # shape (12, 89, 180) float32
    if mask.shape != (12, 89, 180):
        raise ValueError(f"Unexpected mask shape {mask.shape}, expected (12, 89, 180)")

    months, last_year, last_month = load_all_months(args.input_subbox_folder, args.jobs)

    iyr1, mon1 = MIN_YEAR, 1
    iyr2, mon2 = last_year, last_month

    # Build the ordered sequence of (year, month) from iyr1/mon1 to iyr2/mon2,
    # matching MaskRegrid.f's do-loop (m1=mo1 first year, then 1; m2=12 except
    # last year=mo2).
    sequence = []
    y = iyr1
    m1 = mon1
    while y <= iyr2:
        m2 = 12 if y != iyr2 else mon2
        for m in range(m1, m2 + 1):
            sequence.append((y, m))
        m1 = 1
        y += 1
    mnew = len(sequence)
    print(f"Processing {mnew} months: {iyr1}-{mon1:02d} .. {iyr2}-{mon2:02d}")

    # ---- NTRP2EA regrid, one month at a time (natural south->north/west->east order)
    natural_order_series = np.empty((mnew, NUM_SUBBOXES), dtype=np.float32)
    for i, (yr, mo) in enumerate(tqdm(sequence, desc="Regridding to equal-area boxes")):
        data = months.get((yr, mo))
        if data is None:
            raise RuntimeError(f"Missing regridded input data for {yr}-{mo:02d}")
        wta = mask[mo - 1]  # calendar-month mask, shape (89,180)
        natural_order_series[i] = grid.regrid(data, wta, frac=args.frac, skip=float(MISSING_FLAG))

    # ---- Rearrange into Sergei canonical order, full (padded-to-December) length
    NM_padded = 12 * (iyr2 + 1 - MIN_YEAR)
    boxes_full = np.full((NUM_SUBBOXES, NM_padded), MISSING_FLAG, dtype=np.float32)
    latlon_out = np.zeros((NUM_SUBBOXES, 4), dtype=np.int64)

    for n in range(1, NUM_SUBBOXES + 1):
        isbbx = grid.ij_sn_we[n]
        box_idx = n - 1
        boxes_full[box_idx, 0:mnew] = natural_order_series[:, isbbx - 1]
        latlon_out[box_idx] = grid.latlon[isbbx]
    # months mnew..NM_padded-1 are already MISSING_FLAG from np.full()

    title = build_title(iyr2, mon2)

    if args.skip_sbbx_binary and not args.output_sbbx_npz_path:
        raise RuntimeError("--output-sbbx-npz-path is required when --skip-sbbx-binary is set")

    if not args.skip_sbbx_binary:
        write_sbbx_binary(args.output_sbbx_path, boxes_full, latlon_out, NM_padded, title)
        print(f"Wrote {args.output_sbbx_path}")

    if args.output_sbbx_npz_path:
        write_sbbx_npz(args.output_sbbx_npz_path, boxes_full, latlon_out, NM_padded, mon2, title)
        print(f"Wrote {args.output_sbbx_npz_path}")


# ----------------------------------------------------------------------
# trimSBBX.f equivalent -- also used to derive the npz per-box metadata
# ----------------------------------------------------------------------
def compute_trim_records(boxes_full, latlon_out, NM_padded):
    """
    Returns:
      infoo1: int, INFO(1) header value (NM_padded if box[0] has any valid
              data, else 1)
      records: list of 8000 dicts (in Sergei canonical order), each with:
        marker  -- int, the leading integer written for this record (the
                   *next* box's marker, or 0 sentinel for the last box)
        lts, ltn, lnw, lne, nst, nstmn, dmin -- the 7 metadata fields
        has_data -- bool, whether this box has any valid (non-missing) data
        data_idx -- int, index into boxes_full for this box's full time series
    """
    valid_counts = np.count_nonzero(boxes_full != MISSING_FLAG, axis=1)

    infoo1 = int(NM_padded) if valid_counts[0] > 0 else 1

    records = [None] * NUM_SUBBOXES

    LATO = latlon_out[0]
    LATO6 = int(valid_counts[0])
    prev_idx = 0

    for N in range(2, NUM_SUBBOXES + 1):
        cur_idx = N - 1
        LAT = latlon_out[cur_idx]
        LAT6 = int(valid_counts[cur_idx])
        MLN = int(NM_padded) if LAT6 > 0 else 1

        records[prev_idx] = {
            'marker': MLN,
            'lts': int(LATO[0]), 'ltn': int(LATO[1]),
            'lnw': int(LATO[2]), 'lne': int(LATO[3]),
            'nst': 0, 'nstmn': LATO6, 'dmin': 0.0,
            'has_data': LATO6 > 0,
            'data_idx': prev_idx,
        }

        LATO = LAT
        LATO6 = LAT6
        prev_idx = cur_idx

    # final record: box 8000, sentinel marker 0
    records[prev_idx] = {
        'marker': 0,
        'lts': int(LATO[0]), 'ltn': int(LATO[1]),
        'lnw': int(LATO[2]), 'lne': int(LATO[3]),
        'nst': 0, 'nstmn': LATO6, 'dmin': 0.0,
        'has_data': LATO6 > 0,
        'data_idx': prev_idx,
    }

    return infoo1, records


def write_sbbx_binary(output_path, boxes_full, latlon_out, NM_padded, title):
    infoo1, records = compute_trim_records(boxes_full, latlon_out, NM_padded)

    infoo = [infoo1, 1, 6, int(NM_padded), int(NM_padded) + 8, MIN_YEAR, 9999, -9999]

    with open(output_path, 'wb') as f:
        header_payload = np.array(infoo, dtype='>i4').tobytes() + title.encode('ascii')
        write_fortran_record(f, header_payload)

        for rec in tqdm(records, desc="Writing trimmed SBBX records"):
            # marker + LTS,LTN,LNW,LNE,NSt,NstMn = 7 int32 fields ...
            header_ints = np.array(
                [rec['marker'], rec['lts'], rec['ltn'], rec['lnw'], rec['lne'], rec['nst'], rec['nstmn']],
                dtype='>i4',
            )
            # ... followed by Dmin (the 7th LATLON element) as a float32, THEN
            # the data. Dropping this field shifts every subsequent record by
            # 4 bytes -- it must always be present, even though it's always
            # 0.0 for ocean data (no station-distance concept here).
            dmin_bytes = np.array([rec['dmin']], dtype='>f4').tobytes()
            if rec['has_data']:
                data_bytes = boxes_full[rec['data_idx']].astype('>f4').tobytes()
                payload = header_ints.tobytes() + dmin_bytes + data_bytes
            else:
                xbad_bytes = np.array([9999.0], dtype='>f4').tobytes()
                payload = header_ints.tobytes() + dmin_bytes + xbad_bytes
            write_fortran_record(f, payload)


def write_sbbx_npz(output_path, boxes_full, latlon_out, NM_padded, mon2, title):
    """
    Writes the (untrimmed) regridded subbox data as .npz.

    Unlike the binary SBBX file, this does NOT use trimSBBX.f's "next box's
    marker" convention. Each subbox's own meta[0] is simply the length of
    its own data array: NM_padded if it has at least one valid month, or 1
    if every month is missing (in which case the data array itself is
    collapsed to the single value [9999.0] rather than stored at full
    padded length).
    """
    valid_counts = np.count_nonzero(boxes_full != MISSING_FLAG, axis=1)

    arr = np.empty((NUM_SUBBOXES, 2), dtype=object)
    for i in range(NUM_SUBBOXES):
        has_data = valid_counts[i] > 0
        if has_data:
            marker = int(NM_padded)
            data_row = boxes_full[i].astype(np.float64)
        else:
            marker = 1
            data_row = np.array([9999.0], dtype=np.float64)

        meta_row = [
            marker,
            int(latlon_out[i, 0]), int(latlon_out[i, 1]),
            int(latlon_out[i, 2]), int(latlon_out[i, 3]),
            0, int(valid_counts[i]), 0.0,
        ]
        arr[i, 0] = meta_row
        arr[i, 1] = data_row

    meta = np.array(
        [int(NM_padded), 1, mon2, int(NM_padded), int(NM_padded) + 8, MIN_YEAR, 9999, -9999, title.encode('ascii')],
        dtype=object,
    )

    np.savez(output_path, arr=arr, meta=meta)


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def str2bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in ('true', 't', '1', 'yes', 'y'):
        return True
    if value.lower() in ('false', 'f', '0', 'no', 'n'):
        return False
    raise argparse.ArgumentTypeError(f"Expected a boolean value, got {value!r}")


def main():
    parser = argparse.ArgumentParser(
        description="Consolidated, pure-Python ERSSTv6 .nc -> SBBX.ERSSTv6 pipeline "
                    "(replaces MaskRegrid.f + rearrange_ERSST.f + trimSBBX.f)."
    )
    parser.add_argument("-j", "--jobs", type=int, default=1,
                         help="Number of parallel workers for reading .nc files (default: 1)")
    parser.add_argument("--input-subbox-folder", default=SST_NCEI_PATH,
                         help="Folder containing the ERSSTv6 monthly .nc files (default: %(default)s)")
    parser.add_argument("--open-ocean-mask-npz-path", default=str(ERSSTV6_PROGRAM_PATH / "input_files" / "ERSST_open_ocean_mask.npz"),
                         help="Path to the open-ocean mask .npz (default: %(default)s). "
                              "Produce this once with convert_ersst_mask_to_npz.py.")
    parser.add_argument("--output-sbbx-path", default=str(ERSSTV6_PROGRAM_PATH / "SBBX.ERSSTv6"),
                         help="Path to write the trimmed SBBX binary file (default: %(default)s)")
    parser.add_argument("--output-sbbx-npz-path", default=None,
                         help="Optional path to also write the (untrimmed, fixed-length) "
                              "subbox data as .npz (default: not produced)")
    parser.add_argument("--skip-sbbx-binary", type=str2bool, nargs='?', const=True, default=False,
                         help="If true, skip writing the SBBX binary file (requires "
                              "--output-sbbx-npz-path to be set). Default: False")
    parser.add_argument("--frac", type=float, default=0.50,
                         help="Fraction of overlay with good data required to count a box as "
                              "covered (matches MaskRegrid.f's optional 5th argument). (DO NOT CHANGE). Default: 0.50")

    args = parser.parse_args()

    if args.skip_sbbx_binary and not args.output_sbbx_npz_path:
        parser.error("--output-sbbx-npz-path is required when --skip-sbbx-binary is set")

    run(args)


if __name__ == '__main__':
    main()
