#!/usr/bin/env python3
"""
Convert GISTEMP's ERSST_open_ocean_mask.bin (Fortran unformatted,
big-endian, produced for use by MaskRegrid.f) into a .npz file that
ersstv6_to_sbbx.py can read directly, without needing the compiled
Fortran binary or a Fortran-record reader at pipeline run time.

The original .bin file is a single unformatted record containing:
    CHARACTER*80 title
    REAL*4       mask(180,89,12)     ! (lon, lat, month), -88..+88N

Usage:
    python3 convert_ersst_mask_to_npz.py \\
        --input-mask-bin-path input_files/ERSST_open_ocean_mask.bin \\
        --output-mask-npz-path input_files/ERSST_open_ocean_mask.npz
"""
import argparse
import numpy as np

from fortran_io import read_fortran_record

MASK_LON = 180
MASK_LAT = 89
MASK_MONTHS = 12


def convert(input_path, output_path):
    with open(input_path, 'rb') as f:
        record = read_fortran_record(f)
        if record is None:
            raise ValueError(f"{input_path}: no data record found")

    title_bytes = record[:80]
    expected_mask_bytes = MASK_LON * MASK_LAT * MASK_MONTHS * 4
    mask_bytes = record[80:80 + expected_mask_bytes]
    if len(mask_bytes) != expected_mask_bytes:
        raise ValueError(
            f"{input_path}: expected {expected_mask_bytes} bytes of mask "
            f"data after the title, got {len(mask_bytes)}"
        )

    # Fortran mask(180,89,12): first index (lon) varies fastest -> order='F'
    mask_flat = np.frombuffer(mask_bytes, dtype='>f4')
    mask_lon_lat_month = mask_flat.reshape((MASK_LON, MASK_LAT, MASK_MONTHS), order='F')

    # Reorder to (month, lat, lon) -- convenient for per-calendar-month lookup
    # against (89,180) [lat,lon] ERSST arrays.
    mask_month_lat_lon = np.ascontiguousarray(
        mask_lon_lat_month.transpose(2, 1, 0).astype(np.float32)
    )

    title = title_bytes.decode('ascii', errors='replace')

    np.savez(output_path, mask=mask_month_lat_lon, title=title)
    print(f"Read mask title: {title!r}")
    print(f"Wrote {output_path}: mask shape {mask_month_lat_lon.shape} (month, lat, lon)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input-mask-bin-path', default='input_files/ERSST_open_ocean_mask.bin',
                     help='Path to the original Fortran mask binary (default: %(default)s)')
    ap.add_argument('--output-mask-npz-path', default='input_files/ERSST_open_ocean_mask.npz',
                     help='Path to write the converted npz mask (default: %(default)s)')
    args = ap.parse_args()
    convert(args.input_mask_bin_path, args.output_mask_npz_path)


if __name__ == '__main__':
    main()
