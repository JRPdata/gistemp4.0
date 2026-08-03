# ERSSTv6 -> SBBX, consolidated pure-Python pipeline

This replaces the three-Fortran-program pipeline (`MaskRegrid.f` +
`rearrange_ERSST.f` + `trimSBBX.f`, glued together by
`ERSSTv6_to_SBBX.py`) with a single Python script that reads the ERSSTv6
monthly `.nc` files directly and writes `SBBX.ERSSTv6` (and, optionally, an
`.npz` version of the same data) without any intermediate Fortran binaries.

## Files

- **`ea_grid.py`** — the equal-area grid geometry: a line-by-line
  translation of `NTRP2EA0`/`NTRP2EA` (from `MaskRegrid.f`) and
  `def_ea_grid` (from `rearrange_ERSST.f`). This is the part that decides
  which input grid cells feed which of the 8000 output subboxes, and the
  box ordering/boundaries used in the final file.
- **`fortran_io.py`** — tiny reader/writer for gfortran
  `-fconvert=big-endian -frecord-marker=4` sequential unformatted files.
- **`convert_ersst_mask_to_npz.py`** — one-time utility: converts the
  original `ERSST_open_ocean_mask.bin` into a `.npz` file. Run this once;
  the main pipeline never touches the `.bin` file or needs a compiled
  Fortran binary for the mask.
- **`ersstv6_to_sbbx.py`** — the main pipeline (this is what you run).

## Usage

```bash
# One-time: convert the Fortran mask binary to npz
python3 convert_ersst_mask_to_npz.py \
    --input-mask-bin-path input_files/ERSST_open_ocean_mask.bin \
    --output-mask-npz-path input_files/ERSST_open_ocean_mask.npz

# Run the consolidated pipeline
python3 ersstv6_to_sbbx.py \
    --input-subbox-folder ../sst \
    --open-ocean-mask-npz-path input_files/ERSST_open_ocean_mask.npz \
    --output-sbbx-path SBBX.ERSSTv6 \
    -j 4

```

## CLI arguments (`ersstv6_to_sbbx.py`)

| Argument | Default | Notes |
|---|---|---|
| `-j`, `--jobs` | `1` | Parallel workers for reading `.nc` files |
| `--input-subbox-folder` | `../sst` | Folder containing `ersst.v*.nc` |
| `--open-ocean-mask-npz-path` | `input_files/ERSST_open_ocean_mask.npz` | Produced once via `convert_ersst_mask_to_npz.py` |
| `--output-sbbx-path` | `SBBX.ERSSTv6` | Trimmed binary output (same format as the real pipeline's) |
| `--output-sbbx-npz-path` | *(none — not produced)* | Optional `.npz` form of the regridded subbox data |
| `--skip-sbbx-binary` | `False` | If `True`, skip the binary file (requires `--output-sbbx-npz-path`) |
| `--frac` | `0.50` | Same as `MaskRegrid.f`'s optional 5th CLI argument |

## Changelog

- **Fixed a binary record layout bug.** Each subbox record's leading
  metadata is `marker` + `LATLON(7)` = 8 total 4-byte fields (`lts, ltn,
  lnw, lne, nst, nstmn, dmin`) before the data array, per `trimSBBX.f`'s
  `WRITE(11) MLN,LATLON,ARRAY`. An earlier version of this script wrote
  only 7 fields (dropping `dmin`), which shifted every subsequent byte in
  the file by 4 bytes and broke downstream parsing. Fixed.
- **`.nc` files predating 1880 are now skipped** rather than causing a
  contiguity-check failure, so a source folder that happens to contain
  pre-1880 files no longer needs to be pre-filtered by hand.
- **`.npz` semantics changed** (see below): every box's data array is now
  padded to the same December-inclusive length as the binary file
  (`NM_padded`) — except boxes with zero valid months anywhere, which
  collapse to a single `[9999.0]` element. `meta[0]`/`meta[3]` are now
  `NM_padded` rather than the unpadded real-month count. Each box's own
  `arr[i][0][0]` is simply `len(arr[i][1])` (self-referential — not
  `trimSBBX.f`'s "next box's marker" convention, which only applies to the
  binary file). `arr[i][0]` is now a plain Python list rather than an
  object-dtype array.

## Fidelity notes / design decisions

- **Grid translation is literal.** `ea_grid.py` mirrors the Fortran
  `GOTO`-based loops directly (rather than a "cleaner" rewrite) so it can
  be audited line-by-line against `MaskRegrid.f` / `rearrange_ERSST.f`.
  As a sanity check, Sergei-order box `n=1` comes out as
  `lts=6416, ltn=6551, lnw=-18000, lne=-17100` — exactly matching the
  example box you provided.
- **`total_wt` precision.** In `MaskRegrid.f`, `total_wt` is declared
  `REAL*4` while `WEIGHT`/`VALUE` are implicitly `REAL*8`. This script
  reproduces that: `total_wt` is accumulated in `float32`, in the same
  edge order as the Fortran `DO 510` loops, while `WEIGHT`/`VALUE` are
  accumulated in `float64` (vectorized with `numpy.bincount` for speed).
  This only affects the `WEIGHT > total_wt*frac` coverage-threshold gate,
  never the SST value itself — and in practice this gate is rarely near
  its boundary, so it's unlikely to ever cause a mismatch, but it's worth
  knowing about if you ever see a boundary-case discrepancy.
- **Full rebuild only.** Like the original scripts, this only implements
  a from-scratch build starting at 1880-01 (`iyr1 == iyrbeg`). It does not
  implement `rearrange_ERSST.f`'s incremental "update an existing
  `SBBX.ERSSTv6`" code path, since the consolidated in-memory pipeline has
  no reason to need it — every run reprocesses all `.nc` files found.
- **Binary file month count is padded to December**, exactly like
  `rearrange_ERSST.f`/`trimSBBX.f`: `INFO(4) = 12*(iyr2+1-iyrbeg)`. If the
  most recent `.nc` file is e.g. June 2026, the binary's per-box arrays
  are still `NM` months long, with the trailing Jul–Dec months filled
  with `9999`.
- **`.npz` month count matches the binary file's padded length.**
  `meta[0]`/`meta[3]` equal `NM_padded = 12*(iyr2+1-1880)` (padded to
  December of the last year), same as the binary file's `INFO(4)`. Each
  box's `arr[i][1]` is padded to that same length and filled with `9999.`
  for any individual missing month — *except* boxes with zero valid
  months anywhere in their whole history, which are stored as the single
  value `[9999.0]` instead of a full padded array; that box's own
  `arr[i][0][0]` is `1` in that case, and `NM_padded` otherwise (i.e.
  `arr[i][0][0]` always equals `len(arr[i][1])` for that same box — this
  is a simpler, self-referential convention, not `trimSBBX.f`'s "next
  box's marker" trick, which is specific to the compact binary format).
  The rest of `arr[i][0]` (`lts, ltn, lnw, lne, nst, nstmn, dmin`) is a
  plain Python list, not a numpy array.
- **Title.** Same fixed format as `rearrange_ERSST.f`: `"Monthly Sea
  Surface Temperature anom (C) ERSSTv6 01/1880 - mm/yyyy"` (the "anom" is
  a holdover from the original code, kept as-is per your instructions).
- Tested end-to-end against a synthetic dataset (fake `.nc` files +
  synthetic mask) exercising: multi-year ranges, a partial final year
  (June cutoff) vs. a clean December cutoff, fully-missing boxes (compact
  binary encoding + all-`9999` npz rows), parallel vs. serial `.nc`
  reading (byte-identical output), and `--skip-sbbx-binary`.
