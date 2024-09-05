# (UNOFFICIAL) modified gistemp4.0 and utilities to generate approximate ERSSTv5

## EXPERIMENTAL: DO NOT USE! ##

**This is not at all a carbon copy of operational gistemp).**

GISTEMP 4.0 modified to run on later python, numpy (tested on python (pypy) 3.12, numpy 2.0.1).

Experimental!! Not thoroughly tested, so DO NOT USE for anything serious!

I have only examined global monthly means for LOTI (has only 1 month different of official LOTI -- off by 0.01 C).

Has utilities to generate (ersst5_to_sbbx.py), and inspect an approximate ERRSTv5 SBBX (+-0.02 C of monthly global means) from ERSTSv5 netCDF4 files (PSL/NOAA).

SBBX_to_txt.py is a python translation  of the fortran program with the same name from GISS at [https://data.giss.nasa.gov/pub/gistemp/](https://data.giss.nasa.gov/pub/gistemp/).

See images in ersst_calc_comparisons/ for comparisons of different area methods in how it approximates GISS ERSSTv5 analysis (uses global monthly LOTI means for comparison).

Some of the python utilities added to gistemp are hardcoded, so you have to place the files in the correct path.

All utilities in tools/ are meant to be run from the main directory.

## How to run

NOTE: Make sure to delete old parquets in main directory before running on new data. GHCNM data still needs to be downloaded, extracted, and renamed in tmp/input to ghcnm.tavg.qcf.dat and v4.inv respectively from the tar.gz if you are using the most up-to-date version.

Run gistemp (by itself using all GISS data, for reference). results are in tmp/result/ (mixedGLB.Ts.ERSSTV5.GHCN.CL.PA.csv is the gistemp global LOTI)
```
python3 tool/run.py
```

Get the source data sets (the good ERSSTv5):
```
python3 tool/run.py --steps=0
```

Convert a good SBBX to .txt (rename manually the good one to SBBX.ERSSTv5 to .orig for reference)
```
mv tmp/input/SBBX.ERSSTv5 tmp/input/SBBX.ERSSTv5.orig
python3 SBBX_to_txt.py tmp/input/SBBX.ERSSTv5.orig
```

Get the PSL/NOAA (ERSTL) sst monthly means (ERSSTv5):
```
wget -P tmp/input https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc
```

Run the conversion from the PSL/NOAA netCDF4 sst (ERSL) monthly means to get the approximate SBBX (this WILL be different than what GISS provides):
```
python3 tool/ersst5_to_sbbx.py
```

Run gistemp again (it won't overwrite or clobber files already existing such as the newly created SBBX)
```
python3 tool/run.py
```
