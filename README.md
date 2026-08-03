# (UNOFFICIAL) modified gistemp4.0 (for ERSSTv6) and utilities to generate ERSSTv6 SBBX binary

## EXPERIMENTAL: DO NOT USE! ##

**This is not at all a carbon copy of operational gistemp).**

GISTEMP 4.0 modified to run on later python, numpy (tested on python (pypy) 3.12, numpy 2.0.1).

Experimental!! Not thoroughly tested, so DO NOT USE for anything serious!

Has utilities to generate (ersst5_to_sbbx.py), and inspect the ERRSTv6 SBBX from ERSTSv6 monthly netCDF4 files (NCEI).

SBBX_to_txt.py is a python translation of the fortran program with the same name from GISS at [https://data.giss.nasa.gov/pub/gistemp/](https://data.giss.nasa.gov/pub/gistemp/).

'ersstv6' contains the fortran utilities from NASA's GISTEMP to generate the SBBX along with some new helper python tools to convert and generate them
'ersstv6py' is the pure python translation of the above (written by Claude) (use either)
For the above make sure to read the README in their respective subfolders before using (requires some editing for hardcoded paths).

All utilities in tools/ are meant to be run from the main directory.

## How to run

Assuming you are in the main gistemp folder...

You can run the regular gistemp (by itself using all GISS data, for reference). results are in tmp/result/ (mixedGLB.Ts.ERSSTV6.GHCN.CL.PA.csv is the gistemp global LOTI)
```
python3 tool/run.py
```

To get the source data sets only (the official ERSSTv6 SBBX, etc.):
```
python3 tool/run.py --steps=0
```

To get the NCEI (ERSSTv6) sst monthly means to produce the unofficial ERSSTv6 SBBX:

1. First edit tool/toolconfig.py and edit the paths (the GISTEMP_PATH and CURL_PATH).
2. Then run the downloader (it will keep a manifest to only download what changes on successive runs).
```
python3 tool/download_sst.py
```

3. Edit the toolconfig.py in the ersstv6py (pure python) or ersstv6 folder (fortran and python) (depending on which you want to use)
4. Run the conversion utility (either from ersstv6py or ersstv6)
```
python3 ersstv6py/ersst6_to_sbbx.py
```

5. Copy the SBBX to the tmp/input folder (create it if it doesn't exist):
```
mkdir -p tmp/input && cp ersstv6py/SBBX.ERSSTv6 tmp/input/
```

6. Then you can run gistemp again (it won't overwrite or clobber files already existing such as the newly created SBBX)
```
python3 tool/run.py
```
