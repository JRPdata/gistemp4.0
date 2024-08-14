GISTEMP
=======

- Nick Barnes, Climate Code Foundation, 2014-07-11
- David Jones, Climate Code Foundation, 2014-07-11
- Avi Persin, Revision 2016-01-06


### CONTENTS

  1. Introduction
  2. Dependencies
  3. Installation
  4. Input Data
  5. Running
  6. Results

  A. References
  B. Document history
  C. Copyright and license

## 1. INTRODUCTION

This is release 1.0 of the new gistemp code revised by Avi Persin from NASA GISS;
it is based on version 0.6.x of the Clear Climate Code GISTEMP project (ccc-gistemp).

The Clear Climate Code was a reimplementation of GISTEMP (the GISS surface
temperature analysis system), to make it clearer using a single programming
language python instead of a mixture of fortran/python/ksh/c. In addition
there were various bug fixes and improvements to enhance clarity.

The revision consisted mostly of:

    converting the entire code base from Python 2 to Python 3.

    Changes to the algorithm to reflect GISTEMP updates as described here:
    https://data.giss.nasa.gov/gistemp/updates_v3/.

    Changing the internal filesystem structure to utilize numpy instead of
    Fortran-like binaries.

    Providing additional output files such as global seasonal and monthly means.


URLs for further information:

https://data.giss.nasa.gov/gistemp/

http://clearclimatecode.org/ Clear Climate Code website and blog.

https://github.com/ClimateCodeFoundation/ccc-gistemp ccc-gistemp
code repository.


## 2. DEPENDENCIES

You need Python and a machine that can run it, and a network
connection; there are no explicit operating system or CPU architecture
dependencies so "any" operating system or CPU should be okay.

Python comes in several versions.  We recommend Python 3.4, but
gistemp should work on any version of Python from the 3.x branch.

The code should run on OS X, FreeBSD, Windows, and probably a variety of
other Unix-like operating systems.

A network connection is required to download the input files (which
need only be done once). If you use a proxy to access the internet
then Python requires that the "http_proxy" environment variable is set.
The proxy will need to handle both HTTP and FTP requests

Python may already be installed on your machine (for example, it comes
pre-installed on OS X), it may be possible to install it using your
operating system's package manager; for Windows you can download an
installer from http://www.python.org/download/ .  We recommend you use a
stable production release from the Python 3.x series.

To ensure all the necessary Python packages have been installed,
we recommend using Anaconda. Anaconda is an open source scientific
Python distribution. If you are using the standard Python distribution
you will need to install the numpy and scipy modules.


## 3. INSTALLATION

Unpack gistemp-1.0.tar.gz.


## 4. INPUT DATA

gistemp uses input data in the subdirectory tmp/input/.  This input
data includes large files (a few megabytes to a few dozen megabytes)
of temperature records from GHCN, SCAR, sea surface data, and
small files of additional temperature records and station tables from
GISS.  These files are all specified in config/sources, and there is
code in tool/fetch.py to fetch them from the originating organisations
over the internet.

Downloading the input data is a common cause of problems.  Maintaining
the part of the code that does this (which has nothing to do with the
core GISTEMP algorithm) is a significant cost.  If the tools
we provide do not seem to download the input data correctly, you can
download the data "by hand" and install it in the tmp/input/ directory.


## 5. RUNNING

To run gistemp:

    python tool/run.py

That command runs steps 0 through 5.  To run only a single step or a shorter
sequence of steps, use the -s argument.  For instance:

    python tool/run.py -s 3         # Runs just step 3
    python tool/run.py -s 0-3,5     # Runs steps 0,1,2,3,5 (omitting 4)

We use this directory structure:

gistemp-x.x.x    /steps/        Source code for the GISTEMP algorithm only
                 /config/       Configuration files
                 /doc/          Internal developer documentation
                 /tmp/input/    Input data files
                 /tmp/log/      Log files
                 /tool/         Tools - sources other than the GISTEMP algorithm
                 /tmp/work/     Intermediate data files
                 /tmp/result/   Final result files

Running the code should write to the tmp/input/ directory when fetching
input data, but subsequently only write to the tmp/work/ tmp/log/ and tmp/result/
directories.  Before running tool/run.py, these directories can all be
deleted (if you wish, for example, to have a clean run).

In 2016 a complete run takes about 13 minutes on a
Intel(R) Core(TM) i7-4810MQ CPU @ 2.80 GHz running Windows 7 64-bit
operating system.  If you want this to go much faster we recommend that
you run using PyPy (an alternate implementation of
Python http://codespeak.net/pypy/dist/pypy/doc/ ).



## 6. RESULTS

After running run.py, the GISTEMP result files are all in the tmp/result/
directory.


## A. REFERENCES

None.


## B. DOCUMENT HISTORY

Most recent changes first:

2017-01-09 AP  Updated to prepare for 1.0.0.
2010-10-29 DRJ Updated to prepare for 0.6.1.
2010-10-22 DRJ Updated to prepare for 0.6.0.
2010-07-21 DRJ Updated to prepare for 0.5.1.
2010-07-19 DRJ Updated to prepare for 0.5.0.
2010-07-13 DRJ Added note about PyPy.
2010-03-11 DRJ Updated to prepare for 0.4.1.
2010-03-09 DRJ Updated to prepare for 0.4.0.
2010-01-26 NB  Updated to prepare for 0.3.0.
2010-01-25 DRJ Removed PNG result.
2010-01-22 NB  Updated to reflect some code moving to tool/.
2010-01-11 NB  Updated to describe preflight better.
2010-01-06 DRJ Updated for our all-Python status.
2009-12-03 NB  Updated for transfer to GoogleCode project.
2008-09-19 DRJ Added PNG result.
2008-09-13 NB  Updated for CCC 0.1.0.
2008-09-12 NB  Updated for CCC 0.0.3.
2008-09-12 NB  Updated for CCC 0.0.2.
2008-09-11 NB  Updated for CCC 0.0.1.
2008-09-08 NB  Created.
2016-01-06 Avi Revision.
2016-02-25 Michael Michael Changed versioning scheme so that gistemp1.2 is now gistemp4.0


## C. COPYRIGHT AND LICENSE

This document is copyright (C) 2009, 2010 Ravenbrook Limited; and (C)
2010 Climate Code Foundation.  All rights reserved.

Redistribution and use of this document in any form, with or without
modification, is permitted provided that redistributions of this
document retain the above copyright notice, this condition and the
following disclaimer.

THIS DOCUMENT IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDERS AND CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
DOCUMENT, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
