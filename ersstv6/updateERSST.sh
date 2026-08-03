#!/bin/bash

if [ $# -lt 4 ] ; then
echo
echo "Preliminary steps:"
echo
echo "Collect the ERSSTv6 data; they are available as text and as netcdf files."
echo "Convert them to fortran binary files and move them to ./input_files/ersst.v6.yyyymm.bin ."
echo "  Alternatively, modify MaskRegrid.f to read the original text or netcdf files."
echo "Also move the masking file to the input_files directory"
echo
echo "Compile the various <...>.f programs and call the resulting executable files <...>.exe."
echo
echo "To recreate SBBX.ERSSTv6 from scratch, type: ./updateERSST.sh 1880 1 <lastYear> <lastMonth>"
echo "Alternatively, to update an existing file, move it to ./input_files ."
echo "To add e.g. just the latest month YearL monthL, type: ./updateERSST.sh YearL monthL YearL monthL" 
echo
echo "There is a tuning parameter that is currently set to 0.5."
echo "If the data cover only part of an equal area box, it is the fraction of the box"
echo "   that needs to be covered for the data to be accepted."
exit
fi

# Access the arguments
yr1=$1
month1=$2
yr2=$3
month2=$4

# Mask ice-impacted boxes and regrid the modified data to the equal area grid
# compile MaskRegrid.f -> MaskRegrid.exe
#echo "Enter yr1 month1 yr2 month2" ; read yr1 month1 yr2 month2
./MaskRegrid $yr1 $month1 $yr2 $month2 0.5
echo "created:" ; ls -l ERdSST_monthly

# Update SBBX.ERSSTv6 creating SBBX.ERSSTv6.upd
# Compile rearrange.ERSST.f -> rearrange.ERSST.exe and
#         trimSBBX.f        -> trimSBBX.exe
./rearrange.ERSST $yr1 $month1 $yr2 $month2
echo "created:" ; ls -l SBBX.ERSST.upd
./trimSBBX  SBBX.ERSST.upd SBBX.ERSSTv6
echo "created:" ; ls -l SBBX.ERSSTv6
echo "removing the auxiiary files:"
#rm -f ERdSST_monthly SBBX.ERSST.upd ; # clean up
