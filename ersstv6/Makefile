FC = gfortran
#FLAGS = -fconvert=big-endian -frecord-marker=4
FLAGS = -fconvert=big-endian -frecord-marker=4 -O2 -march=x86-64 -mtune=generic -ffp-contract=off

all: MaskRegrid rearrangeERSST trimSBBX

MaskRegrid: MaskRegrid.f
	$(FC) $(FLAGS) MaskRegrid.f -o MaskRegrid

rearrangeERSST: rearrange.ERSST.f
	$(FC) $(FLAGS) rearrange.ERSST.f -o rearrange.ERSST

trimSBBX: trimSBBX.f
	$(FC) ${FLAGS} trimSBBX.f -o trimSBBX

clean:
	rm -f MaskRegrid rearrangeERSST trimSBBX
