!** gfortran -fconvert=big-endian -frecord-marker=4 test.f -o test
      implicit none
      character*80 title
      character*31 :: fileer='input_files/ersst.v6.188001.bin'
      real er(180,89)
      open(1,file=fileer,form='unformatted',status='old',err=991)
      read(1) title,er
      close(1)
      print *, title
      print *, er
      stop
991   print *, 'Error reading file'
      stop
      end
