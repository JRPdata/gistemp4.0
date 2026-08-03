!** gfortran -fconvert=big-endian -frecord-marker=4 test.f -o test
      implicit none
      real :: er(180,89)
      integer :: i, j
      character*80 :: title
      character*31 :: fileer='input_files/ersst.v6.188001.bin'
      title = 'Repeating seq'
      do i=1,180
         do j=1,89
            er(i,j) = mod(i-1,5) * 2**(mod(j-1,4))
         enddo
      enddo

      open(1,file=fileer,form='unformatted',status='replace',err=991)
      write(1) title,er
      close(1)
      print *, 'Data written to file:'
      print *, er
      stop
991   print *, 'Error writing file'
      stop
      end
