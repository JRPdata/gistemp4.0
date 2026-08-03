!! gfortran -fconvert=big-endian -frecord-marker=4 rearrange.ERSST.f
c     reformats ERdSST: time series of maps to sequence of time series
      parameter (im=8000)
      parameter (iyrbeg=1880)
      real*4, allocatable :: sst(:,:)
c
c     sst is the input SST-data real*4 array (C)
c
      integer info(8),latlon(4,8000),ij_SN_EW(8000)
c
c     land/sea array not needed - built into climatology
c
      real,allocatable :: to(:) !!! to(monx)
      character*80 title,line
      title(1:40)= 'Monthly Sea Surface Temperature anom (C)'
      title(41:80)=' ERSSTv6 01/1880 - mm/yyyy'
      info(1)=1
      info(2)=1
      info(3)=6
      info(4)=0  ! will be set later
      info(5)=0  ! will be set later
      info(6)=iyrbeg
      info(7)=9999
      info(8)=-9999

c     command line arguments: new_data range
      if (iargc()<4) then
        write(0,*) 'Usage: rearrange.ERSST y1 mo1 y2 mo2'
        stop
      end if
      call getarg(1,line) ; read(line,*) iyr1
      call getarg(2,line) ; read(line,*) mon1
      if(iyr1.lt.iyrbeg) then
        write(0,*) 'iyr1=',iyr1,' cannot be less than',iyrbeg
        stop
      end if
      call getarg(3,line) ; read(line,*) iyr2
      call getarg(4,line) ; read(line,*) mon2
      if(iyr2.lt.iyr1) then
        write(*,*) 'iyr2=',iyr2,' must be > ',iyr1
        stop
      end if

c     Using iyr1 decide whether this is an update
      if (iyr1 > iyrbeg) then !               this is an update
        open(1,file='input_files/SBBX.ERSSTv6',form='unformatted')
        read(1) info,line
        mnow=info(1)
        info(1)=1  ! input is trimmed, output is NOT trimmed/reorganized
        if(info(6).ne.iyrbeg) then
          write(0,*) 'iyrbeg-new',iyrbeg,' iyrbeg-old',info(6)
          stop
        end if
        write(*,*) 'SBBX.ERSSTv6 opened'
      end if
      bad=info(7)
c
c   Read in the SST data for recent years (unit 11)
c
      open(11,file='ERdSST_monthly',form='unformatted')
      if(info(4).lt.12*(iyr2+1-iyrbeg)) info(4)=12*(iyr2+1-iyrbeg)
      info(5)=info(4)+4  ! output has old, not trimmed format
      monx = info(4) ; allocate (to(monx))
      mnew = (iyr2-iyr1)*12 + mon2 - mon1 +1 ; allocate (sst(im,mnew))
c**** Collect a slab of new data (skip over data already collected)
      mon=mon1 ; iyr=iyr1
      do 100 m=1,mnew
   80 read(11,end=110) line,(sst(i,m),i=1,8000)
      read(line(40:80),*) iyre,moe
      if(iyre<iyr1) go to 80
      if(iyre == iyr1 .and. moe < mon1) go to 80
      if(mon.ne.moe .or. iyr.ne.iyre) then
        write(*,*) 'positioning problem: month year=',mon,moe,iyr,iyre
        stop 'positioning problem'
      end if
      mon=mon+1
      if (mon==13) then
        mon = 1
        iyr = iyr+1
      end if
  100 continue

c****
c**** Rearrange series of maps into records of time series
c****
  110 open(2,file='SBBX.ERSST.upd',form='unformatted')
      write(title(60:66),'(I2.2,''/'',I4)') moe,iyre
      write(2) info,title
      call def_ea_grid (latlon,ij_SN_EW)
      moff=(iyr1-iyrbeg)*12 + mon1 - 1

      do 500 n=1,8000
      to = bad
      if (moff>0) then
        call sread(1,to,mnow,lts,ltn,lnw,lne,mnext)
        mnow=mnext
      end if
      isbbx=ij_SN_EW(n)
      do 220 m=1,12*(iyre-iyr1)+(moe-mon1)+1
      to(m+moff)=sst(isbbx,m)
  220 continue
      do 230 m=(iyre-iyrbeg)*12+moe+1,info(4)
  230 to(m)=bad

cc       if(iw.lt.1) then
cc          write(*,*) n,'th subbox',lts,ltn,lnw,lne,'->',js,jn,iw,ie
cc          stop
cc       end if
      write(2) (to(i),i=1,info(4)),(latlon(i,isbbx),i=1,4)  ! untrimmed
  500 continue
      stop
      end

      subroutine sread(in,a,na,lts,ltn,lnw,lne,mnext)
      real a(na)
      integer info(7)
      read(in) mnext,lts,ltn,lnw,lne,x,x,x,a
      return
      end

      subroutine def_ea_grid (latlon,ij_SN_WE)
!
!**** define corners of the grid boxes in their SN->WE ordering
!**** south->north (outer) and west->east (inner) starting at date line
!                   loop                   loop
!**** Also find algorithm to write out this information in Sergei's
!**** ordering: ij_SN_WE(ij_sb) = SN_WE index for ij_sb
!
!**** boxes:                    1-4   north , west->east ...
!****                            .      to
!****                          76-80  south , west->east
!**** Order of subboxes:        1-10  south , west->east ...
!****                            .      to
!****                         91-100  north , west->east
!****
!**** Grid constants for latitude zones
!!!!! REAL BANDS(9)/-90.,-64.2,-44.4,-23.6,0.,23.6,44.4,64.2,90./


      integer latlon(4,8000),nbefor(81),LatSb(81)
      integer ij_SN_WE(8000)

      ! Set  LatSb(j), the j-th Southern latitude edge in .01 degrees
      ! and IBEFOR(j), the number of grid boxes South of LATS(j)
      xbypi = 9000d0/asin(1d0)        ! 100 * 180/pi
      j=1 ; nbefor(1)=0 ; sband = -1. ! sine of southern edge of band
      do iband=1,8
        iz = iband ; if(iz>4) iz = 9-iband
        dsband = .1d0 * iz
        do jzs=1,10
          nbefor(j+1) = nbefor(j) + 40*iz
          LatSb(J) = nint(xbypi*asin(sband + .1d0*(jzs-1)*dsband))
          j = j + 1
        end do
        sband = sband + dsband
      end do
      LatSb(81)  = -LatSb(1)   ! 9000

      do j=1,80
        latsouth = LatSb(j)
        latnorth = LatSb(j+1)
        im = nbefor(j+1) - nbefor(j)
        do i=1,im
           lonwest = -18000 + nint( (i-1)*36000d0/im )
           loneast = -18000 + nint(     i*36000d0/im )
           latlon(1,nbefor(j)+i) = latsouth
           latlon(2,nbefor(j)+i) = latnorth
           latlon(3,nbefor(j)+i) = lonwest
           latlon(4,nbefor(j)+i) = loneast
        end do
      end do

      ij=0
      do jband=1,8
        jup = jband*10 + 1
        iz = jband ; if(iz>4) iz = 9-jband
        nbefore = 8000 - nbefor(jup) ; ibefore=0
        do j=1,10
        do nb = 1, 4*iz
          nbeforj = nbefore + (nb-1)*100
          do i=1,10
            ij = ij+1
            ij_SN_WE( nbeforj + (j-1)*10 + i ) = ij
!!          write(*,*) ij,ij_SN_WE(ij)
          end do
        end do
        end do
      end do

!     do ij=1,8000
!       write(*,*) ij,latlon(1:4,ij_SN_WE(ij))
!     end do

      return
      end
