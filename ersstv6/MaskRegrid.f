!** gfortran -fconvert=big-endian -frecord-marker=4 MaskRegrid.f -o MaskRegrid.exe
      real er(180,89),mask(180,89,12)

      real er_ea(8000)

      character*80 title
      character*31 :: fileer='input_files/ersst.v6.yyyymm.bin'

!**** ERSST open ocean mask (2x2 deg grid - -88->+88N)
      open(1,file='input_files/ERSST_open_ocean_mask.bin',
     *     form='unformatted')
      read(1) title,mask
      close (1)

!**** get time range and other parameters
      if (iargc() < 4) then
        write(0,*) 'usage: mkEA_SST yr1 mo1 yr2 mo2 [ fr_ea'
        stop
      end if

      call getarg(1,title) ; read(title,*) iy1
      call getarg(2,title) ; read(title,*) mo1
      call getarg(3,title) ; read(title,*) iy2
      call getarg(4,title) ; read(title,*) mo2

      frac=.50 ! fraction of overlay with good data to be counted
      if(iargc()>4) then
        call getarg(5,title) ; read(title,*) frac
      end if

!**** prepare interpolation to equal area grid (for years before 1982)
      call ntrp2ea0 (180,89, 89.5,90., 9999.)

!**** Final output file on equal-area grid
      open(3,file='ERdSST_monthly',form='unformatted')

!**** Read requested monthly ERSST records
      m1 = mo1 ; m2 = 12

      do iy = iy1,iy2
        if(iy == iy2) m2=mo2
        do m=m1,m2
          write(fileer(22:27),'(i4,i2.2)') iy,m
          open(1,file=fileer,form='unformatted',status='old',err=991)
          read(1) title,er
          close (1)

!**** interpolate to equal area grid
          call ntrp2ea(mask(1,1,m),er,er_ea,frac)
          write(3) title,er_ea
        end do
        m1 = 1
      end do

      stop
  991 write(0,*) 'could not find file ',fileer
      stop
      end

C**** NTRP2EA.f   Horizontal Interpolation Program REAL*4   2009/02/03
C****
      SUBROUTINE NTRP2EA0 (INA,JNA,OFFIA,DIVJA, SKIB)

C****
C**** NTRPEA performs a horizontal interpolation of per unit area or per
C**** unit mass quantities defined on grid A, calculating the quantity
C**** on grid B.  B grid values that cannot be calculated because the
C**** covering A grid boxes have WTA = 0, are set to the value of SKIP.
C**** The area weighted integral of the quantity is conserved.
C**** In this version the output grid (B) is Sergei's equal area grid
C****
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (TWOPI=6.283185307179586477d0)
      REAL*4 WTA(*),A(*),B(8000),
     *       OFFIA,DIVJA, SKIB,SKIP,frac,total_wt
      COMMON /NTRPEAcb/ SINA(0:361),SINB(0:80),
     *       FMIN(160,4),FMAX(160,4),GMIN(81),GMAX(81),
     *       IMIN(160,4),IMAX(160,4),JMIN(81),JMAX(81),nbefor(81)
      DATA IMA,JMA/2*0/, SKIP/0/
      save ima,jma,jmb
C****
      IMA = INA
      JMA = JNA
C     IMB = 40*izone(iband) izone=1,2,3,4,4,3,2,1
      JMB = 80
      SKIP = SKIB

      IF(IMA.lt.1 .or. IMA.gt.720 .or. JMA.lt.1 .or. JMA.gt.361)
     *   GO TO 400
C****
C**** Partitions in the I direction
C**** RIA = longitude in degrees of right edge of grid box IA on grid A
C**** RIB = longitude in degrees of right edge of grid box IB of grid B
C**** IMIN(IB) = box on grid A containing left edge of box IB on B
C**** IMAX(IB) = box on grid A containing right edge of box IB on B
C**** FMIN(IB) = fraction of box IMIN(IB) on A that is left of box IB
C**** FMAX(IB) = fraction of box IMAX(IB) on A that is right of box IB
C****
      DIA = 360d0/IMA
      do izone=1,4
      IMB = 40*izone
      DIB = 360d0/IMB
      IA  = 1
      RIA = (IA+OFFIA)*DIA - 360
      IB  = IMB
      DO 150 IBP1=1,IMB
      RIB = (IBP1-1)*DIB   ! offIB=0. for EA-grid, divJB irrelevant
  110 IF(RIA-RIB)  120,130,140
  120 IA  = IA  + 1
      RIA = RIA + DIA
      GO TO 110
C**** Right edge of A box IA and right edge of B box IB coincide
  130 IMAX(IB,izone) = IA
      FMAX(IB,izone) = 0.
      IA  = IA  + 1
      RIA = RIA + DIA
      IMIN(IBP1,izone) = IA
      FMIN(IBP1,izone) = 0.
      GO TO 150
C**** A box IA contains right edge of B box IB
  140 IMAX(IB,izone) = IA
      FMAX(IB,izone) = (RIA-RIB)/DIA
      IMIN(IBP1,izone) = IA
      FMIN(IBP1,izone) = 1.-FMAX(IB,izone)
  150 IB = IBP1
      IMAX(IMB,izone) = IMAX(IMB,izone) + IMA
!       WRITE (0,*) 'Zone=',izone,' IMA=',IMA
!       WRITE (0,*) 'IMIN=',(IMIN(I,izone),I=1,IMB)
!       WRITE (0,*) 'IMAX=',(IMAX(I,izone),I=1,IMB)
!       WRITE (0,*) 'FMIN=',(FMIN(I,izone),I=1,IMB)
!       WRITE (0,*) 'FMAX=',(FMAX(I,izone),I=1,IMB)
      end do
C****
C**** Partitions in the J direction
C****
C**** RJA = latitude in radians at top edge of box JA on grid A
C**** SINA(JA) = sine of latitude of top edge of box JA on grid A
      OFFJA = (DIVJA-JMA)/2.
      DJA   = .5*TWOPI/DIVJA
      DO 220 JA=1,JMA-1
      RJA = (JA+OFFJA)*DJA - .25*TWOPI
  220 SINA(JA) = DSIN(RJA)
      SINA(0)  = -1.
      SINA(JMA)=  1.
C**** RJB = latitude in radians at top edge of box JB on grid B
C**** SINB(JB) = sine of latitude of top edge of box JB on grid B
      SINB(0)  = -1. ; sband = SINB(0) ! sine of northern edge of band
      do izone=1,4
        dsband = .1d0 * izone
        do jzs=1,10
          jb = jzs + (izone-1)*10
          SINB(JB) = sband + .1d0*jzs*dsband
        end do
        sband = sband + dsband
      end do
      SINB(40)  = 0.
      do jb=41,80
        sinb(jb) = -sinb(80-jb)
      end do
      jb=1
      nbefor(1) = 0
      do izone=1,8
        iz = izone ; if(iz>4) iz = 9-izone
        do jzs=1,10
          nbefor(jb+1) = nbefor(jb) + 40*iz
          jb=jb+1
        end do
      end do
!       WRITE (0,*) 'JB, sinb(jb)           ',0,sinb(0)
!       do jb=1,80
!         WRITE (0,*) 'Jb, sinb(jb), nbef(jb)=',jb,sinb(jb),nbefor(jb)
!       end do
!       WRITE (0,*) 'Jb,           nbef(jb)=',81,nbefor(81)
C****
C**** JMIN(JB) = index of box of A that contains bottom edge of box JB
C**** JMAX(JB) = index of box of A that contains top edge of box JB
C**** GMIN(JB) = fraction of box JMIN(JB) on A grid that is below box JB
C**** GMAX(JB) = fraction of box JMAX(JB) on A grid that is above box JB
C****
      JMIN(1) = 1
      GMIN(1) = 0.
      JA = 1
      DO 350 JB=1,JMB-1
  310 IF(SINA(JA)-SINB(JB))  320,330,340
  320 JA = JA + 1
      GO TO 310
C**** Top edge of A box JA and top edge of B box JB coincide
  330 JMAX(JB) = JA
      GMAX(JB) = 0.
      JA = JA + 1
      JMIN(JB+1) = JA
      GMIN(JB+1) = 0.
      GO TO 350
C**** A box JA contains top edge of B box JB
  340 JMAX(JB) = JA
      GMAX(JB) = SINA(JA)-SINB(JB)
      JMIN(JB+1) = JA
      GMIN(JB+1) = SINB(JB)-SINA(JA-1)
  350 CONTINUE
      JMAX(JMB) = JMA
      GMAX(JMB) = 0.
!       WRITE (0,*) 'JMIN=',(JMIN(J),J=1,JMB)
!       WRITE (0,*) 'JMAX=',(JMAX(J),J=1,JMB)
!       WRITE (0,*) 'GMIN=',(GMIN(J),J=1,JMB)
!       WRITE (0,*) 'GMAX=',(GMAX(J),J=1,JMB)
      RETURN
C****
C**** Invalid parameters or dimensions out of range
C****
  400 WRITE (0,940) IMA,JMA,OFFIA,DIVJA, SKIP
      STOP 400
  940 FORMAT ('0Arguments received by NTRPEA0 in order:'/
     *   2I12,' = IMA,JMA = array dimensions for A grid'/
     *  E24.8,' = OFFIA   = fractional number of grid boxes from',
     *                    ' IDL to left edge of grid box I=1'/
     *  E24.8,' = DIVJA   = number of whole grid boxes from SP to NP'/
     *  E24.8,' = SKIP    = value to be put in B array when B',
     *  ' grid box is subset of A grid boxes with WTA = 0'/
     *  '0These arguments are invalid or out of range.')
C****

      ENTRY NTRP2EA (WTA,A,B,frac)
C****
C**** NTRPEA performs the horizontal interpolation
C**** Input: WTA = weighting array for values on the A grid
C****          A = per unit area or per unit mass quantity
C**** Output:  B = horizontally interpolated quantity on B grid
C****
C****
C**** Interpolate the A grid onto the B grid
C****
      if(frac > .99) frac=.99 ! insure a little leeway !
      if(frac < 0.) frac=0.
      DO 520 JB=1,JMB
      izone = (jb + 9)/10
      if(izone > 4) izone=9-izone
      imb = izone*40
      JAMIN = JMIN(JB)
      JAMAX = JMAX(JB)
      DO 520 IB=1,IMB
      IJB  = IB + nbefor(JB)
      WEIGHT= 0. ; total_wt=0.
      VALUE = 0.
      IAMIN = IMIN(IB,izone)
      IAMAX = IMAX(IB,izone)
      DO 510 JA=JAMIN,JAMAX
      G = SINA(JA)-SINA(JA-1)
      IF(JA.eq.JAMIN)  G = G - GMIN(JB)
      IF(JA.eq.JAMAX)  G = G - GMAX(JB)
      DO 510 IAREV=IAMIN,IAMAX
      IA  = 1+MOD(IAREV-1,IMA)
      IJA = IA + IMA*(JA-1)
      F   = 1.
      IF(IAREV.eq.IAMIN)  F = F - FMIN(IB,izone)
      IF(IAREV.eq.IAMAX)  F = F - FMAX(IB,izone)

      IF(A(IJA) /= SKIP)  WEIGHT = WEIGHT + F*G*WTA(IJA)
      total_wt = total_wt + F*G

  510 IF(A(IJA) /= SKIP)  VALUE  = VALUE  + F*G*WTA(IJA)*A(IJA)
      B(IJB) = SKIP
      IF(WEIGHT > total_wt*frac)  B(IJB) = VALUE/WEIGHT
  520 continue
      RETURN
      END
